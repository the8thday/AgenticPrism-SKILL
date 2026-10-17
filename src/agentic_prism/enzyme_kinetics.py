"""Reversible enzyme inhibition: global fits of initial rates over substrate and inhibitor (0.14.2).

v = Vmax S / (Km (1 + I/Ki) + S (1 + I/(alpha Ki)))   (mixed; Kiu = alpha Ki)
  competitive    alpha -> infinity: v = Vmax S / (Km (1 + I/Ki) + S)
  uncompetitive  only the S term:   v = Vmax S / (Km + S (1 + I/Ki))
  noncompetitive alpha = 1:         v = Vmax S / ((Km + S)(1 + I/Ki))
Unweighted least squares on log-parameterized Vmax, Km, Ki (and alpha), one global fit per
declared independent experiment (shared parameters across inhibitor concentrations, as
GraphPad's global fit). Profile-F intervals (SS <= SSmin (1 + F(1, df)/df)); an open profile
withholds that parameter. The declared mechanism gives the reported Ki; the other models are
fitted only to describe compatibility (extra-sum-of-squares F against mixed, AICc), never to
choose the reported mechanism after seeing the data. Gates: declared initial-rate conditions;
tight binding (Ki below 10 x total enzyme withholds Ki: this model ignores depletion; use tight_binding);
substrate and inhibitor ranges. Summaries across experiments: geometric t intervals.

0.14.3 adds the opt-in `tight_binding` block (Morrison 1969): inhibitor titrations at one substrate
concentration per experiment, with inhibitor depletion by enzyme accounted for,
  v = v0 (1 - [EI]/Et),  [EI] = ((Et + I + Ki_app) - sqrt((Et + I + Ki_app)^2 - 4 Et I))/2,
evaluated in the cancellation-free form v = v0 (1 - 2 I/(a + sqrt((Et - I + Ki_app)^2 + 4 I Ki_app))),
a = Et + I + Ki_app. Et is the declared total enzyme (fixed) or estimated as active enzyme. Ki_app
has profile-F intervals; Ki follows from the declared mechanism (competitive Ki_app/(1 + S/Km),
uncompetitive Ki_app/(1 + Km/S), noncompetitive Ki_app; mixed reports Ki_app only) with Km declared
from a separate experiment and treated as known.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'enzyme_inhibition'}
MECHANISMS = ('competitive', 'noncompetitive', 'uncompetitive', 'mixed')
PARAMS = {'competitive': ['Vmax', 'Km', 'Ki'], 'noncompetitive': ['Vmax', 'Km', 'Ki'], 'uncompetitive': ['Vmax', 'Km', 'Ki'],
          'mixed': ['Vmax', 'Km', 'Ki', 'alpha']}
MOLAR = {'M': 1., 'mM': 1e-3, 'uM': 1e-6, 'µM': 1e-6, 'nM': 1e-9, 'pM': 1e-12}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'enzyme_inhibition', 'input': None, 'source': None,
            'design': {'enzyme': None, 'substrate_unit': None, 'inhibitor_unit': None, 'rate_unit': None, 'initial_rate_rationale': None,
                       'independent_unit': None, 'independence_source': None},
            'enzyme_concentration': {'value': None, 'unit': None, 'rationale': None},
            'mechanism': {'declared': None, 'rationale': None},
            'uncertainty': {'level': .95}, 'tight_binding': None, 'report': {'plot_style': 'prism_like'}}
TIGHT = 10.
TB_KEYS = {'model', 'enzyme', 'rationale', 'km'}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported enzyme-inhibition configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    tb = c.pop('tight_binding')
    if tb is not None:
        if not isinstance(tb, dict) or not {'model', 'enzyme', 'rationale'} <= set(tb) <= TB_KEYS: raise ValueError('tight_binding needs model, enzyme and rationale (optional km)')
        if tb['model'] != 'morrison': raise ValueError('tight_binding.model must be morrison')
        if tb['enzyme'] not in ('fixed', 'estimate'): raise ValueError('tight_binding.enzyme must be fixed (declared total enzyme is active) or estimate (fit active enzyme)')
        if not _text(tb['rationale']): raise ValueError('Declare tight_binding.rationale')
        km = tb.get('km')
        if km is not None and (not isinstance(km, dict) or set(km) != {'value', 'source'} or type(km['value']) not in (int, float) or km['value'] <= 0 or not _text(km['source'])):
            raise ValueError('tight_binding.km needs a positive value (substrate unit) and its source')
        if c['enzyme_concentration']['value'] is None: raise ValueError('Morrison analysis needs the total enzyme concentration (value and molar unit)')
        if c['mechanism']['declared'] in ('competitive', 'uncompetitive') and km is None:
            raise ValueError('Converting Ki_app to Ki for a competitive or uncompetitive inhibitor needs tight_binding.km with its source; declare mechanism mixed to report Ki_app only')
        c['tight_binding'] = tb
    if c['schema_version'] != 1 or c['analysis_type'] != 'enzyme_inhibition': raise ValueError('Unsupported enzyme-inhibition schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    m = c['mechanism']
    if m['declared'] not in MECHANISMS or not _text(m['rationale']):
        raise ValueError('Declare mechanism.declared (competitive, noncompetitive, uncompetitive or mixed) with a rationale fixed before the results; use mixed when the mechanism is not known')
    e = c['enzyme_concentration']
    if e['value'] is None:
        if not _text(e['rationale']): raise ValueError('Declare the total enzyme concentration, or why it is unknown (tight binding then cannot be excluded)')
    else:
        if type(e['value']) not in (int, float) or e['value'] <= 0 or e['unit'] not in MOLAR: raise ValueError('enzyme_concentration needs a positive value and a molar unit (M, mM, uM, nM, pM)')
        if c['design']['inhibitor_unit'] not in MOLAR: raise ValueError('With a declared enzyme concentration the inhibitor unit must be molar (M, mM, uM, nM, pM)')
    lv = c['uncertainty']['level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'experiment': str})
    need = {'experiment', 'substrate', 'inhibitor', 'rate'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded rate needs an exclusion_reason')
    else: d['exclude'] = False
    for col in ('substrate', 'inhibitor', 'rate'):
        v = pd.to_numeric(d[col], errors='coerce')
        if v.isna().any() or not np.isfinite(v).all(): raise ValueError(f'{col} must be numeric')
        d[col] = v
    if (d.substrate < 0).any() or (d.inhibitor < 0).any(): raise ValueError('Substrate and inhibitor concentrations cannot be negative')
    if d.experiment.isna().any(): raise ValueError('experiment is required on every row')
    if c.get('tight_binding'):
        for e, part in d[~d.exclude].groupby('experiment'):
            if part.substrate.nunique() != 1 or part.substrate.iloc[0] <= 0: raise ValueError(f'Experiment {e}: a Morrison titration uses one nonzero substrate concentration per experiment')
            if not (part.inhibitor == 0).any(): raise ValueError(f'Experiment {e} has no uninhibited (inhibitor 0) rates')
            if part.inhibitor.nunique() < 7: raise ValueError(f'Experiment {e} needs at least 6 nonzero inhibitor concentrations plus 0 for a Morrison fit')
        return d
    for e, part in d[~d.exclude].groupby('experiment'):
        if not (part.inhibitor == 0).any(): raise ValueError(f'Experiment {e} has no uninhibited (inhibitor 0) rates')
        if part.inhibitor.nunique() < 3: raise ValueError(f'Experiment {e} needs at least 3 inhibitor concentrations including 0')
        if part.loc[part.substrate > 0, 'substrate'].nunique() < 4: raise ValueError(f'Experiment {e} needs at least 4 nonzero substrate concentrations')
    return d


def rate(mech, z, s, i):
    vmax, km, ki = np.exp(z[0]), np.exp(z[1]), np.exp(z[2])
    if mech == 'competitive': return vmax*s/(km*(1 + i/ki) + s)
    if mech == 'uncompetitive': return vmax*s/(km + s*(1 + i/ki))
    if mech == 'noncompetitive': return vmax*s/((km + s)*(1 + i/ki))
    return vmax*s/(km*(1 + i/ki) + s*(1 + i/(np.exp(z[3])*ki)))


def _fit(mech, s, i, y, fixed=None, start=None):
    """Least squares on log parameters; `fixed` = (index, value) holds one parameter for profiling."""
    k = len(PARAMS[mech]); free = [j for j in range(k) if not fixed or j != fixed[0]]

    def full(zf):
        z = np.empty(k); z[free] = zf
        if fixed: z[fixed[0]] = fixed[1]
        return z

    def resid(zf): return rate(mech, full(zf), s, i) - y
    # log-parameter bounds keep exp() finite; alpha within 1e-6..1e6 (beyond that the S or Km term is negligible)
    lb = np.full(k, -60.); ub = np.full(k, 60.)
    if mech == 'mixed': lb[3], ub[3] = np.log(1e-6), np.log(1e6)
    pos = s > 0; ymax = max(float(y.max()), 1e-12); smid = float(np.median(s[pos])); imid = float(np.median(i[i > 0]))
    starts = [start] if start is not None else []
    if start is not None and mech == 'mixed':  # alpha can run off to the boundary; restart it on both sides
        starts += [np.r_[start[:3], 0.], np.r_[start[:3], np.log(1e3)]]
    for km in (() if start is not None else (smid/4, smid, smid*4)):  # profile steps warm-start from the previous step only
        for ki in (imid/4, imid, imid*4):
            z = [np.log(ymax*1.2), np.log(km), np.log(ki)] + ([0.] if mech == 'mixed' else [])
            starts.append(np.array(z))
    best = None
    for z0 in starts:
        try:
            r = optimize.least_squares(resid, np.clip(np.asarray(z0, float)[free], lb[free] + 1e-9, ub[free] - 1e-9), bounds=(lb[free], ub[free]),
                                       method='trf', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if np.isfinite(r.cost) and (best is None or r.cost < best.cost): best = r
    if best is None: return None
    return {'z': full(best.x), 'ss': float(2*best.cost), 'jac': best.jac}


def _profile(mech, s, i, y, fit, j, df, level):
    target = fit['ss']*(1 + stats.f.ppf(level, 1, df)/df); zhat = fit['z'][j]; out = []
    for direction in (-1, 1):
        prev, prev_ss, found = zhat, fit['ss'], None; step = .02; start = fit['z']
        for _ in range(400):
            cur = prev + direction*step
            f = _fit(mech, s, i, y, fixed=(j, cur), start=start)
            if f is None: break
            if f['ss'] >= target:
                lo, hi = (prev, cur) if direction > 0 else (cur, prev)
                for _ in range(60):
                    mid = (lo + hi)/2; g = _fit(mech, s, i, y, fixed=(j, mid), start=f['z'])
                    inside = g is not None and g['ss'] < target
                    if direction > 0: lo, hi = (mid, hi) if inside else (lo, mid)
                    else: lo, hi = (lo, mid) if inside else (mid, hi)
                found = (lo + hi)/2; break
            if abs(cur - zhat) > np.log(1e6): break
            prev, prev_ss, start = cur, f['ss'], f['z']; step = min(step*1.5, .5)
        out.append(None if found is None else float(np.exp(found)))
    return out


def fit_experiment(part, mech, level):
    s, i, y = part.substrate.to_numpy(float), part.inhibitor.to_numpy(float), part.rate.to_numpy(float)
    names = PARAMS[mech]; n, k = len(y), len(names); df = n - k
    f = _fit(mech, s, i, y)
    if f is None or df < 1: return {'status': 'failed', 'reason': 'fit_failed_or_no_residual_df', 'n_points': n}
    jac_rank = int(np.linalg.matrix_rank(f['jac'])) if f['jac'] is not None else 0
    pars = {}
    for j, nm in enumerate(names):
        ci = _profile(mech, s, i, y, f, j, df, level)
        pars[nm] = {'estimate': float(np.exp(f['z'][j])), 'ci': ci, 'ci_status': 'two_sided_profile_f' if None not in ci else 'open'}
    if mech == 'mixed':
        kiu = pars['Ki']['estimate']*pars['alpha']['estimate']
        pars['Kiu_derived'] = {'estimate': float(kiu), 'note': 'alpha x Ki; no separate interval'}
    return {'status': 'estimated', 'n_points': n, 'df': df, 'ss': f['ss'], 'sigma': float(np.sqrt(f['ss']/df)), 'jacobian_rank': jac_rank,
            'parameters': pars, 'z': f['z'].tolist()}


def compare_models(part):
    s, i, y = part.substrate.to_numpy(float), part.inhibitor.to_numpy(float), part.rate.to_numpy(float); n = len(y); out = {}
    fits = {m: _fit(m, s, i, y) for m in MECHANISMS}
    full = fits['mixed']
    for m, f in fits.items():
        if f is None: out[m] = {'status': 'failed'}; continue
        k = len(PARAMS[m]); aicc = n*np.log(f['ss']/n) + 2*k + 2*k*(k + 1)/max(n - k - 1, 1)
        r = {'ss': f['ss'], 'k': k, 'aicc': float(aicc)}
        if m != 'mixed' and full is not None and full['ss'] > 0:
            fs = max(f['ss'] - full['ss'], 0.)/(full['ss']/(n - 4)); r.update(f_vs_mixed=float(fs), p_vs_mixed=float(stats.f.sf(fs, 1, n - 4)))
        out[m] = r
    ok = [m for m in out if 'aicc' in out[m]]; best = min(out[m]['aicc'] for m in ok) if ok else None
    for m in ok: out[m]['delta_aicc'] = float(out[m]['aicc'] - best)
    return out


def _to_molar(v, unit): return v*MOLAR[unit]


def compute(d, c):
    if c.get('tight_binding'): return compute_morrison(d, c)
    mech = c['mechanism']['declared']; lv = c['uncertainty']['level']; u = d[~d.exclude]
    e = c['enzyme_concentration']; iu = c['design']['inhibitor_unit']; fits, failing = [], []
    for ex, part in u.groupby('experiment', sort=False):
        f = fit_experiment(part, mech, lv); f['experiment'] = str(ex); diag = []
        if f['status'] != 'estimated':
            f['reportable'] = False; fits.append(f); failing.append({'experiment': str(ex), 'reason': f['reason']}); continue
        p = f['parameters']; km, ki = p['Km']['estimate'], p['Ki']['estimate']
        smax, smin = float(part.substrate.max()), float(part.loc[part.substrate > 0, 'substrate'].min()); imax = float(part.inhibitor.max())
        if smax < 2*km: diag.append('substrate_range_below_2_Km')
        if smin > km: diag.append('lowest_substrate_above_Km')
        if imax < ki: diag.append('inhibitor_range_below_Ki')
        if f['jacobian_rank'] < len(PARAMS[mech]): diag.append('parameters_not_identifiable')
        cmp = compare_models(part); f['model_comparison'] = cmp
        declared_rejected = mech != 'mixed' and cmp.get(mech, {}).get('p_vs_mixed', 1.) < .05
        if declared_rejected: diag.append('declared_mechanism_rejected_against_mixed')
        tight = None
        if e['value'] is not None:
            tight = _to_molar(ki, iu) < TIGHT*_to_molar(e['value'], e['unit']); f['ki_over_enzyme'] = float(_to_molar(ki, iu)/_to_molar(e['value'], e['unit']))
            if tight: diag.append('tight_binding_Ki_below_10x_enzyme')
        for nm, pr in p.items():
            if nm == 'Kiu_derived': pr['reportable'] = False; continue
            pr['reportable'] = bool(pr.get('ci_status') == 'two_sided_profile_f' and 'parameters_not_identifiable' not in diag
                                    and not (nm in ('Ki', 'alpha') and (declared_rejected or tight)))
        f['diagnostics'] = diag; f['reportable'] = p['Ki']['reportable']
        if not f['reportable']:
            why = ','.join(x for x in diag if x in ('declared_mechanism_rejected_against_mixed', 'tight_binding_Ki_below_10x_enzyme', 'parameters_not_identifiable')) or 'profile_interval_open'
            failing.append({'experiment': str(ex), 'reason': 'Ki_withheld:' + why})
        fits.append(f)
    summaries = {}
    est = [f for f in fits if f.get('status') == 'estimated']
    for nm in PARAMS[mech]:
        vals = [f['parameters'][nm] for f in est]
        if len(fits) < 2: summaries[nm] = {'status': 'withheld', 'reason': 'fewer_than_two_experiments'}; continue
        if len(vals) < len(fits) or not all(v['reportable'] for v in vals): summaries[nm] = {'status': 'withheld', 'reason': 'a_per_experiment_estimate_was_withheld'}; continue
        lg = np.log10([v['estimate'] for v in vals]); m = float(lg.mean()); h = float(stats.t.ppf((1+lv)/2, len(lg)-1)*lg.std(ddof=1)/np.sqrt(len(lg)))
        summaries[nm] = {'status': 'estimated', 'n_experiments': len(lg), 'geometric_mean': 10**m, 'ci': [10**(m-h), 10**(m+h)]}
    units = {'Vmax': c['design']['rate_unit'], 'Km': c['design']['substrate_unit'], 'Ki': iu, 'alpha': 'dimensionless'}
    notes = [f"Mechanism declared before the results: {mech} ({c['mechanism']['rationale']}). Ki comes from that model; the other mechanisms were fitted only to describe compatibility (F against the mixed model, AICc) and were not used to choose the reported mechanism.",
             'Global unweighted least squares across inhibitor concentrations within each independent experiment (' + c['design']['independent_unit'] + '; ' + c['design']['independence_source'] + '); profile-F intervals assume independent, equal-variance Gaussian errors in the rates.',
             'Rates are assumed to be initial rates (' + c['design']['initial_rate_rationale'] + ') with reversible, rapid-equilibrium inhibition; time-dependent or irreversible inhibition needs a different analysis.',
             'Not rejecting a mechanism against the mixed model does not establish it; competitive and uncompetitive models sit at the boundary of the mixed model, so their F tests are approximate.']
    if e['value'] is None:
        notes.append('Total enzyme concentration was not declared (' + e['rationale'] + '); tight binding (Ki near the enzyme concentration, where free inhibitor is depleted) cannot be excluded.')
    if any('tight_binding_Ki_below_10x_enzyme' in f.get('diagnostics', []) for f in fits):
        notes.append('Ki is below 10 x the total enzyme concentration in at least one experiment: inhibitor depletion makes the standard equations invalid; Ki was withheld there. Analyse a titration at one substrate concentration with the tight_binding (Morrison) option instead.')
    if any('declared_mechanism_rejected_against_mixed' in f.get('diagnostics', []) for f in fits):
        notes.append('The declared mechanism fitted significantly worse than the mixed model (F test, p < 0.05) in at least one experiment; Ki from the declared model was withheld there. Do not switch to the best-fitting model and report it as if declared.')
    for flag, text in (('substrate_range_below_2_Km', 'The highest substrate concentration is below 2 x Km in some experiment; Vmax and Km are poorly determined.'),
                       ('lowest_substrate_above_Km', 'The lowest nonzero substrate concentration is above Km in some experiment; Km is poorly determined.'),
                       ('inhibitor_range_below_Ki', 'The highest inhibitor concentration is below Ki in some experiment; Ki is an extrapolation.')):
        if any(flag in f.get('diagnostics', []) for f in fits): notes.append(text)
    if failing: notes.append('Withheld: ' + '; '.join(f"{x['experiment']} ({x['reason']})" for x in failing) + '.')
    primary = {'status': 'estimated' if any(f.get('reportable') for f in fits) else 'withheld', 'reportable': any(f.get('reportable') for f in fits),
               'mechanism': mech, 'confidence_level': lv, 'units': units, 'fits': fits, 'summaries': summaries,
               'enzyme_concentration': e}
    return {'analysis_type': 'enzyme_inhibition', 'primary': primary, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['Reversible rapid-equilibrium inhibition of a single-substrate enzyme by one inhibitor; no tight binding in this model (see the tight_binding Morrison option), and no time-dependent, irreversible, partial or multi-substrate mechanisms.',
                            'Unweighted least squares only; no substrate inhibition or progress-curve analysis.',
                            'Model comparison is descriptive; mechanism is a declaration. No Prism equivalence beyond the recorded benchmark.']}


def morrison(z, i, et_fixed=None):
    """Rate at inhibitor i; z = log(v0), log(Ki_app)[, log(Et)]. Cancellation-free bound fraction."""
    v0, k = np.exp(z[0]), np.exp(z[1]); et = np.exp(z[2]) if et_fixed is None else et_fixed
    a = et + i + k; return v0*(1 - 2*i/(a + np.sqrt((et - i + k)**2 + 4*i*k)))


def _mfit(i, y, et, estimate_et, fixed=None, start=None):
    k = 3 if estimate_et else 2; free = [j for j in range(k) if not fixed or j != fixed[0]]

    def full(zf):
        z = np.empty(k); z[free] = zf
        if fixed: z[fixed[0]] = fixed[1]
        return z

    def resid(zf): return morrison(full(zf), i, None if estimate_et else et) - y
    lb, ub = np.full(k, -60.), np.full(k, 60.)
    if start is not None: starts = [start]
    else:
        v0 = max(float(y[i == 0].mean()), 1e-12); pos = i[i > 0]; starts = []
        for kk in np.geomspace(max(pos.min(), et*1e-3), max(pos.max(), et), 6):
            starts.append(np.array([np.log(v0), np.log(kk)] + ([np.log(et)] if estimate_et else [])))
            if estimate_et: starts += [np.r_[np.log(v0), np.log(kk), np.log(et/3)], np.r_[np.log(v0), np.log(kk), np.log(et*3)]]
    best = None
    for z0 in starts:
        try:
            r = optimize.least_squares(resid, np.clip(np.asarray(z0, float)[free], lb[free] + 1e-9, ub[free] - 1e-9), bounds=(lb[free], ub[free]),
                                       method='trf', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if np.isfinite(r.cost) and (best is None or r.cost < best.cost): best = r
    if best is None: return None
    return {'z': full(best.x), 'ss': float(2*best.cost), 'jac': best.jac}


def _mprofile(i, y, et, est, fit, j, df, level):
    target = fit['ss']*(1 + stats.f.ppf(level, 1, df)/df); zhat = fit['z'][j]; out = []
    for direction in (-1, 1):
        prev, found, step, start = zhat, None, .02, fit['z']
        for _ in range(400):
            cur = prev + direction*step; f = _mfit(i, y, et, est, fixed=(j, cur), start=start)
            if f is None: break
            if f['ss'] >= target:
                lo, hi = (prev, cur) if direction > 0 else (cur, prev)
                for _ in range(60):
                    mid = (lo + hi)/2; g = _mfit(i, y, et, est, fixed=(j, mid), start=f['z'])
                    inside = g is not None and g['ss'] < target
                    if direction > 0: lo, hi = (mid, hi) if inside else (lo, mid)
                    else: lo, hi = (lo, mid) if inside else (mid, hi)
                found = (lo + hi)/2; break
            if abs(cur - zhat) > np.log(1e6): break
            prev, start = cur, f['z']; step = min(step*1.5, .5)
        out.append(None if found is None else float(np.exp(found)))
    return out


def ki_factor(mech, s, km):
    """Ki = Ki_app / factor; None when Ki is not identifiable from one substrate concentration."""
    if mech == 'noncompetitive': return 1.
    if mech == 'competitive': return 1 + s/km
    if mech == 'uncompetitive': return 1 + km/s
    return None


def compute_morrison(d, c):
    mech = c['mechanism']['declared']; lv = c['uncertainty']['level']; u = d[~d.exclude]; tb = c['tight_binding']; est = tb['enzyme'] == 'estimate'
    e = c['enzyme_concentration']; iu = c['design']['inhibitor_unit']; et = _to_molar(e['value'], e['unit'])/MOLAR[iu]
    km = tb['km']['value'] if tb.get('km') else None; fits, failing = [], []
    for ex, part in u.groupby('experiment', sort=False):
        i, y = part.inhibitor.to_numpy(float), part.rate.to_numpy(float); s = float(part.substrate.iloc[0]); n = len(y); k = 3 if est else 2; df = n - k
        f = _mfit(i, y, et, est)
        if f is None or df < 1:
            fits.append({'experiment': str(ex), 'status': 'failed', 'reason': 'fit_failed_or_no_residual_df', 'reportable': False, 'n_points': n})
            failing.append({'experiment': str(ex), 'reason': 'fit_failed_or_no_residual_df'}); continue
        names = ['v0', 'Ki_app'] + (['Et_active'] if est else []); diag = []
        rank = int(np.linalg.matrix_rank(f['jac']))
        if rank < k: diag.append('parameters_not_identifiable')
        pars = {}
        for j, nm in enumerate(names):
            ci = _mprofile(i, y, et, est, f, j, df, lv)
            pars[nm] = {'estimate': float(np.exp(f['z'][j])), 'ci': ci, 'ci_status': 'two_sided_profile_f' if None not in ci else 'open'}
            pars[nm]['reportable'] = bool(pars[nm]['ci_status'] == 'two_sided_profile_f' and not diag)
        etu = pars['Et_active']['estimate'] if est else et; kapp = pars['Ki_app']['estimate']; ic50 = kapp + etu/2
        pos = i[i > 0]
        if pos.max() < 2*ic50: diag.append('inhibitor_range_below_2x_IC50')
        if pos.min() > ic50: diag.append('lowest_inhibitor_above_IC50')
        if kapp > 100*etu: diag.append('Ki_app_above_100x_enzyme_depletion_negligible')
        if kapp < etu/100: diag.append('Ki_app_below_enzyme_over_100_titration_regime')
        if est and not .5 <= etu/et <= 2: diag.append('active_enzyme_differs_from_declared_by_more_than_2x')
        fac = ki_factor(mech, s, km)
        if fac is not None:
            ci = pars['Ki_app']['ci']
            pars['Ki'] = {'estimate': kapp/fac, 'ci': [None if x is None else x/fac for x in ci], 'ci_status': pars['Ki_app']['ci_status'],
                          'reportable': pars['Ki_app']['reportable'], 'conversion': {'competitive': 'Ki_app/(1 + S/Km)', 'uncompetitive': 'Ki_app/(1 + Km/S)', 'noncompetitive': 'Ki = Ki_app'}[mech]}
        f_out = {'experiment': str(ex), 'status': 'estimated', 'n_points': n, 'df': df, 'ss': f['ss'], 'sigma': float(np.sqrt(f['ss']/df)), 'jacobian_rank': rank,
                 'substrate': s, 'enzyme_total_declared': et, 'enzyme_treatment': tb['enzyme'], 'parameters': pars, 'IC50_derived': {'estimate': ic50, 'note': 'Ki_app + Et/2; no interval'},
                 'z': f['z'].tolist(), 'diagnostics': diag}
        f_out['reportable'] = pars['Ki_app']['reportable']
        if not f_out['reportable']:
            failing.append({'experiment': str(ex), 'reason': 'Ki_app_withheld:' + ('parameters_not_identifiable' if 'parameters_not_identifiable' in diag else 'profile_interval_open')})
        fits.append(f_out)
    summaries = {}; est_fits = [f for f in fits if f.get('status') == 'estimated']
    keys = ['Ki_app'] + (['Ki'] if ki_factor(mech, 1., km or 1.) is not None else []) + (['Et_active'] if est else [])
    for nm in keys:
        if len(fits) < 2: summaries[nm] = {'status': 'withheld', 'reason': 'fewer_than_two_experiments'}; continue
        if nm == 'Ki_app' and len({f.get('substrate') for f in est_fits}) > 1: summaries[nm] = {'status': 'withheld', 'reason': 'Ki_app_depends_on_substrate_and_experiments_differ'}; continue
        vals = [f['parameters'][nm] for f in est_fits]
        if len(vals) < len(fits) or not all(v['reportable'] for v in vals): summaries[nm] = {'status': 'withheld', 'reason': 'a_per_experiment_estimate_was_withheld'}; continue
        lg = np.log10([v['estimate'] for v in vals]); m = float(lg.mean()); h = float(stats.t.ppf((1+lv)/2, len(lg)-1)*lg.std(ddof=1)/np.sqrt(len(lg)))
        summaries[nm] = {'status': 'estimated', 'n_experiments': len(lg), 'geometric_mean': 10**m, 'ci': [10**(m-h), 10**(m+h)]}
    units = {'v0': c['design']['rate_unit'], 'Ki_app': iu, 'Ki': iu, 'Et_active': iu, 'Km': c['design']['substrate_unit']}
    notes = [f"Morrison tight-binding analysis ({tb['rationale']}): the fit accounts for inhibitor bound to enzyme, so it remains valid when Ki_app is near or below the enzyme concentration. Total enzyme {e['value']:g} {e['unit']} "
             + ("was estimated as active enzyme (Et_active), together with v0 and Ki_app." if est else "was fixed as fully active; if part of the enzyme is inactive, Ki_app is biased (estimate active enzyme or titrate it)."),
             f"Mechanism declared before the results: {mech} ({c['mechanism']['rationale']}). "
             + ({'competitive': f"Ki = Ki_app/(1 + S/Km) with Km {km:g} {c['design']['substrate_unit']} ({tb['km']['source']}) treated as known; its uncertainty is not propagated." if km else '',
                 'uncompetitive': f"Ki = Ki_app/(1 + Km/S) with Km {km:g} {c['design']['substrate_unit']} ({tb['km']['source']}) treated as known; its uncertainty is not propagated." if km else '',
                 'noncompetitive': 'For a noncompetitive inhibitor Ki equals Ki_app at any substrate concentration.',
                 'mixed': 'A mixed inhibitor\'s Ki and alpha are not identifiable from one substrate concentration; only Ki_app (specific to this substrate concentration) is reported.'}[mech]),
             'Profile-F intervals assume independent, equal-variance Gaussian rate errors, rapid-equilibrium reversible binding and initial rates (' + c['design']['initial_rate_rationale'] + '). A single substrate concentration cannot establish the mechanism.',
             'IC50 in a tight-binding titration depends on the enzyme concentration (IC50 = Ki_app + Et/2); it is shown for orientation only.']
    for flag, text in (('inhibitor_range_below_2x_IC50', 'The highest inhibitor concentration is below twice the IC50 in some experiment; the lower part of the curve is not observed.'),
                       ('lowest_inhibitor_above_IC50', 'The lowest nonzero inhibitor concentration is above the IC50 in some experiment; the upper part of the curve is not observed.'),
                       ('Ki_app_above_100x_enzyme_depletion_negligible', 'Ki_app is more than 100 x the enzyme concentration in some experiment; depletion is negligible there and the Morrison fit reduces to the classical one.'),
                       ('Ki_app_below_enzyme_over_100_titration_regime', 'Ki_app is below 1/100 of the enzyme concentration in some experiment (titration regime); the curve mainly reports the enzyme concentration and Ki_app is poorly determined. Lower the enzyme concentration if possible.'),
                       ('active_enzyme_differs_from_declared_by_more_than_2x', 'The estimated active enzyme differs from the declared total by more than 2-fold in some experiment; check enzyme activity and the declared concentration.')):
        if any(flag in f.get('diagnostics', []) for f in fits): notes.append(text)
    if failing: notes.append('Withheld: ' + '; '.join(f"{x['experiment']} ({x['reason']})" for x in failing) + '.')
    primary = {'status': 'estimated' if any(f.get('reportable') for f in fits) else 'withheld', 'reportable': any(f.get('reportable') for f in fits),
               'model': 'morrison', 'mechanism': mech, 'confidence_level': lv, 'units': units, 'fits': fits, 'summaries': summaries, 'enzyme_concentration': e, 'tight_binding': tb}
    return {'analysis_type': 'enzyme_inhibition', 'primary': primary, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['Morrison quadratic for one inhibitor binding one enzyme site at rapid equilibrium; no slow-binding (time-dependent) onset, irreversible inhibition or multiple binding sites.',
                            'One substrate concentration per experiment: mechanism is a declaration and Km is an external input without propagated uncertainty.',
                            'Unweighted least squares; no Prism equivalence beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    d = pd.read_csv(Path(run)/'input.csv', dtype={'experiment': str})
    if 'exclude' in d: d = d[~d.exclude.astype(str).str.lower().eq('true')]
    fits = [f for f in p['fits'] if f.get('status') == 'estimated'][:6]
    if not fits: return []
    if p.get('model') == 'morrison':
        with plt.rc_context(pstyle.rc(token, *font_setup())):
            fig, axes = plt.subplots(1, len(fits), figsize=(3.4*len(fits), 3), layout='constrained', squeeze=False)
            for ax, f in zip(axes.ravel(), fits):
                part = d[d.experiment == f['experiment']]; pos = part[part.inhibitor > 0]; col = pstyle.color(token, 0)
                ax.scatter(pos.inhibitor, pos.rate, s=12, color=col); ax.set_xscale('log')
                v0 = part[part.inhibitor == 0]; ax.axhline(v0.rate.mean(), ls=':', color=pstyle.MUTED, lw=.8)
                ii = np.geomspace(pos.inhibitor.min(), pos.inhibitor.max(), 200)
                ax.plot(ii, morrison(np.array(f['z']), ii, None if f['enzyme_treatment'] == 'estimate' else f['enzyme_total_declared']), color=pstyle.color(token, 1), lw=1.2)
                ax.set(xlabel=f"Inhibitor ({p['units']['Ki_app']})", ylabel=f"Rate ({p['units']['v0']})", title=f"{f['experiment']}: Morrison (S = {f['substrate']:g})")
            pstyle.save(fig, out/'enzyme_inhibition', token, 300)
        return ['enzyme_inhibition']
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, axes = plt.subplots(1, len(fits), figsize=(3.4*len(fits), 3), layout='constrained', squeeze=False)
        for ax, f in zip(axes.ravel(), fits):
            part = d[d.experiment == f['experiment']]
            for k, (iv, g) in enumerate(part.groupby('inhibitor')):
                col = pstyle.color(token, k % 8); ax.scatter(g.substrate, g.rate, s=12, color=col, label=f"I = {iv:g}")
                ss = np.linspace(0, part.substrate.max(), 200); ax.plot(ss, rate(p['mechanism'], np.array(f['z']), ss, iv), color=col, lw=1)
            ax.set(xlabel=f"Substrate ({p['units']['Km']})", ylabel=f"Rate ({p['units']['Vmax']})", title=f"{f['experiment']}: {p['mechanism']} (global fit)")
            ax.legend(frameon=False, fontsize=7, title=f"Inhibitor ({p['units']['Ki']})", title_fontsize=7)
        pstyle.save(fig, out/'enzyme_inhibition', token, 300)
    return ['enzyme_inhibition']
