"""Population PK with linear, Michaelis-Menten or QSS-TMDD elimination, estimated by R nlmixr2 (0.15.1).

nlmixr2 FOCEi estimates every value through the fixed bridge script `r_scripts_pop/pop_runner.R`
(model templates, whitelisted parameter names). This module checks the declared contract, builds
the event table, applies identifiability gates to nlmixr2's output and writes the report. The
structural model, random effects, residual model, fixed parameters and starting values are
declared before fitting; comparison models are fitted only when declared and are descriptive.

Gates (fixed here, not tuned per study):
- covariance step: CIs need nlmixr2's "r,s" sandwich; "r" or "s" alone is a flagged fallback;
  none withholds every interval (estimates remain audit values);
- |correlation| > 0.95 between two estimated fixed effects withholds both;
- nonlinearity: Km (Michaelis-Menten) or KSS (QSS) outside the observed concentration range
  withholds the nonlinear parameters (the data do not show saturation);
- ETA shrinkage > 30% withholds individual (EBE-based) parameter claims for that effect;
- natural-scale RSE (100 sqrt(exp(SE_log^2) - 1)) > 100% withholds a parameter as not identified;
  50-100% flags it as imprecise; a withheld nonlinear parameter withholds every nonlinear one.
"""
from copy import deepcopy
import math
import numpy as np
import pandas as pd

TYPES = {'population_pk'}
MODELS = {'linear_2cmt': ('cl', 'v1', 'q', 'v2'),
          'linear_mm_2cmt': ('cl', 'v1', 'q', 'v2', 'vmax', 'km'),
          'qss_tmdd_2cmt': ('cl', 'v1', 'q', 'v2', 'kss', 'kint', 'kdeg', 'r0')}
NONLINEAR = {'linear_mm_2cmt': ('vmax', 'km'), 'qss_tmdd_2cmt': ('kss', 'kint', 'kdeg', 'r0')}
ROUTES = ('iv_bolus', 'iv_infusion', 'sc')
CORR_LIMIT, SHRINK_LIMIT, RSE_FLAG, RSE_WITHHOLD = .95, 30., 50., 100.
DEFAULTS = {'schema_version': 1, 'analysis_type': 'population_pk', 'input': None, 'source': None,
            'design': {'analyte': None, 'matrix': None, 'species': None, 'subject_unit': None, 'route': None,
                       'dose_unit': None, 'concentration_unit': None, 'time_unit': None},
            'model': {'structural': None, 'rationale': None, 'observed_species': None, 'initial_estimates': None,
                      'initial_estimates_source': None, 'fixed': None, 'fixed_source': None, 'bioavailability': None,
                      'random_effects': None, 'random_effect_initial_variance': .1, 'residual': None, 'residual_initial': None},
            'estimation': {'method': 'focei', 'seed': None, 'max_outer_iterations': 5000},
            'blq': None, 'comparison_models': None, 'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def _pos(v):
    return type(v) in (int, float) and math.isfinite(v) and v > 0


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported population_pk configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict) and k in ('design', 'model', 'estimation', 'report'):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'population_pk': raise ValueError('Unsupported population_pk schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + (name if name in ('input', 'source') else 'design.' + name))
    route = c['design']['route']
    if route not in ROUTES: raise ValueError('design.route must be iv_bolus, iv_infusion or sc')
    m = c['model']; s = m['structural']
    if s not in MODELS: raise ValueError('model.structural must be linear_2cmt, linear_mm_2cmt or qss_tmdd_2cmt')
    if not _text(m['rationale']): raise ValueError('Declare model.rationale: the structural model is chosen before fitting, not by fit')
    if s == 'qss_tmdd_2cmt':
        if m['observed_species'] not in ('total', 'free'): raise ValueError('qss_tmdd_2cmt needs model.observed_species total or free (what the assay measures)')
    elif m['observed_species'] is not None: raise ValueError('observed_species applies only to qss_tmdd_2cmt')
    params = MODELS[s] + (('ka',) if route == 'sc' else ())
    fixed = m['fixed'] or {}
    if not isinstance(fixed, dict) or set(fixed) - set(params) or not all(_pos(v) for v in fixed.values()):
        raise ValueError('model.fixed must map model parameters to positive values')
    if fixed and not _text(m['fixed_source']): raise ValueError('Declare model.fixed_source for fixed parameters')
    init = m['initial_estimates']
    if not isinstance(init, dict) or set(init) != set(params) - set(fixed) or not all(_pos(v) for v in init.values()):
        raise ValueError('model.initial_estimates must give a positive starting value for every estimated parameter: ' + ', '.join(sorted(set(params) - set(fixed))))
    if not _text(m['initial_estimates_source']): raise ValueError('Declare model.initial_estimates_source')
    if route == 'sc':
        b = m['bioavailability']
        if not isinstance(b, dict) or set(b) != {'value', 'source'} or not _pos(b['value']) or b['value'] > 1 or not _text(b['source']):
            raise ValueError('SC dosing needs model.bioavailability {value in (0, 1], source}; F is not identifiable from SC data alone (use 1 to estimate apparent CL/F, V/F)')
    elif m['bioavailability'] is not None: raise ValueError('bioavailability applies only to SC dosing')
    re = m['random_effects']
    if not isinstance(re, list) or not re or set(re) - set(params) or set(re) & set(fixed) or len(set(re)) != len(re):
        raise ValueError('model.random_effects must list distinct estimated parameters')
    if not _pos(m['random_effect_initial_variance']): raise ValueError('random_effect_initial_variance must be positive')
    if m['residual'] not in ('proportional', 'combined'): raise ValueError('model.residual must be proportional or combined')
    ri = m['residual_initial']; need = {'proportional'} | ({'additive'} if m['residual'] == 'combined' else set())
    if not isinstance(ri, dict) or set(ri) != need or not all(_pos(v) for v in ri.values()): raise ValueError('residual_initial needs ' + ' and '.join(sorted(need)))
    e = c['estimation']
    if e['method'] != 'focei': raise ValueError('estimation.method must be focei')
    if type(e['seed']) is not int: raise ValueError('Declare estimation.seed (integer)')
    if type(e['max_outer_iterations']) is not int or e['max_outer_iterations'] < 100: raise ValueError('max_outer_iterations must be an integer of at least 100')
    if c['blq'] is not None and (not isinstance(c['blq'], dict) or set(c['blq']) != {'rule', 'rationale'} or c['blq']['rule'] != 'drop' or not _text(c['blq']['rationale'])):
        raise ValueError('blq must be {rule: drop, rationale}; censored-likelihood (M3) BLQ handling is not available in this release')
    cm = c['comparison_models']
    if cm is None: c.pop('comparison_models')
    else:
        if not isinstance(cm, list) or not cm: raise ValueError('comparison_models must be a non-empty list or omitted')
        for x in cm:
            if not isinstance(x, dict) or set(x) != {'structural', 'additional_initial_estimates'} or x['structural'] not in MODELS or x['structural'] == s:
                raise ValueError('Each comparison model needs a different structural model and additional_initial_estimates')
            extra = set(MODELS[x['structural']]) - set(MODELS[s])
            a = x['additional_initial_estimates']
            if not isinstance(a, dict) or set(a) != extra or not all(_pos(v) for v in a.values()):
                raise ValueError(f"Comparison {x['structural']} needs additional_initial_estimates for: " + (', '.join(sorted(extra)) or 'none ({})'))
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def _bool(d, col):
    v = d[col].astype(str).str.strip().str.lower()
    if not v.isin(['true', 'false']).all(): raise ValueError(f'{col} must be true or false')
    return v.eq('true')


def load_data(path, c):
    d = pd.read_csv(path, dtype={'subject': str, 'group': str})
    need = {'subject', 'time', 'concentration', 'dose'} | ({'infusion_duration'} if c['design']['route'] == 'iv_infusion' else set())
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'group', 'blq', 'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    d['time'] = pd.to_numeric(d.time, errors='coerce')
    if d.time.isna().any() or (d.time < 0).any(): raise ValueError('time must be numeric and non-negative')
    d['dose'] = pd.to_numeric(d.dose, errors='coerce'); d['concentration'] = pd.to_numeric(d.concentration, errors='coerce')
    d['blq'] = _bool(d, 'blq') if 'blq' in d else False
    d['exclude'] = _bool(d, 'exclude') if 'exclude' in d else False
    if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
        raise ValueError('Every excluded row needs an exclusion_reason')
    is_dose = d.dose.notna()
    if (is_dose & (d.concentration.notna() | d.blq)).any(): raise ValueError('A row is either a dose (dose filled) or an observation (concentration or blq), not both')
    if (~is_dose & d.concentration.isna() & ~d.blq).any(): raise ValueError('Observation rows need a concentration or blq = true')
    if (d.blq & d.concentration.notna()).any(): raise ValueError('A BLQ row must have an empty concentration')
    if (d.dose[is_dose] <= 0).any() or (d.concentration.dropna() < 0).any(): raise ValueError('Doses must be positive and concentrations non-negative')
    if 'infusion_duration' in d:
        dur = pd.to_numeric(d.infusion_duration, errors='coerce')
        if (is_dose & ~(dur > 0)).any(): raise ValueError('Every infusion dose row needs a positive infusion_duration')
        d['infusion_duration'] = dur
    if d[~is_dose].duplicated(['subject', 'time']).any(): raise ValueError('Duplicate observation times within a subject')
    first = d[is_dose].groupby('subject').time.min()
    if set(d.subject) - set(first.index): raise ValueError('Every subject needs at least one dose row')
    obs = d[~is_dose & ~d.exclude & ~d.blq]
    if (obs.groupby('subject').time.min() <= first.reindex(obs.subject.unique())).any(): raise ValueError('Observations at or before the first dose belong to no fitted profile; mark predose samples excluded with a reason')
    if d.blq.any() and c['blq'] is None: raise ValueError('BLQ observations present: declare blq {rule: drop, rationale}')
    return d


def _events(d, c):
    route = c['design']['route']; use = d[~d.exclude].copy(); ids = {s: i + 1 for i, s in enumerate(sorted(use.subject.unique()))}
    rows = []; dropped = int(use.blq.sum())
    for _, x in use.sort_values(['subject', 'time']).iterrows():
        if pd.notna(x.dose):
            rate = float(x.dose/x.infusion_duration) if route == 'iv_infusion' else 0.
            rows.append((ids[x.subject], float(x.time), 0., float(x.dose), 1, rate, 'depot' if route == 'sc' else 'central'))
        elif not x.blq:
            rows.append((ids[x.subject], float(x.time), float(x.concentration), 0., 0, 0., 'central'))
    ev = pd.DataFrame(rows, columns=['id', 'time', 'dv', 'amt', 'evid', 'rate', 'cmt']).sort_values(['id', 'time', 'evid'], ascending=[True, True, False])
    return ev, ids, dropped


def _payload(ev, c, structural, init, eta):
    m = c['model']; fixed = m['fixed'] or {}
    return {'model': structural, 'route': 'sc' if c['design']['route'] == 'sc' else 'iv', 'species': m['observed_species'],
            'init': {**init, **fixed}, 'fixed': sorted(fixed), 'eta': eta, 'eta_init': m['random_effect_initial_variance'],
            'residual': m['residual'], 'residual_init': m['residual_initial'],
            'bioavailability': (m['bioavailability'] or {}).get('value', 1.), 'seed': c['estimation']['seed'],
            'max_outer_iterations': c['estimation']['max_outer_iterations'], 'data': ev.to_dict('list')}


def _units(dz):
    dose, conc, t = dz['dose_unit'], dz['concentration_unit'], dz['time_unit']
    vol = f'({dose})/({conc})'
    return {'cl': f'{vol}/{t}', 'q': f'{vol}/{t}', 'v1': vol, 'v2': vol, 'vmax': f'{dose}/{t}', 'km': conc, 'kss': conc, 'r0': conc,
            'kint': f'1/{t}', 'kdeg': f'1/{t}', 'ka': f'1/{t}'}


def _gate(fit, c, structural, obs_range):
    th = {r: i for i, r in enumerate(fit['theta']['name'])}; T = fit['theta']
    cov = fit['cov_method'] or ''; kind = cov.split(' ')[0]  # nlmixr2 reports e.g. "r,s (full)"
    finite = all(se is not None and math.isfinite(se) for se, fx in zip(T['se'], T['fixed']) if not fx)
    cov_status = ('sandwich' if kind == 'r,s' else 'fallback' if kind in ('r', 's') else 'failed') if finite else 'failed'
    corr_pairs = []
    if fit.get('theta_correlation'):
        names = fit['theta_correlation']['names']; M = np.array(fit['theta_correlation']['matrix'], float)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                if abs(M[i, j]) > CORR_LIMIT: corr_pairs.append([names[i], names[j], float(M[i, j])])
    params = MODELS[structural] + (('ka',) if c['design']['route'] == 'sc' else ())
    fixed = c['model']['fixed'] or {}; out = {}; units = _units(c['design'])
    nonlin_bad = None
    if structural in NONLINEAR:
        key = 'km' if structural == 'linear_mm_2cmt' else 'kss'
        if key not in fixed:
            v = float(np.exp(T['estimate'][th['l' + key]]))
            if not obs_range[0] <= v <= obs_range[1]: nonlin_bad = f'{key} {v:.3g} outside the observed concentration range {obs_range[0]:.3g}-{obs_range[1]:.3g}'
    for p in params:
        i = th.get('l' + p)
        if p in fixed:
            out[p] = {'value': fixed[p], 'status': 'fixed', 'unit': units[p], 'reasons': []}; continue
        est, se = T['estimate'][i], T['se'][i]
        # nlmixr2's %RSE divides by the log-scale estimate (meaningless near log = 0); use the
        # natural-scale approximation 100 sqrt(exp(se^2) - 1) for the log-normal parameter instead.
        # A huge log-scale SE (unidentified parameter) would overflow exp; it is withheld below.
        huge = se is not None and math.isfinite(se) and se*se > 700
        rse = float(100*math.sqrt(math.expm1(se*se))) if se is not None and math.isfinite(se) and not huge else None
        rec = {'estimate': float(np.exp(est)), 'log_estimate': est, 'log_se': se, 'rse_percent': rse,
               'ci95': [T['ci_low'][i], T['ci_high'][i]], 'unit': units[p], 'reasons': [], 'flags': []}
        if not (est is not None and math.isfinite(est)): rec['reasons'].append('non-finite estimate')
        if cov_status == 'failed' or se is None or not math.isfinite(se): rec['reasons'].append('no covariance (standard error unavailable)')
        if any('l' + p in pr[:2] for pr in corr_pairs): rec['reasons'].append('|correlation| > 0.95 with another fixed effect')
        if nonlin_bad and p in NONLINEAR[structural]: rec['reasons'].append('nonlinearity not supported: ' + nonlin_bad)
        if huge: rec['reasons'].append(f'not identified: log-scale SE {se:.3g}')
        elif rec['rse_percent'] is not None and rec['rse_percent'] > RSE_WITHHOLD: rec['reasons'].append(f"not identified: natural-scale RSE {rec['rse_percent']:.3g}% > {RSE_WITHHOLD:g}%")
        elif rec['rse_percent'] is not None and rec['rse_percent'] > RSE_FLAG: rec['flags'].append('rse_above_50_percent')
        if cov_status == 'fallback': rec['flags'].append('covariance_fallback_' + cov.replace(',', ''))
        out[p] = rec
    if structural in NONLINEAR:
        bad = [p for p in NONLINEAR[structural] if p in out and out[p].get('reasons')]
        for p in NONLINEAR[structural]:
            if bad and p in out and 'estimate' in out[p] and p not in bad:
                out[p]['reasons'].append('withheld with the other nonlinear parameters (' + ', '.join(bad) + ') because they are interpreted jointly')
    for p, rec in out.items():
        if 'estimate' not in rec: continue
        rec['status'] = 'withheld' if rec['reasons'] else 'reportable'
        if rec['status'] == 'withheld': rec['ci95'] = None
    iiv = {}
    for k, w in fit['omega_diagonal'].items():
        p = k.replace('eta_', ''); sh = fit['shrinkage_sd_percent'].get(k)
        iiv[p] = {'omega_variance': w, 'cv_percent': float(100*math.sqrt(math.expm1(w))) if w is not None and 0 <= w < 700 else None, 'shrinkage_percent': sh,
                  'individual_estimates': 'withheld' if sh is None or sh > SHRINK_LIMIT else 'usable_with_caution'}
    resid = {r: float(T['estimate'][th[r]]) for r in ('prop_sd', 'add_sd') if r in th}
    return {'parameters': out, 'iiv': iiv, 'residual': resid, 'covariance': {'method': cov, 'status': cov_status},
            'high_correlation_pairs': corr_pairs, 'nonlinearity_check': nonlin_bad}


def compute(d, c):
    from .r_bridge import call_pop
    ev, ids, dropped = _events(d, c); s = c['model']['structural']
    pos = ev[(ev.evid == 0) & (ev.dv > 0)].dv
    obs_range = [float(pos.min()), float(pos.max())]
    fit = call_pop('nlmixr2_fit', _payload(ev, c, s, c['model']['initial_estimates'], c['model']['random_effects']))
    g = _gate(fit, c, s, obs_range)
    comps = []
    for x in c.get('comparison_models', []) or []:
        keep = set(MODELS[x['structural']]) | ({'ka'} if c['design']['route'] == 'sc' else set())
        init = {k: v for k, v in c['model']['initial_estimates'].items() if k in keep} | x['additional_initial_estimates']
        eta = [k for k in c['model']['random_effects'] if k in keep]
        cc = deepcopy(c); cc['model']['fixed'] = {k: v for k, v in (c['model']['fixed'] or {}).items() if k in keep}
        cc['model']['observed_species'] = c['model']['observed_species'] if x['structural'] == 'qss_tmdd_2cmt' else None
        try:
            f2 = call_pop('nlmixr2_fit', _payload(ev, cc, x['structural'], init, eta))
            comps.append({'structural': x['structural'], 'objective': f2['objective'], 'aic': f2['aic'], 'bic': f2['bic'],
                          'delta_aic_vs_primary': f2['aic'] - fit['aic'], 'covariance_method': f2['cov_method'], 'status': 'fitted'})
        except (ArithmeticError, RuntimeError) as exc:
            comps.append({'structural': x['structural'], 'status': 'failed', 'reason': str(exc)[:300]})
    failing = [{'parameter': p, 'status': 'withheld', 'reason': '; '.join(r['reasons'])} for p, r in g['parameters'].items() if r['status'] == 'withheld']
    failing += [{'parameter': 'eta_' + p, 'status': 'individual_estimates_withheld', 'reason': f"shrinkage {v['shrinkage_percent']}% > {SHRINK_LIMIT:g}%"}
                for p, v in g['iiv'].items() if v['individual_estimates'] == 'withheld']
    tab = pd.DataFrame(fit['table']); inv = {v: k for k, v in ids.items()}
    if 'ID' in tab: tab['subject'] = tab.ID.astype(int).map(inv)
    primary = {'structural_model': s, 'route': c['design']['route'], 'observed_species': c['model']['observed_species'],
               **g, 'fit': {k: fit[k] for k in ('objective', 'aic', 'bic', 'n_obs', 'n_subjects', 'message')},
               'observed_concentration_range': obs_range, 'blq_observations_dropped': dropped,
               'excluded_rows': d[d.exclude][['subject', 'time', 'exclusion_reason']].to_dict('records') if d.exclude.any() else [],
               'comparisons': comps, 'nlmixr2': {'version': fit['r_environment']['packages'].get('nlmixr2est'),
                                                'rxode2': fit['r_environment']['packages'].get('rxode2'), 'warnings': fit['warnings']},
               'model_code': fit['model_code']}
    return {'analysis_type': 'population_pk', 'primary': primary, 'fit_table': tab.to_dict('list'), 'r_environment': fit['r_environment'],
            'must_mention': _notes(c, primary), 'failing_items': failing,
            'limitations': ['Population model of declared structure only: no covariate models, interoccasion variability, PD, full (kon/koff) TMDD or censored-likelihood (M3) BLQ.',
                            'FOCEi Wald intervals on the log scale; with few animals they can be too narrow (see the registered calibration). No VPC or bootstrap is run.',
                            'QSS target parameters (KSS, R0) are in the drug concentration unit as drug-binding capacity, assuming 1:1 binding; they are not measured target levels.',
                            'Estimates by R nlmixr2; no NONMEM, Monolix, Phoenix or Prism equivalence beyond the recorded checks.']}


def _notes(c, p):
    m = c['model']; dz = c['design']; notes = []
    notes.append(f"Population PK by R nlmixr2 {p['nlmixr2']['version']} (FOCEi) for {dz['analyte']} in {dz['matrix']} ({dz['species']}), route {dz['route']}: "
                 f"{p['structural_model']} declared before fitting ({m['rationale']}); {p['fit']['n_subjects']} subjects, {p['fit']['n_obs']} observations.")
    if p['structural_model'] == 'qss_tmdd_2cmt':
        notes.append(f"QSS TMDD fitted to {m['observed_species']} drug concentrations; KSS and R0 are in drug concentration units (drug-binding capacity, 1:1).")
    if m['fixed']: notes.append('Fixed (not estimated): ' + ', '.join(f'{k} = {v:g}' for k, v in m['fixed'].items()) + f" ({m['fixed_source']}).")
    if dz['route'] == 'sc':
        b = m['bioavailability']; notes.append(f"SC bioavailability fixed at {b['value']:g} ({b['source']})" + ('; CL and V are therefore apparent (CL/F, V/F).' if b['value'] == 1 else '.'))
    w = [k for k, r in p['parameters'].items() if r['status'] == 'withheld']
    notes.append('Covariance: ' + {'sandwich': 'r,s (sandwich) available.', 'fallback': f"fallback {p['covariance']['method']} only; intervals are less robust.",
                                   'failed': 'failed; every interval is withheld and estimates are audit values only.'}[p['covariance']['status']])
    if w: notes.append('Withheld parameters: ' + '; '.join(f"{k} ({'; '.join(p['parameters'][k]['reasons'])})" for k in w) + '.')
    if p['nonlinearity_check']: notes.append('Saturation is not supported by these data (' + p['nonlinearity_check'] + '); a nonlinear model fitted to them does not demonstrate TMDD.')
    hs = [k for k, v in p['iiv'].items() if v['individual_estimates'] == 'withheld']
    if hs: notes.append('Individual (EBE) parameter claims withheld for ' + ', '.join(hs) + f' (shrinkage > {SHRINK_LIMIT:g}%).')
    if p['blq_observations_dropped']: notes.append(f"{p['blq_observations_dropped']} BLQ observation(s) dropped from the likelihood (M1, declared: {c['blq']['rationale']}).")
    if p['comparisons']:
        notes.append('Comparison models are descriptive and were declared before fitting; the reported model is the declared one: ' +
                     '; '.join(f"{x['structural']} delta AIC {x['delta_aic_vs_primary']:+.1f}" if x['status'] == 'fitted' else f"{x['structural']} failed" for x in p['comparisons']) + '.')
    if p['nlmixr2']['warnings']: notes.append('nlmixr2 warnings retained: ' + '; '.join(p['nlmixr2']['warnings'][:3]) + '.')
    return notes


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    t = pd.DataFrame(result['fit_table']); t = t[t.DV > 0] if 'DV' in t else t
    if t.empty or 'IPRED' not in t: return []
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.3), layout='constrained')
        lo, hi = t.DV.min(), t.DV.max()
        for a, col in zip(ax[:2], ('PRED', 'IPRED')):
            a.scatter(t[col], t.DV, s=8, color=pstyle.color(token, 0)); a.plot([lo, hi], [lo, hi], color=pstyle.MUTED, lw=.8)
            a.set_xscale('log'); a.set_yscale('log'); a.set(xlabel=col, ylabel='Observed', title=f'Observed vs {col}')
        if 'IWRES' in t: ax[2].scatter(t.TIME, t.IWRES, s=8, color=pstyle.color(token, 1)); ax[2].axhline(0, color=pstyle.MUTED, lw=.8)
        ax[2].set(xlabel='Time', ylabel='IWRES', title='Individual weighted residuals')
        pstyle.save(fig, out/'goodness_of_fit', token, 300)
        subs = list(dict.fromkeys(t.subject))[:12] if 'subject' in t else []
        if subs:
            n = len(subs); cols = 4; rows = math.ceil(n/cols)
            fig, axes = plt.subplots(rows, cols, figsize=(11, 2.4*rows), layout='constrained', squeeze=False)
            for a, s in zip(axes.flat, subs):
                q = t[t.subject == s].sort_values('TIME')
                a.scatter(q.TIME, q.DV, s=10, color=pstyle.color(token, 0)); a.plot(q.TIME, q.IPRED, color=pstyle.color(token, 1), lw=1)
                a.plot(q.TIME, q.PRED, color=pstyle.MUTED, lw=.8, ls='--'); a.set_yscale('log'); a.set_title(str(s), fontsize=8)
            for a in list(axes.flat)[n:]: a.axis('off')
            pstyle.save(fig, out/'individual_fits', token, 300)
            return ['goodness_of_fit', 'individual_fits']
    return ['goodness_of_fit']
