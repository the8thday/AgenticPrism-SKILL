"""Growth-rate inhibition (GR) metrics for cell proliferation assays (0.14.0).

GR(c) = 2^(log2(x(c)/x0) / log2(x_ctrl/x0)) - 1 (Hafner et al., Nat Methods 2016),
with x0 the time-zero count and x_ctrl the untreated count at the endpoint. Each
curve is fitted as in the authors' reference implementations (datarail gr_metrics,
R GRmetrics / Python gr50): GR(c) = GRinf + (1 - GRinf)/(1 + (c/GEC50)^h) by least
squares with bounds GRinf in [-1, 1], log10 GEC50 in [log10(min c/100), log10(max c*100)],
h in [0.1, 5]; an F test against a flat fit (the mean; alpha 0.05, denominator df n - 3 + 1 as
in the reference code) selects the flat curve when not significant. GR50 inverts the
fit at 0.5; GRmax is the minimum observed GR at the two highest concentrations;
GR_AOC is the trapezoid area of 1 - GR over log10 concentration divided by its range.
GR50 outside the tested range is withheld (no extrapolation). Summaries across
declared independent experiments use t intervals on log10 GR50 and on GRmax/GR_AOC.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'growth_rate_inhibition'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'growth_rate_inhibition', 'input': None, 'source': None,
            'design': {'cell_line_rationale': None, 'readout': None, 'concentration_unit': None, 'time_unit': None,
                       'independent_unit': None, 'independence_source': None},
            'fit': {'alpha': .05}, 'summary': {'confidence_level': .95},
            'report': {'plot_style': 'prism_like'}}
KEYS = ['cell_line', 'agent', 'time', 'experiment']


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported GR configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'growth_rate_inhibition': raise ValueError('Unsupported GR schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    if c['fit']['alpha'] != .05: raise ValueError('The flat-fit test uses the published alpha 0.05')
    lv = c['summary']['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'cell_line': str, 'agent': str, 'experiment': str})
    need = set(KEYS) | {'concentration', 'cell_count', 'cell_count__ctrl', 'cell_count__time0'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded well needs an exclusion_reason')
    else: d['exclude'] = False
    for col in ('time', 'concentration', 'cell_count', 'cell_count__ctrl', 'cell_count__time0'):
        v = pd.to_numeric(d[col], errors='coerce')
        if v.isna().any() or not np.isfinite(v).all(): raise ValueError(f'{col} must be numeric')
        d[col] = v
    if (d.concentration <= 0).any(): raise ValueError('Concentrations must be positive; untreated controls go in cell_count__ctrl, not as concentration 0')
    if (d[['cell_count', 'cell_count__ctrl', 'cell_count__time0']] < 0).any().any(): raise ValueError('Counts cannot be negative')
    if (d.cell_count__time0 <= 0).any(): raise ValueError('Time-zero counts must be positive')
    if d.duplicated(KEYS + ['concentration']).any(): raise ValueError('One value per curve and concentration; average technical wells first by a declared rule')
    return d


def gr_values(x, xc, x0):
    a = np.log2(np.maximum(x/x0, 1e-6)); b = np.log2(np.maximum(xc/x0, 1e-6))
    return 2**(a/b) - 1


def logistic(c, p):
    einf, lmid, h = p
    return (1. - einf)/(1 + (c/10**lmid)**h) + einf


def _fit(fn, c, y, prior, bounds):
    rss = lambda p: float(((y - fn(c, p))**2).sum())
    return optimize.minimize(rss, prior, bounds=bounds, method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-12, 'maxiter': 20000})


def fit_curve(conc, gr, alpha=.05):
    order = np.argsort(conc, kind='mergesort'); c, y = np.asarray(conc, float)[order], np.asarray(gr, float)[order]; n = len(c)
    bounds = [(-1, 1), tuple(np.log10([c.min()/100, c.max()*100])), (.1, 5)]
    prior = np.array([.1, np.log10(np.median(c)), 2.])
    best = None
    for start in (prior, np.array([y.min(), np.log10(np.median(c)), 1.]), np.array([-.5, np.log10(c.min()*10), 1.])):
        r = _fit(logistic, c, y, np.clip(start, [b[0] for b in bounds], [b[1] for b in bounds]), bounds)
        if best is None or r.fun < best.fun: best = r
    flat_level = float(y.mean())  # least-squares flat fit, unbounded as in R GRmetrics
    rss2, rss1 = float(best.fun), float(((y - flat_level)**2).sum()); df1, df2 = 2, n - 3 + 1
    f = ((rss1 - rss2)/df1)/(rss2/df2) if rss2 > 0 else np.inf
    p = float(stats.f.sf(f, df1, df2)) if np.isfinite(f) else 0.
    tss = float(((y - y.mean())**2).sum())
    out = {'n_concentrations': int(n), 'pval': p, 'tested_range': [float(c.min()), float(c.max())],
           'GRmax': float(min(y[-2:])), 'GR_AOC': float(np.trapezoid(1 - y, np.log10(c))/(np.log10(c.max()) - np.log10(c.min())))}
    # Optimizer status flags are not used: they reported converged bounded fits as failures.
    if p > alpha or not np.isfinite(rss2):
        out.update(fit='flat', GRinf=flat_level, GEC50=None, h_GR=None, r2=1 - rss1/tss if tss > 0 else None,
                   GR50_audit=float('inf') if flat_level > .5 else float('-inf'))
    else:
        einf, lmid, h = best.x
        gr50 = 10**lmid*((.5 - 1)/(einf - .5))**(1/h) if min(1, einf) <= .5 <= max(1, einf) else float('inf')
        out.update(fit='sigmoid', GRinf=float(einf), GEC50=float(10**lmid), h_GR=float(h), r2=1 - rss2/tss if tss > 0 else None,
                   GR50_audit=float(gr50))
    g = out['GR50_audit']; lo, hi = out['tested_range']
    if out['fit'] == 'flat':
        out.update(GR50=None, GR50_status='flat_fit_no_concentration_response')
    elif not np.isfinite(g):
        out.update(GR50=None, GR50_status='not_reached_GRinf_above_0.5')
    elif g > hi or g < lo:
        out.update(GR50=None, GR50_status='outside_tested_range')
    else:
        out.update(GR50=float(g), GR50_status='estimated')
    if out['GEC50'] is not None: out['GEC50_within_tested_range'] = bool(lo <= out['GEC50'] <= hi)
    return out


def compute(d, c):
    u = d[~d.exclude].copy(); lv = c['summary']['confidence_level']
    u['GR'] = gr_values(u.cell_count.to_numpy(float), u.cell_count__ctrl.to_numpy(float), u.cell_count__time0.to_numpy(float))
    fits, gvals, failing = [], [], []
    for key, part in u.groupby(KEYS, sort=False):
        rec = dict(zip(KEYS, [str(k) if i != 2 else float(k) for i, k in enumerate(key)])); diag = []
        growth = part.cell_count__ctrl.iloc[0]/part.cell_count__time0.iloc[0]
        if part.cell_count__ctrl.nunique() > 1 or part.cell_count__time0.nunique() > 1: diag.append('multiple_control_or_time0_values_in_curve')
        # Post-freeze gate (0.14.0 live-agent finding): below half a doubling the GR denominator amplifies count noise (GR up to +29 observed).
        if growth < 2**.5:
            rec.update(status='withheld', reportable=False, reason='untreated_control_grew_less_than_half_a_doubling', control_fold_growth=float(growth)); fits.append(rec)
            failing.append({**{k: rec[k] for k in KEYS}, 'reason': rec['reason']}); continue
        if growth < 2: diag.append('control_less_than_one_doubling')
        if part.concentration.nunique() < 5:
            rec.update(status='withheld', reportable=False, reason='fewer_than_5_concentrations'); fits.append(rec)
            failing.append({**{k: rec[k] for k in KEYS}, 'reason': rec['reason']}); continue
        f = fit_curve(part.concentration.to_numpy(float), part.GR.to_numpy(float), c['fit']['alpha'])
        if not np.isfinite(f['GR50_audit']): f['GR50_audit'] = 'inf' if f['GR50_audit'] > 0 else '-inf'  # JSON has no infinity (post-freeze format fix)
        rec.update(f, status='estimated', reportable=True, control_fold_growth=float(growth), diagnostics=diag)
        fits.append(rec)
        gvals += [{**{k: rec[k] for k in KEYS}, 'concentration': float(cc), 'GR': float(g)} for cc, g in zip(part.concentration, part.GR)]
    summaries = []
    for key, grp in pd.DataFrame([f for f in fits]).groupby(['cell_line', 'agent', 'time'], sort=False):
        rows = grp.to_dict('records'); n = len(rows); s = {'cell_line': key[0], 'agent': key[1], 'time': float(key[2]), 'n_experiments': n}
        ok = [r for r in rows if r['reportable']]
        def tci(v):
            v = np.asarray(v, float); m = float(v.mean())
            if len(v) < 2: return {'mean': m, 'ci': None}
            h = float(stats.t.ppf((1+lv)/2, len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))); return {'mean': m, 'ci': [m-h, m+h]}
        if n < 2 or len(ok) < n:
            s.update(status='withheld', reason='fewer_than_two_experiments' if n < 2 else 'a_curve_was_withheld')
        else:
            s.update(status='estimated', GRmax=tci([r['GRmax'] for r in ok]), GR_AOC=tci([r['GR_AOC'] for r in ok]))
            if all(r['GR50_status'] == 'estimated' for r in ok):
                t = tci(np.log10([r['GR50'] for r in ok])); s['GR50_geometric'] = {'mean': 10**t['mean'], 'ci': [10**x for x in t['ci']]}
            else:
                s['GR50_geometric'] = None; s['GR50_summary_status'] = 'withheld_some_GR50_not_estimated:' + ','.join(sorted({r['GR50_status'] for r in ok if r['GR50_status'] != 'estimated'}))
        summaries.append(s)
    notes = ['GR values follow Hafner et al. (2016): growth relative to untreated controls from time-zero, treated and control counts, which corrects relative-viability IC50s for differences in division rate and assay length.',
             'Curve fitting follows the authors\' reference implementations (bounded least squares, flat-fit F test at alpha 0.05 with their denominator df n - 3 + 1). GRmax and GR_AOC are computed from observed values; GR50 outside the tested concentrations is withheld, never extrapolated.',
             'Independent unit for summaries: ' + c['design']['independent_unit'] + ' (' + c['design']['independence_source'] + '). Replicate wells within one plate are not experiments.',
             'GR values near -1 indicate cell death; GR between 0 and 1 is partial growth inhibition (cytostatic). GRinf and GRmax describe efficacy; GEC50/GR50 describe potency.']
    if any('control_less_than_one_doubling' in f.get('diagnostics', []) for f in fits): notes.append('Some untreated controls grew less than one doubling; GR values are noisy when control growth is small.')
    if any(f.get('fit') == 'flat' for f in fits): notes.append('Curves with a non-significant concentration response were fitted flat; their GEC50 and GR50 are not defined.')
    if failing: notes.append('Withheld curves: ' + '; '.join(f"{x['cell_line']}/{x['agent']}/{x['experiment']} ({x['reason']})" for x in failing) + '.')
    primary = {'status': 'estimated' if any(f['reportable'] for f in fits) else 'withheld', 'reportable': any(f['reportable'] for f in fits),
               'confidence_level': lv, 'fits': fits, 'summaries': summaries, 'concentration_unit': c['design']['concentration_unit']}
    return {'analysis_type': 'growth_rate_inhibition', 'primary': primary, 'gr_values': gvals, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['Static GR metrics from end-point counts only; no time-course GR, toxic/static decomposition or dead-cell readouts.',
                            'No interval for single-curve GEC50/GR50; uncertainty comes only from independent experiments.',
                            'No Prism equivalence; agreement with the authors\' published toy-example outputs is the benchmark.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    g = pd.DataFrame(result['gr_values']); names = []
    if g.empty: return names
    groups = list(g.groupby(['cell_line', 'agent', 'time'], sort=False))[:12]
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        cols = min(3, len(groups)); rows = int(np.ceil(len(groups)/cols))
        fig, axes = plt.subplots(rows, cols, figsize=(3.2*cols, 2.7*rows), layout='constrained', squeeze=False)
        for ax, ((cl, ag, tm), part) in zip(axes.ravel(), groups):
            for i, (ex, e) in enumerate(part.groupby('experiment', sort=False)):
                ax.scatter(e.concentration, e.GR, s=12, color=pstyle.color(token, i % 8))
                f = next((x for x in p['fits'] if x['cell_line'] == cl and x['agent'] == ag and x['time'] == tm and x['experiment'] == ex and x['reportable']), None)
                if f:
                    cc = np.geomspace(e.concentration.min(), e.concentration.max(), 100)
                    yy = np.full_like(cc, f['GRinf']) if f['fit'] == 'flat' else logistic(cc, (f['GRinf'], np.log10(f['GEC50']), f['h_GR']))
                    ax.plot(cc, yy, color=pstyle.color(token, i % 8), lw=1)
            for yref in (0, .5, -1): ax.axhline(yref, ls=':', color=pstyle.MUTED, lw=.8)
            ax.set_xscale('log'); ax.set(ylim=(-1.05, 1.25), title=f'{cl} / {ag} / {tm:g}', xlabel=f"Concentration ({p['concentration_unit']})", ylabel='GR value')
        for ax in axes.ravel()[len(groups):]: ax.axis('off')
        pstyle.save(fig, out/'gr_curves', token, 300); names.append('gr_curves')
    return names
