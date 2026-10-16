"""Prespecified outlier flagging for one sample per declared group (0.13.6).

Flag only: nothing is removed from any other analysis. Methods:
- Grubbs (1969) single-outlier test, two-sided or one-sided (max/min), normal.
- Generalized ESD (Rosner 1983) for up to a declared number r of outliers.
- ROUT (Motulsky & Brown 2006) for a location model: Lorentzian robust fit with
  RSDR (68.27th percentile of |residuals| x sqrt(N/(N-1))) recomputed each
  iteration, FDR-style threshold alpha_i = Q (N - i + 1) / N over the 30% largest
  residuals, t reference with N - 1 df. This follows the published description;
  the robust optimizer is ours, so Prism equivalence is not claimed.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'outlier_screen'}
METHODS = ('grubbs', 'generalized_esd', 'rout')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'outlier_screen', 'input': None, 'source': None, 'method': None,
            'design': {'unit': None, 'rationale': None, 'outcome': None, 'outcome_unit': None,
                       'purpose': None, 'prespecified_source': None},
            'test': {'alpha': None, 'q': None, 'max_outliers': None, 'sides': None, 'scale': 'original', 'scale_rationale': None},
            'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported outlier configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'outlier_screen' or c['method'] not in METHODS:
        raise ValueError('Declare method grubbs, generalized_esd or rout')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    t = c['test']; m = c['method']
    if m in ('grubbs', 'generalized_esd'):
        if type(t['alpha']) not in (int, float) or not 0 < t['alpha'] < .5: raise ValueError('Declare alpha (0 to 0.5) before screening')
        if t['q'] is not None: raise ValueError('q applies only to rout')
    else:
        if type(t['q']) not in (int, float) or not 0 < t['q'] < .5: raise ValueError('Declare the ROUT false discovery level q (for example 0.01)')
        if t['alpha'] is not None: raise ValueError('alpha does not apply to rout; declare q')
    if m == 'grubbs':
        if t['sides'] not in ('two_sided', 'max', 'min'): raise ValueError('Declare Grubbs sides: two_sided, max or min')
    elif t['sides'] is not None: raise ValueError('sides applies only to grubbs')
    if m == 'generalized_esd':
        if type(t['max_outliers']) is not int or t['max_outliers'] < 1: raise ValueError('Declare max_outliers (r) before screening')
    elif t['max_outliers'] is not None: raise ValueError('max_outliers applies only to generalized_esd')
    if t['scale'] not in ('original', 'log'): raise ValueError('scale must be original or log')
    if t['scale'] == 'log' and not _text(t['scale_rationale']): raise ValueError('Declare why the log scale is appropriate')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'unit_id': str, 'group': str})
    need = {'unit_id', 'group', 'value'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every pre-existing exclusion needs a documented reason')
    else: d['exclude'] = False
    v = pd.to_numeric(d.value, errors='coerce')
    if d.unit_id.isna().any() or d.group.isna().any() or v.isna().any() or not np.isfinite(v).all():
        raise ValueError('Missing or non-numeric values; no silent deletion, imputation or BLQ substitution')
    d['value'] = v
    if d.duplicated(['group', 'unit_id']).any(): raise ValueError('Each unit appears once per group')
    if c['test']['scale'] == 'log' and (d.loc[~d.exclude, 'value'] <= 0).any(): raise ValueError('Log scale requires positive values; no pseudocounts')
    return d


def grubbs(x, alpha, sides):
    x = np.asarray(x, float); n = len(x); mean, sd = x.mean(), x.std(ddof=1)
    if sd <= 0: return {'status': 'withheld', 'reason': 'zero_variance'}
    dev = x - mean
    j = int(np.argmax(np.abs(dev) if sides == 'two_sided' else dev if sides == 'max' else -dev))
    G = float(abs(dev[j]) if sides == 'two_sided' else dev[j] if sides == 'max' else -dev[j])/sd
    a = alpha/(2*n) if sides == 'two_sided' else alpha/n
    t = stats.t.ppf(1-a, n-2); crit = float((n-1)/np.sqrt(n)*np.sqrt(t*t/(n-2+t*t)))
    # Bonferroni-type p value (as in the common implementations), capped at 1.
    inner = n*(n-2)*G*G/((n-1)**2 - n*G*G) if (n-1)**2 - n*G*G > 0 else np.inf
    tt = np.sqrt(inner) if np.isfinite(inner) else np.inf
    p = float(min(1., (2*n if sides == 'two_sided' else n)*stats.t.sf(tt, n-2))) if np.isfinite(tt) else 0.
    return {'status': 'tested', 'G': float(G), 'critical_value': crit, 'p_value': p, 'candidate_index': j, 'flagged': [j] if G > crit else []}


def generalized_esd(x, alpha, r):
    x = np.asarray(x, float); n = len(x); idx = list(range(n)); rows = []; last = 0
    if n - r < 3: return {'status': 'withheld', 'reason': 'max_outliers_leaves_fewer_than_three_values'}
    for i in range(1, r+1):
        z = x[idx]; sd = z.std(ddof=1)
        if sd <= 0: rows.append({'step': i, 'status': 'stopped_zero_variance'}); break
        d = np.abs(z - z.mean()); k = int(np.argmax(d)); R = float(d[k]/sd)
        p = 1 - alpha/(2*(n-i+1)); t = stats.t.ppf(p, n-i-1); lam = float((n-i)*t/np.sqrt((n-i-1+t*t)*(n-i+1)))
        rows.append({'step': i, 'candidate_index': idx[k], 'value': float(x[idx[k]]), 'mean': float(z.mean()), 'sd': float(sd),
                     'R': R, 'critical_value': lam, 'exceeds': R > lam})
        if R > lam: last = i
        idx.pop(k)
    return {'status': 'tested', 'steps': rows, 'n_outliers': last, 'flagged': [r_['candidate_index'] for r_ in rows[:last]]}


def _rsdr(res, k=1):
    a = np.sort(np.abs(res)); n = len(a); pos = .6827*n
    lo = int(np.floor(pos)); hi = min(lo+1, n)
    p68 = a[max(lo-1, 0)] + (pos-lo)*(a[hi-1]-a[max(lo-1, 0)]) if lo >= 1 else a[0]*pos
    return float(p68*np.sqrt(n/(n-k)))


def rout(x, q, max_iter=200):
    x = np.asarray(x, float); n = len(x); mu = float(np.median(x))
    for _ in range(max_iter):
        s = _rsdr(x - mu)
        if s <= 0: return {'status': 'withheld', 'reason': 'zero_robust_scale'}
        span = max(np.ptp(x), s)
        new = optimize.minimize_scalar(lambda m: np.log1p(((x-m)/s)**2).sum(), bounds=(mu-span, mu+span), method='bounded',
                                       options={'xatol': 1e-12*max(1., abs(mu))}).x
        if abs(new - mu) <= 1e-10*max(1., abs(mu), s): mu = float(new); break
        mu = float(new)
    else: return {'status': 'withheld', 'reason': 'robust_fit_not_converged'}
    res = x - mu; s = _rsdr(res); order = np.argsort(np.abs(res)); flagged = []; tested = []
    for i in range(int(.7*n), n+1):
        if i == 0: continue
        j = int(order[i-1]); t = abs(res[j])/s; p = float(2*stats.t.sf(t, n-1)); a_i = q*(n-(i-1))/n
        tested.append({'rank': i, 'index': j, 'value': float(x[j]), 't_ratio': float(t), 'p_value': p, 'threshold': float(a_i)})
        if p < a_i:
            flagged = [int(v) for v in order[i-1:]]; break
    return {'status': 'tested', 'robust_location': mu, 'rsdr': s, 'tested': tested, 'flagged': flagged}


def compute(d, c):
    t = c['test']; m = c['method']; u = d[~d.exclude]; groups = []
    for g, part in u.groupby('group', sort=False):
        part = part.reset_index(drop=True); raw = part.value.to_numpy(float)
        x = np.log(raw) if t['scale'] == 'log' else raw; n = len(x)
        row = {'group': g, 'n': n, 'scale': t['scale'], 'diagnostics': []}
        minimum = {'grubbs': 3, 'generalized_esd': 3 + (t['max_outliers'] or 0), 'rout': 4}[m]
        if n < minimum:
            row.update(status='withheld', reportable=False, reason=f'fewer_than_{minimum}_values', flagged=[]); groups.append(row); continue
        r = grubbs(x, t['alpha'], t['sides']) if m == 'grubbs' else generalized_esd(x, t['alpha'], t['max_outliers']) if m == 'generalized_esd' else rout(x, t['q'])
        if r['status'] == 'withheld':
            row.update(status='withheld', reportable=False, reason=r['reason'], flagged=[]); groups.append(row); continue
        flag = r.pop('flagged')
        for key in ('candidate_index',):
            if key in r: r['candidate_unit_id'] = part.unit_id[r.pop(key)]
        for s in r.get('steps', []) + r.get('tested', []):
            k = s.pop('candidate_index', s.pop('index', None))
            if k is not None: s['unit_id'] = part.unit_id[k]
        keep = np.ones(n, bool); keep[flag] = False
        def summary(v):
            return {'n': int(len(v)), 'mean': float(v.mean()), 'sd': float(v.std(ddof=1)) if len(v) > 1 else None,
                    'median': float(np.median(v))}
        if n < 10: row['diagnostics'].append('fewer_than_10_values_low_power')
        row.update(status='tested', reportable=True, test=r, flagged=[{'unit_id': part.unit_id[i], 'value': float(raw[i])} for i in flag],
                   n_flagged=len(flag), summary_all=summary(raw), summary_without_flagged=summary(raw[keep]) if keep.sum() >= 2 else None)
        groups.append(row)
    flagged_total = sum(g.get('n_flagged', 0) for g in groups)
    notes = ['Flagging only: no value was removed from any analysis. Excluding a flagged value requires a documented assignable cause (for example a recorded pipetting or instrument error), not the test result alone.',
             f"Method {m} with {'alpha ' + str(t['alpha']) if m != 'rout' else 'Q ' + str(t['q'])}" + (f", up to {t['max_outliers']} outliers" if m == 'generalized_esd' else '') +
             (f", {t['sides']}" if m == 'grubbs' else '') + f", declared before screening ({c['design']['prespecified_source']}). Purpose: {c['design']['purpose']}.",
             'Each group is screened separately as one homogeneous sample; values from different conditions, doses or time points are not one sample.',
             'Grubbs and generalized ESD assume an approximately normal sample apart from the outliers; skewed or heavy-tailed data produce false flags. ' +
             ('Values were screened on the declared log scale.' if t['scale'] == 'log' else 'Values were screened on the original scale.'),
             'Report results with all values, and any sensitivity analysis without flagged values, side by side.']
    if m == 'grubbs': notes.append('Grubbs tests a single outlier; repeating it after deleting a point inflates the error rate and suffers masking. Use generalized ESD with a declared maximum for several suspects.')
    if m == 'rout': notes.append('ROUT follows the Motulsky & Brown (2006) description for a location model; it is not a Prism-equivalent implementation. Q is a false discovery level for flagged points, not a per-test alpha.')
    if any('fewer_than_10_values_low_power' in g['diagnostics'] for g in groups): notes.append('Groups with fewer than 10 values have little power; a non-flag does not show the absence of outliers.')
    withheld = [g for g in groups if g['status'] == 'withheld']
    for g in withheld: notes.append(f"Screening withheld for group {g['group']}: {g['reason']}.")
    primary = {'status': 'tested' if not withheld else 'limited', 'reportable': not withheld, 'method': m, 'n_groups': len(groups),
               'n_flagged': flagged_total, 'groups': groups}
    return {'analysis_type': 'outlier_screen', 'primary': primary, 'must_mention': notes,
            'preprocessing': {'pre_existing_exclusions': d[d.exclude][['unit_id', 'group'] + (['exclusion_reason'] if 'exclusion_reason' in d else [])].to_dict('records')},
            'failing_items': [{'group': g['group'], 'reason': g['reason']} for g in withheld],
            'limitations': ['Univariate screening within declared groups only; no regression/curve residual outliers, multivariate or time-series outliers.',
                            'No automatic exclusion and no Prism-equivalence claim.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    d = pd.read_csv(Path(run)/'input.csv', dtype={'unit_id': str, 'group': str})
    if 'exclude' in d: d = d[~d.exclude.astype(str).str.lower().eq('true')]
    groups = [g['group'] for g in result['primary']['groups']]
    flagged = {(g['group'], f['unit_id']) for g in result['primary']['groups'] for f in g.get('flagged', [])}
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(1.8+1.1*len(groups), 3.6), layout='constrained')
        labelled = False
        for i, g in enumerate(groups):
            part = d[d.group == g]; jitter = np.linspace(-.15, .15, len(part)) if len(part) > 1 else np.zeros(1)
            hit = np.array([(g, u) in flagged for u in part.unit_id], bool)
            ax.scatter(i + jitter[~hit], part.value[~hit], s=18, color=pstyle.color(token, 0))
            if hit.any():
                ax.scatter(i + jitter[hit], part.value[hit], s=36, facecolors='none', edgecolors=pstyle.color(token, 1), linewidths=1.4,
                           label=None if labelled else 'flagged (not removed)')
                labelled = True
        if result['primary']['groups'] and result['primary']['groups'][0]['scale'] == 'log': ax.set_yscale('log')
        ax.set(xticks=range(len(groups)), xticklabels=groups, ylabel='Value', title='Outlier screen (' + result['primary']['method'].replace('_', ' ') + ')')
        if flagged: ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'outlier_screen', token, 300)
    return ['outlier_screen']
