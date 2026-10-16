"""Empirical ROC analysis against a declared binary reference standard (0.13.6).

One row per independent unit. AUC is the Mann-Whitney probability with ties
counted one half, oriented by a declared marker direction that is never flipped
after seeing the data. DeLong (1988) placement-value variance; the AUC interval
is Wald on the logit scale, the paired AUC difference interval Wald on the
original scale. Sensitivity/specificity at declared cut-offs use Clopper-Pearson
intervals. The Youden cut-off is optional, in-sample and exploratory only.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'roc_analysis'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'roc_analysis', 'input': None, 'source': None,
            'design': {'unit': None, 'rationale': None, 'reference_standard': None, 'reference_source': None},
            'markers': None, 'comparison': None, 'thresholds': [], 'exploratory_youden': False,
            'confidence_level': .95, 'report': {'plot_style': 'prism_like'}}
MARKER_KEYS = {'column', 'label', 'unit', 'direction', 'direction_source'}
RESERVED = {'unit_id', 'condition', 'exclude', 'exclusion_reason'}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported ROC configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'roc_analysis': raise ValueError('Unsupported ROC schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    m = c['markers']
    if not isinstance(m, list) or not 1 <= len(m) <= 2: raise ValueError('Declare one or two markers')
    for marker in m:
        if not isinstance(marker, dict) or set(marker) != MARKER_KEYS: raise ValueError('Each marker needs ' + ', '.join(sorted(MARKER_KEYS)))
        if not all(_text(marker[k]) for k in MARKER_KEYS): raise ValueError('Marker declarations must be non-empty text')
        if marker['column'] in RESERVED: raise ValueError('Marker column name is reserved')
        if marker['direction'] not in ('higher_indicates_condition', 'lower_indicates_condition'):
            raise ValueError('Marker direction must be declared before analysis; it is never chosen from the data')
    if len({x['column'] for x in m}) != len(m): raise ValueError('Markers must be different columns')
    if c['comparison'] not in (None, 'paired_delong'): raise ValueError('comparison must be null or paired_delong')
    if c['comparison'] and len(m) != 2: raise ValueError('paired_delong needs exactly two markers measured on the same units')
    if not isinstance(c['thresholds'], list): raise ValueError('thresholds must be a list')
    columns = {x['column'] for x in m}
    for t in c['thresholds']:
        if not isinstance(t, dict) or set(t) != {'marker', 'cutoff', 'source'} or t['marker'] not in columns or not _text(t['source']):
            raise ValueError('Each threshold declares marker, cutoff and its source')
        if type(t['cutoff']) not in (int, float) or not np.isfinite(t['cutoff']): raise ValueError('Threshold cutoff must be finite')
    if not isinstance(c['exploratory_youden'], bool): raise ValueError('exploratory_youden must be true or false')
    lv = c['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'unit_id': str})
    need = {'unit_id', 'condition'} | {m['column'] for m in c['markers']}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if d.unit_id.isna().any() or d.unit_id.duplicated().any(): raise ValueError('One row per independent unit; unit_id must be unique')
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded unit needs an exclusion_reason')
    else: d['exclude'] = False
    used = d[~d.exclude]
    if not used.condition.isin([0, 1]).all(): raise ValueError('condition must be 1 (present by the reference standard) or 0')
    for m in c['markers']:
        v = pd.to_numeric(used[m['column']], errors='coerce')
        if v.isna().any() or not np.isfinite(v).all():
            raise ValueError(f"Marker {m['column']} has missing or non-numeric values; no silent deletion or imputation")
    if min((used.condition == 1).sum(), (used.condition == 0).sum()) < 2:
        raise ValueError('At least two units with and two without the condition are required')
    return d


def _psi(pos, neg):
    diff = pos[:, None] - neg[None, :]
    return (diff > 0) + .5*(diff == 0)


def delong(pos, neg):
    """AUC and DeLong placement values (structural components) for one oriented marker."""
    k = _psi(pos, neg)
    return float(k.mean()), k.mean(axis=1), k.mean(axis=0)


def _oriented(values, direction):
    v = np.asarray(values, float)
    return v if direction == 'higher_indicates_condition' else -v


def clopper_pearson(x, n, level):
    a = 1 - level
    lo = 0. if x == 0 else float(stats.beta.ppf(a/2, x, n-x+1))
    hi = 1. if x == n else float(stats.beta.ppf(1-a/2, x+1, n-x))
    return [lo, hi]


def marker_result(u, marker, level):
    s = _oriented(u[marker['column']], marker['direction']); y = u.condition.to_numpy(int)
    pos, neg = s[y == 1], s[y == 0]; m, n = len(pos), len(neg)
    auc, v10, v01 = delong(pos, neg)
    var = float(v10.var(ddof=1)/m + v01.var(ddof=1)/n); se = float(np.sqrt(max(var, 0.)))
    z = float(stats.norm.ppf((1+level)/2))
    tied = float((pos[:, None] == neg[None, :]).mean())
    mw = stats.mannwhitneyu(pos, neg, alternative='two-sided', method='asymptotic')
    r = {'marker': marker['column'], 'label': marker['label'], 'unit': marker['unit'], 'direction': marker['direction'],
         'direction_source': marker['direction_source'], 'n_condition': m, 'n_no_condition': n, 'auc': auc,
         'delong_se': se, 'tied_pair_fraction': tied, 'mann_whitney_p_two_sided': float(mw.pvalue),
         'interval_method': 'DeLong variance, Wald interval on logit(AUC)', 'diagnostics': []}
    if min(m, n) < 10: r['diagnostics'].append('small_class_size')
    if auc < .5: r['diagnostics'].append('auc_below_half_under_declared_direction')
    if tied > .2: r['diagnostics'].append('many_tied_pairs')
    if auc in (0., 1.):
        r.update(status='limited', reportable=False, ci=None, reason='complete_separation_no_delong_interval')
    elif se <= 0:
        r.update(status='limited', reportable=False, ci=None, reason='zero_delong_variance')
    else:
        lg = np.log(auc/(1-auc)); h = z*se/(auc*(1-auc))
        r.update(status='estimated', reportable=True, ci=[float(1/(1+np.exp(-(lg-h)))), float(1/(1+np.exp(-(lg+h))))])
    return r, (pos, neg, v10, v01)


def paired(parts, results, level):
    (p1, n1, a10, a01), (p2, n2, b10, b01) = parts
    m, n = len(a10), len(a01)
    s10 = np.cov(np.vstack([a10, b10]), ddof=1); s01 = np.cov(np.vstack([a01, b01]), ddof=1)
    var = float((s10[0, 0]+s10[1, 1]-2*s10[0, 1])/m + (s01[0, 0]+s01[1, 1]-2*s01[0, 1])/n)
    diff = results[1]['auc'] - results[0]['auc']
    r = {'method': 'paired_delong', 'estimand': 'AUC(marker 2) - AUC(marker 1) on the same units',
         'marker_1': results[0]['marker'], 'marker_2': results[1]['marker'], 'difference': float(diff),
         'correlation_condition': float(np.corrcoef(a10, b10)[0, 1]) if a10.std() > 0 and b10.std() > 0 else None,
         'correlation_no_condition': float(np.corrcoef(a01, b01)[0, 1]) if a01.std() > 0 and b01.std() > 0 else None}
    if var <= 0 or not np.isfinite(var):
        r.update(status='withheld', reportable=False, reason='nonpositive_difference_variance'); return r
    se = float(np.sqrt(var)); z = float(stats.norm.ppf((1+level)/2)); stat = diff/se
    r.update(status='estimated', reportable=True, standard_error=se, ci=[float(diff-z*se), float(diff+z*se)],
             z_statistic=float(stat), p_two_sided=float(2*stats.norm.sf(abs(stat))))
    return r


def roc_points(u, marker):
    s = _oriented(u[marker['column']], marker['direction']); y = u.condition.to_numpy(int)
    m, n = (y == 1).sum(), (y == 0).sum(); rows = [{'marker': marker['column'], 'threshold': None, 'rule': 'none positive', 'sensitivity': 0., 'specificity': 1.}]
    for t in np.unique(s)[::-1]:
        call = s >= t
        cutoff = float(t) if marker['direction'] == 'higher_indicates_condition' else float(-t)
        rows.append({'marker': marker['column'], 'threshold': cutoff,
                     'rule': ('>= ' if marker['direction'] == 'higher_indicates_condition' else '<= ') + repr(cutoff),
                     'sensitivity': float(call[y == 1].sum()/m), 'specificity': float((~call)[y == 0].sum()/n)})
    return rows


def compute(d, c):
    u = d[~d.exclude].reset_index(drop=True); level = c['confidence_level']
    results, parts = [], []
    for marker in c['markers']:
        r, p = marker_result(u, marker, level); results.append(r); parts.append(p)
    curve = [row for marker in c['markers'] for row in roc_points(u, marker)]
    thresholds = []
    for t in c['thresholds']:
        marker = next(x for x in c['markers'] if x['column'] == t['marker'])
        v = u[marker['column']].to_numpy(float); y = u.condition.to_numpy(int)
        call = v >= t['cutoff'] if marker['direction'] == 'higher_indicates_condition' else v <= t['cutoff']
        tp, fn = int((call & (y == 1)).sum()), int((~call & (y == 1)).sum())
        tn, fp = int((~call & (y == 0)).sum()), int((call & (y == 0)).sum())
        thresholds.append({'marker': marker['column'], 'cutoff': t['cutoff'], 'source': t['source'],
                           'rule': ('>= ' if marker['direction'] == 'higher_indicates_condition' else '<= ') + repr(t['cutoff']) + ' called positive',
                           'true_positive': tp, 'false_negative': fn, 'true_negative': tn, 'false_positive': fp,
                           'sensitivity': tp/(tp+fn), 'sensitivity_ci': clopper_pearson(tp, tp+fn, level),
                           'specificity': tn/(tn+fp), 'specificity_ci': clopper_pearson(tn, tn+fp, level),
                           'interval_method': 'Clopper-Pearson exact (conservative)'})
    youden = []
    if c['exploratory_youden']:
        for marker in c['markers']:
            pts = [p for p in curve if p['marker'] == marker['column'] and p['threshold'] is not None]
            best = max(pts, key=lambda p: p['sensitivity']+p['specificity']-1)
            youden.append({'marker': marker['column'], 'status': 'exploratory', 'reportable': False, 'cutoff': best['threshold'],
                           'rule': best['rule'], 'in_sample_sensitivity': best['sensitivity'], 'in_sample_specificity': best['specificity'],
                           'youden_j': best['sensitivity']+best['specificity']-1,
                           'reason': 'data_selected_cutoff_optimistic_requires_independent_validation'})
    comparison = paired(parts, results, level) if c['comparison'] else None
    notes = ['Independent unit: ' + c['design']['unit'] + '. Reference standard: ' + c['design']['reference_standard'] +
             ' (' + c['design']['reference_source'] + '). Accuracy is relative to that reference, which is assumed error-free.',
             'Marker directions were declared before analysis and are never flipped; AUC below 0.5 under the declared direction is reported as observed.',
             'AUC is a rank-based discrimination summary; it does not choose a clinical cut-off, establish calibration or clinical utility, or transfer to another population.',
             'Sensitivity and specificity are conditional on the observed condition groups. Predictive values depend on prevalence and are not reported.']
    for r in results:
        if not r['reportable']: notes.append(f"AUC interval for {r['marker']} withheld: {r['reason']}; the observed AUC is descriptive only.")
        if 'small_class_size' in r['diagnostics']: notes.append(f"{r['marker']}: fewer than 10 units in a condition group; DeLong/Wald intervals are approximate.")
        if 'auc_below_half_under_declared_direction' in r['diagnostics']: notes.append(f"{r['marker']}: AUC below 0.5 under the declared direction; do not reverse the direction post hoc.")
    if youden: notes.append('Youden cut-offs are selected from these data and are optimistic; they are exploratory, carry no interval and need independent validation.')
    if comparison and not comparison['reportable']: notes.append('Paired AUC comparison withheld: ' + comparison['reason'])
    if comparison: notes.append('A non-significant AUC difference does not establish equivalence of the markers.')
    failing = [{'marker': r['marker'], 'reason': r['reason']} for r in results if not r['reportable']]
    if comparison and not comparison['reportable']: failing.append({'comparison': 'paired_delong', 'reason': comparison['reason']})
    primary = {'status': 'estimated' if not failing else 'limited', 'reportable': not failing, 'confidence_level': level,
               'n_units': len(u), 'n_excluded': int(d.exclude.sum()), 'markers': results, 'comparison': comparison,
               'thresholds': thresholds, 'exploratory_youden': youden}
    return {'analysis_type': 'roc_analysis', 'primary': primary, 'roc_curve': curve,
            'preprocessing': {'excluded_units': d[d.exclude][['unit_id'] + (['exclusion_reason'] if 'exclusion_reason' in d else [])].to_dict('records')},
            'must_mention': notes, 'failing_items': failing,
            'limitations': ['Empirical (nonparametric) ROC only; no binormal/smoothed curves, partial AUC, covariate adjustment, verification-bias correction or multi-reader designs.',
                            'Independent units with one reference-standard classification each; clustered or repeated measurements are outside scope.',
                            'No Prism, pROC or software-equivalence claim beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(4.2, 4.2), layout='constrained')
        ax.plot([0, 1], [0, 1], ls='--', color=pstyle.MUTED, lw=1)
        for i, r in enumerate(result['primary']['markers']):
            # Successive operating points joined by straight segments (tied values give diagonal segments).
            pts = [p for p in result['roc_curve'] if p['marker'] == r['marker']]
            ax.plot([1-p['specificity'] for p in pts], [p['sensitivity'] for p in pts], color=pstyle.color(token, i),
                    marker='o', ms=2.5, label=f"{r['label']} (AUC {r['auc']:.3f})")
        ax.set(xlim=(0, 1), ylim=(0, 1.01), xlabel='1 - specificity', ylabel='Sensitivity', title='Empirical ROC')
        ax.legend(loc='lower right', frameon=False)
        pstyle.save(fig, out/'roc_curve', token, 300)
    return ['roc_curve']
