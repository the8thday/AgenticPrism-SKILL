"""Area under time curves per independent unit, compared between groups (0.13.2).

Linear trapezoid over a declared interval [t_start, t_end] for each unit, after the declared baseline
handling: none, the unit's value at t_start, or a declared constant (net area: segments below the
baseline count negative). A unit must be observed exactly at t_start and t_end; there is no
extrapolation. Units that end early follow the declared policy: withhold the unit, or shorten every
unit to the common last time. Group comparisons use Welch t on the unit AUCs with Holm p values and
Bonferroni intervals over the declared family.
"""
from copy import deepcopy
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'curve_auc'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'curve_auc', 'input': None, 'source': None,
            'design': {'unit': None, 'unit_rationale': None, 'outcome': None, 'outcome_unit': None, 'time_unit': None},
            'auc': {'interval': None, 'baseline': None, 'baseline_value': None, 'incomplete_policy': None, 'policy_rationale': None},
            'comparison': {'groups': None, 'post_hoc': 'vs_control', 'control_group': None, 'confidence_level': .95},
            'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported curve-AUC configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'curve_auc' or c['schema_version'] != 1:
        raise ValueError('Unsupported curve-AUC schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    for k in ('unit', 'unit_rationale', 'outcome', 'outcome_unit', 'time_unit'):
        _text(c['design'][k], 'design.' + k)
    a = c['auc']; iv = a['interval']
    if not isinstance(iv, list) or len(iv) != 2 or any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in iv) or not iv[0] < iv[1]:
        raise ValueError('Declare auc.interval [t_start, t_end] before seeing the curves')
    if a['baseline'] not in ('none', 'first_value', 'declared_constant'):
        raise ValueError('auc.baseline must be none, first_value or declared_constant')
    if (a['baseline'] == 'declared_constant') != (a['baseline_value'] is not None):
        raise ValueError('baseline_value is required for declared_constant and refused otherwise')
    if a['incomplete_policy'] not in ('withhold_unit', 'common_interval'):
        raise ValueError('Declare auc.incomplete_policy: withhold_unit or common_interval')
    _text(a['policy_rationale'], 'auc.policy_rationale (why units may end early and why the policy is unbiased enough)')
    q = c['comparison']; g = q['groups']
    if not isinstance(g, list) or len(g) < 2 or len(set(g)) != len(g):
        raise ValueError('Declare at least two distinct groups')
    if q['post_hoc'] not in ('vs_control', 'all_pairs'):
        raise ValueError('post_hoc must be vs_control or all_pairs')
    if q['post_hoc'] == 'vs_control' and q['control_group'] not in g:
        raise ValueError('Declare control_group as one of the groups')
    lv = q['confidence_level']
    if isinstance(lv, bool) or not isinstance(lv, (int, float)) or not .5 < lv < 1:
        raise ValueError('confidence_level must lie in (0.5, 1)')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'unit_id', 'group', 'time', 'value'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    for k in ('time', 'value'):
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'Nonfinite {k}; leave missing time points out instead of encoding them')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    if (d.groupby('unit_id').group.nunique() > 1).any():
        raise ValueError('A unit appears in more than one group')
    if d.duplicated(['unit_id', 'time']).any():
        raise ValueError('Duplicate time within a unit; average technical replicates upstream only if that is the declared unit value')
    if set(d.group) != set(cfg['comparison']['groups']):
        raise ValueError('Input groups do not match declared groups')
    return d


def trapezoid(t, y):
    t, y = np.asarray(t, float), np.asarray(y, float)
    return float(np.sum(np.diff(t) * (y[:-1] + y[1:]) / 2))


def _holm(p):
    order = np.argsort(p); m = len(p); out = np.empty(m); running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., (m - rank) * p[i])); out[i] = running
    return out.tolist()


def compute(d, cfg):
    a, q = cfg['auc'], cfg['comparison']; groups = q['groups']; t0, t1 = a['interval']; lv = q['confidence_level']
    used = d[~d.exclude].sort_values(['unit_id', 'time'])
    last = used.groupby('unit_id').time.max(); first = used.groupby('unit_id').time.min()
    must, diag = [], []
    if a['incomplete_policy'] == 'common_interval':
        end = float(min(t1, last.min()))
        if end < t1:
            must.append(f'Units ended early, so every AUC uses the common interval [{t0:g}, {end:g}] instead of the declared [{t0:g}, {t1:g}].')
        if end <= t0:
            raise ValueError('No common interval remains after the earliest dropout')
    else:
        end = float(t1)
    units = []
    for uid, u in used.groupby('unit_id', sort=False):
        g = u.group.iloc[0]; row = {'unit_id': uid, 'group': g, 'auc': None}
        times = u.time.to_numpy(float)
        if not np.any(np.isclose(times, t0)):
            row['status'] = 'no_observation_at_interval_start'
        elif not np.any(np.isclose(times, end)):
            row['status'] = 'ended_before_interval_end' if times.max() < end else 'no_observation_at_interval_end'
        else:
            w = u[(u.time >= t0 - 1e-12) & (u.time <= end + 1e-12)]; y = w.value.to_numpy(float)
            base = 0. if a['baseline'] == 'none' else float(y[0]) if a['baseline'] == 'first_value' else float(a['baseline_value'])
            row.update(status='computed', auc=trapezoid(w.time, y - base), baseline=base, n_points=int(len(w)), interval=[float(t0), end])
        units.append(row)
    ut = pd.DataFrame(units)
    withheld = ut[ut.status != 'computed']
    if len(withheld):
        diag.append('units_without_auc')
        by = withheld.groupby('group').size().reindex(groups, fill_value=0).to_dict()
        must.append(f'{len(withheld)} unit(s) have no AUC on the interval ({by}); excluding them can bias the comparison when dropout relates to the outcome.')
        # groups with no withheld unit must count as 0, not NaN (NaN silently disabled this check)
        frac = withheld.groupby('group').size().reindex(groups, fill_value=0) / ut.groupby('group').size().reindex(groups)
        if frac.max() - frac.min() > .2:
            diag.append('unequal_withholding_between_groups')
            must.append('Units without AUC differ by more than 20 percentage points between groups (for example humane-endpoint dropout); AUC comparisons are then biased. Consider tumor-growth (mixed model) or time-to-event analysis.')
    ok = ut[ut.status == 'computed']
    summaries = []
    for g in groups:
        v = ok[ok.group == g].auc.to_numpy(float); row = {'group': g, 'n': int(len(v)), 'mean_auc': float(v.mean()) if len(v) else None,
                                                          'sd': float(v.std(ddof=1)) if len(v) > 1 else None, 'ci': [None, None]}
        if len(v) > 1:
            h = stats.t.ppf((1 + lv) / 2, len(v) - 1) * row['sd'] / np.sqrt(len(v)); row['ci'] = [row['mean_auc'] - h, row['mean_auc'] + h]
        summaries.append(row)
    pairs = [(q['control_group'], g) for g in groups if g != q['control_group']] if q['post_hoc'] == 'vs_control' else list(combinations(groups, 2))
    m = len(pairs); comps = []
    for a_, b_ in pairs:
        x, y = ok[ok.group == a_].auc.to_numpy(float), ok[ok.group == b_].auc.to_numpy(float)
        row = {'contrast': f'{b_} - {a_}', 'estimate': None, 'ci': [None, None], 'p_unadjusted': None}
        if len(x) < 2 or len(y) < 2:
            row['status'] = 'withheld_fewer_than_two_units'; diag.append('comparison_withheld')
        else:
            vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y); se = np.sqrt(vx + vy)
            df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)); est = float(y.mean() - x.mean())
            qn = stats.t.ppf(1 - (1 - lv) / (2 * m), df)
            row.update(status='estimated', estimate=est, standard_error=float(se), df=float(df), ci=[est - qn * se, est + qn * se],
                       p_unadjusted=float(2 * stats.t.sf(abs(est / se), df)), method='Welch t on unit AUCs')
        comps.append(row)
    est = [i for i, r in enumerate(comps) if r['status'] == 'estimated']
    for i, p in zip(est, _holm(np.array([comps[i]['p_unadjusted'] for i in est]))):
        comps[i]['p_holm'] = p
    must.insert(0, f"AUC is computed per {cfg['design']['unit']} by the linear trapezoid with baseline '{a['baseline']}'; areas below the baseline count negative. Group comparisons use one AUC per unit.")
    return {'analysis_type': 'curve_auc', 'primary': {'interval_used': [float(t0), end], 'declared_interval': [float(t0), float(t1)], 'baseline': a['baseline'],
            'incomplete_policy': a['incomplete_policy'], 'group_summaries': summaries, 'comparisons': comps, 'units': units, 'design': cfg['design'],
            'confidence_level': lv, 'interval_adjustment': 'Bonferroni across the declared family' if m > 1 else 'none (single contrast)'},
            'must_mention': must, 'failing_items': [r for r in comps if r['status'] != 'estimated'],
            'limitations': ['Linear trapezoid between observed times; no curve fitting, smoothing or extrapolation.',
                            'AUC summarises the whole interval; it does not show when groups diverge. Equal sampling schedules across groups are assumed.']}


def render(run, result, style):
    """Group mean curves (from saved inputs) on the used interval."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); p = result['primary']
    d = pd.read_csv(Path(run)/'input.csv', dtype={'unit_id': str, 'group': str}, keep_default_na=False)
    if 'exclude' in d:
        d = d[~d.exclude.astype(str).str.lower().eq('true')]
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.8, 3.8), layout='constrained')
        for i, s in enumerate(p['group_summaries']):
            g = d[d.group == s['group']]; col = pstyle.color(token, i)
            for _, u in g.groupby('unit_id'):
                ax.plot(u.time, u.value, '-', color=col, lw=.6, alpha=.35)
            mean = g.groupby('time').value.mean(); ax.plot(mean.index, mean.values, '-', color=col, lw=2, label=f"{s['group']} (n = {s['n']})")
        for x in p['interval_used']:
            ax.axvline(x, color=pstyle.MUTED, lw=.8, ls=(0, (3, 2)))
        dz = p['design']; ax.set(xlabel=f"Time ({dz['time_unit']})", ylabel=f"{dz['outcome']} ({dz['outcome_unit']})", title='Unit curves and group means (dashed: AUC interval)')
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'auc_curves', token, 300)
    return ['auc_curves']
