"""Nested comparisons: technical replicates within independent units (0.13.2).

Model: y_gur = mu_g + b_u + e_gur, b_u ~ N(0, tau), e ~ N(0, sigma^2), each unit in one group.
Fixed group means use cell-means coding; REML fit and Satterthwaite t/F as lmerTest
(`repeated.fit_random_intercept`, `mixed_inference.satterthwaite`). Contrasts form one
predeclared family: Holm-adjusted p values and Bonferroni-adjusted intervals. The unit-means
analysis (average replicates within each unit, then Welch t for two groups or classical
one-way ANOVA contrasts) is always reported; it is the primary result when the between-unit
variance is estimated at the boundary (singular fit), where Satterthwaite df are undefined.
"""
from copy import deepcopy
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
from .repeated import fit_random_intercept
from .mixed_inference import RandomInterceptModel, satterthwaite

TYPES = {'nested_comparison'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'nested_comparison', 'input': None, 'source': None,
            'design': {'unit': None, 'unit_rationale': None, 'replicate': None, 'replicate_rationale': None,
                       'outcome': None, 'outcome_unit': None},
            'comparison': {'groups': None, 'post_hoc': 'vs_control', 'control_group': None, 'confidence_level': .95},
            'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported nested-comparison configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'nested_comparison' or c['schema_version'] != 1:
        raise ValueError('Unsupported nested-comparison schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    for k in ('unit', 'unit_rationale', 'replicate', 'replicate_rationale', 'outcome', 'outcome_unit'):
        _text(c['design'][k], 'design.' + k + ' (unit_rationale: why units are independent; replicate_rationale: what the technical replicates are)')
    q = c['comparison']; g = q['groups']
    if not isinstance(g, list) or len(g) < 2 or len(set(g)) != len(g) or any(not isinstance(x, str) or not x.strip() for x in g):
        raise ValueError('Declare at least two distinct groups')
    if q['post_hoc'] not in ('vs_control', 'all_pairs', 'none'):
        raise ValueError('post_hoc must be vs_control, all_pairs or none')
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
    req = {'group', 'unit_id', 'value'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    for k in ('group', 'unit_id'):
        if d[k].str.strip().eq('').any():
            raise ValueError(f'Missing {k}')
    d['value'] = pd.to_numeric(d.value, errors='raise')
    if not np.isfinite(d.value).all():
        raise ValueError('Nonfinite value; do not encode missing replicates as numbers')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    if (d.groupby('unit_id').group.nunique() > 1).any():
        raise ValueError('A unit appears in more than one group: that is a repeated-measures design (use repeated-measures)')
    groups = cfg['comparison']['groups']
    if set(d.group) != set(groups):
        raise ValueError('Input groups do not match declared groups')
    used = d[~d.exclude]
    units = used.groupby('group').unit_id.nunique().reindex(groups).fillna(0)
    if units.min() < 2:
        raise ValueError('Every group needs at least two independent units')
    if used.groupby('unit_id').size().max() < 2:
        raise ValueError('No unit has replicates: this is not a nested design (use group-comparison on the unit values)')
    return d


def _contrasts(groups, q):
    k = len(groups); rows = []
    pairs = [(q['control_group'], g) for g in groups if g != q['control_group']] if q['post_hoc'] == 'vs_control' else \
        list(combinations(groups, 2)) if q['post_hoc'] == 'all_pairs' else []
    for a, b in pairs:  # estimate b - a
        l = np.zeros(k); l[groups.index(b)] = 1; l[groups.index(a)] = -1
        rows.append((f'{b} - {a}', l))
    return rows


def _holm(p):
    order = np.argsort(p); m = len(p); out = np.empty(m); running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., (m - rank) * p[i])); out[i] = running
    return out.tolist()


def unit_means_analysis(means, groups, contrasts, level):
    """Average replicates within units, then two-group Welch or classical one-way ANOVA contrasts."""
    vals = [means[means.group == g].value.to_numpy() for g in groups]; m = max(len(contrasts), 1); rows = []
    if len(groups) == 2:
        a, b = vals; va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b); se = np.sqrt(va + vb)
        df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        est = b.mean() - a.mean(); q = stats.t.ppf(1 - (1 - level) / 2, df)
        rows.append({'contrast': f'{groups[1]} - {groups[0]}', 'estimate': float(est), 'standard_error': float(se), 'df': float(df),
                     'p_unadjusted': float(2 * stats.t.sf(abs(est / se), df)), 'ci_low': float(est - q * se), 'ci_high': float(est + q * se), 'method': 'Welch t on unit means'})
        omnibus = {'method': 'Welch t on unit means', 'p_value': rows[0]['p_unadjusted']}
    else:
        n = np.array([len(v) for v in vals]); dfe = int(n.sum() - len(groups))
        mse = sum(((v - v.mean()) ** 2).sum() for v in vals) / dfe; mu = np.array([v.mean() for v in vals])
        f = stats.f_oneway(*vals)
        omnibus = {'method': 'One-way ANOVA on unit means', 'f_statistic': float(f.statistic), 'df_numerator': len(groups) - 1, 'df_denominator': dfe, 'p_value': float(f.pvalue)}
        for name, l in contrasts:
            est = float(l @ mu); se = float(np.sqrt(mse * np.sum(l ** 2 / n))); q = stats.t.ppf(1 - (1 - level) / (2 * m), dfe)
            rows.append({'contrast': name, 'estimate': est, 'standard_error': se, 'df': dfe, 'p_unadjusted': float(2 * stats.t.sf(abs(est / se), dfe)),
                         'ci_low': est - q * se, 'ci_high': est + q * se, 'method': 'pooled-variance t on unit means'})
    if rows:
        for r, p in zip(rows, _holm(np.array([r['p_unadjusted'] for r in rows]))):
            r['p_holm'] = p
    return {'omnibus': omnibus, 'contrasts': rows,
            'interval_adjustment': 'none (single contrast)' if len(rows) <= 1 else 'Bonferroni across the declared family'}


def compute(d, cfg):
    q = cfg['comparison']; groups = q['groups']; lv = q['confidence_level']
    used = d[~d.exclude].copy()
    used['group'] = pd.Categorical(used.group, groups)
    used = used.sort_values(['group', 'unit_id'], kind='stable').reset_index(drop=True)
    y = used.value.to_numpy(float); units = used.unit_id.to_numpy()
    x = np.column_stack([(used.group == g).to_numpy(float) for g in groups])
    contrasts = _contrasts(groups, q)
    reps = used.groupby('unit_id').size(); means = used.groupby(['group', 'unit_id'], observed=True).value.mean().reset_index()
    summaries = []
    for g in groups:
        s = used[used.group == g]; um = means[means.group == g].value
        summaries.append({'group': g, 'units': int(s.unit_id.nunique()), 'observations': int(len(s)),
                          'replicates_per_unit': sorted(set(int(v) for v in reps[s.unit_id.unique()])),
                          'mean_of_unit_means': float(um.mean()), 'sd_of_unit_means': float(um.std(ddof=1)) if len(um) > 1 else None})
    balanced = reps.nunique() == 1
    unit_means = unit_means_analysis(means, groups, contrasts, lv)
    fit = fit_random_intercept(y, x, units)
    diag, must = [], []
    mixed = {'status': 'failed'}
    if fit['ok']:
        tau, sigma = fit['tau'], fit['sigma']; icc = tau / (tau + sigma); m_bar = float(reps.mean())
        mixed = {'status': 'estimated', 'between_unit_variance': tau, 'within_unit_variance': sigma, 'icc': float(icc),
                 'design_effect': float(1 + (m_bar - 1) * icc), 'effective_n_total': float(len(y) / (1 + (m_bar - 1) * icc)),
                 'group_means': [{'group': g, 'estimate': float(b), 'standard_error': float(np.sqrt(fit['covariance'][i, i]))} for i, (g, b) in enumerate(zip(groups, fit['beta']))],
                 'boundary': bool(fit['boundary'])}
        if fit['boundary']:
            mixed['status'] = 'boundary'; diag.append('between_unit_variance_at_boundary')
            must.append('The between-unit variance was estimated at zero (singular fit); Satterthwaite df are undefined, so the unit-means analysis is the primary result.')
        else:
            k = len(groups); joint = np.column_stack([-np.ones(k - 1), np.eye(k - 1)])
            try:
                rows, omnibus, _ = satterthwaite(RandomInterceptModel(y, x, units), (tau, sigma), fit['beta'],
                                                 [l for _, l in contrasts], joint=joint, level=lv, n_family=max(len(contrasts), 1))
                for (name, _), r in zip(contrasts, rows):
                    r['contrast'] = name
                if rows:
                    for r, p in zip(rows, _holm(np.array([r['p_unadjusted'] for r in rows]))):
                        r['p_holm'] = p
                mixed.update(omnibus=omnibus, contrasts=rows, inference='REML, Satterthwaite t/F (lmerTest)',
                             interval_adjustment='none (single contrast)' if len(rows) <= 1 else 'Bonferroni across the declared family')
            except ArithmeticError as e:
                mixed['status'] = 'satterthwaite_failed'; diag.append('satterthwaite_failed'); must.append(f'Satterthwaite inference failed ({e}); the unit-means analysis is the primary result.')
    else:
        diag.append('mixed_model_failed'); must.append('The mixed model did not converge; the unit-means analysis is the primary result.')
    primary = 'mixed_model' if mixed['status'] == 'estimated' else 'unit_means'
    few = [s['group'] for s in summaries if s['units'] < 3]
    if few:
        diag.append('fewer_than_three_units_in_a_group'); must.append(f"Groups {few} have fewer than three independent units; intervals are wide and variance estimates unstable.")
    must.insert(0, f"The sample size is the number of {cfg['design']['unit']}s ({int(used.unit_id.nunique())}), not the number of {cfg['design']['replicate']}s ({len(used)}); technical replicates only reduce within-unit noise.")
    if balanced and primary == 'mixed_model':
        must.append('With equal replicates per unit, the mixed model and the unit-means analysis give the same group comparison.')
    return {'analysis_type': 'nested_comparison', 'primary': {'primary_inference': primary, 'balanced': bool(balanced), 'groups': summaries,
            'mixed_model': mixed, 'unit_means_analysis': unit_means, 'design': cfg['design'], 'post_hoc': q['post_hoc'], 'confidence_level': lv},
            'unit_means': means.assign(group=means.group.astype(str)).to_dict('records'),
            'must_mention': must, 'failing_items': [{'item': x} for x in diag if x in ('mixed_model_failed', 'satterthwaite_failed')],
            'limitations': ['One level of nesting (replicates within units within groups); deeper hierarchies and crossed factors need variance-components or repeated-measures.',
                            'Gaussian errors with equal within-unit variance across groups; no heteroscedastic mixed model.']}


def render(run, result, style):
    """Unit means and individual replicates per group from saved values only."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    data = pd.read_csv(Path(run)/'input.csv', dtype={'unit_id': str, 'group': str}, keep_default_na=False)
    data = data[~data.get('exclude', pd.Series('false', index=data.index)).astype(str).str.lower().eq('true')]
    groups = [s['group'] for s in result['primary']['groups']]; means = pd.DataFrame(result['unit_means'])
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(1.6 + 1.3 * len(groups), 3.8), layout='constrained')
        rng = np.random.default_rng(0)
        for i, g in enumerate(groups):
            r = data[data.group == g]; m = means[means.group == g]
            ax.plot(i + rng.uniform(-.18, .18, len(r)), pd.to_numeric(r.value), '.', ms=3, color=pstyle.MUTED, alpha=.6)
            ax.plot(np.full(len(m), i), m.value, 'o', ms=6, color=pstyle.color(token, i), mfc='white', mew=1.4)
            ax.hlines(m.value.mean(), i - .3, i + .3, color=pstyle.color(token, i), lw=2)
        d = result['primary']['design']
        ax.set(xticks=range(len(groups)), xticklabels=groups, ylabel=f"{d['outcome']} ({d['outcome_unit']})",
               title=f"Open circles: {d['unit']} means; dots: {d['replicate']}s")
        pstyle.save(fig, out/'nested_units', token, 300)
    return ['nested_units']
