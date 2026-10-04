"""ANCOVA: comparison of independent groups adjusted for prespecified pre-treatment covariates (0.13.2).

Model y = mu_g + sum_k beta_k (x_k - xbar_k) + e by ordinary least squares, cell-means coding, common
slopes. Adjusted means are the group means at the overall covariate means (emmeans convention).
Homogeneity of slopes is tested beforehand by the extra-sum-of-squares F for all group x covariate
interactions at a declared alpha; if it is rejected, adjusted comparisons are withheld because the group
difference then depends on the covariate value. Contrasts form one declared family: Holm p values and
Bonferroni intervals. The unadjusted comparison is reported alongside.
"""
from copy import deepcopy
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'ancova'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'ancova', 'input': None, 'source': None,
            'design': {'unit': None, 'unit_rationale': None, 'outcome': None, 'outcome_unit': None,
                       'covariates': None, 'covariate_timing': None, 'covariate_rationale': None},
            'comparison': {'groups': None, 'post_hoc': 'vs_control', 'control_group': None, 'confidence_level': .95,
                           'slope_homogeneity_alpha': .05},
            'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported ANCOVA configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'ancova' or c['schema_version'] != 1:
        raise ValueError('Unsupported ANCOVA schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    d = c['design']
    for k in ('unit', 'unit_rationale', 'outcome', 'outcome_unit', 'covariate_rationale'):
        _text(d[k], 'design.' + k)
    cov = d['covariates']
    if not isinstance(cov, list) or not 1 <= len(cov) <= 3 or len(set(cov)) != len(cov) or any(not isinstance(x, str) or not x.strip() for x in cov) \
            or set(cov) & {'unit_id', 'group', 'outcome'}:
        raise ValueError('Declare 1-3 distinct covariate column names')
    if d['covariate_timing'] != 'before_treatment':
        raise ValueError('Covariates must be declared covariate_timing: before_treatment (measured before randomization or dosing); '
                         'a post-treatment covariate can absorb the treatment effect')
    q = c['comparison']; g = q['groups']
    if not isinstance(g, list) or len(g) < 2 or len(set(g)) != len(g) or any(not isinstance(x, str) or not x.strip() for x in g):
        raise ValueError('Declare at least two distinct groups')
    if q['post_hoc'] not in ('vs_control', 'all_pairs', 'none'):
        raise ValueError('post_hoc must be vs_control, all_pairs or none')
    if q['post_hoc'] == 'vs_control' and q['control_group'] not in g:
        raise ValueError('Declare control_group as one of the groups')
    for k, lo, hi in (('confidence_level', .5, 1), ('slope_homogeneity_alpha', 0, .5)):
        v = q[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not lo < v < hi:
            raise ValueError(f'{k} out of range')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    cov = cfg['design']['covariates']; req = {'unit_id', 'group', 'outcome', *cov}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    if d.unit_id.duplicated().any():
        raise ValueError('One row per independent unit; repeated measurements need repeated-measures or MMRM')
    for k in ('outcome', *cov):
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'{k} must be finite for every unit (no silent case deletion)')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    if set(d.group) != set(cfg['comparison']['groups']):
        raise ValueError('Input groups do not match declared groups')
    if d[~d.exclude].groupby('group').size().min() < 2:
        raise ValueError('Every group needs at least two units')
    return d


def _ols(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None); resid = y - X @ beta
    df = len(y) - np.linalg.matrix_rank(X); sse = float(resid @ resid)
    return beta, sse, df


def _holm(p):
    order = np.argsort(p); m = len(p); out = np.empty(m); running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., (m - rank) * p[i])); out[i] = running
    return out.tolist()


def _pairs(groups, q):
    if q['post_hoc'] == 'vs_control':
        return [(q['control_group'], g) for g in groups if g != q['control_group']]
    return list(combinations(groups, 2)) if q['post_hoc'] == 'all_pairs' else []


def _contrast_rows(beta, cov, df, groups, pairs, level, k_cols):
    m = max(len(pairs), 1); rows = []
    for a, b in pairs:
        l = np.zeros(k_cols); l[groups.index(b)] = 1; l[groups.index(a)] = -1
        est = float(l @ beta); se = float(np.sqrt(l @ cov @ l)); q = stats.t.ppf(1 - (1 - level) / (2 * m), df)
        rows.append({'contrast': f'{b} - {a}', 'estimate': est, 'standard_error': se, 'df': int(df), 'statistic': est / se,
                     'p_unadjusted': float(2 * stats.t.sf(abs(est / se), df)), 'ci_low': est - q * se, 'ci_high': est + q * se})
    if rows:
        for r, p in zip(rows, _holm(np.array([r['p_unadjusted'] for r in rows]))):
            r['p_holm'] = p
    return rows


def compute(d, cfg):
    q = cfg['comparison']; groups = q['groups']; cov = cfg['design']['covariates']; lv = q['confidence_level']
    used = d[~d.exclude].reset_index(drop=True); y = used.outcome.to_numpy(float); n = len(y); k = len(groups)
    G = np.column_stack([(used.group == g).to_numpy(float) for g in groups])
    xbar = used[cov].mean().to_numpy(float); Xc = used[cov].to_numpy(float) - xbar
    X = np.column_stack([G, Xc]); beta, sse, df = _ols(y, X)
    if df < 2 or np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError('Too few units or a covariate collinear with the groups')
    s2 = sse / df; covb = s2 * np.linalg.inv(X.T @ X)
    # homogeneity of slopes: add group x covariate interactions
    inter = np.column_stack([G[:, i] * Xc[:, j] for j in range(len(cov)) for i in range(1, k)])
    beta_h, sse_h, df_h = _ols(y, np.column_stack([X, inter]))
    num = len(cov) * (k - 1); f_h = ((sse - sse_h) / num) / (sse_h / df_h) if df_h > 0 else np.nan
    p_h = float(stats.f.sf(f_h, num, df_h)) if df_h > 0 else None
    homogeneous = p_h is None or p_h >= q['slope_homogeneity_alpha']
    pairs = _pairs(groups, q); diag, must = [], []
    adjusted = _contrast_rows(beta, covb, df, groups, pairs, lv, X.shape[1])
    # omnibus adjusted group F (equal adjusted means)
    L = np.column_stack([-np.ones(k - 1), np.eye(k - 1), np.zeros((k - 1, len(cov)))]); lb = L @ beta
    f_g = float(lb @ np.linalg.solve(L @ covb @ L.T, lb) / (k - 1))
    tq = stats.t.ppf((1 + lv) / 2, df)
    emm = [{'group': g, 'adjusted_mean': float(beta[i]), 'standard_error': float(np.sqrt(covb[i, i])), 'ci': [float(beta[i] - tq * np.sqrt(covb[i, i])), float(beta[i] + tq * np.sqrt(covb[i, i]))],
            'unadjusted_mean': float(used.outcome[used.group == g].mean()), 'n': int((used.group == g).sum()),
            'covariate_means': {c: float(used[c][used.group == g].mean()) for c in cov}} for i, g in enumerate(groups)]
    slopes = [{'covariate': c, 'slope': float(beta[k + j]), 'standard_error': float(np.sqrt(covb[k + j, k + j])),
               'ci': [float(beta[k + j] - tq * np.sqrt(covb[k + j, k + j])), float(beta[k + j] + tq * np.sqrt(covb[k + j, k + j]))],
               'p_value': float(2 * stats.t.sf(abs(beta[k + j] / np.sqrt(covb[k + j, k + j])), df))} for j, c in enumerate(cov)]
    # unadjusted comparison (one-way ANOVA, pooled variance) for context
    beta_u, sse_u, df_u = _ols(y, G); covu = sse_u / df_u * np.linalg.inv(G.T @ G)
    unadjusted = _contrast_rows(beta_u, covu, df_u, groups, pairs, lv, k)
    # covariate balance and support
    for c in cov:
        groups_c = [used[c][used.group == g] for g in groups]
        p_bal = float(stats.f_oneway(*groups_c).pvalue)
        if p_bal < .05:
            diag.append(f'covariate_imbalance:{c}')
            must.append(f'{c} differs between groups (one-way ANOVA p = {p_bal:.3g}); the adjustment then relies on the common-slope model over a region some groups barely cover.')
        outside = [g for g, s in zip(groups, groups_c) if not s.min() <= xbar[cov.index(c)] <= s.max()]
        if outside:
            diag.append(f'adjustment_point_outside_group_range:{c}')
            must.append(f'The adjustment point {c} = {xbar[cov.index(c)]:.4g} lies outside the observed range of {outside}; their adjusted means are extrapolated.')
    reportable = homogeneous and not any(x.startswith('adjustment_point_outside') for x in diag)
    if not homogeneous:
        diag.append('slopes_differ_between_groups')
        must.append(f'Slopes differ between groups (interaction F p = {p_h:.3g} < {q["slope_homogeneity_alpha"]}); adjusted group differences depend on the covariate value and are withheld.')
    must.insert(0, 'Covariates were declared as measured before treatment; adjusted differences are conditional on them and assume a common linear slope.')
    r2 = 1 - sse / float(((y - y.mean()) ** 2).sum())
    return {'analysis_type': 'ancova', 'primary': {'reportable': bool(reportable), 'adjusted_means': emm, 'adjusted_contrasts': adjusted if reportable else [],
            'adjusted_contrasts_audit': [] if reportable else adjusted,
            'adjusted_group_test': {'f_statistic': f_g, 'df_numerator': k - 1, 'df_denominator': int(df), 'p_value': float(stats.f.sf(f_g, k - 1, df))} if reportable else None,
            'covariate_slopes': slopes, 'slope_homogeneity': {'f_statistic': None if p_h is None else float(f_h), 'df': [num, int(df_h)], 'p_value': p_h, 'alpha': q['slope_homogeneity_alpha']},
            'adjustment_point': dict(zip(cov, map(float, xbar))), 'residual_sd': float(np.sqrt(s2)), 'residual_df': int(df), 'r_squared': float(r2),
            'unadjusted_contrasts': unadjusted, 'post_hoc': q['post_hoc'], 'confidence_level': lv, 'design': cfg['design']},
            'must_mention': must, 'failing_items': [] if reportable else [{'item': 'adjusted_comparison', 'diagnostics': diag}],
            'limitations': ['Linear covariate effects with a common slope; nonlinear or group-specific slopes are not modelled.',
                            'Independent units with equal residual variance; repeated outcomes need MMRM with a baseline covariate.']}


def render(run, result, style):
    """Outcome against the first covariate with the fitted common-slope lines, from saved values."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); p = result['primary']
    d = pd.read_csv(Path(run)/'input.csv', dtype={'group': str}, keep_default_na=False)
    if 'exclude' in d:
        d = d[~d.exclude.astype(str).str.lower().eq('true')]
    c = p['covariate_slopes'][0]['covariate']; slope = p['covariate_slopes'][0]['slope']; x0 = p['adjustment_point'][c]
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.6, 4.), layout='constrained')
        for i, m in enumerate(p['adjusted_means']):
            s = d[d.group == m['group']]; xs = pd.to_numeric(s[c]); col = pstyle.color(token, i)
            ax.plot(xs, pd.to_numeric(s.outcome), 'o', ms=5, color=col, label=m['group'])
            grid = np.linspace(xs.min(), xs.max(), 20)
            ax.plot(grid, m['adjusted_mean'] + slope * (grid - x0), '-', color=col, lw=1.3)  # other covariates held at their means
        ax.axvline(x0, color=pstyle.MUTED, lw=.8, ls=(0, (3, 2)))
        ax.set(xlabel=c, ylabel=f"{p['design']['outcome']} ({p['design']['outcome_unit']})", title='Common-slope ANCOVA (dashed: adjustment point)')
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'ancova', token, 300)
    return ['ancova']
