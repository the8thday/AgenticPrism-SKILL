"""Inter-rater agreement for independent rated units (0.13.6).

Cohen's kappa (unweighted, or weighted with declared linear/quadratic or
sourced agreement weights) for two raters and nominal/ordinal categories, with
Fleiss, Cohen & Everitt (1969) large-sample variances. Intraclass correlation in
the McGraw & Wong (1996) forms (one-way random; two-way consistency or absolute
agreement; single or average of k ratings) with exact F or approximate
(agreement) confidence limits; average-measure limits are the Spearman-Brown
transform of the single-measure limits. Complete balanced designs only.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'rater_agreement'}
METHODS = ('cohen_kappa', 'icc')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'rater_agreement', 'input': None, 'source': None, 'method': None,
            'design': {'unit': None, 'rationale': None, 'raters': None, 'rater_rationale': None, 'outcome': None, 'outcome_unit': None},
            'kappa': {'categories': None, 'ordered': None, 'weights': None, 'weights_source': None},
            'icc': {'model': None, 'type': None, 'unit': None, 'rationale': None},
            'confidence_level': .95, 'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported agreement configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'rater_agreement' or c['method'] not in METHODS:
        raise ValueError('Declare method cohen_kappa or icc')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    lv = c['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    if c['method'] == 'cohen_kappa':
        q = c['kappa']
        if any(v is not None for v in c['icc'].values()): raise ValueError('icc settings do not apply to cohen_kappa')
        cats = q['categories']
        if not isinstance(cats, list) or len(cats) < 2 or not all(_text(x) for x in cats) or len(set(cats)) != len(cats):
            raise ValueError('Declare at least two distinct category labels, in order if ordinal')
        if not isinstance(q['ordered'], bool): raise ValueError('Declare whether categories are ordered')
        w = q['weights']
        if w in (None, 'none'): q['weights'] = 'none'
        elif w in ('linear', 'quadratic'):
            if not q['ordered']: raise ValueError('Linear/quadratic weights need ordered categories')
        elif isinstance(w, list):
            a = np.asarray(w, float)
            if a.shape != (len(cats), len(cats)) or not np.isfinite(a).all() or not np.allclose(np.diag(a), 1) or (a < 0).any() or (a > 1).any():
                raise ValueError('Custom agreement weights: a square matrix in category order, diagonal 1, entries in [0, 1]')
        else: raise ValueError('weights must be none, linear, quadratic or a sourced agreement-weight matrix')
        if q['weights'] != 'none' and not _text(q['weights_source']): raise ValueError('Declare weights_source for weighted kappa')
        if q['weights'] == 'none' and q['weights_source'] is not None: raise ValueError('weights_source applies only to weighted kappa')
    else:
        q = c['icc']
        if any(v is not None for v in c['kappa'].values()): raise ValueError('kappa settings do not apply to icc')
        if q['model'] not in ('one_way_random', 'two_way_random', 'two_way_mixed'): raise ValueError('Declare icc.model')
        if q['model'] == 'one_way_random':
            if q['type'] not in (None, 'agreement'): raise ValueError('One-way ICC measures absolute agreement only')
            q['type'] = 'agreement'
        elif q['type'] not in ('agreement', 'consistency'): raise ValueError('Declare icc.type agreement or consistency')
        if q['unit'] not in ('single', 'average'): raise ValueError('Declare icc.unit single or average')
        if not _text(q['rationale']): raise ValueError('Declare icc.rationale (why this form answers the question)')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'unit_id': str, 'rater': str, 'rating': str})
    need = {'unit_id', 'rater', 'rating'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if d[list(need)].isna().any().any(): raise ValueError('Missing unit, rater or rating; no imputation or silent deletion')
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded row needs an exclusion_reason')
        if d.groupby('unit_id').exclude.nunique().max() > 1: raise ValueError('Exclude whole units, not single ratings')
    else: d['exclude'] = False
    u = d[~d.exclude]
    if u.duplicated(['unit_id', 'rater']).any(): raise ValueError('Each rater rates each unit once; average technical repeats first by a declared rule')
    counts = u.groupby('unit_id').size()
    if counts.nunique() != 1 or counts.iloc[0] < 2: raise ValueError('Complete balanced design required: every unit has the same number (>= 2) of ratings')
    two_way = c['method'] == 'cohen_kappa' or c['icc']['model'] != 'one_way_random'
    if two_way:
        raters = sorted(u.rater.unique())
        if any(sorted(g.rater) != raters for _, g in u.groupby('unit_id')):
            raise ValueError('Two-way designs need the same raters on every unit; use the one-way model only if raters differ by design')
    if c['method'] == 'cohen_kappa':
        if u.rater.nunique() != 2: raise ValueError("Cohen's kappa needs exactly two raters")
        bad = set(u.rating) - set(c['kappa']['categories'])
        if bad: raise ValueError('Ratings outside the declared categories: ' + ', '.join(sorted(bad)))
    else:
        v = pd.to_numeric(u.rating, errors='coerce')
        if v.isna().any() or not np.isfinite(v).all(): raise ValueError('ICC needs numeric ratings')
    if u.unit_id.nunique() < 3: raise ValueError('At least three rated units are required')
    return d


def kappa_weights(kind, k, matrix=None):
    if kind == 'none': return np.eye(k)
    if kind == 'linear': return 1 - np.abs(np.subtract.outer(np.arange(k), np.arange(k)))/(k-1)
    if kind == 'quadratic': return 1 - (np.subtract.outer(np.arange(k), np.arange(k))/(k-1))**2
    return np.asarray(matrix, float)


def kappa_from_table(table, w):
    """Kappa and Fleiss-Cohen-Everitt (1969) non-null and null variances from a count table."""
    table = np.asarray(table, float); n = table.sum(); p = table/n; r, s = p.sum(1), p.sum(0)
    po = float((w*p).sum()); pe = float((w*np.outer(r, s)).sum())
    if pe >= 1 - 1e-12: return {'n': int(n), 'observed_agreement': po, 'expected_agreement': pe, 'kappa': None}
    kap = (po-pe)/(1-pe); wr, wc = w @ s, w.T @ r
    var = ((p*(w - np.add.outer(wr, wc)*(1-kap))**2).sum() - (kap - pe*(1-kap))**2)/(n*(1-pe)**2)
    var0 = ((np.outer(r, s)*(w - np.add.outer(wr, wc))**2).sum() - pe**2)/(n*(1-pe)**2)
    return {'n': int(n), 'observed_agreement': po, 'expected_agreement': pe, 'kappa': float(kap),
            'variance': float(var), 'null_variance': float(var0)}


def kappa(u, c):
    q = c['kappa']; cats = q['categories']; k = len(cats); level = c['confidence_level']
    raters = sorted(u.rater.unique())
    wide = u.pivot(index='unit_id', columns='rater', values='rating')
    table = pd.crosstab(pd.Categorical(wide[raters[0]], cats), pd.Categorical(wide[raters[1]], cats), dropna=False).to_numpy()
    w = kappa_weights(q['weights'], k, q['weights'] if isinstance(q['weights'], list) else None)
    r = kappa_from_table(table, w)
    out = {'method': 'cohen_kappa', 'weights': q['weights'] if not isinstance(q['weights'], list) else 'custom',
           'weight_matrix': w.tolist(), 'weights_source': q['weights_source'], 'categories': cats,
           'rater_rows': raters[0], 'rater_columns': raters[1], 'table': table.tolist(), 'n_units': r['n'],
           'observed_agreement': r['observed_agreement'], 'expected_agreement': r['expected_agreement'],
           'confidence_level': level, 'interval_method': 'Fleiss-Cohen-Everitt (1969) large-sample Wald; null-variance z test',
           'diagnostics': []}
    rs, cs = table.sum(1)/r['n'], table.sum(0)/r['n']
    if q['weights'] == 'none':
        pmax = float(np.minimum(rs, cs).sum())
        out['maximum_kappa_given_marginals'] = (pmax - r['expected_agreement'])/(1 - r['expected_agreement']) if r['kappa'] is not None else None
    if k == 2 and q['weights'] == 'none':
        a, b, cc, dd = table.ravel()/r['n']
        out.update(prevalence_index=float(a-dd), bias_index=float(b-cc), pabak=float(2*r['observed_agreement']-1))
        if abs(a-dd) > .5: out['diagnostics'].append('high_prevalence_index')
    if (table.sum(1) == 0).any() or (table.sum(0) == 0).any(): out['diagnostics'].append('declared_category_unused_by_a_rater')
    if r['n'] < 30: out['diagnostics'].append('small_sample_wald_interval')
    if r['kappa'] is None:
        out.update(status='withheld', reportable=False, reason='expected_agreement_one_kappa_undefined'); return out
    se = float(np.sqrt(max(r['variance'], 0.))); se0 = float(np.sqrt(max(r['null_variance'], 0.)))
    z = float(stats.norm.ppf((1+level)/2)); est = r['kappa']
    lo, hi = est - z*se, est + z*se
    if lo < -1 or hi > 1: out['diagnostics'].append('interval_truncated_to_kappa_range')
    out.update(status='estimated' if se > 0 else 'limited', reportable=se > 0, estimate=est, standard_error=se,
               null_standard_error=se0, ci=[float(max(lo, -1.)), float(min(hi, 1.))] if se > 0 else None)
    if se0 > 0:
        stat = est/se0; out.update(z_statistic=float(stat), p_two_sided_vs_zero=float(2*stats.norm.sf(abs(stat))))
    if se <= 0: out['reason'] = 'zero_large_sample_variance'
    return out


def mean_squares(y):
    """Balanced n x k table: one-way and two-way ANOVA mean squares."""
    n, k = y.shape; g = y.mean(); rm, cm = y.mean(1), y.mean(0)
    ssr = k*((rm-g)**2).sum(); ssc = n*((cm-g)**2).sum(); sst = ((y-g)**2).sum()
    sse = sst - ssr - ssc; ssw = sst - ssr
    return {'n': n, 'k': k, 'MSR': ssr/(n-1), 'MSC': ssc/(k-1), 'MSE': sse/((n-1)*(k-1)), 'MSW': ssw/(n*(k-1)),
            'df_rows': n-1, 'df_columns': k-1, 'df_error': (n-1)*(k-1), 'df_within': n*(k-1)}


def icc_from_table(y, model, kind, unit, level):
    m = mean_squares(np.asarray(y, float)); n, k = m['n'], m['k']; a = (1-level)/2
    MSR, MSC, MSE, MSW = m['MSR'], m['MSC'], m['MSE'], m['MSW']
    if model == 'one_way_random':
        resid, df2 = MSW, m['df_within']
    else: resid, df2 = MSE, m['df_error']
    if resid <= 0 or MSR <= 0:
        return {'mean_squares': m, 'estimate': None, 'reason': 'zero_residual_or_between_unit_mean_square'}
    F0 = MSR/resid; p = float(stats.f.sf(F0, m['df_rows'], df2))
    if model == 'one_way_random' or kind == 'consistency':
        est = (MSR-resid)/(MSR+(k-1)*resid)
        FL = F0/stats.f.ppf(1-a, m['df_rows'], df2); FU = F0*stats.f.ppf(1-a, df2, m['df_rows'])
        lo, hi = (FL-1)/(FL+k-1), (FU-1)/(FU+k-1); approx = False; v = None
    else:
        est = (MSR-MSE)/(MSR+(k-1)*MSE+k*(MSC-MSE)/n)
        aa = k*est/(n*(1-est)); bb = 1 + k*est*(n-1)/(n*(1-est))
        v = (aa*MSC+bb*MSE)**2/((aa*MSC)**2/(k-1)+(bb*MSE)**2/((n-1)*(k-1)))
        Fs = stats.f.ppf(1-a, n-1, v); Fu = stats.f.ppf(1-a, v, n-1)
        lo = n*(MSR-Fs*MSE)/(Fs*(k*MSC+(k*n-k-n)*MSE)+n*MSR)
        hi = n*(Fu*MSR-MSE)/(k*MSC+(k*n-k-n)*MSE+n*Fu*MSR); approx = True
    if unit == 'average':
        sb = lambda r: k*r/(1+(k-1)*r)
        est, lo, hi = sb(est), sb(lo), sb(hi)
    return {'mean_squares': m, 'estimate': float(est), 'ci': [float(lo), float(hi)], 'F': float(F0),
            'F_df': [m['df_rows'], df2], 'p_vs_zero': p, 'approximate_df': None if v is None else float(v), 'approximate_interval': approx}


def icc(u, c):
    q = c['icc']; level = c['confidence_level']
    if q['model'] == 'one_way_random':
        u = u.assign(slot=u.groupby('unit_id').cumcount())
        y = u.pivot(index='unit_id', columns='slot', values='rating').astype(float)
    else: y = u.pivot(index='unit_id', columns='rater', values='rating').astype(float)
    r = icc_from_table(y.to_numpy(), q['model'], q['type'], q['unit'], level)
    label = {'one_way_random': 'ICC(1', 'two_way_random': 'ICC(A' if q['type'] == 'agreement' else 'ICC(C',
             'two_way_mixed': 'ICC(A' if q['type'] == 'agreement' else 'ICC(C'}[q['model']] + (',1)' if q['unit'] == 'single' else f',k={y.shape[1]})')
    out = {'method': 'icc', 'form': label, 'model': q['model'], 'type': q['type'], 'unit': q['unit'], 'rationale': q['rationale'],
           'n_units': int(y.shape[0]), 'k_ratings': int(y.shape[1]), 'confidence_level': level, 'mean_squares': r['mean_squares'],
           'interval_method': ('McGraw-Wong approximate F (absolute agreement)' if r.get('approximate_interval') else 'Exact F') +
                              (', Spearman-Brown transform for the average of k ratings' if q['unit'] == 'average' else ''),
           'rater_means': {str(k): float(v) for k, v in y.mean(axis=0).items()} if q['model'] != 'one_way_random' else None,
           'diagnostics': []}
    if r['estimate'] is None:
        out.update(status='withheld', reportable=False, reason=r['reason']); return out
    if y.shape[0] < 10: out['diagnostics'].append('fewer_than_10_units_wide_interval')
    if r['estimate'] < 0: out['diagnostics'].append('negative_icc_estimate')
    if q['model'] != 'one_way_random' and r['mean_squares']['MSC'] > 4*r['mean_squares']['MSE'] and q['type'] == 'consistency':
        out['diagnostics'].append('systematic_rater_differences_ignored_by_consistency')
    out.update(status='estimated', reportable=True, estimate=r['estimate'], ci=r['ci'], F=r['F'], F_df=r['F_df'],
               p_vs_zero=r['p_vs_zero'], approximate_df=r['approximate_df'])
    return out


def compute(d, c):
    u = d[~d.exclude].copy(); res = kappa(u, c) if c['method'] == 'cohen_kappa' else icc(u, c)
    res.update(unit=c['design']['unit'], outcome=c['design']['outcome'], outcome_unit=c['design']['outcome_unit'])
    notes = ['Independent rated unit: ' + c['design']['unit'] + '. Raters: ' + c['design']['raters'] + ' (' + c['design']['rater_rationale'] + ').',
             'Agreement statistics depend on the spread of the rated units and the category prevalence in this sample; they do not transfer to other populations or rater pools without assumption.',
             'Reliability is not accuracy: agreement between raters does not show that either rater is correct.']
    if c['method'] == 'cohen_kappa':
        notes.append('Kappa was computed with ' + ('no weights' if res['weights'] == 'none' else res['weights'] + ' agreement weights (' + str(c['kappa']['weights_source']) + ')') +
                     '. Category labels and weights were declared before analysis; they are not chosen to maximize kappa.')
        notes.append('Benchmark labels such as "substantial" are conventions without inferential meaning; report kappa with its interval, the observed agreement and the table.')
        if 'high_prevalence_index' in res['diagnostics']: notes.append('Skewed category prevalence: kappa can be low despite high observed agreement (prevalence index, PABAK and maximum kappa are descriptive context, not substitutes).')
        if 'small_sample_wald_interval' in res['diagnostics']: notes.append('Fewer than 30 units: the large-sample kappa interval is approximate.')
        if 'interval_truncated_to_kappa_range' in res['diagnostics']: notes.append('The Wald interval exceeded the kappa range and was truncated; it is not reliable near the boundary.')
    else:
        notes.append(f"ICC form {res['form']} ({res['model']}, {res['type']}, {res['unit']}) was declared before analysis: {res['rationale']}. Forms answer different questions and are not chosen by value.")
        if res['type'] == 'consistency': notes.append('Consistency ICC ignores systematic differences between raters; it does not measure absolute agreement.')
        if res['unit'] == 'average': notes.append('An average-measure ICC describes the mean of k ratings; it does not apply to a single rater.')
        for flag, text in (('fewer_than_10_units_wide_interval', 'Fewer than 10 rated units: the ICC interval is wide and the estimate unstable.'),
                           ('negative_icc_estimate', 'Negative ICC estimate: within-unit variation exceeds between-unit variation; reported as observed, not truncated to zero.'),
                           ('systematic_rater_differences_ignored_by_consistency', 'Rater mean squares are large relative to error; consistency ICC ignores these systematic rater differences.')):
            if flag in res['diagnostics']: notes.append(text)
    if not res['reportable']: notes.append('Agreement inference withheld: ' + res['reason'])
    return {'analysis_type': 'rater_agreement', 'primary': res, 'must_mention': notes,
            'preprocessing': {'excluded_units': sorted(d[d.exclude].unit_id.unique().tolist())},
            'failing_items': [] if res['reportable'] else [{'reason': res['reason']}],
            'limitations': ['Complete balanced designs only: no missing ratings, unbalanced or multi-rater kappa (Fleiss/Light), Krippendorff alpha, or mixed-model ICC with missing data.',
                            'Agreement on a continuous scale between two methods is method comparison (Bland-Altman), not ICC.',
                            'No Prism, psych or software-equivalence claim beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    d = pd.read_csv(Path(run)/'input.csv', dtype={'unit_id': str, 'rater': str, 'rating': str})
    if 'exclude' in d: d = d[~d.exclude.astype(str).str.lower().eq('true')]
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        if p['method'] == 'cohen_kappa':
            t = np.asarray(p['table']); fig, ax = plt.subplots(figsize=(1.6+.7*len(t), 1.4+.7*len(t)), layout='constrained')
            ax.imshow(t, cmap='Greys'); cats = p['categories']
            for i in range(len(t)):
                for j in range(len(t)): ax.text(j, i, int(t[i, j]), ha='center', va='center', color='white' if t[i, j] > t.max()/2 else 'black')
            ax.set(xticks=range(len(cats)), yticks=range(len(cats)), xticklabels=cats, yticklabels=cats,
                   xlabel=p['rater_columns'], ylabel=p['rater_rows'], title='Rating table')
            pstyle.save(fig, out/'agreement_table', token, 300); return ['agreement_table']
        fig, ax = plt.subplots(figsize=(5.5, 3.6), layout='constrained')
        units = list(dict.fromkeys(d.unit_id)); pos = {x: i for i, x in enumerate(units)}
        for j, (rater, g) in enumerate(d.groupby('rater')):
            ax.scatter([pos[x] for x in g.unit_id], g.rating.astype(float), s=16, color=pstyle.color(token, j % 8), label=rater if p['model'] != 'one_way_random' else None)
        ax.set(xticks=range(len(units)), xticklabels=units if len(units) <= 30 else [], xlabel='Rated unit', ylabel=p['outcome'], title=p['form'])
        if p['model'] != 'one_way_random' and d.rater.nunique() <= 8: ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'agreement_ratings', token, 300); return ['agreement_ratings']
