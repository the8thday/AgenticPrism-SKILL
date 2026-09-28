"""ICH Q1E single-attribute, fixed-batch linear regression.

Slopes first, then intercepts (alpha .25); selected-model pointwise mean
confidence limits, not prediction/tolerance bounds for future batches.
"""
import numpy as np
from scipy import stats, optimize
from .cmc_common import merge, common, text, number, table

DEFAULTS = {'analysis_type': 'stability', 'schema_version': 1, 'input': None, 'source': None,
    'assay': {'attribute': None, 'unit': None, 'time_unit': 'month', 'direction': None,
              'direction_source': None, 'storage_condition': None, 'storage_source': None,
              'storage_class': None, 'linear_model_appropriate': None,
              'independent_errors': None, 'common_variance': None, 'rationale': None},
    'specification': {'lower': None, 'upper': None, 'source': None},
    'extrapolation': {'requested': False, 'accelerated_change': None,
                      'intermediate_no_significant_change': None, 'supporting_data': None,
                      'change_pattern_continues': None, 'commitment_study': None, 'source': None},
    'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    c = merge(raw, DEFAULTS); common(c, 'stability')
    a, s, e = c['assay'], c['specification'], c['extrapolation']
    for k in ('attribute', 'unit', 'direction_source', 'storage_condition', 'storage_source', 'rationale'):
        text(a[k], 'assay.'+k)
    for k in ('linear_model_appropriate', 'independent_errors', 'common_variance'):
        if a[k] is not True:
            raise ValueError('Declare literal true with rationale: assay.'+k)
    if a['time_unit'] != 'month' or a['direction'] not in ('decreasing', 'increasing', 'two_sided'):
        raise ValueError('Use months and declare decreasing, increasing or two_sided direction')
    if a['storage_class'] not in ('room', 'refrigerated', 'frozen'):
        raise ValueError('Declare storage_class room, refrigerated or frozen')
    text(s['source'], 'specification.source')
    needed = {'decreasing': ('lower',), 'increasing': ('upper',), 'two_sided': ('lower', 'upper')}[a['direction']]
    for k in ('lower', 'upper'):
        if k in needed: number(s[k], 'specification.'+k)
        elif s[k] is not None: raise ValueError('Use two_sided to apply both specification limits')
    if len(needed) == 2 and s['lower'] >= s['upper']:
        raise ValueError('Lower specification must be below upper')
    if type(e['requested']) is not bool: raise ValueError('extrapolation.requested must be boolean')
    if e['accelerated_change'] not in (None, 'none_6_months', 'within_3_months', 'between_3_6_months'):
        raise ValueError('Unsupported accelerated_change')
    for k in ('intermediate_no_significant_change', 'supporting_data', 'change_pattern_continues', 'commitment_study'):
        if e[k] is not None and type(e[k]) is not bool: raise ValueError(k+' must be a boolean declaration')
    if e['requested']:
        text(e['source'], 'extrapolation.source')
    return c


def load_data(path, cfg):
    d = table(path, ('observation_id', 'batch', 'time', 'value'), ('time', 'value'))
    if (d.time < 0).any(): raise ValueError('Time must be nonnegative months')
    if d.batch.nunique() < 3: raise ValueError('Q1E analysis requires at least three batches')
    for batch, g in d.groupby('batch'):
        if g.time.nunique() < 3 or g.time.min() != 0 or len(g) < 4:
            raise ValueError(f'Batch {batch} needs baseline zero, three distinct times and at least four observations')
    if 'exclude' in d and (d.exclude != 'false').any():
        raise ValueError('This stability scope does not exclude observations; resolve documented invalid data upstream')
    return d


def ols(y, x):
    y, x = np.asarray(y, float), np.asarray(x, float)
    if np.linalg.matrix_rank(x) < x.shape[1] or len(y) <= x.shape[1]:
        raise ValueError('Regression needs full rank and positive residual degrees of freedom')
    b = np.linalg.lstsq(x, y, rcond=None)[0]; resid = y-x@b
    sse = float(resid@resid); df = len(y)-x.shape[1]
    if sse <= np.finfo(float).eps*max(1., float(np.var(y))):
        raise ValueError('Residual variation is numerically degenerate')
    return {'coefficients': b.tolist(), 'covariance': (sse/df*np.linalg.inv(x.T@x)).tolist(),
            'sse': sse, 'df': df, 'residual_variance': sse/df, 'residuals': resid.tolist()}


def nested_test(small, large):
    df = small['df']-large['df']
    f = max(0., (small['sse']-large['sse'])/df/large['residual_variance'])
    return {'f': f, 'df_numerator': df, 'df_denominator': large['df'],
            'p_value': float(stats.f.sf(f, df, large['df'])), 'alpha': .25}


def regressions(d):
    batches = sorted(d.batch.unique()); t = d.time.to_numpy(float); y = d.value.to_numpy(float)
    z = np.column_stack([(d.batch == b).to_numpy(float) for b in batches])
    full = ols(y, np.c_[z, z*t[:, None]])
    partial = ols(y, np.c_[z, t]); pooled = ols(y, np.c_[np.ones(len(t)), t])
    slopes = nested_test(partial, full); intercepts = None
    if slopes['p_value'] < .25:
        model, selected = 'separate_batches', full
    else:
        intercepts = nested_test(pooled, partial)
        model, selected = ('common_slope', partial) if intercepts['p_value'] < .25 else ('pooled', pooled)
    lines = []
    for i, batch in enumerate(batches):
        if model == 'pooled': idx = [0, 1]
        elif model == 'common_slope': idx = [i, len(batches)]
        else: idx = [i, i+len(batches)]
        lines.append({'batch': batch, 'coefficients': np.asarray(selected['coefficients'])[idx].tolist(),
                      'covariance': np.asarray(selected['covariance'])[np.ix_(idx, idx)].tolist(), 'df': selected['df']})
    return {'model': model, 'slope_test': slopes, 'intercept_test': intercepts,
            'intercept_test_status': 'not_tested_slopes_differ' if intercepts is None else 'tested_after_slopes',
            'selected_fit': selected, 'lines': lines,
            'individual_batch_fits': {b: ols(g.value, np.c_[np.ones(len(g)), g.time]) for b, g in d.groupby('batch')}}


def mean_bound(line, time, side='lower', two_sided=False):
    v = np.array([1., time]); mean = float(v@line['coefficients'])
    se = np.sqrt(max(0., float(v@np.asarray(line['covariance'])@v)))
    q = stats.t.ppf(.975 if two_sided else .95, line['df'])
    return float(mean + (-1 if side == 'lower' else 1)*q*se)


def crossing(line, limit, side, two_sided=False):
    """Earliest nonnegative confidence-bound crossing; None = never crosses.

    Solve the quadratic then reject roots introduced by squaring. No arbitrary
    search horizon is promoted to a shelf-life estimate.
    """
    b = np.asarray(line['coefficients']); c = np.asarray(line['covariance'])
    sign = 1 if side == 'lower' else -1
    a, slope = sign*(b[0]-limit), sign*b[1]
    q = stats.t.ppf(.975 if two_sided else .95, line['df'])
    if a-q*np.sqrt(c[0, 0]) <= 0: return 0.
    poly = [slope*slope-q*q*c[1, 1], 2*a*slope-2*q*q*c[0, 1], a*a-q*q*c[0, 0]]
    roots = np.roots(np.trim_zeros(poly, 'f'))
    candidates = [float(r.real) for r in roots if abs(r.imag) < 1e-7 and r.real >= 0
                  and a+slope*r.real >= -1e-8
                  and abs(mean_bound(line, r.real, side, two_sided)-limit) < 1e-6*max(1., abs(limit))]
    return min(candidates) if candidates else None


def horizon(observed, cfg):
    a, e = cfg['assay'], cfg['extrapolation']; reason = 'No extrapolation requested'
    cap = observed
    supported = all(e[k] is True for k in ('supporting_data', 'change_pattern_continues', 'commitment_study'))
    if e['requested']:
        if not supported: reason = 'Supporting data, continuation of change pattern and commitment study must all be declared'
        elif a['storage_class'] == 'frozen': reason = 'Frozen storage: long-term observations only (Q1E 2.5.2-3)'
        elif e['accelerated_change'] == 'within_3_months': reason = 'Early accelerated change: no extrapolation; assess excursions and shorter shelf life'
        elif e['accelerated_change'] == 'none_6_months':
            cap = min(2*observed, observed+12) if a['storage_class'] == 'room' else min(1.5*observed, observed+6)
            reason = 'Q1E 2.4.1.2 statistical-analysis branch' if a['storage_class'] == 'room' else 'Q1E 2.5.1.1 statistical-analysis branch'
        elif a['storage_class'] == 'room' and e['accelerated_change'] == 'between_3_6_months' and e['intermediate_no_significant_change'] is True:
            cap = min(1.5*observed, observed+6); reason = 'Q1E 2.4.2.1 statistical-analysis branch'
        else: reason = 'Accelerated/intermediate conditions do not support extrapolation'
    return {'observed_common_duration_months': observed, 'allowed_horizon_months': cap,
            'extrapolation_allowed': cap > observed, 'decision': reason, 'declarations': e}


def compute(d, cfg):
    r = regressions(d); a, spec = cfg['assay'], cfg['specification']
    duration = float(d.groupby('batch').time.max().min()); h = horizon(duration, cfg)
    notes = [h['decision'], 'Poolability tested in Q1E order: slopes before intercepts at alpha 0.25.']
    rows = []
    for line in r['lines']:
        times = {side: crossing(line, spec[side], side, a['direction'] == 'two_sided') for side in ('lower','upper') if spec[side] is not None}
        finite = [v for v in times.values() if v is not None]
        raw = min(finite) if finite else None
        rows.append({**line, 'bound_crossings_months': times, 'statistical_crossing_months': raw,
                     'supported_months': min(raw, h['allowed_horizon_months']) if raw is not None else h['allowed_horizon_months']})
    shelf = min(row['supported_months'] for row in rows)
    failures = []
    for side in ('lower', 'upper'):
        if spec[side] is not None:
            mask = d.value < spec[side] if side == 'lower' else d.value > spec[side]
            failures.extend({'observation_id': row.observation_id, 'reason': 'observed_'+side+'_specification_failure'} for row in d[mask].itertuples())
    if failures: notes.append('Observed values outside specification exist; inspect every flagged observation.')
    # Declared linearity is still accompanied by a response-based diagnostic.
    t = d.time.to_numpy(float); z = np.column_stack([(d.batch == b).to_numpy(float) for b in sorted(d.batch.unique())])
    curvature = None
    if len(d) > 3*z.shape[1] and all(g.time.nunique() >= 4 for _,g in d.groupby('batch')):
        linear = ols(d.value, np.c_[z,z*t[:,None]])
        quadratic = ols(d.value, np.c_[z,z*t[:,None],z*t[:,None]**2])
        curvature = nested_test(linear, quadratic); curvature['alpha'] = .05
    reportable = not (curvature and curvature['p_value'] < .05)
    if not reportable:
        notes.append('Batch-specific curvature diagnostic failed at 0.05; linear shelf life withheld.')
        failures.append({'reason': 'linear_model_curvature'})
    if shelf == 0: failures.append({'reason': 'confidence_bound_fails_at_baseline'})
    extrapolated = shelf > duration
    notes.append('Shelf life is extrapolated within declared Q1E limits; verify with commitment data.' if extrapolated else 'Supported shelf life is not extrapolated beyond the common observed duration.')
    return {'analysis_type': 'stability', 'schema_version': 1, 'poolability': r,
        'batch_estimates': rows, 'specification': spec, 'assay': a, 'decision_tree': h, 'curvature_diagnostic': curvature,
        'primary': {'status': 'estimated' if reportable else 'withheld', 'shelf_life_months': shelf if reportable else None,
                    'audit_supported_months': shelf, 'extrapolated': extrapolated if reportable else False,
                    'model': r['model'], 'confidence_level': .95, 'bound': a['direction'],
                    'governing_batches': [row['batch'] for row in rows if row['supported_months'] == shelf]},
        'must_mention': notes, 'failing_items': failures,
        'limitations': ['Unweighted linear Gaussian fixed-batch regression with common residual variance; one attribute and storage condition per analysis.',
            'Confidence bounds concern the mean of studied batches, not individual units or a population tolerance bound for future batches.',
            'Intervals condition on the selected poolability model; model selection can reduce coverage.',
            'The shortest supported attribute and batch governs the overall proposal; regulatory approval is outside this computation.']}
