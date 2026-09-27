"""Bioanalytical method validation summaries (ICH M10-style experiments).

Inputs are back-calculated concentrations; calibration-curve fitting is upstream
(elisa-quantification or the lab's system). Acceptance criteria are always
declared by the user with a source; nothing is defaulted from a guideline name.

Accuracy and precision per QC level use the one-way random-run ANOVA
(classical unbalanced n0 estimator, as VCA::anovaVCA): within-run CV from the
within-run mean square, between-run (intermediate) CV from the between-run
component plus within-run variance, total error = |bias%| + between-run CV%.
Supplements that go beyond the guideline's point-estimate rules are labelled:
t interval for the mean bias, MLS interval for between-run precision and the
beta-expectation tolerance interval of Mee (1984) used in accuracy profiles.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats
from .variance_intervals import mls_limits

EXPERIMENTS = ('accuracy_precision', 'dilution_linearity', 'parallelism', 'selectivity', 'specificity', 'stability')
CRITERIA = ('accuracy_percent', 'accuracy_percent_edge', 'precision_cv_percent', 'precision_cv_percent_edge',
            'total_error_percent', 'total_error_percent_edge', 'required_pass_fraction', 'blank_pass_fraction',
            'parallelism_cv_percent')
NEEDED = {'accuracy_precision': ('accuracy_percent', 'accuracy_percent_edge', 'precision_cv_percent', 'precision_cv_percent_edge',
                                 'total_error_percent', 'total_error_percent_edge'),
          'dilution_linearity': ('accuracy_percent', 'precision_cv_percent'),
          'parallelism': ('parallelism_cv_percent',),
          'selectivity': ('accuracy_percent', 'accuracy_percent_edge', 'required_pass_fraction', 'blank_pass_fraction'),
          'specificity': ('accuracy_percent_edge', 'required_pass_fraction', 'blank_pass_fraction'),
          'stability': ('accuracy_percent',)}
COLUMNS = {'accuracy_precision': ('level', 'role', 'run_id', 'nominal'),
           'dilution_linearity': ('series_id', 'dilution_factor', 'nominal'),
           'parallelism': ('sample_id', 'dilution_factor'),
           'selectivity': ('source_id', 'role', 'nominal'),
           'specificity': ('source_id', 'role', 'nominal', 'interferent'),
           'stability': ('condition', 'level', 'nominal')}
DEFAULTS = {'analysis_type': 'method_validation', 'schema_version': 1, 'input': None, 'source': None, 'experiment': None,
    'assay': {'platform': None, 'matrix': None, 'concentration_unit': None, 'lloq': None, 'uloq': None,
              'concentrations_back_calculated': None, 'independent_runs': None, 'rationale': None},
    'criteria': {'source': None, **{k: None for k in CRITERIA}},
    'statistics': {'confidence_level': .90, 'cv_denominator': None, 'tolerance_beta': None, 'profile_limit_percent': None,
                   'slope_margin': None},
    'report': {'plot_style': 'prism_like'}}


def _positive(v, name, upper=None):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) or v <= 0 or (upper is not None and v > upper):
        raise ValueError(f'{name} must be a finite positive number' + (f' not above {upper}' if upper else ''))


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported method-validation configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'method_validation' or type(c['schema_version']) is not int or c['schema_version'] != 1:
        raise ValueError('Unsupported analysis/schema')
    e = c['experiment']
    if e not in EXPERIMENTS:
        raise ValueError('experiment must be one of ' + ', '.join(EXPERIMENTS))
    a, q, s = c['assay'], c['criteria'], c['statistics']
    for obj, keys in ((c, ('input', 'source')), (a, ('platform', 'matrix', 'concentration_unit', 'rationale')), (q, ('source',))):
        if any(not isinstance(obj[k], str) or not obj[k].strip() for k in keys):
            raise ValueError('input, source, platform, matrix, concentration_unit, rationale and criteria.source are required strings')
    if a['concentrations_back_calculated'] is not True or a['independent_runs'] is not True:
        raise ValueError('Declare literal true: concentrations_back_calculated and independent_runs (with evidence in rationale)')
    for k in ('lloq', 'uloq'):
        if a[k] is not None:
            _positive(a[k], 'assay.' + k)
    if a['lloq'] is not None and a['uloq'] is not None and a['lloq'] >= a['uloq']:
        raise ValueError('lloq must be below uloq')
    if e in ('dilution_linearity', 'selectivity', 'specificity') and a['lloq'] is None:
        raise ValueError(f'{e} requires assay.lloq')
    if e == 'dilution_linearity' and a['uloq'] is None:
        raise ValueError('dilution_linearity requires assay.uloq for the hook-effect check')
    needed = NEEDED[e]
    for k in CRITERIA:
        if k in needed:
            if k.endswith('fraction'):
                _positive(q[k], 'criteria.' + k, 1)
            else:
                _positive(q[k], 'criteria.' + k, 100)
        elif q[k] is not None:
            raise ValueError(f'criteria.{k} does not apply to {e}')
    level = s['confidence_level']
    if isinstance(level, bool) or not isinstance(level, (int, float)) or not .5 < level < 1:
        raise ValueError('confidence_level must lie between .5 and 1')
    if e == 'accuracy_precision':
        if s['cv_denominator'] not in ('observed_mean', 'nominal'):
            raise ValueError('Declare statistics.cv_denominator: observed_mean or nominal')
        if (s['tolerance_beta'] is None) != (s['profile_limit_percent'] is None):
            raise ValueError('tolerance_beta and profile_limit_percent are declared together')
        if s['tolerance_beta'] is not None:
            _positive(s['tolerance_beta'], 'tolerance_beta', .999)
            _positive(s['profile_limit_percent'], 'profile_limit_percent', 100)
    elif s['cv_denominator'] is not None or s['tolerance_beta'] is not None or s['profile_limit_percent'] is not None:
        raise ValueError('cv_denominator/tolerance settings apply only to accuracy_precision')
    if s['slope_margin'] is not None:
        if e != 'parallelism':
            raise ValueError('slope_margin applies only to parallelism')
        _positive(s['slope_margin'], 'slope_margin', 1)
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Unsupported plot style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    e = cfg['experiment']
    required = {'observation_id', 'value', 'status', 'exclude', 'exclusion_reason', *COLUMNS[e]}
    if not required <= set(d.columns) or not len(d):
        raise ValueError('Missing columns: ' + ', '.join(sorted(required - set(d.columns))) if not required <= set(d.columns) else 'Empty data')
    if d.observation_id.duplicated().any():
        raise ValueError('Duplicate observation_id')
    for k in required - {'value', 'exclusion_reason', 'nominal'}:
        if (d[k].str.strip() == '').any():
            raise ValueError(f'Empty {k}')
    if set(d.exclude) - {'true', 'false'}:
        raise ValueError('exclude must be literal true/false')
    d['exclude'] = d.exclude == 'true'
    if (d.exclude & (d.exclusion_reason.str.strip() == '')).any():
        raise ValueError('Excluded rows require a documented reason (M10: only obvious, documented errors)')
    if set(d.status) - {'quantified', 'below_lloq', 'above_uloq'}:
        raise ValueError('status must be quantified, below_lloq or above_uloq')
    q = d.status == 'quantified'
    if (d.value[q].str.strip() == '').any() or (d.value[~q].str.strip() != '').any():
        raise ValueError('value is required for quantified rows and must be empty otherwise; do not substitute LLOQ/ULOQ numbers')
    d['value'] = pd.to_numeric(d.value.where(q, None), errors='raise')
    if not np.isfinite(d.value[q]).all() or (d.value[q] < 0).any():
        raise ValueError('Quantified values must be finite and nonnegative')
    if 'nominal' in d:
        blank = d.nominal.str.strip() == ''
        allowed = e in ('selectivity', 'specificity') and (d.role[blank] == 'blank').all()
        if blank.any() and not allowed:
            raise ValueError('nominal may be empty only for selectivity/specificity blank rows')
        d['nominal'] = pd.to_numeric(d.nominal.where(~blank, None), errors='raise')
        if (d.nominal[~blank] <= 0).any():
            raise ValueError('Nominal concentrations must be positive')
    if 'dilution_factor' in d:
        d['dilution_factor'] = pd.to_numeric(d.dilution_factor, errors='raise')
        if not np.isfinite(d.dilution_factor).all() or (d.dilution_factor < 1).any():
            raise ValueError('dilution_factor must be finite and >= 1')
    if e == 'accuracy_precision':
        if set(d.role) - {'lloq', 'qc', 'uloq'}:
            raise ValueError('role must be lloq, qc or uloq')
        if (d.groupby('level').nominal.nunique() != 1).any() or (d.groupby('level').role.nunique() != 1).any():
            raise ValueError('Each level needs one nominal concentration and one role')
    if e in ('selectivity', 'specificity'):
        roles = {'blank', 'lloq', 'high'} if e == 'selectivity' else {'blank', 'lloq', 'uloq'}
        if set(d.role) - roles:
            raise ValueError(f'{e} roles must be ' + ', '.join(sorted(roles)))
    if e == 'stability' and (d.groupby(['condition', 'level']).nominal.nunique() != 1).any():
        raise ValueError('Each condition/level needs one nominal concentration')
    return d


def one_way(groups):
    """Classical one-way random ANOVA; returns components and mean squares (n0 for unbalanced)."""
    y = [np.asarray(g, float) for g in groups]
    n = np.array([len(g) for g in y]); p = len(y); N = int(n.sum())
    if p < 2 or N - p < 1:
        raise ValueError('At least two runs and one within-run degree of freedom required')
    means = np.array([g.mean() for g in y]); grand = float(np.concatenate(y).mean())
    msb = float(np.sum(n*(means-grand)**2)/(p-1))
    msw = float(sum(((g-g.mean())**2).sum() for g in y)/(N-p))
    n0 = float((N - (n**2).sum()/N)/(p-1))
    between = max(0., (msb-msw)/n0)
    return {'runs': p, 'n': N, 'n0': n0, 'balanced': bool(np.all(n == n[0])), 'grand_mean': grand, 'run_means': means.tolist(),
            'ms_between': msb, 'ms_within': msw, 'df_between': p-1, 'df_within': N-p,
            'var_between_raw': (msb-msw)/n0, 'var_between': between, 'var_within': msw, 'var_intermediate': between+msw,
            'mean_variance': float(((n**2).sum()/N**2)*between + msw/N)}


def mee_tolerance(ow, beta):
    """Mee (1984) beta-expectation interval for a new result in a new run (Satterthwaite df)."""
    p, n0, sb, sw = ow['runs'], ow['n0'], ow['var_between'], ow['var_within']
    r = sb/sw if sw > 0 else np.inf
    b2 = (r+1)/(n0*r+1)
    nu = (r+1)**2/((r+1/n0)**2/(p-1) + (1-1/n0)/(p*n0))
    half = float(stats.t.ppf((1+beta)/2, nu)*np.sqrt(sb+sw)*np.sqrt(1+1/(p*n0*b2)))
    return {'beta': beta, 'df': float(nu), 'variance_ratio': float(r), 'lower': ow['grand_mean']-half, 'upper': ow['grand_mean']+half,
            'method': 'Mee (1984) beta-expectation tolerance interval, one-way random runs, Satterthwaite df; n0 for unbalanced runs'}


def _limit(q, key, edge):
    return q[key+'_edge'] if edge and q.get(key+'_edge') is not None else q[key]


def accuracy_precision(d, cfg):
    q, s = cfg['criteria'], cfg['statistics']; level = s['confidence_level']
    rows, runs, diag = [], [], []
    used = d[~d.exclude]
    if (used.status != 'quantified').any():
        raise ValueError('Accuracy/precision QCs must all be quantified; a QC outside the range is a failed result that needs investigation, not a missing value')
    for lv, g in used.groupby('level', sort=False):
        nominal, role = float(g.nominal.iloc[0]), g.role.iloc[0]
        edge = role in ('lloq', 'uloq')
        by_run = [r.value.to_numpy() for _, r in g.groupby('run_id', sort=True)]
        ow = one_way(by_run)
        denom = ow['grand_mean'] if s['cv_denominator'] == 'observed_mean' else nominal
        bias = 100*(ow['grand_mean']/nominal-1)
        within_cv = 100*np.sqrt(ow['var_within'])/denom
        between_cv = 100*np.sqrt(ow['var_intermediate'])/denom
        te = abs(bias)+between_cv
        t = stats.t.ppf((1+level)/2, ow['df_between'])
        se = np.sqrt(ow['mean_variance'])
        # MLS interval for intermediate variance = MS_B/n0 + (1-1/n0) MS_W (positive coefficients).
        mls = mls_limits([ow['ms_between'], ow['ms_within']], [ow['df_between'], ow['df_within']], [1/ow['n0'], 1-1/ow['n0']], level)
        lim_acc, lim_cv, lim_te = (_limit(q, k, edge) for k in ('accuracy_percent', 'precision_cv_percent', 'total_error_percent'))
        checks = {'between_run_accuracy': bool(abs(bias) <= lim_acc), 'within_run_precision': bool(within_cv <= lim_cv),
                  'between_run_precision': bool(between_cv <= lim_cv), 'total_error': bool(te <= lim_te)}
        row = {'level': lv, 'role': role, 'nominal': nominal, 'runs': ow['runs'], 'n': ow['n'], 'balanced': ow['balanced'],
               'mean': ow['grand_mean'], 'bias_percent': float(bias),
               'bias_interval_percent': [float(100*((ow['grand_mean']-t*se)/nominal-1)), float(100*((ow['grand_mean']+t*se)/nominal-1))],
               'within_run_cv_percent': float(within_cv), 'between_run_cv_percent': float(between_cv),
               'between_run_cv_interval_percent': (100*np.sqrt(np.array(mls['variance_interval']))/denom).tolist(),
               'total_error_percent': float(te), 'limits': {'accuracy': lim_acc, 'precision_cv': lim_cv, 'total_error': lim_te},
               'checks': checks, 'passes_criteria': all(checks.values()), 'anova': ow, 'cv_denominator': denom}
        if ow['var_between_raw'] < 0:
            row['note'] = 'Negative raw between-run component set to zero; between-run CV equals within-run CV.'
        if s['tolerance_beta'] is not None:
            ti = mee_tolerance(ow, s['tolerance_beta'])
            ti['lower_percent'], ti['upper_percent'] = 100*(ti['lower']/nominal-1), 100*(ti['upper']/nominal-1)
            ti['within_profile_limits'] = bool(max(abs(ti['lower_percent']), abs(ti['upper_percent'])) <= s['profile_limit_percent'])
            row['tolerance_interval'] = ti
        rows.append(row)
        for run, r in g.groupby('run_id', sort=True):
            v = r.value.to_numpy()
            runs.append({'level': lv, 'run_id': run, 'n': len(v), 'mean': float(v.mean()), 'bias_percent': float(100*(v.mean()/nominal-1)),
                         'cv_percent': float(100*v.std(ddof=1)/v.mean()) if len(v) > 1 and v.mean() > 0 else None,
                         'within_run_accuracy_pass': bool(abs(100*(v.mean()/nominal-1)) <= lim_acc),
                         'within_run_precision_pass': bool(len(v) > 1 and 100*v.std(ddof=1)/v.mean() <= lim_cv)})
        if ow['runs'] < 6:
            diag.append(f'{lv}: {ow["runs"]} runs; ICH M10 asks for at least 6 runs over 2 or more days.')
        if min(len(v) for v in by_run) < 3:
            diag.append(f'{lv}: fewer than 3 replicates in some run; ICH M10 asks for at least 3 per run and level.')
    rows.sort(key=lambda r: r['nominal'])
    roles = [r['role'] for r in rows]
    if roles.count('lloq') != 1 or roles.count('uloq') != 1 or len(rows) < 5:
        diag.append('ICH M10 LBA A&P uses 5 levels: LLOQ, low, medium, high and ULOQ; this design differs.')
    passing = [r['nominal'] for r in rows if r['passes_criteria']]
    ends = {r['role']: r for r in rows if r['role'] in ('lloq', 'uloq')}
    summary = {'all_levels_pass': all(r['passes_criteria'] for r in rows),
               'lloq_level_passes': ends['lloq']['passes_criteria'] if 'lloq' in ends else None,
               'uloq_level_passes': ends['uloq']['passes_criteria'] if 'uloq' in ends else None,
               'passing_levels': passing}
    if s['tolerance_beta'] is not None:
        summary['accuracy_profile'] = profile_range(rows, s['profile_limit_percent'])
    return {'levels': rows, 'runs': runs, 'summary': summary, 'design_diagnostics': diag}


def profile_range(rows, limit):
    """Concentrations where the beta-expectation interval stays inside +/-limit (log-linear interpolation)."""
    c = np.log([r['nominal'] for r in rows])
    g = np.array([max(abs(r['tolerance_interval']['lower_percent']), abs(r['tolerance_interval']['upper_percent'])) - limit for r in rows])
    inside = g <= 0
    if not inside.any():
        return {'status': 'no_level_inside', 'lower': None, 'upper': None, 'note': 'No tested level has its tolerance interval inside the declared limits.'}
    first, last = int(np.argmax(inside)), int(len(inside)-1-np.argmax(inside[::-1]))
    if not inside[first:last+1].all():
        return {'status': 'not_contiguous', 'lower': None, 'upper': None, 'note': 'Levels inside the limits are not contiguous; no single range is reported.'}
    def cross(i, j):
        w = g[i]/(g[i]-g[j])
        return float(np.exp(c[i]+w*(c[j]-c[i])))
    lower = float(rows[first]['nominal']) if first == 0 else cross(first-1, first)
    upper = float(rows[last]['nominal']) if last == len(c)-1 else cross(last, last+1)
    return {'status': 'estimated', 'lower': lower, 'upper': upper, 'lower_is_lowest_tested': first == 0, 'upper_is_highest_tested': last == len(c)-1,
            'note': 'Accuracy-profile limits interpolate linearly on log concentration between tested levels; not extrapolated. Exploratory supplement to the ICH M10 level-based rule.'}


def dilution_linearity(d, cfg):
    q, a = cfg['criteria'], cfg['assay']; used = d[~d.exclude].copy(); rows, diag = [], []
    used['expected_in_well'] = used.nominal/used.dilution_factor
    used['corrected'] = used.value*used.dilution_factor
    hook = used[(used.expected_in_well > a['uloq']) & (used.status != 'above_uloq')]
    for df, g in used.groupby('dilution_factor', sort=True):
        quantified = g[g.status == 'quantified']
        row = {'dilution_factor': float(df), 'n': len(g), 'n_quantified': len(quantified), 'n_series': g.series_id.nunique(),
               'expected_in_well': float(g.expected_in_well.mean()), 'statuses': g.status.value_counts().to_dict()}
        if len(quantified) == len(g) and len(g) >= 2:
            nominal = float(g.nominal.iloc[0]); m = float(quantified.corrected.mean())
            row.update(mean_corrected=m, accuracy_percent=100*(m/nominal-1), cv_percent=float(100*quantified.corrected.std(ddof=1)/m))
            row['passes_criteria'] = bool(abs(row['accuracy_percent']) <= q['accuracy_percent'] and row['cv_percent'] <= q['precision_cv_percent'])
            row['in_range'] = True
        else:
            row['in_range'] = False; row['passes_criteria'] = None
        rows.append(row)
    in_range = [r for r in rows if r['in_range']]
    if len(in_range) < 3:
        diag.append('Fewer than 3 dilution factors within the range; ICH M10 asks for at least 3.')
    if used.series_id.nunique() < 3:
        diag.append('Fewer than 3 independent dilution series; ICH M10 asks for at least 3.')
    top = used.dilution_factor.min()
    if (used[used.dilution_factor == top].expected_in_well <= a['uloq']).all():
        diag.append('No undiluted or above-ULOQ QC was tested, so the hook effect was not assessed.')
    return {'dilutions': rows, 'hook_effect_suspected': bool(len(hook)),
            'hook_rows': hook[['observation_id', 'series_id', 'dilution_factor', 'status']].to_dict('records'),
            'summary': {'all_in_range_dilutions_pass': bool(in_range and all(r['passes_criteria'] for r in in_range)),
                        'passes_without_hook': bool(in_range and all(r['passes_criteria'] for r in in_range) and not len(hook)),
                        'dilution_factor_range_validated': [min(r['dilution_factor'] for r in in_range), max(r['dilution_factor'] for r in in_range)] if in_range else None,
                        'hook_effect_suspected': bool(len(hook))},
            'design_diagnostics': diag}


def parallelism(d, cfg):
    q, s = cfg['criteria'], cfg['statistics']; level = s['confidence_level']
    used = d[~d.exclude & (d.status == 'quantified')].copy(); rows, diag = [], []
    used['corrected'] = used.value*used.dilution_factor
    fit = used[used.groupby('sample_id').dilution_factor.transform('nunique') >= 3]
    for sid, g in used.groupby('sample_id', sort=True):
        row = {'sample_id': sid, 'n_quantified_dilutions': int(g.dilution_factor.nunique())}
        if row['n_quantified_dilutions'] >= 3:
            v = g.corrected.to_numpy(); cv = float(100*v.std(ddof=1)/v.mean())
            x, yl = np.log(g.dilution_factor.to_numpy()), np.log(v)
            slope = float(np.polyfit(x, yl, 1)[0])
            row.update(mean_corrected=float(v.mean()), cv_percent=cv, passes_criteria=bool(cv <= q['parallelism_cv_percent']),
                       log_slope=slope, fold_change_over_range=float(np.exp(slope*(x.max()-x.min()))))
            if s['slope_margin'] is not None and abs(slope) > s['slope_margin']:
                row['trend_exceeds_margin'] = True
                diag.append(f'{sid}: log slope {slope:.3f} exceeds the declared margin {s["slope_margin"]} (fold change {row["fold_change_over_range"]:.2f} over the dilution range) even though the CV criterion is {"met" if row["passes_criteria"] else "not met"}.')
        else:
            row['passes_criteria'] = None
            diag.append(f'{sid}: fewer than 3 quantified dilutions; ICH M10 asks for at least three.')
        rows.append(row)
    trend = None
    if fit.sample_id.nunique() >= 1 and len(fit) > fit.sample_id.nunique()+1:
        # Common within-sample slope of log(corrected) on log(dilution factor): 0 means parallel.
        x = np.log(fit.dilution_factor.to_numpy()); y = np.log(fit.corrected.to_numpy()); g = pd.factorize(fit.sample_id)[0]
        xm = x-pd.Series(x).groupby(g).transform('mean').to_numpy(); ym = y-pd.Series(y).groupby(g).transform('mean').to_numpy()
        b = float(xm@ym/(xm@xm)); df = len(x)-g.max()-2
        if df >= 1:
            se = float(np.sqrt(((ym-b*xm)**2).sum()/df/(xm@xm))); t = stats.t.ppf((1+level)/2, df)
            trend = {'common_log_slope': b, 'interval': [b-t*se, b+t*se], 'df': int(df), 'confidence_level': level,
                     'p_value_slope_zero': float(2*stats.t.sf(abs(b/se), df)),
                     'method': 'OLS on log corrected concentration vs log dilution factor with sample intercepts'}
            if s['slope_margin'] is not None:
                trend['equivalent_within_margin'] = bool(-s['slope_margin'] <= b-t*se and b+t*se <= s['slope_margin'])
                trend['margin'] = s['slope_margin']
    evaluable = [r for r in rows if r['passes_criteria'] is not None]
    return {'samples': rows, 'trend': trend,
            'summary': {'all_samples_pass_cv': bool(evaluable and all(r['passes_criteria'] for r in evaluable)), 'n_evaluable_samples': len(evaluable)},
            'design_diagnostics': diag}


def clopper_pearson(k, n, level):
    a = 1-level
    lo = 0. if k == 0 else float(stats.beta.ppf(a/2, k, n-k+1))
    hi = 1. if k == n else float(stats.beta.ppf(1-a/2, k+1, n-k))
    return [lo, hi]


def pass_rates(d, cfg):
    e, q, a, s = cfg['experiment'], cfg['criteria'], cfg['assay'], cfg['statistics']
    used = d[~d.exclude]; keys = ['source_id', 'role'] + (['interferent'] if e == 'specificity' else []); per, diag = [], []
    for key, g in used.groupby(keys, sort=True):
        info = dict(zip(keys, key)); role = info['role']
        if role == 'blank':
            ok = bool(((g.status == 'below_lloq') | ((g.status == 'quantified') & (g.value < a['lloq']))).all())
            per.append({**info, 'n': len(g), 'passes': ok, 'criterion': 'response below LLOQ'})
            continue
        if (g.status != 'quantified').any():
            per.append({**info, 'n': len(g), 'passes': False, 'criterion': 'not quantified'}); continue
        nominal = float(g.nominal.iloc[0]); m = float(g.value.mean()); acc = 100*(m/nominal-1)
        limit = q['accuracy_percent_edge'] if role in ('lloq', 'uloq') else q['accuracy_percent']
        per.append({**info, 'n': len(g), 'mean': m, 'accuracy_percent': acc, 'limit_percent': limit, 'passes': bool(abs(acc) <= limit)})
    table = pd.DataFrame(per); groups = []
    for key, g in table.groupby(['role'] + (['interferent'] if e == 'specificity' else []), sort=True):
        key = key if isinstance(key, tuple) else (key,)
        k, n = int(g.passes.sum()), len(g); role = key[0]
        need = q['blank_pass_fraction'] if role == 'blank' else q['required_pass_fraction']
        groups.append({'role': role, **({'interferent': key[1]} if e == 'specificity' else {}), 'sources': n, 'passing': k,
                       'pass_fraction': k/n, 'required_fraction': need, 'passes_criteria': bool(k/n >= need-1e-12),
                       'pass_fraction_interval': clopper_pearson(k, n, s['confidence_level']),
                       'interval_method': 'Clopper-Pearson exact binomial (supplement; the guideline rule uses the observed fraction)'})
    if e == 'selectivity' and table.source_id.nunique() < 10:
        diag.append('Fewer than 10 individual matrix sources; ICH M10 asks for at least 10 (fewer only for rare matrices).')
    return {'sources': per, 'groups': groups, 'summary': {'all_groups_pass': all(r['passes_criteria'] for r in groups)}, 'design_diagnostics': diag}


def stability(d, cfg):
    q, s = cfg['criteria'], cfg['statistics']; level = s['confidence_level']; rows, diag = [], []
    used = d[~d.exclude]
    for (cond, lv), g in used.groupby(['condition', 'level'], sort=True):
        nominal = float(g.nominal.iloc[0]); row = {'condition': cond, 'level': lv, 'nominal': nominal, 'n': len(g)}
        if (g.status != 'quantified').any():
            row.update(passes_criteria=False, note='Some aliquots were not quantified.'); rows.append(row); continue
        v = g.value.to_numpy(); m = float(v.mean()); acc = 100*(m/nominal-1)
        row.update(mean=m, accuracy_percent=acc, passes_criteria=bool(abs(acc) <= q['accuracy_percent']))
        if len(v) >= 2:
            t = stats.t.ppf((1+level)/2, len(v)-1); se = v.std(ddof=1)/np.sqrt(len(v))
            ci = [100*((m-t*se)/nominal-1), 100*((m+t*se)/nominal-1)]
            row.update(accuracy_interval_percent=ci, interval_within_limits=bool(-q['accuracy_percent'] <= ci[0] and ci[1] <= q['accuracy_percent']),
                       interval_note='Supplement: the interval inside +/-limit is a stricter equivalence-style view; the M10 rule uses the mean.')
        if len(v) < 3:
            diag.append(f'{cond}/{lv}: fewer than 3 aliquots; ICH M10 asks for at least three.')
        rows.append(row)
    return {'conditions': rows, 'summary': {'all_conditions_pass': all(r['passes_criteria'] for r in rows)}, 'design_diagnostics': diag}


def compute(d, cfg):
    e = cfg['experiment']
    result = {'accuracy_precision': accuracy_precision, 'dilution_linearity': dilution_linearity, 'parallelism': parallelism,
              'selectivity': pass_rates, 'specificity': pass_rates, 'stability': stability}[e](d, cfg)
    result.update(schema_version=1, analysis_type='method_validation', experiment=e,
                  criteria=cfg['criteria'], excluded=d.loc[d.exclude, ['observation_id', 'exclusion_reason']].to_dict('records'))
    result['limitations'] = [
        'Acceptance criteria are the user-declared values and source; the software does not decide which guideline applies.',
        'Inputs are back-calculated concentrations; calibration-curve fitting and run acceptance are upstream.',
        'Passing these experiments is part of a validation, not a complete validation or regulatory acceptance.']
    return result
