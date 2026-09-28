"""Replicated log-RP random-run combination and nominal-level validation.

Unweighted determinations, common residual variance; technical wells are not
independent RP determinations. Per-run 4PL fitting remains in dose_fit and
parallelism equivalence in equivalence.py, via verified dose-response runs.
"""
import numpy as np
from scipy import stats
from .cmc_common import merge, common, text, number, table
from .variance_components import fit_components
from .variance_intervals import quadratic_design, precision_intervals

DEFAULTS = {'analysis_type': 'potency_assay', 'schema_version': 1, 'input': None, 'source': None,
    'mode': 'combine', 'input_mode': 'declared_estimates',
    'assay': {'independent_runs': None, 'independent_determinations': None, 'common_log_variance': None,
              'same_material_per_level': None, 'rationale': None},
    'criteria': {'source': None, 'max_abs_log_bias': None, 'max_gcv_percent': None,
                 'linearity_equivalence': False, 'slope_limits': None, 'intercept_limits': None,
                 'margins_prespecified': None},
    'suitability': {'source': None, 'require_parallelism_equivalence': True},
    'outliers': {'method': 'none', 'alpha': None, 'prespecified': None, 'source': None},
    'statistics': {'confidence_level': .95}, 'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    c = merge(raw, DEFAULTS); common(c, 'potency_assay')
    if c['mode'] not in ('combine', 'validate') or c['input_mode'] not in ('declared_estimates', 'dose_runs'):
        raise ValueError('Unsupported potency mode/input_mode')
    for k in ('independent_runs','independent_determinations','common_log_variance','same_material_per_level'):
        if c['assay'][k] is not True: raise ValueError('Declare literal true with rationale: assay.'+k)
    text(c['assay']['rationale'], 'assay.rationale'); text(c['suitability']['source'], 'suitability.source')
    if c['suitability']['require_parallelism_equivalence'] is not True:
        raise ValueError('Potency combination requires demonstrated prespecified equivalence parallelism')
    conf = c['statistics']['confidence_level']; number(conf, 'confidence_level')
    if not .5 < conf < 1: raise ValueError('confidence_level must lie between .5 and 1')
    q = c['criteria']
    if type(q['linearity_equivalence']) is not bool: raise ValueError('linearity_equivalence must be boolean')
    if c['mode'] == 'validate':
        text(q['source'], 'criteria.source')
        for k in ('max_abs_log_bias', 'max_gcv_percent'): number(q[k], k, positive=True)
    elif any(q[k] is not None for k in ('max_abs_log_bias','max_gcv_percent','slope_limits','intercept_limits')) or q['linearity_equivalence']:
        raise ValueError('Validation criteria require mode=validate')
    if q['linearity_equivalence']:
        if q['margins_prespecified'] is not True: raise ValueError('Equivalence margins must be prespecified, never selected from these data')
        for k, target in (('slope_limits',1.),('intercept_limits',0.)):
            v = q[k]
            if not isinstance(v,list) or len(v)!=2: raise ValueError(k+' requires two declared limits')
            for edge in v: number(edge,k)
            if not v[0] < target < v[1]: raise ValueError(k+' must bracket its target')
    elif q['slope_limits'] is not None or q['intercept_limits'] is not None:
        raise ValueError('Declare linearity_equivalence to use margins')
    o = c['outliers']
    if o['method'] not in ('none','grubbs_run_means'): raise ValueError('Only prespecified single Grubbs screening of run means is supported')
    if o['method'] != 'none':
        if c['mode'] != 'combine': raise ValueError('Run-mean Grubbs screening applies only to single-level combination')
        if o['prespecified'] is not True: raise ValueError('Outlier screening must be prespecified')
        text(o['source'],'outliers.source'); number(o['alpha'],'outliers.alpha')
        if not 0 < o['alpha'] < .5: raise ValueError('outliers.alpha must lie between 0 and .5')
    elif any(o[k] is not None for k in ('alpha','prespecified','source')):
        raise ValueError('Outlier settings supplied with method=none')
    return c


def load_data(path, cfg):
    required = ['observation_id','run_id','nominal_rp']
    if cfg['input_mode'] == 'dose_runs':
        return table(path, required+['dose_run','comparison_id'], ['nominal_rp'])
    d = table(path, required+['relative_potency','reference_quality','parallelism','fit_source'], ['nominal_rp','relative_potency'])
    return validate_data(d, cfg)


def validate_data(d, cfg):
    if (d.nominal_rp <= 0).any() or (d.relative_potency <= 0).any():
        raise ValueError('RP and nominal RP must be positive')
    for k in ('reference_quality','parallelism'):
        if set(d[k]) - {'pass','fail','not_evaluable'}: raise ValueError(k+' requires pass/fail/not_evaluable')
    if 'exclude' in d and (d.exclude != 'false').any():
        raise ValueError('No automatic or manual run dropping in this scope; retain failing runs and withhold combination')
    if d.run_id.nunique() < 3: raise ValueError('At least three independent runs required')
    levels = sorted(d.nominal_rp.unique())
    if cfg['mode'] == 'combine' and len(levels) != 1: raise ValueError('Combine one nominal RP level at a time')
    if cfg['mode'] == 'validate' and len(levels) < 3: raise ValueError('Validation needs at least three nominal RP levels')
    counts = d.groupby(['run_id','nominal_rp']).size().unstack(fill_value=0)
    if (counts.to_numpy() < 2).any():
        raise ValueError('Each run needs at least two independent RP determinations at every level to separate run and residual variance')
    return d


def random_run(y, runs, x=None, confidence=.95):
    y = np.asarray(y,float); runs = np.asarray(runs)
    x = np.ones((len(y),1)) if x is None else np.asarray(x,float)
    kernel = (runs[:,None] == runs[None,:]).astype(float)
    fit = fit_components(y, x, [kernel], ['run'], confidence)
    design = quadratic_design(x, [kernel])
    v = [fit['components'][k]['variance'] for k in ('run','residual')]
    centered = y-x@np.linalg.lstsq(x,y,rcond=None)[0]
    interval = precision_intervals(centered, design, v, confidence)['sums']['total']
    fit['intermediate_precision_mls'] = interval
    variance = fit['intermediate_precision']['variance']
    fit['gcv_percent'] = float(100*np.sqrt(np.expm1(variance)))
    fit['gcv_interval_percent'] = [float(100*np.sqrt(np.expm1(edge))) for edge in interval['variance_interval']]
    df = len(np.unique(runs))-1; q = stats.t.ppf((1+confidence)/2,df)
    beta = np.array(fit['fixed_coefficients']); se = np.sqrt(np.diag(fit['fixed_covariance']))
    fit['fixed_intervals'] = np.c_[beta-q*se,beta+q*se].tolist()
    fit['fixed_interval_df'] = df
    fit['fixed_interval_method'] = 'GLS covariance with t quantile, independent runs minus one df; approximate outside balanced intercept-only design'
    return fit


def grubbs(d, cfg):
    o = cfg['outliers']
    if o['method'] == 'none': return {'method':'none','flagged_runs':[], 'dropped_runs':[]}
    means = d.assign(log_rp=np.log(d.relative_potency)).groupby('run_id').log_rp.mean()
    counts = d.groupby('run_id').size()
    if counts.nunique() != 1: raise ValueError('Grubbs screening requires equal run sizes (exchangeable run means)')
    n=len(means); sd=float(means.std(ddof=1)); dev=abs(means-means.mean())
    g=float(dev.max()/sd) if sd>0 else 0.
    t=stats.t.ppf(1-o['alpha']/(2*n),n-2); critical=(n-1)/np.sqrt(n)*np.sqrt(t*t/(n-2+t*t))
    flagged=[str(means.index[int(np.argmax(dev))])] if g>critical else []
    return {'method':o['method'],'alpha':o['alpha'],'source':o['source'], 'statistic':g,'critical':float(critical),
            'flagged_runs':flagged,'dropped_runs':[], 'note':'Single two-sided screen under Gaussian run means; flags do not remove observations.'}


def compute(d, cfg):
    validate_data(d,cfg); conf=cfg['statistics']['confidence_level']; q=cfg['criteria']
    screening=grubbs(d,cfg)
    per_run=[]
    for run,g in d.groupby('run_id',sort=True):
        failed=g[(g.reference_quality!='pass') | (g.parallelism!='pass')]
        per_run.append({'run_id':str(run),'determinations':len(g),'suitable':not len(failed),
            'flagged_outlier':run in screening['flagged_runs'],
            'determination_diagnostics':g[['observation_id','reference_quality','parallelism','fit_source']].to_dict('records')})
    failures=[{'run_id':r['run_id'],'reason':'system_suitability_failed_or_not_evaluable'} for r in per_run if not r['suitable']]
    notes=['All supplied runs are retained; a failing or unevaluable run withholds the combined result.',
           'Natural-log RP model; run random intercept plus independent determination residual, common log variance.',
           'Intermediate precision is for one future determination; it includes run and residual variance.']
    if screening['flagged_runs']: notes.append('Prespecified Grubbs flags: '+', '.join(screening['flagged_runs'])+'; every flagged run remains included.')
    result={'analysis_type':'potency_assay','schema_version':1,'mode':cfg['mode'],'per_run':per_run,'outlier_screening':screening,
        'failing_items':failures,'must_mention':notes,'criteria':q,
        'limitations':['Replicated RP estimates must be independent within run (no shared reference fits or technical-well pseudoreplication).',
            'Unweighted estimates with common log residual variance; heterogeneous known-SE meta-analysis is not implemented.',
            'Fixed-effect t intervals use runs minus one df; MLS/MOVER and unbalanced inference are approximate.',
            'Passing declared criteria does not establish USP compliance or full assay validation.']}
    if failures:
        result.update(primary={'status':'withheld','reason':'system_suitability','combined_rp':None},levels=[],linearity=None)
        return result
    levels=[]
    for nominal,g in d.groupby('nominal_rp',sort=True):
        fit=random_run(np.log(g.relative_potency),g.run_id,confidence=conf)
        mu=fit['fixed_coefficients'][0]; interval=fit['fixed_intervals'][0]; bias=mu-np.log(nominal)
        row={'nominal_rp':float(nominal),'combined_rp':float(np.exp(mu)), 'combined_rp_interval':np.exp(interval).tolist(),
             'log_bias':float(bias),'log_bias_interval':(np.asarray(interval)-np.log(nominal)).tolist(),
             'relative_bias_percent':float(100*np.expm1(bias)),'fit':fit}
        if cfg['mode']=='validate':
            row['checks']={'log_bias':bool(abs(bias)<=q['max_abs_log_bias']), 'intermediate_gcv':bool(fit['gcv_percent']<=q['max_gcv_percent'])}
            row['passes']=all(row['checks'].values())
            if not row['passes']: failures.append({'nominal_rp':float(nominal),'reason':'declared_level_criteria','checks':row['checks']})
        levels.append(row)
        if fit['boundary_components']: notes.append(f'Nominal RP {nominal}: boundary components {fit["boundary_components"]}; component intervals may be unavailable.')
    result['levels']=levels
    if cfg['mode']=='combine':
        result['primary']={'status':'estimated','combined_rp':levels[0]['combined_rp'],'interval':levels[0]['combined_rp_interval'],
                           'confidence_level':conf,'gcv_percent':levels[0]['fit']['gcv_percent'],'gcv_interval_percent':levels[0]['fit']['gcv_interval_percent']}
        result['linearity']=None
    else:
        x=np.c_[np.ones(len(d)),np.log(d.nominal_rp)]
        fit=random_run(np.log(d.relative_potency),d.run_id,x,conf)
        lin={'intercept':fit['fixed_coefficients'][0],'slope':fit['fixed_coefficients'][1],
             'intervals':fit['fixed_intervals'],'confidence_level':conf,'fit':fit,'equivalence':None}
        if q['linearity_equivalence']:
            checks={k:bool(lim[0]<ci[0] and ci[1]<lim[1]) for k,lim,ci in zip(('intercept','slope'),(q['intercept_limits'],q['slope_limits']),fit['fixed_intervals'])}
            lin['equivalence']={'checks':checks,'equivalent':all(checks.values()),'alpha_each_side':(1-conf)/2,
                                 'source':q['source'],'limits':{'intercept':q['intercept_limits'],'slope':q['slope_limits']}}
            if not all(checks.values()): failures.append({'reason':'linearity_equivalence_not_demonstrated','checks':checks})
        # Range contains contiguous tested passing levels only; no fitted extension through failed levels.
        blocks=[]; block=[]
        for row in levels:
            if row['passes']: block.append(row['nominal_rp'])
            elif block: blocks.append(block); block=[]
        if block: blocks.append(block)
        linear_pass=lin['equivalence'] is not None and lin['equivalence']['equivalent']
        validated=[{'lower':b[0],'upper':b[-1],'tested_levels':b} for b in blocks if len(b)>=3] if linear_pass else []
        if not q['linearity_equivalence']: notes.append('Linearity intervals are statistical supplements; no declared equivalence margins, so the validated range is withheld.')
        notes.append('Bias and precision acceptance use declared point-estimate criteria; their intervals are labelled statistical supplements.')
        result['linearity']=lin
        result['primary']={'status':'evaluated','all_levels_pass':all(row['passes'] for row in levels),
            'validated_ranges':validated,'range_status':'estimated_on_tested_levels' if validated else 'withheld',
            'range_rule':'At least three contiguous passing tested levels plus declared global linearity equivalence; no extrapolation or interpolation claim.'}
    return result
