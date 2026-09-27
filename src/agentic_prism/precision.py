"""Declared nested/crossed random-intercept precision models."""
from copy import deepcopy
import numpy as np
import pandas as pd
from .variance_components import fit_components
DEFAULTS={'analysis_type':'variance_components','schema_version':1,'input':None,'source':None,'response_unit':None,
 'design':{'random_terms':[],'fixed_numeric':[],'fixed_categorical':[],
           'independent_random_effects':False,'gaussian_homoscedastic_errors':False,
           'representative_levels':False,'rationale':None},
 'inference':{'confidence_level':.95,'cv_reference':None,'cv_rationale':None,'sum_interval':'mls','mean_square_intervals_applicable':None,'interval_rationale':None,'positive_sums':{}},
 'report':{'plot_style':'prism_like'}}


def resolve_config(raw):
    def merge(dst,src):
        if not isinstance(src,dict):raise ValueError('Config sections must be objects')
        for k,v in src.items():
            if k not in dst:raise ValueError('Unknown config key: '+k)
            if isinstance(dst[k],dict) and k!='positive_sums':merge(dst[k],v)
            else:dst[k]=v
    c=deepcopy(DEFAULTS);merge(c,raw)
    if c['analysis_type']!='variance_components' or c['schema_version']!=1:raise ValueError('Unsupported analysis/schema')
    for k in ('input','source','response_unit'):
        if not isinstance(c[k],str) or not c[k].strip():raise ValueError(k+' required')
    d=c['design']
    for k in ('independent_random_effects','gaussian_homoscedastic_errors','representative_levels'):
        if d[k] is not True:raise ValueError(k+' must be literal true with design evidence')
    if not isinstance(d['rationale'],str) or not d['rationale'].strip():raise ValueError('Design rationale required')
    terms=d['random_terms']
    if not isinstance(terms,list) or not terms or len(terms)>8:raise ValueError('Declare 1 to 8 random terms')
    for t in terms:
        if not isinstance(t,list) or not t or any(not isinstance(k,str) for k in t) or len(set(t))!=len(t):raise ValueError('Each random term is a list of distinct factor columns')
    for cols in [*terms,d['fixed_numeric'],d['fixed_categorical']]:
        if not isinstance(cols,list) or any(not isinstance(s,str) or not s.strip() or s in ('value','observation_id') for s in cols):raise ValueError('Invalid design columns')
        if len(set(cols))!=len(cols):raise ValueError('Duplicate design columns')
    if len({tuple(sorted(t)) for t in terms})!=len(terms):raise ValueError('Duplicate random terms')
    if set(d['fixed_numeric']) & set(d['fixed_categorical']):raise ValueError('A fixed column cannot be numeric and categorical')
    inf=c['inference'];level=inf['confidence_level']
    if isinstance(level,bool) or not isinstance(level,(int,float)) or not .5<level<1:raise ValueError('confidence_level must lie between .5 and 1')
    if inf['sum_interval'] not in ('satterthwaite','mls'):raise ValueError('sum_interval must be satterthwaite or mls')
    # 0.9.3: MLS/MOVER is the default for sums (Satterthwaite undercovered, 83-85%, in
    # the 0.9.1 calibration). Mean-square applicability is checked from the design
    # (quadratic_design refuses confounded strata); an explicit declaration is optional.
    if not any(inf['mean_square_intervals_applicable'] is v for v in (None,True,False)):raise ValueError('mean_square_intervals_applicable must be true, false or null')
    if inf['interval_rationale'] is not None and (not isinstance(inf['interval_rationale'],str) or not inf['interval_rationale'].strip()):raise ValueError('interval_rationale must be a nonempty string when given')
    if inf['sum_interval']=='mls' and inf['mean_square_intervals_applicable'] is False:raise ValueError('MLS was selected but mean-square intervals were declared inapplicable; choose satterthwaite explicitly')
    if inf['sum_interval']=='satterthwaite' and inf['mean_square_intervals_applicable'] is True:raise ValueError('Mean-square applicability was declared but satterthwaite was selected')
    sums=inf['positive_sums'];labels={':'.join(t) for t in terms}|{'residual'}
    if not isinstance(sums,dict):raise ValueError('positive_sums must be an object')
    for name,weights in sums.items():
        if not isinstance(name,str) or not name.strip() or name=='total' or not isinstance(weights,dict) or not weights or set(weights)-labels:raise ValueError('Invalid positive sum name/components')
        if any(type(v) not in (int,float) or not np.isfinite(v) or v<=0 for v in weights.values()):raise ValueError('Sum weights must be finite and positive')
    if sums and inf['sum_interval']!='mls':raise ValueError('Additional sums currently require mls')
    ref=inf['cv_reference']
    if ref is not None:
        if isinstance(ref,bool) or not isinstance(ref,(int,float)) or not np.isfinite(ref) or ref<=0:raise ValueError('cv_reference must be a finite positive declared reference')
        if not isinstance(inf['cv_rationale'],str) or not inf['cv_rationale'].strip():raise ValueError('CV reference rationale required')
    if c['report']['plot_style'] not in ('prism_like','standard'):raise ValueError('Unsupported plot style')
    return c


def load_data(path,cfg):
    d=pd.read_csv(path,dtype=str,keep_default_na=False)
    des=cfg['design']; required={'observation_id','value',*des['fixed_numeric'],*des['fixed_categorical'],*[s for t in des['random_terms'] for s in t]}
    if not required.issubset(d.columns):raise ValueError('Missing columns: '+str(required-set(d.columns)))
    if d[list(required)].eq('').any().any():raise ValueError('Missing values require explicit upstream review; no silent deletion')
    if d.observation_id.duplicated().any():raise ValueError('Duplicate observation_id')
    for k in ['value',*des['fixed_numeric']]:
        d[k]=pd.to_numeric(d[k],errors='raise')
        if not np.isfinite(d[k]).all():raise ValueError('Non-finite '+k)
    return d.sort_values('observation_id').reset_index(drop=True)


def matrices(d,cfg):
    des=cfg['design'];n=len(d)
    if n>1000:raise ValueError('Dense REML limited to 1000 observations; do not subsample silently')
    xs=[np.ones(n)];xn=['intercept']
    for col in des['fixed_numeric']:xs.append(d[col].to_numpy(float));xn.append(col)
    for col in des['fixed_categorical']:
        levels=sorted(d[col].unique())
        if len(levels)<2:raise ValueError('Fixed factor must have at least two levels')
        for level in levels[1:]:xs.append((d[col]==level).to_numpy(float));xn.append(col+'='+level)
    ks=[];names=[];counts={}
    for i,term in enumerate(des['random_terms']):
        codes=pd.factorize(pd.MultiIndex.from_frame(d[term]),sort=True)[0]
        levels=len(set(codes))
        if levels<3 or levels>=n:raise ValueError('Random terms need >=3 levels and replicated levels')
        ks.append((codes[:,None]==codes[None,:]).astype(float));name=':'.join(term)
        names.append(name);counts[name]=levels
    return np.asarray(xs).T,ks,names,xn,counts


def compute(d,cfg):
    x,ks,names,xn,counts=matrices(d,cfg);inf=cfg['inference']
    result=fit_components(d.value.to_numpy(float),x,ks,names,inf['confidence_level'],inf['cv_reference'])
    if inf['sum_interval']=='mls':
        from .variance_intervals import quadratic_design,precision_intervals
        labels=result['component_order'];weights={'total':np.ones(len(labels)),**{name:np.array([w.get(k,0.) for k in labels]) for name,w in inf['positive_sums'].items()}}
        variances=[result['components'][k]['variance'] for k in labels]
        residual=d.value.to_numpy(float)-x@np.linalg.lstsq(x,d.value.to_numpy(float),rcond=None)[0]
        try:design=quadratic_design(x,ks)
        except ValueError as exc:raise ValueError(f'MLS/MOVER intervals are not available for this design ({exc}); set inference.sum_interval to satterthwaite explicitly') from exc
        detail=precision_intervals(residual,design,variances,inf['confidence_level'],weights)
        result['method']='REML point estimates; MLS/MOVER sum intervals; Satterthwaite individual two-sided intervals'
        result['mean_square_inference']=detail;result['positive_sums']={}
        result['limitations']=[t for t in result['limitations'] if not t.startswith('Satterthwaite intervals') and not t.startswith('Boundary component intervals')]
        result['limitations'].append('Individual two-sided intervals retain Satterthwaite and may be unavailable at a boundary; separate one-sided upper bounds use mean squares.')
        for name,interval in detail['sums'].items():
            summary=result['intermediate_precision'] if name=='total' else {'variance':float(weights[name]@variances)}
            summary['sd']=float(np.sqrt(summary['variance']));summary['cv_percent']=100*summary['sd']/inf['cv_reference'] if inf['cv_reference'] else None
            if name=='total':summary['satterthwaite_interval']=summary['variance_interval'];summary['satterthwaite_df']=summary.pop('df')
            summary.update(interval_method=interval['method'],interval_center=interval['center'],variance_interval=interval['variance_interval'],interval_status='estimated')
            summary['sd_interval']=np.sqrt(interval['variance_interval']).tolist();summary['cv_interval']=(100*np.sqrt(interval['variance_interval'])/inf['cv_reference']).tolist() if inf['cv_reference'] else None
            if name!='total':result['positive_sums'][name]=summary
        for k,u in zip(labels,detail['component_upper_bounds']):result['components'][k]['one_sided_upper_bound']=u
        result['limitations'].append(detail['limitation'])
        result['limitations'].append('REML point estimates are retained; MLS/MOVER intervals are centered on unconstrained ANOVA quadratic estimates, which can differ at boundaries and in unbalanced designs. One-sided component upper bounds apply to every component, not only selected zero estimates.')
    result.update(response_unit=cfg['response_unit'],analysis_type='variance_components',fixed_coefficient_names=xn,random_levels=counts,
                  status='boundary_caution' if result['boundary_components'] else 'estimated')
    if inf['cv_reference'] is None:result['limitations'].append('CV unavailable: no positive reference and rationale declared.')
    return result
