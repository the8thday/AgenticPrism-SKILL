"""One-sample and ratio comparisons on independent experimental units (0.13.4)."""
from copy import deepcopy
import numpy as np
from scipy import stats
from .groups import read_unit_table

TYPES = {'location_test'}
METHODS = ('one_sample_t', 'one_sample_ratio_t', 'paired_ratio_t', 'independent_ratio_welch')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'location_test', 'input': None, 'source': None, 'method': None,
            'design': {'unit': None, 'rationale': None, 'outcome': None, 'outcome_unit': None},
            'comparison': {'group_a': None, 'group_b': None, 'null_value': None, 'null_source': None, 'confidence_level': .95},
            'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw)-set(DEFAULTS): raise ValueError('Unsupported location-test configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v)-set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'location_test' or c['method'] not in METHODS:
        raise ValueError('Declare a supported location-test method')
    q = c['comparison']
    for name, v in [('input',c['input']),('source',c['source']),*c['design'].items(),('group_a',q['group_a']),('null_source',q['null_source'])]:
        if not isinstance(v,str) or not v.strip(): raise ValueError('Declare '+name)
    one = c['method'].startswith('one_sample')
    if one and q['group_b'] is not None: raise ValueError('One-sample tests have one group; group_b must be null')
    if not one and (not isinstance(q['group_b'],str) or not q['group_b'].strip() or q['group_a']==q['group_b']):
        raise ValueError('Declare two different groups for ratio comparison')
    null = q['null_value']
    if type(null) not in (int,float) or not np.isfinite(null): raise ValueError('Declare a finite null_value with its source')
    if c['method'] != 'one_sample_t' and null <= 0: raise ValueError('Ratio tests require a positive null_value')
    lv = q['confidence_level']
    if type(lv) not in (float,int) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard','prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = read_unit_table(path); q = c['comparison']
    groups = {q['group_a']} if c['method'].startswith('one_sample') else {q['group_a'],q['group_b']}
    if set(d.group)!=groups or set(d.outcome)!={c['design']['outcome']} or set(d.unit)!={c['design']['outcome_unit']}:
        raise ValueError('Input group/outcome/unit does not match the declared comparison')
    used = d[~d.exclude]
    if used.duplicated(['independent_unit_id','group']).any(): raise ValueError('Technical replicates are not independent units')
    if c['method']=='paired_ratio_t':
        if set(used[used.group==q['group_a']].independent_unit_id)!=set(used[used.group==q['group_b']].independent_unit_id):
            raise ValueError('Paired ratio tests require complete matching unit IDs; no silent pair deletion')
    elif used.independent_unit_id.duplicated().any(): raise ValueError('Independent groups share unit IDs')
    if c['method']!='one_sample_t' and (used.value<=0).any(): raise ValueError('Ratio tests require strictly positive observations; no pseudocounts')
    return d


def t_summary(values, null=0., level=.95):
    """Student interval/test for the mean on the supplied analysis scale."""
    v = np.asarray(values,float)
    if len(v)<3: return {'status':'withheld','reportable':False,'reason':'fewer_than_three_units_or_pairs'}
    estimate = float(v.mean()); se = float(v.std(ddof=1)/np.sqrt(len(v))); df = len(v)-1
    if se <= np.finfo(float).eps*max(1.,float(np.max(np.abs(v)))):
        return {'status':'withheld','reportable':False,'reason':'zero_or_numerically_constant_variance'}
    t = (estimate-null)/se; h = float(stats.t.ppf((1+level)/2,df)*se)
    return {'status':'estimated','reportable':True,'estimate':estimate,'ci':[estimate-h,estimate+h],
            'standard_error':se,'df':float(df),'t_statistic':float(t),'p_two_sided':float(2*stats.t.sf(abs(t),df))}


def compute(d, c):
    q = c['comparison']; m = c['method']; lv = q['confidence_level']; u = d[~d.exclude]
    a = u[u.group==q['group_a']].set_index('independent_unit_id').value
    b = u[u.group==q['group_b']].set_index('independent_unit_id').value
    ratio = m!='one_sample_t'; null = float(q['null_value']); analysis_null = float(np.log(null)) if ratio else null
    if ratio and (u.value<=0).any(): raise ValueError('Ratio tests require strictly positive values')
    if u.duplicated(['independent_unit_id','group']).any(): raise ValueError('Duplicate independent units')
    if m=='paired_ratio_t':
        if set(a.index)!=set(b.index): raise ValueError('Complete matching pairs required')
        b = b.reindex(a.index)
        z = t_summary(np.log(b.to_numpy(float))-np.log(a.to_numpy(float)),analysis_null,lv)
    elif m.startswith('one_sample'):
        z = t_summary(np.log(a.to_numpy(float)) if ratio else a.to_numpy(float),analysis_null,lv)
    else:
        if set(a.index)&set(b.index): raise ValueError('Independent groups share unit IDs')
        if min(len(a),len(b))<3: z={'status':'withheld','reportable':False,'reason':'fewer_than_three_units_per_group'}
        else:
            x,y=np.log(a.to_numpy(float)),np.log(b.to_numpy(float))
            va,vb=x.var(ddof=1)/len(x),y.var(ddof=1)/len(y);se=float(np.sqrt(va+vb));est=float(y.mean()-x.mean())
            if se <= np.finfo(float).eps*max(1.,float(np.max(np.abs(np.r_[x,y])))):
                z={'status':'withheld','reportable':False,'reason':'zero_or_numerically_constant_variance'}
            else:
                df=float((va+vb)**2/(va**2/(len(x)-1)+vb**2/(len(y)-1)));h=float(stats.t.ppf((1+lv)/2,df)*se);t=(est-analysis_null)/se
                z={'status':'estimated','reportable':True,'estimate':est,'ci':[est-h,est+h],'standard_error':se,'df':df,'t_statistic':t,'p_two_sided':float(2*stats.t.sf(abs(t),df))}
    z.update(method=m,group_a=q['group_a'],group_b=q['group_b'],n_a=len(a),n_b=len(b),null_value=null,null_source=q['null_source'],
             confidence_level=lv,analysis_scale='natural_log' if ratio else 'original',
             estimand='arithmetic_mean' if not ratio else 'geometric_mean' if m=='one_sample_ratio_t' else 'geometric_mean_ratio_B_over_A',
             outcome=c['design']['outcome'],unit='ratio' if ratio and not m.startswith('one_sample') else c['design']['outcome_unit'])
    if m=='paired_ratio_t': z['n_pairs']=len(a)
    summaries=[]
    for g in [q['group_a']]+([q['group_b']] if q['group_b'] else []):
        v=u[u.group==g].value.to_numpy(float)
        summaries.append({'group':g,'n':len(v),'arithmetic_mean':float(v.mean()) if len(v) else None,
                          'geometric_mean':float(np.exp(np.log(v).mean())) if ratio and len(v) else None})
    z['group_summaries']=summaries
    if z['reportable'] and ratio:
        z['log_estimate']=z['estimate'];z['log_ci']=list(z['ci']);z['log_standard_error']=z.pop('standard_error')
        with np.errstate(over='ignore',under='ignore'):
            transformed=np.exp([z['estimate'],*z['ci']])
        if not np.isfinite(transformed).all() or (transformed<=0).any():
            z.update(status='withheld',reportable=False,reason='backtransform_outside_numeric_range')
            z.pop('estimate');z.pop('ci');z.pop('p_two_sided')
        else: z.update(estimate=float(transformed[0]),ci=transformed[1:].tolist())
    notes=['The independent unit is '+c['design']['unit']+'. No automatic exclusions, normalization or pseudocounts.',
           'A nonsignificant test against the declared null does not establish equivalence or assay acceptance.']
    if ratio: notes.append('Ratio inference is on natural logarithms and back-transformed; it estimates geometric means/ratios, not arithmetic mean ratios. Normality concerns the analysis scale (paired log ratios for a paired test).')
    if m=='paired_ratio_t':notes.append('Pairing is by independent_unit_id. The numerator is group B and denominator group A.')
    if not z['reportable']:notes.append('Inference withheld: '+z['reason'])
    return {'analysis_type':'location_test','primary':z,'must_mention':notes,'failing_items':[] if z['reportable'] else [{'reason':z['reason']}],
            'limitations':['Two-sided Student/Welch inference under the declared unit and distribution assumptions; no robustness or equivalence claim.']}


def render(run, result, style):
    """Display saved unit observations and the saved effect interval; no inference."""
    from pathlib import Path
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p=result['primary'];d=pd.read_csv(Path(run)/'input.csv',dtype={'independent_unit_id':str,'group':str})
    if 'exclude' in d:d=d[~d.exclude.astype(str).str.lower().eq('true')]
    token=pstyle.THEMES[style];out=Path(run)/'figures';out.mkdir(exist_ok=True)
    groups=[p['group_a']]+([p['group_b']] if p['group_b'] else [])
    with plt.rc_context(pstyle.rc(token,*font_setup())):
        fig,ax=plt.subplots(figsize=(5,3.5),layout='constrained')
        if p['method']=='paired_ratio_t':
            w=d.pivot(index='independent_unit_id',columns='group',values='value')
            for _,r in w.iterrows():ax.plot([0,1],[r[groups[0]],r[groups[1]]],color=pstyle.MUTED,alpha=.35,lw=.7)
        for i,g in enumerate(groups):
            v=d[d.group==g].value.to_numpy(float);ax.scatter(np.full(len(v),i),v,s=22,color=pstyle.color(token,i))
        if p['analysis_scale']=='natural_log':ax.set_yscale('log')
        if len(groups)==1:ax.axhline(p['null_value'],ls='--',color=pstyle.MUTED,lw=1)
        ax.set(xticks=range(len(groups)),xticklabels=groups,ylabel=p['outcome'],title=p['method'].replace('_',' '))
        pstyle.save(fig,out/'location_test',token,300)
    return ['location_test']
