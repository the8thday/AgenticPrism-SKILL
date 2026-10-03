"""Declared HTS plate quality, median-polish B scores and exploratory FDR hits."""
import zlib
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests
from .cmc_common import merge,number,text
TYPES={'hts_qc'}
DEFAULTS={'schema_version':1,'analysis_type':'hts_qc','input':None,'source':None,
 'plate_shape':None,'readout':None,'response_unit':None,'direction':None,'independent_wells':False,'independence_source':None,
 'randomized_layout':False,'layout_source':None,'majority_inactive':False,'majority_inactive_source':None,
 'correction':'median_polish','mad_scale':1.4826,'mad_scale_source':None,
 'median_polish':{'epsilon':.01,'max_iterations':10},
 'criteria':{'predeclared':False,'source':None,'zprime_min':.5,'ssmd_abs_min':3.,'fdr':.05},
 'hit_reference':{'method':'predictive_t','simulations':999,'seed':None,'source':None},
 'report':{'plot_style':'prism_like'}}

def resolve_config(raw):
 c=merge(raw,DEFAULTS)
 if c['analysis_type']!='hts_qc' or c['schema_version']!=1:raise ValueError('Unsupported HTS schema')
 for k in ('input','source','readout','response_unit','independence_source','layout_source','majority_inactive_source','mad_scale_source'):text(c[k],k)
 if c['plate_shape'] not in ([8,12],[16,24]):raise ValueError('Declare plate_shape as [8,12] or [16,24] before analysis')
 if c['direction'] not in ('increasing','decreasing'):raise ValueError('Declare hit direction')
 if c['independent_wells'] is not True or c['randomized_layout'] is not True:raise ValueError('Exploratory well inference requires independent wells and a declared randomized layout')
 if c['majority_inactive'] is not True:raise ValueError('Median polish requires a justified majority of inactive wells')
 if c['correction']!='median_polish':raise ValueError('Hit calling requires prespecified plate-effect correction')
 number(c['mad_scale'],'MAD scale',True);m=c['median_polish'];number(m['epsilon'],'median-polish epsilon',True)
 if type(m['max_iterations'])is not int or not 1<=m['max_iterations']<=1000:raise ValueError('Invalid median-polish iteration limit')
 a=c['criteria'];text(a['source'],'criteria source')
 if a['predeclared'] is not True:raise ValueError('QC and FDR criteria must be declared before analysis')
 for k in ('zprime_min','ssmd_abs_min','fdr'):number(a[k],k)
 if not 0<a['fdr']<1 or not 0<=a['zprime_min']<1 or a['ssmd_abs_min']<=0:raise ValueError('Invalid QC criteria')
 if c['report']['plot_style'] not in ('prism_like','standard'):raise ValueError('Invalid style')
 h=c['hit_reference']
 if h['method']=='predictive_t':
  if h!=DEFAULTS['hit_reference']:raise ValueError('Simulation settings apply only to hit_reference.method=layout_simulation')
  c.pop('hit_reference')  # Preserve legacy resolved config bytes when opt-in is absent.
 elif h['method']=='layout_simulation':
  text(h['source'],'hit_reference source')
  if type(h['simulations']) is not int or not 199<=h['simulations']<=9999:raise ValueError('hit_reference.simulations must be an integer in [199, 9999]')
  if type(h['seed']) is not int or h['seed']<0:raise ValueError('Declare a non-negative integer hit_reference.seed')
 else:raise ValueError('hit_reference.method must be predictive_t or layout_simulation')
 return c


def load_data(path,cfg):
 d=pd.read_csv(path,dtype={'plate_id':str,'well_id':str,'compound_id':str,'independent_unit_id':str})
 req={'plate_id','well_id','compound_id','independent_unit_id','row','column','role','value'}
 if not req<=set(d) or d.empty:raise ValueError('Require canonical plate/well/compound/unit identities, row/column, role and value')
 for k in ('plate_id','well_id','compound_id','independent_unit_id'):
  if d[k].isna().any() or d[k].str.strip().eq('').any():raise ValueError('Missing well identity')
 if d.duplicated(['plate_id','well_id']).any() or d.duplicated(['plate_id','row','column']).any():raise ValueError('Duplicate well')
 if d.independent_unit_id.duplicated().any():raise ValueError('Repeated unit is not an independent well')
 if not d.role.isin(['positive','negative','sample']).all():raise ValueError('Explicit positive/negative/sample roles required')
 for k in ('row','column','value'):
  d[k]=pd.to_numeric(d[k],errors='raise')
  if not np.isfinite(d[k]).all():raise ValueError('Nonfinite plate data')
 for k in ('row','column'):
  if (d[k]<1).any() or (d[k]%1!=0).any():raise ValueError('Positive integer grid positions required')
 for _,g in d.groupby('plate_id'):
  nr,nc=cfg['plate_shape']
  if set(g.row)!=set(range(1,nr+1)) or set(g.column)!=set(range(1,nc+1)) or len(g)!=nr*nc:raise ValueError('Missing well, entire row or column relative to declared plate_shape')
  if len(g)!=g.row.nunique()*g.column.nunique():raise ValueError('Complete rectangular plate required; do not impute missing wells')
  for role in ('positive','negative'):
   h=g[g.role==role]
   if len(h)<4 or h.row.nunique()<2 or h.column.nunique()<2:raise ValueError('Each control needs >=4 wells spanning >=2 rows and columns')
 return d


def median_polish(x,epsilon=.01,max_iterations=10):
 z=np.asarray(x,float).copy();nr,nc=z.shape;r=np.zeros(nr);c=np.zeros(nc);overall=0.;oldsum=0.;converged=False
 for it in range(max_iterations):
  delta=np.median(z,axis=1);z-=delta[:,None];r+=delta;delta=np.median(c);c-=delta;overall+=delta
  delta=np.median(z,axis=0);z-=delta[None,:];c+=delta;delta=np.median(r);r-=delta;overall+=delta
  newsum=float(np.sum(abs(z)));converged=newsum==0 or abs(newsum-oldsum)<epsilon*newsum
  if converged:break
  oldsum=newsum
 return dict(overall=float(overall),row=r,column=c,residuals=z,iterations=it+1,converged=bool(converged))


def quality(positive,negative):
 p,n=np.asarray(positive),np.asarray(negative);difference=float(p.mean()-n.mean());sp=float(p.std(ddof=1));sn=float(n.std(ddof=1));den=np.sqrt(sp*sp+sn*sn)
 return {'zprime':float(1-3*(sp+sn)/abs(difference)) if difference else None,'ssmd':float(difference/den) if den else None,'positive_mean':float(p.mean()),'negative_mean':float(n.mean()),'positive_sd':sp,'negative_sd':sn}


def _statistic(x,negative,sign,cfg):
 z=median_polish(x,cfg['median_polish']['epsilon'],cfg['median_polish']['max_iterations'])['residuals'];n=z[negative]
 return sign*(z-n.mean())/(n.std(ddof=1)*np.sqrt(1+1/n.size))


def layout_reference(x,roles,sign,cfg,plate_id):
 """Layout-conditional parametric null for the predictive statistic of each sample well.

 Null plates keep the observed additive plate pattern, the positive-control offset and
 the control/sample layout; every non-positive well gets independent Gaussian noise with
 the negative-control residual SD. Each simulated plate is median-polished and scored
 exactly like the observed plate. Statistics are standardized per well by their null
 mean and SD, then referred to the pooled standardized null across wells and plates."""
 h=cfg['hit_reference'];B=h['simulations'];rng=np.random.default_rng(np.random.SeedSequence([h['seed'],zlib.crc32(str(plate_id).encode())]))
 m=median_polish(x,cfg['median_polish']['epsilon'],cfg['median_polish']['max_iterations']);z=m['residuals'];negative=roles=='negative';positive=roles=='positive';sample=roles=='sample'
 noise=float(z[negative].std(ddof=1));offset=float(np.median(z[positive]));fit=m['overall']+m['row'][:,None]+m['column'][None,:]+offset*positive
 null=np.empty((B,int(sample.sum())))
 for b in range(B):null[b]=_statistic(fit+rng.normal(0.,noise,x.shape),negative,sign,cfg)[sample]
 mean=null.mean(axis=0);sd=null.std(axis=0,ddof=1);pool=np.sort(((null-mean)/sd).ravel());observed=(_statistic(x,negative,sign,cfg)[sample]-mean)/sd
 p=(1+pool.size-np.searchsorted(pool,observed,side='left'))/(1+pool.size)
 return [float(v) for v in p],{'method':'layout_simulation','simulations':B,'seed':h['seed'],'plate_stream':'SeedSequence([seed, crc32(plate_id)])','noise_sd':noise,'positive_offset':offset,'minimum_attainable_p':float(1/(1+pool.size)),'source':h['source']}


def compute(d,cfg):
 plates=[];audit_tests=[];notes=['Z-prime and SSMD describe control separation; they do not validate a biological hit.','B scores use median-polish residuals divided by the declared scaled MAD.','BH-adjusted predictive t tests are exploratory well-level screens; confirm hits in independent experiments.','Plate correction can induce dependence; reported FDR is conditional on assumptions and recorded calibration.'];fail=[]
 if 'hit_reference' in cfg:notes.append('Hit p-values use a layout-conditional simulated null (Gaussian additive plate model, declared seed and simulation count); they are Monte Carlo p-values, not a validated biological activity claim.')
 for pid,g in d.groupby('plate_id',sort=False):
  g=g.copy();rs=sorted(g.row.unique());cs=sorted(g.column.unique());x=g.pivot(index='row',columns='column',values='value').loc[rs,cs].to_numpy();m=median_polish(x,cfg['median_polish']['epsilon'],cfg['median_polish']['max_iterations']);z=m['residuals'];median=float(np.median(z));mad=float(cfg['mad_scale']*np.median(abs(z-median)));q=quality(g[g.role=='positive'].value,g[g.role=='negative'].value);reasons=[]
  if not m['converged']:reasons.append('median_polish_did_not_converge')
  if mad<=0:reasons.append('zero_residual_MAD')
  if q['zprime'] is None or q['zprime']<cfg['criteria']['zprime_min']:reasons.append('zprime_below_declared_criterion')
  if q['ssmd'] is None or abs(q['ssmd'])<cfg['criteria']['ssmd_abs_min']:reasons.append('ssmd_below_declared_criterion')
  sign=1 if cfg['direction']=='increasing' else -1
  if sign*(q['positive_mean']-q['negative_mean'])<=0:reasons.append('control_direction_conflicts_with_declared_hit_direction')
  residual={(r,c):float(z[i,j]) for i,r in enumerate(rs) for j,c in enumerate(cs)};g['residual']=[residual[r,c] for r,c in zip(g.row,g.column)];g['b_score']=(g.residual-median)/mad if mad>0 else np.nan
  neg=g[g.role=='negative'].residual;sd=float(neg.std(ddof=1));samples=g[g.role=='sample'];pv=[];reference=None
  if sd<=0:reasons.append('zero_negative_control_variance')
  if 'hit_reference' in cfg and sd>0:
   roles=g.pivot(index='row',columns='column',values='role').loc[rs,cs].to_numpy();order={rc:k for k,rc in enumerate([(r,c) for i,r in enumerate(rs) for j,c in enumerate(cs) if roles[i,j]=='sample'])}
   simulated,reference=layout_reference(x,roles,sign,cfg,pid);pv=[simulated[order[r,c]] for r,c in zip(samples.row,samples.column)]
  else:
   for v in samples.residual:pv.append(float(student_t.sf(sign*(v-neg.mean())/(sd*np.sqrt(1+1/len(neg))),len(neg)-1)) if sd>0 else 1.)
  adjusted=multipletests(pv,method='fdr_bh')[1] if pv else [];wells=[]
  for (_,w),p,adj in zip(samples.iterrows(),pv,adjusted):wells.append({'well_id':w.well_id,'compound_id':w.compound_id,'value':float(w.value),'residual':float(w.residual),'b_score':float(w.b_score) if mad>0 else None,'p_value':p if not reasons else None,'adjusted_p':float(adj) if not reasons else None,'hit':bool(adj<=cfg['criteria']['fdr']) if not reasons else None,'reportable':not reasons})
  audit_tests.append({'plate_id':str(pid),'tests':[{'well_id':str(w.well_id),'p_value':float(p),'adjusted_p':float(adj)} for (_,w),p,adj in zip(samples.iterrows(),pv,adjusted)]})
  edge=g.row.isin([rs[0],rs[-1]])|g.column.isin([cs[0],cs[-1]]);edge_diff=float(g.loc[edge,'value'].median()-g.loc[~edge,'value'].median()) if (~edge).any() else None
  row={'plate_id':str(pid),'status':'limited' if reasons else 'estimated','reportable':not reasons,'quality':q,'correction':{'overall':m['overall'],'row_effects':dict(zip(map(str,rs),map(float,m['row']))),'column_effects':dict(zip(map(str,cs),map(float,m['column']))),'iterations':m['iterations'],'converged':m['converged'],'residual_MAD':mad,'raw_edge_minus_interior_median':edge_diff},**({'hit_reference':reference} if reference else {}),'wells':wells,'diagnostics':reasons};plates.append(row)
  if reasons:fail.append({'plate_id':str(pid),'reasons':reasons})
 return {'analysis_type':'hts_qc','primary':{'status':'limited' if fail else 'estimated','plates':plates,'assay_context':{k:cfg[k] for k in ('plate_shape','readout','response_unit','direction','criteria','independence_source','majority_inactive_source','layout_source')}},'audit_tests':audit_tests,'must_mention':notes,'failing_items':fail,'limitations':['Z-prime is a control-only metric; SSMD is the sample plug-in effect size.','Median polish assumes most wells inactive and no systematic biological row/column layout.','One-sided predictive t p-values use negative controls and normal equal-variance assumptions.','A well-level hit is not a validated compound effect; no automatic pooling across plates.']}
