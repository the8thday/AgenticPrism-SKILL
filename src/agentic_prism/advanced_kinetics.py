"""Opt-in surface mechanisms with predeclared models and withheld weak inference."""
from copy import deepcopy
import hashlib
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import f as f_dist
from .binding_ode import response
from .cmc_common import merge,number,text
MODELS=('one_to_one_drift','heterogeneous_ligand','bivalent_analyte','mass_transport','off_rate_screening')
DEFAULTS={'predeclared':False,'mechanistic_rationale':None,'reference_valid':False,
 'response_is_bound_analyte':False,'response_scales_comparable':False,'valency':None,
 'drift_source':None,'surface_conversion_nM_per_response':None,'surface_conversion_source':None,
 'second_rebinding':None,'comparison_rule':'diagnostic_only','aicc_improvement':10.,
 'profile_log_span':1.,'max_parameter_correlation':.995,'max_condition_number':1e10,
 'ode_rtol':1e-7,'ode_atol':1e-9,'profile':True,
 'ka2_bounds':[1e-7,.1],'kd2_bounds':[1e-6,1.],'km_bounds':[1e-5,10.],
 'fraction_bounds':[.01,.99],'drift_bounds':[-1.,1.],'max_nfev':250}


def validate(cfg):
 a=cfg['advanced'];model=cfg['model']
 for k in ('predeclared','reference_valid','response_is_bound_analyte','response_scales_comparable','profile'):
  if a[k] is not True:raise ValueError(k+' must be declared true for advanced kinetics')
 text(a['mechanistic_rationale'],'mechanistic_rationale')
 if a['valency'] not in ('monovalent','bivalent'):raise ValueError('Declare analyte valency')
 if model=='bivalent_analyte':
  if a['valency']!='bivalent' or type(a['second_rebinding']) is not bool:raise ValueError('Bivalent model requires valency and second_rebinding declaration')
 elif model!='off_rate_screening' and a['valency']!='monovalent':raise ValueError('This mechanism assumes monovalent analyte; use a justified bivalent model')
 if cfg['assay']['injection_design']!='multi_cycle':raise ValueError('Advanced models currently require independently regenerated cycles')
 if model=='one_to_one_drift':text(a['drift_source'],'drift_source')
 if model=='mass_transport':number(a['surface_conversion_nM_per_response'],'surface conversion',True);text(a['surface_conversion_source'],'surface conversion source')
 if model!='off_rate_screening' and cfg['fit']['rmax']!='shared':raise ValueError('Advanced global models require a declared common surface response scale and shared Rmax')
 if cfg.get('steady_state',{}).get('enabled'):raise ValueError('Steady-state diagnostic currently requires default reportable 1:1 kinetics')
 if a['comparison_rule'] not in ('diagnostic_only','aicc_and_identifiability'):raise ValueError('Invalid predeclared comparison rule')
 for k in ('aicc_improvement','profile_log_span','max_parameter_correlation','max_condition_number','ode_rtol','ode_atol'):number(a[k],k,True)
 if not .5<a['max_parameter_correlation']<1 or not .5<=a['profile_log_span']<=3:raise ValueError('Invalid identifiability gate')
 for k in ('ka2_bounds','kd2_bounds','km_bounds','fraction_bounds','drift_bounds'):
  b=a[k]
  if not isinstance(b,list) or len(b)!=2 or not np.isfinite(b).all() or not b[0]<b[1]:raise ValueError('Invalid '+k)
  if k!='drift_bounds' and b[0]<=0:raise ValueError('Positive '+k+' required')
 if a['fraction_bounds'][1]>=1:raise ValueError('Fractions must lie inside (0,1)')
 if type(a['max_nfev']) is not int or a['max_nfev']<30:raise ValueError('max_nfev must be >=30')
 if cfg['uncertainty']['method']!='block_bootstrap':raise ValueError('Advanced models require segment-wise block bootstrap')


def setup(curves,cfg,model=None):
 model=model or cfg['model'];a=cfg['advanced'];amp=max(max(c['y']) for c in curves);amp=max(amp,1e-3)
 names=['ka','kd','rmax'];lo=[np.log10(cfg['fit']['kon_bounds'][0]*1e-9),np.log10(cfg['fit']['koff_bounds'][0]),np.log10(amp*.01)];hi=[np.log10(cfg['fit']['kon_bounds'][1]*1e-9),np.log10(cfg['fit']['koff_bounds'][1]),np.log10(amp*100)]
 start=[-3.,-2.,np.log10(amp*1.2)]
 if model in ('heterogeneous_ligand','bivalent_analyte'):
  names+=['ka2','kd2'];lo+=list(np.log10([a['ka2_bounds'][0],a['kd2_bounds'][0]]));hi+=list(np.log10([a['ka2_bounds'][1],a['kd2_bounds'][1]]));start += [-4.,-1.3 if model=='heterogeneous_ligand' else -3.]
 if model=='heterogeneous_ligand':names+=['fraction'];lo+=[a['fraction_bounds'][0]];hi+=[a['fraction_bounds'][1]];start+=[.5]
 if model=='mass_transport':names+=['km'];lo+=[np.log10(a['km_bounds'][0])];hi+=[np.log10(a['km_bounds'][1])];start+=[-2.]
 if model=='one_to_one_drift':
  for i in range(len(curves)):names+=[f'drift_{i}'];lo+=[a['drift_bounds'][0]];hi+=[a['drift_bounds'][1]];start+=[0.]
 if cfg['fit']['offset']=='fitted_per_curve':
  for i in range(len(curves)):names+=[f'offset_{i}'];lo+=[-amp*2];hi+=[amp*2];start+=[0.]
 return names,np.array(lo),np.array(hi),np.clip(start,np.array(lo)+1e-6,np.array(hi)-1e-6)


def unpack(z,names):return {k:float(v if k=='fraction' or k.startswith(('drift_','offset_')) else 10**v) for k,v in zip(names,z)}

def predict(curves,z,names,cfg,model=None,tight=False):
 p=unpack(z,names);a=cfg['advanced'];out=[]
 if (model or cfg['model'])=='heterogeneous_ligand' and p['kd']/p['ka']>=p['kd2']/p['ka2']:raise ValueError('Heterogeneous sites use ordered KD1 < KD2')
 for i,c in enumerate(curves):
  pp={**p,'drift':p.get(f'drift_{i}',0.)}
  y=response(model or cfg['model'],c['t'],c['c']*1e9,c['duration'],pp,gamma=a['surface_conversion_nM_per_response'] or .1,second_rebinding=a['second_rebinding'] is not False,rtol=a['ode_rtol']/(10 if tight else 1),atol=a['ode_atol']/(10 if tight else 1))
  out.append(y+p.get(f'offset_{i}',0.))
 return out


def solve(curves,cfg,ys=None,start=None,model=None,fixed=None):
 names,lo,hi,initial=setup(curves,cfg,model);ys=[c['y'] for c in curves] if ys is None else ys;scale=max(np.std(np.concatenate(ys)),1e-8);fixed=fixed or {};free=[i for i in range(len(names)) if i not in fixed]
 def expand(z):
  zz=np.array(initial);zz[free]=z
  for i,v in fixed.items():zz[i]=v
  return zz
 def fun(z):
  try:return np.concatenate([p-y for p,y in zip(predict(curves,expand(z),names,cfg,model),ys)])/scale
  except (ValueError,ArithmeticError,FloatingPointError):return np.full(sum(map(len,ys)),1e8)
 starts=[np.clip(start,lo+1e-8,hi-1e-8)] if start is not None else [initial]
 if start is None:
  for shift in (-.8,.8):
   z=initial.copy();z[:2]+=shift
   if 'ka2' in names:z[names.index('ka2')]-=shift
   starts.append(np.clip(z,lo+1e-8,hi-1e-8))
 fits=[least_squares(fun,z[free],bounds=(lo[free],hi[free]),max_nfev=cfg['advanced']['max_nfev'],ftol=1e-8,xtol=1e-8,gtol=1e-8,diff_step=1e-4) for z in starts]
 fit=min(fits,key=lambda s:float(s.fun@s.fun));z=expand(fit.x);sse=float(fit.fun@fit.fun)*scale**2
 return {'z':z,'names':names,'lo':lo,'hi':hi,'sse':sse,'jac':fit.jac,'success':bool(fit.success),'predictions':predict(curves,z,names,cfg,model),'parameters':unpack(z,names)}


def _aicc(sse,n,k):return n*np.log(max(sse/n,1e-300))+2*k+2*k*(k+1)/(n-k-1) if n>k+1 else None

def fit_group(group,cfg,*,_selection_only=False):
 from .kinetics_fit import _arrays,_resample_blocks
 curves=_arrays(group);gid=str(group.fit_group_id.iloc[0]);a=cfg['advanced'];u=cfg['uncertainty'];notes=[]
 if len({c['c'] for c in curves})<3:raise ValueError('Advanced kinetic model requires >=3 concentrations spanning association curvature')
 if any(min(np.sum(c['phase']==p) for p in ('association','dissociation'))<max(8,2*u['block_length']) for c in curves):raise ValueError('Too few observations per declared phase for block bootstrap')
 sol=solve(curves,cfg);z=sol['z'];names=sol['names'];n=sum(len(c['y']) for c in curves);k=len(z);df=n-k;ss=sol['sse'];pred=sol['predictions'];rates=[i for i,v in enumerate(names) if v in ('ka','kd','ka2','kd2','km')]
 if not sol['success']:notes.append('optimizer_failed')
 if np.any((z-sol['lo']<1e-4)|(sol['hi']-z<1e-4)):notes.append('numerical_boundary_hit')
 J=sol['jac'];cov=np.linalg.pinv(J.T@J);sd=np.sqrt(np.maximum(np.diag(cov),1e-300));corr=cov/np.outer(sd,sd);maxcorr=float(np.max(abs(corr-np.eye(k))));condition=float(np.linalg.cond(J))
 if maxcorr>a['max_parameter_correlation'] or condition>a['max_condition_number']:notes.append('parameters_correlated_or_ill_conditioned')
 profiles=[];threshold=ss*(1+f_dist.ppf(u['level'],1,df)/df)
 for j in rates:
  endpoints=[]
  for direction in (-1,1):
   value=float(np.clip(z[j]+direction*a['profile_log_span'],sol['lo'][j]+1e-6,sol['hi'][j]-1e-6));p=solve(curves,cfg,start=z,fixed={j:value});endpoints.append({'log10_value':value,'objective_sse':p['sse'],'above_threshold':bool(p['sse']>threshold),'converged':p['success']})
  flat=any(not e['above_threshold'] or not e['converged'] for e in endpoints)
  profiles.append({'parameter':names[j],'endpoints':endpoints,'flat_or_open':flat,'threshold_sse':threshold})
  if flat:notes.append('profile_flat_or_open_'+names[j])
 tight=predict(curves,z,names,cfg,tight=True);tol_error=max(float(np.max(abs(x-y))) for x,y in zip(tight,pred))
 if tol_error>1e-5*max(1.,max(float(max(abs(p))) for p in pred)):notes.append('ode_tolerance_sensitive')
 base=solve(curves,cfg,model='one_to_one');base_aicc=_aicc(base['sse'],n,len(base['z']));aicc=_aicc(ss,n,k);improvement=base_aicc-aicc
 if a['comparison_rule']=='aicc_and_identifiability' and improvement<a['aicc_improvement']:notes.append('predeclared_comparison_rule_not_met')
 curve_diagnostics=[]
 for i,(c,p) in enumerate(zip(curves,pred)):
  acfs=[]
  for phase in ('association','dissociation'):
   mask=c['phase']==phase;res=c['y'][mask]-p[mask];ac=float(np.corrcoef(res[:-1],res[1:])[0,1]) if np.std(res)>1e-12 else 0.;acfs.append(ac if np.isfinite(ac) else 0.)
   dt=np.diff(c['t'][mask])
   if len(dt)>1 and np.std(dt)>.05*np.mean(dt):notes.append('irregular_sampling_block_duration_varies')
  diss=p[c['phase']=='dissociation'];decay=float((diss[0]-diss[-1])/max(abs(diss[0]),1e-12))
  if decay<.02:notes.append('little_dissociation_in_observation_window')
  if max(map(abs,acfs))>.5:notes.append('serially_correlated_residuals_review_model_and_block_length')
  curve_diagnostics.append({'curve_id':c['curve_id'],'association_residual_lag1':acfs[0],'dissociation_residual_lag1':acfs[1],'fractional_decay':decay})
 windows=[]
 for frac in (.75,.5):
  shorter=[]
  for c in curves:
   td=c['t'][c['phase']=='dissociation'];mask=(c['phase']!='dissociation')|(c['t']<=td.min()+frac*np.ptp(td));shorter.append({**c,**{key:c[key][mask] for key in ('t','y','phase','segment','indices','final_dissociation')}})
  sw=solve(shorter,cfg,start=z);ratios={names[j]:sw['parameters'][names[j]]/sol['parameters'][names[j]] for j in rates};bad=not sw['success'] or max(max(v,1/v) for v in ratios.values())>2
  if bad:notes.append('fit_window_sensitive')
  windows.append({'fraction':frac,'ratios':ratios,'passed':not bad})
 rng=np.random.default_rng(np.random.SeedSequence([u['seed'],int(hashlib.sha256(gid.encode()).hexdigest()[:8],16)]));draws=[]
 # Validation-only short circuit: bootstrap can only append rejection notes.
 # It cannot rescue a fit already rejected by a pre-bootstrap reliability gate.
 skip_bootstrap=bool(_selection_only and notes)
 for it in range(0 if skip_bootstrap else u['replicates']):
  ys=[]
  for c,p in zip(curves,pred):
   y=p.copy()
   for seg in np.unique(c['segment']):
    mask=c['segment']==seg;y[mask]+=_resample_blocks(c['y'][mask]-p[mask],u['block_length'],rng)
   ys.append(y)
  try:
   sb=solve(curves,cfg,ys,start=z);accepted=sb['success'] and not np.any((sb['z']-sb['lo']<1e-4)|(sb['hi']-sb['z']<1e-4));draws.append({'iteration':it,'accepted':bool(accepted),'parameters':sb['parameters']})
  except (ValueError,ArithmeticError):draws.append({'iteration':it,'accepted':False,'parameters':{}})
 accepted=[b for b in draws if b['accepted']];intervals={}
 if skip_bootstrap:notes.append('bootstrap_skipped_after_decisive_rejection')
 elif len(accepted)>=.9*u['replicates']:
  for key in names:intervals[key]=list(map(float,np.quantile([b['parameters'][key] for b in accepted],[(1-u['level'])/2,(1+u['level'])/2])))
 else:notes.append('insufficient_valid_bootstrap_refits')
 reportable=not notes;pars=sol['parameters'];audit=deepcopy(pars);audit['ka_M_inv_s_inv']=pars['ka']*1e9
 if cfg['model']=='heterogeneous_ligand':
  audit['ka2_M_inv_s_inv']=pars['ka2']*1e9;audit['surface_KDs_M']=[pars['kd']/pars['ka']*1e-9,pars['kd2']/pars['ka2']*1e-9]
 elif cfg['model']!='bivalent_analyte':audit['kd_M']=pars['kd']/pars['ka']*1e-9
 return {'fit_group_id':gid,'model':cfg['model'],'n_curves':len(curves),'status':'estimated' if reportable else 'limited','reportable':reportable,
  'parameters':audit if reportable else None,'audit_parameters':audit,'intervals':intervals if reportable else {},'audit_intervals':intervals,
  'parameter_units':{'ka':'nM^-1 s^-1','ka_M_inv_s_inv':'M^-1 s^-1','kd':'s^-1','ka2':'nM^-1 s^-1 for heterogeneous; response^-1 s^-1 for bivalent','kd2':'s^-1','km':'s^-1','rmax':'response units','drift':'response/s'},
  'objective_sse':ss,'aicc':aicc,'comparison':{'one_to_one_aicc':base_aicc,'aicc_improvement':improvement,'rule':a['comparison_rule'],'threshold':a['aicc_improvement'],'same_windows':True},
  'identifiability':{'max_parameter_correlation':maxcorr,'jacobian_condition':condition,'profiles':profiles,'ode_tolerance_absolute_difference':tol_error},
  'window_sensitivity':windows,'curve_parameters':curve_diagnostics,'bootstrap_successes':len(accepted),'bootstrap':draws,'diagnostics':list(dict.fromkeys(notes)),
  'mass_transport_conclusion':('supported transport model; km is identified' if reportable else 'withheld: km or other parameters/reliability not supported') if cfg['model']=='mass_transport' else None,
  'predictions':[{'curve_id':c['curve_id'],'time_s':float(t),'observed':float(y),'predicted':float(p),'residual':float(y-p)} for c,prediction in zip(curves,pred) for t,y,p in zip(c['t'],c['y'],prediction)]}


def off_rate(group,cfg):
 """Single-concentration dissociation screening; no affinity conversion."""
 from .kinetics_fit import _resample_blocks
 u=cfg['uncertainty'];results=[];predictions=[]
 for cid,g in group.groupby('curve_id',sort=False):
  if g.concentration_M.nunique()!=1:raise ValueError('Off-rate screening requires one concentration per curve')
  h=g[g.used_in_fit & (g.phase=='dissociation')];t=h.time_s.to_numpy();y=h.response_processed.to_numpy()
  if len(t)<max(8,2*u['block_length']):raise ValueError('Insufficient dissociation observations')
  t=t-t[0];offset=cfg['fit']['offset']=='fitted_per_curve';bounds=np.log10(cfg['fit']['koff_bounds'])
  def fit(t,y,start=-2):
   def project(z):
    e=np.exp(-10**float(z[0])*t);X=np.column_stack([e,np.ones(len(t))]) if offset else e[:,None];beta=np.linalg.lstsq(X,y,rcond=None)[0];return X@beta,beta
   f=least_squares(lambda z:project(z)[0]-y,[np.clip(start,*bounds)],bounds=([bounds[0]],[bounds[1]]));pr,beta=project(f.x);return float(10**f.x[0]),pr,beta,bool(f.success)
  rate,pred,beta,success=fit(t,y);notes=[]
  predictions.extend({'curve_id':str(cid),'time_s':float(tt),'observed':float(yy),'predicted':float(pp),'residual':float(yy-pp)} for tt,yy,pp in zip(t,y,pred))
  if not success or beta[0]<=0:notes.append('nondecaying_or_failed_fit')
  if -np.expm1(-rate*np.ptp(t))<.02:notes.append('little_dissociation_in_observation_window')
  res=y-pred;acf=float(np.corrcoef(res[:-1],res[1:])[0,1]) if np.std(res)>1e-12 else 0
  if abs(acf)>.5:notes.append('serially_correlated_residuals_review_model_and_block_length')
  dt=np.diff(t)
  if np.std(dt)>.05*np.mean(dt):notes.append('irregular_sampling_block_duration_varies')
  windows=[]
  for f in (.75,.5):
   mask=t<=f*t.max();r,_,_,ok=fit(t[mask],y[mask],np.log10(rate));windows.append({'fraction':f,'ratio':r/rate})
   if not ok or max(r/rate,rate/r)>2:notes.append('fit_window_sensitive')
  rng=np.random.default_rng(np.random.SeedSequence([u['seed'],int(hashlib.sha256(str(cid).encode()).hexdigest()[:8],16)]));draws=[]
  for it in range(u['replicates']):
   rr,_,_,ok=fit(t,pred+_resample_blocks(res,u['block_length'],rng),np.log10(rate));draws.append({'iteration':it,'koff_s_inv':rr,'accepted':ok})
  accepted=[b['koff_s_inv'] for b in draws if b['accepted']];ci=list(map(float,np.quantile(accepted,[(1-u['level'])/2,(1+u['level'])/2]))) if len(accepted)>=.9*u['replicates'] else None
  if ci is None:notes.append('insufficient_valid_bootstrap_refits')
  if not bounds[0]+1e-4<np.log10(rate)<bounds[1]-1e-4:notes.append('numerical_boundary_hit')
  results.append({'curve_id':str(cid),'reportable':not notes,'status':'limited' if notes else 'estimated','koff_s_inv':rate if not notes else None,'audit_koff_s_inv':rate,'interval':ci if not notes else None,'audit_interval':ci,'diagnostics':list(dict.fromkeys(notes)),'window_sensitivity':windows,'bootstrap':draws})
 ranked=sorted([r for r in results if r['reportable']],key=lambda x:x['koff_s_inv'])
 for i,r in enumerate(ranked):r['rank_slowest_first']=i+1
 return {'fit_group_id':str(group.fit_group_id.iloc[0]),'model':'off_rate_screening','n_curves':len(results),'status':'estimated' if len(ranked)==len(results) else 'limited','reportable':bool(ranked),'curves':results,'diagnostics':['Dissociation-only apparent koff screening; no kon or KD is estimated.'],'predictions':predictions}
