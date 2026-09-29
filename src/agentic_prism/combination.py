"""Named two-agent reference models with independent-matrix uncertainty."""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import t as student_t
from .cmc_common import merge,text,number
from .dose_schema import DEFAULTS as DOSE_DEFAULTS
from .dose_fit import _evaluate,_best,response as logistic
TYPES={'drug_combination'}
MODEL_ASSUMPTIONS={'Bliss':'Independent fractional inhibition probabilities on a common 0 to 100 scale; uses observed single-agent responses.','Loewe':'Dose equivalence between invertible increasing single-agent curves with compatible maximal and baseline effects.','HSA':'Combination compared with the strongest observed single-agent effect at the same component doses.','ZIP':'No potency interaction relative to fitted single-agent independence; average conditional 4PL predictions in both directions.'}
DEFAULTS={'schema_version':1,'analysis_type':'drug_combination','input':None,'source':None,
 'agent_a':None,'agent_b':None,'dose_unit_a':None,'dose_unit_b':None,'condition_id':None,
 'readout':'percent_inhibition','readout_description':None,'normalization_source':None,'same_response_scale':False,
 'independent_experiments':None,'independence_supported':False,'independence_source':None,
 'models':['Bliss','Loewe','HSA','ZIP'],'loewe':{'compatible_maximal_effects':False,'source':None,'maximum_plateau_difference':10.},
 'fit':{'fixed_bottom':None,'fixed_top':None,'hill_bounds':[.05,8.],'multistart':9,'max_nfev':2000},
 'uncertainty':{'method':'independent_matrix_t','level':.95},'report':{'plot_style':'prism_like'}}

def resolve_config(raw):
 c=merge(raw,DEFAULTS)
 if c['analysis_type']!='drug_combination' or c['schema_version']!=1:raise ValueError('Unsupported combination schema')
 for k in ('input','source','agent_a','agent_b','condition_id','readout_description','normalization_source','independence_source'):text(c[k],k)
 if c['readout']!='percent_inhibition' or c['same_response_scale'] is not True:raise ValueError('Declare same-scale percent inhibition with normalization source; no implicit viability conversion')
 if c['dose_unit_a'] not in ('M','mM','uM','nM','pM') or c['dose_unit_b'] not in ('M','mM','uM','nM','pM'):raise ValueError('Declare molar dose units for each agent')
 ids=c['independent_experiments']
 if not isinstance(ids,list) or not ids or len(set(ids))!=len(ids) or any(not isinstance(x,str) or not x.strip() for x in ids):raise ValueError('Declare matrix experiment IDs; wells are not experiments')
 if type(c['independence_supported'])is not bool:raise ValueError('Declare whether matrix experiments are independent')
 if not isinstance(c['models'],list) or not c['models'] or len(set(c['models']))!=len(c['models']) or set(c['models'])-{'Bliss','Loewe','HSA','ZIP'}:raise ValueError('Select named reference models')
 if 'Loewe' in c['models']:
  if c['loewe']['compatible_maximal_effects'] is not True:raise ValueError('Loewe requires declared compatible maximal effects; choose other models otherwise')
  text(c['loewe']['source'],'Loewe applicability source');number(c['loewe']['maximum_plateau_difference'],'plateau tolerance',True)
 f=c['fit']
 for k in ('fixed_bottom','fixed_top'):
  if f[k] is not None:number(f[k],k)
 if f['fixed_bottom'] is not None and f['fixed_top'] is not None and f['fixed_bottom']>=f['fixed_top']:raise ValueError('Top must exceed bottom')
 b=f['hill_bounds']
 if not isinstance(b,list) or len(b)!=2 or not 0<b[0]<b[1]:raise ValueError('Invalid Hill bounds')
 if type(f['multistart'])is not int or not 1<=f['multistart']<=25 or type(f['max_nfev'])is not int or f['max_nfev']<50:raise ValueError('Invalid fitter settings')
 if c['uncertainty']['method']!='independent_matrix_t' or not .5<c['uncertainty']['level']<1:raise ValueError('Unsupported uncertainty')
 if c['report']['plot_style'] not in ('prism_like','standard'):raise ValueError('Invalid style')
 return c


def load_data(path,cfg):
 d=pd.read_csv(path,dtype={'experiment_id':str,'condition_id':str});req={'experiment_id','condition_id','dose_a','dose_b','response'}
 if not req<=set(d) or d.empty:raise ValueError('Require matrix experiment, condition, two doses and response')
 if set(d.condition_id)!={cfg['condition_id']}:raise ValueError('Do not pool cell lines or conditions in a combination fit')
 if set(d.experiment_id)!=set(cfg['independent_experiments']):raise ValueError('Matrix IDs must match the declaration')
 if d.duplicated(['experiment_id','dose_a','dose_b']).any():raise ValueError('Technical wells must be summarized explicitly; duplicate matrix cell is not independent replication')
 for k in ('dose_a','dose_b','response'):
  d[k]=pd.to_numeric(d[k],errors='raise')
  if not np.isfinite(d[k]).all():raise ValueError('Nonfinite matrix')
 if (d[['dose_a','dose_b']]<0).any().any():raise ValueError('Nonnegative doses required')
 grid=None
 for _,g in d.groupby('experiment_id'):
  aa=sorted(g.dose_a.unique());bb=sorted(g.dose_b.unique());pairs=set(zip(g.dose_a,g.dose_b))
  if aa[0]!=0 or bb[0]!=0 or len(aa)<5 or len(bb)<5 or len(pairs)!=len(aa)*len(bb):raise ValueError('Complete matrix with zero-dose axes and >=4 positive doses per agent required')
  if grid is not None and pairs!=grid:raise ValueError('Independent matrices must use the same prespecified dose grid')
  grid=pairs
 return d


def fit_single(x,y,cfg,fixed_bottom=None,conditional=False):
 """Reuse the existing separable 4PL fitter and bounded linear plateau solver."""
 x=np.asarray(x,float);y=np.asarray(y,float);dc=deepcopy(DOSE_DEFAULTS);dc['fit'].update(cfg['fit']);dc['assay']['direction']='increasing'
 if conditional:dc['fit']['fixed_bottom']=float(fixed_bottom);dc['fit']['fixed_top']=None
 positive=x[x>0];lo,hi=np.log10(positive.min()),np.log10(positive.max());lower=np.array([lo-3,np.log10(dc['fit']['hill_bounds'][0])]);upper=np.array([hi+3,np.log10(dc['fit']['hill_bounds'][1])]);m=int(np.ceil(np.sqrt(dc['fit']['multistart'])));starts=[(a,b) for a in np.linspace(lo,hi,m) for b in np.linspace(np.log10(.5),np.log10(2),m)][:dc['fit']['multistart']]
 def evaluate(z):return _evaluate([(x,y)],[z[0]],[z[1]],[[0]],dc)
 fit=_best(lambda z:evaluate(z)[2]/max(np.ptp(y),1.),starts,lower,upper,dc['fit']['max_nfev'])
 if fit is None:return {'reportable':False,'reason':'4PL optimization failed'},None
 pred,plateaus,res=evaluate(fit.x);bottom,top=plateaus[0];hill=float(10**fit.x[1]);mid=float(10**fit.x[0]);reason=[]
 if np.any(fit.x-lower<1e-4) or np.any(upper-fit.x<1e-4):reason.append('4PL numerical boundary')
 if not positive.min()<=mid<=positive.max():reason.append('4PL midpoint outside tested positive doses')
 if top-bottom<1e-6:reason.append('4PL flat response')
 if np.linalg.matrix_rank(fit.jac)<2:reason.append('4PL shape not identifiable')
 pars={'bottom':float(bottom),'top':float(top),'midpoint':mid,'hill':hill,'reportable':not reason,'reason':'; '.join(reason),'sse':float(res@res)}
 return pars,lambda dose:logistic(dose,bottom,top,float(fit.x[0]),hill,'increasing')


def matrix_scores(g,cfg):
 aa=np.sort(g.dose_a.unique());bb=np.sort(g.dose_b.unique());y=g.pivot(index='dose_a',columns='dose_b',values='response').loc[aa,bb].to_numpy();pa,fa=fit_single(aa,y[:,0],cfg);pb,fb=fit_single(bb,y[0,:],cfg);scores={m:np.full(y.shape,np.nan) for m in cfg['models']};reasons=[];conditional=[]
 if 'Bliss' in scores:scores['Bliss']=y-(y[:,0,None]+y[None,0,:]-y[:,0,None]*y[None,0,:]/100)
 if 'HSA' in scores:scores['HSA']=y-np.maximum(y[:,0,None],y[None,0,:])
 usable=pa['reportable'] and pb['reportable'];available=fa is not None and fb is not None;model_reportable={m:(usable if m in ('Loewe','ZIP') else True) for m in cfg['models']}
 if not usable:reasons.append('Single-agent 4PL unsupported; Loewe and ZIP withheld.')
 if 'Loewe' in scores and available:
  compatible=abs(pa['top']-pb['top'])<=cfg['loewe']['maximum_plateau_difference'] and abs(pa['bottom']-pb['bottom'])<=cfg['loewe']['maximum_plateau_difference']
  if not compatible:
   reasons.append('Fitted single-agent plateaus violate declared Loewe compatibility tolerance.');model_reportable['Loewe']=False
  if min(pa['top'],pb['top'])>max(pa['bottom'],pb['bottom']):
   low=max(pa['bottom'],pb['bottom'])+1e-8;high=min(pa['top'],pb['top'])-1e-8
   def dose(effect,p):return p['midpoint']*((effect-p['bottom'])/(p['top']-effect))**(1/p['hill'])
   for i,a in enumerate(aa[1:],1):
    for j,b in enumerate(bb[1:],1):
     try:ref=brentq(lambda e:a/dose(e,pa)+b/dose(e,pb)-1,low,high);scores['Loewe'][i,j]=y[i,j]-ref
     except (ValueError,ZeroDivisionError,OverflowError):reasons.append('Loewe equal-effect inverse unavailable in common response range.')
 if 'ZIP' in scores and available:
  ref=fa(aa)[:,None]+fb(bb)[None,:]-fa(aa)[:,None]*fb(bb)[None,:]/100;pred_a=np.full_like(y,np.nan);pred_b=np.full_like(y,np.nan)
  for j in range(len(bb)):
   pp,fn=fit_single(aa,y[:,j],cfg,float(y[0,j]),True)
   conditional.append({'varied_agent':'agent_a','fixed_agent_dose':float(bb[j]),'fit':pp})
   if fn is not None:pred_a[:,j]=fn(aa)
   if not pp['reportable']:model_reportable['ZIP']=False
  for i in range(len(aa)):
   pp,fn=fit_single(bb,y[i,:],cfg,float(y[i,0]),True)
   conditional.append({'varied_agent':'agent_b','fixed_agent_dose':float(aa[i]),'fit':pp})
   if fn is not None:pred_b[i,:]=fn(bb)
   if not pp['reportable']:model_reportable['ZIP']=False
  scores['ZIP']=(pred_a+pred_b)/2-ref
  if not model_reportable['ZIP'] or not np.isfinite(scores['ZIP'][1:,1:]).all():reasons.append('Conditional 4PL unsupported in some ZIP slices; complete-matrix ZIP summary withheld.')
 for score in scores.values():score[0,:]=0.;score[:,0]=0.
 return aa,bb,scores,{'agent_a':pa,'agent_b':pb,'conditional_fits':conditional,'model_reportable':model_reportable},list(dict.fromkeys(reasons))


def compute(d,cfg):
 runs=[];notes=['Scores depend on the named reference model; Bliss, Loewe, HSA and ZIP are not interchangeable.','Scores are percentage-point differences on the declared inhibition scale.','ZIP uses conditional fitted slices and fitted single-agent Bliss reference.','Technical wells and repeated rows do not establish independent experiments.'];fail=[]
 for exp,g in d.groupby('experiment_id',sort=False):
  aa,bb,scores,pars,reasons=matrix_scores(g,cfg);runs.append({'experiment_id':str(exp),'scores':scores,'single_agent_fits':pars,'diagnostics':reasons})
  if reasons:fail.append({'experiment_id':str(exp),'reasons':reasons})
 summaries=[];audit_summaries=[];cells=[];level=cfg['uncertainty']['level']
 def summarize(values,allowed=True):
  valid=all(v is not None and np.isfinite(v) for v in values);n=len(values);mean=float(np.mean(values)) if valid else None;ci=None
  if valid and n>=2 and cfg['independence_supported']:
   half=float(student_t.ppf((1+level)/2,n-1)*np.std(values,ddof=1)/np.sqrt(n));ci=[mean-half,mean+half]
  return {'mean':mean if allowed else None,'interval':ci if allowed else None,'audit_mean':mean,'audit_interval':ci,'n_matrices':n,'reportable':bool(valid and allowed),'synergy_claim':bool(allowed and ci and ci[0]>0),'status':'withheld' if not allowed or not valid else 'estimated' if ci else 'descriptive'}
 for model in cfg['models']:
  values=[float(np.mean(r['scores'][model][1:,1:])) if np.isfinite(r['scores'][model][1:,1:]).all() else None for r in runs];allowed=all(r['single_agent_fits']['model_reportable'][model] for r in runs);summary={'model':model,**summarize(values,allowed)};audit_summaries.append(summary);summaries.append({k:v for k,v in summary.items() if not k.startswith('audit_')})
  for i,a in enumerate(aa):
   for j,b in enumerate(bb):cells.append({'model':model,'dose_a':float(a),'dose_b':float(b),**{k:v for k,v in summarize([float(r['scores'][model][i,j]) if np.isfinite(r['scores'][model][i,j]) else None for r in runs],allowed).items() if not k.startswith('audit_')}})
 if len(cfg['models'])<4:notes.append('Only the requested reference models were assessed; their scores do not establish agreement with unassessed models.')
 if len(runs)<2 or not cfg['independence_supported']:notes.append('No supported independent matrix replication: no uncertainty interval or synergy claim, even if scores are positive.')
 signs={np.sign(s['mean']) for s in summaries if s['mean'] is not None}
 if -1 in signs and 1 in signs:notes.append('Reference models disagree in score direction; report all requested models and do not select the favorable one.')
 if (d.response<0).any() or (d.response>100).any():notes.append('Measured inhibition outside [0,100] is retained without clipping; independence interpretations require review of normalization/noise.')
 for r in runs:
  r['scores']={m:[[None if not np.isfinite(v) else float(v) for v in row] for row in arr] for m,arr in r['scores'].items()}
 return {'analysis_type':'drug_combination','primary':{'status':'limited' if fail else 'descriptive' if len(runs)<2 or not cfg['independence_supported'] else 'estimated','summaries':summaries,'weighting':'unweighted 4PL fits; equal weights across independent matrices','model_assumptions':{m:MODEL_ASSUMPTIONS[m] for m in cfg['models']},'cells':cells,'assay_context':{k:cfg[k] for k in ('agent_a','agent_b','dose_unit_a','dose_unit_b','condition_id','readout','readout_description','normalization_source','independent_experiments','independence_supported','independence_source','models','loewe','fit','uncertainty')}},'matrix_results':runs,'audit_summaries':audit_summaries,'must_mention':notes,'failing_items':fail,'limitations':['Equal weights across declared independent matrices; Student t intervals conditional on that design.','No pooling of cell lines or response scales.','Loewe requires compatible effect ranges and invertible increasing single-agent curves.','No uncertainty from an unreplicated matrix.','Combination scores are assay- and dose-grid-specific, not intrinsic drug mechanisms.']}
