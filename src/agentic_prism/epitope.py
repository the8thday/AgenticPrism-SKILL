"""Directed blocking matrices and conditional panel clustering stability."""
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy.sparse.csgraph import connected_components
from scipy.stats import t as student_t
from .cmc_common import merge,number,text
TYPES={'epitope_binning'}
DEFAULTS={'schema_version':1,'analysis_type':'epitope_binning','input':None,'source':None,
 'format':None,'assay_rationale':None,'valency':None,'readout':None,
 'controls_share_response_scale':False,'independent_experiments':None,'independence_source':None,
 'background_handling':None,'background_source':None,
 'thresholds_predeclared':False,'threshold_source':None,'blocked_min':.8,'nonblocked_max':.2,
 'minimum_control_span':1.,'control_span_source':None,'self_block_min':.8,
 'asymmetry_difference':.3,'graph_rule':'reciprocal_blocked_components',
 'clustering':{'method':'average','distance':'euclidean','cut_height':.5,'cut_source':None,'bootstrap_replicates':1000,'seed':121260930},
 'report':{'plot_style':'prism_like'}}

def resolve_config(raw):
 c=merge(raw,DEFAULTS)
 if c['analysis_type']!='epitope_binning' or c['schema_version']!=1:raise ValueError('Unsupported epitope schema')
 for k in ('input','source','assay_rationale','readout','independence_source','background_source','threshold_source','control_span_source'):text(c[k],k)
 if c['format'] not in ('tandem','premix','sandwich'):raise ValueError('Declare assay format')
 expected={'tandem':'additional_binding_response','sandwich':'second_antibody_binding_response','premix':'available_antigen_binding_response'}
 if c['readout']!=expected[c['format']]:raise ValueError('Readout must match the declared binning format')
 if c['valency'] not in ('monovalent','bivalent','mixed_declared'):raise ValueError('Declare valency')
 if c['controls_share_response_scale'] is not True or c['thresholds_predeclared'] is not True:raise ValueError('Controls must share a justified response scale and thresholds must be prespecified')
 if c['background_handling'] not in ('subtract_declared','already_corrected'):raise ValueError('Declare explicit background handling')
 ids=c['independent_experiments']
 if not isinstance(ids,list) or not ids or len(set(ids))!=len(ids) or any(not isinstance(x,str) or not x.strip() for x in ids):raise ValueError('List unique independent experiment IDs; wells are not experiments')
 for k in ('blocked_min','nonblocked_max','minimum_control_span','self_block_min','asymmetry_difference'):number(c[k],k)
 if not 0<=c['nonblocked_max']<c['blocked_min']<=1 or c['minimum_control_span']<=0 or not 0<c['self_block_min']<=1 or not 0<c['asymmetry_difference']<=1:raise ValueError('Invalid prespecified thresholds')
 if c['graph_rule']!='reciprocal_blocked_components':raise ValueError('Only explicitly reciprocal blocking graph is supported')
 cl=c['clustering'];text(cl['cut_source'],'clustering cut source');number(cl['cut_height'],'cut height',True)
 if cl['method']!='average' or cl['distance']!='euclidean':raise ValueError('Only average linkage and Euclidean directed profiles supported')
 if type(cl['bootstrap_replicates'])is not int or not 100<=cl['bootstrap_replicates']<=10000 or type(cl['seed'])is not int or cl['seed']<0:raise ValueError('Invalid bootstrap count/seed')
 if c['report']['plot_style'] not in ('prism_like','standard'):raise ValueError('Invalid style')
 return c


def load_data(path,cfg):
 d=pd.read_csv(path,dtype={'experiment_id':str,'first':str,'second':str});required={'experiment_id','first','second','response','reference_response','self_response','background'}
 if not required<=set(d) or d.empty:raise ValueError('Require ordered antibodies, experiment and matched controls')
 for k in ('experiment_id','first','second'):
  if d[k].isna().any() or d[k].str.strip().eq('').any():raise ValueError('Missing identity')
 if d.duplicated(['experiment_id','first','second']).any():raise ValueError('One measurement per ordered pair per independent experiment; summarize technical wells explicitly before import')
 if set(d.experiment_id)!=set(cfg['independent_experiments']):raise ValueError('Input must match declared independent experiment IDs')
 for k in required-{'experiment_id','first','second'}:
  d[k]=pd.to_numeric(d[k],errors='raise')
  if not np.isfinite(d[k]).all():raise ValueError('Nonfinite response/control')
 if cfg['background_handling']=='already_corrected' and (d.background!=0).any():raise ValueError('Already-corrected data must have zero background column')
 return d


def branches(z,labels):
 sets={i:frozenset([labels[i]]) for i in range(len(labels))};out=[]
 for i,row in enumerate(z):
  s=sets[int(row[0])]|sets[int(row[1])];sets[len(labels)+i]=s
  if len(s)<len(labels):out.append(s)
 return out


def cluster_profiles(matrix,labels,settings,indices=None):
 x=np.asarray(matrix,float)
 if x.ndim!=2 or len(labels)!=x.shape[0] or not np.isfinite(x).all():raise ValueError('Complete finite directed profiles required')
 if len(labels)<3:raise ValueError('At least three antibodies required for clustering')
 z=linkage(pdist(x),method='average');base=branches(z,labels);rng=np.random.default_rng(settings['seed']);counts={s:0 for s in base};n=settings['bootstrap_replicates']
 if indices is None:indices=rng.integers(0,x.shape[1],size=(n,x.shape[1]))
 for sample in indices:
  found=set(branches(linkage(pdist(x[:,sample]),method='average'),labels))
  for s in counts:counts[s]+=s in found
 n=len(indices)
 return {'linkage':z.tolist(),'labels':labels,'clusters':dict(zip(labels,map(int,fcluster(z,settings['cut_height'],criterion='distance')))),
  'branches':[{'members':sorted(s),'bootstrap_probability':counts[s]/n,'mcse':float(np.sqrt(counts[s]/n*(1-counts[s]/n)/n)),'replicates':n} for s in base],
  'resampling_unit':'profile feature (panel antibody), not independent experiment','approximately_unbiased_support':None}


def compute(d,cfg):
 d=d.copy();notes=['Blocking does not establish an identical structural epitope.','Reciprocal directions are retained; asymmetric pairs are never silently averaged.','Panel bootstrap stability is conditional on observed antibody features, not biological confidence.'];fail=[]
 if cfg['background_handling']=='subtract_declared':
  for column in ('response','reference_response','self_response'):d[column]=d[column]-d.background
 d['control_span']=d.reference_response-d.self_response
 d['blocking']=(d.reference_response-d.response)/d.control_span.where(d.control_span>=cfg['minimum_control_span'])
 d['valid']=np.isfinite(d.blocking)
 ids=sorted(set(d['first'])|set(d['second']));exps=cfg['independent_experiments'];pairs=[];matrix=np.full((len(ids),len(ids)),np.nan)
 # Each experiment must demonstrate its own self-block controls for every antibody.
 control_ok={}
 for exp in exps:
  h=d[(d.experiment_id==exp)&(d['first']==d['second'])]
  for ab in ids:
   row=h[h['first']==ab];control_ok[exp,ab]=len(row)==1 and bool(row.valid.iloc[0]) and row.blocking.iloc[0]>=cfg['self_block_min'] and row.reference_response.iloc[0]>0 and row.self_response.iloc[0]/row.reference_response.iloc[0]<=1-cfg['self_block_min']
 for i,a in enumerate(ids):
  for j,b in enumerate(ids):
   h=d[(d['first']==a)&(d['second']==b)];reason=[]
   if set(h.experiment_id)!=set(exps):reason.append('missing_declared_experiment_or_pair')
   if not h.valid.all():reason.append('insufficient_reference_self_control_span')
   if not all(control_ok[e,a] and control_ok[e,b] for e in exps):reason.append('self_block_control_failed_or_missing')
   value=None;interval=None;state='missing' if h.empty else 'withheld'
   if not reason:
    value=float(h.blocking.mean());matrix[i,j]=value;state='blocked' if value>=cfg['blocked_min'] else 'nonblocked' if value<=cfg['nonblocked_max'] else 'ambiguous'
    if len(h)>=2:
     half=float(student_t.ppf(.975,len(h)-1)*h.blocking.std(ddof=1)/np.sqrt(len(h)));interval=[value-half,value+half]
    if value<0 or value>1:notes.append('Normalized blocking outside [0,1] retained without clipping; review controls and readout.')
   else:fail.append({'first':a,'second':b,'reasons':reason})
   pairs.append(dict(first=a,second=b,n_experiments=len(h),blocking=value,interval=interval,status=state,reportable=not reason,reasons=reason))
 index={(p['first'],p['second']):p for p in pairs};asymmetric=[];adj=np.zeros((len(ids),len(ids)))
 for i,a in enumerate(ids):
  for j in range(i+1,len(ids)):
   b=ids[j];ab=index[a,b];ba=index[b,a]
   if ab['blocking'] is not None and ba['blocking'] is not None:
    if ab['status']!=ba['status'] or abs(ab['blocking']-ba['blocking'])>=cfg['asymmetry_difference']:asymmetric.append({'first':a,'second':b,'forward_status':ab['status'],'reverse_status':ba['status'],'difference':ab['blocking']-ba['blocking']})
    if ab['status']==ba['status']=='blocked':adj[i,j]=adj[j,i]=1
 components,labels=connected_components(adj,directed=False);communities=[]
 for n in range(components):
  ix=np.flatnonzero(labels==n);members=[ids[i] for i in ix];clique=bool(np.sum(adj[np.ix_(ix,ix)])==len(ix)*(len(ix)-1));communities.append({'members':members,'all_pairs_reciprocally_blocked':clique,'reportable':not fail,'status':'provisional_incomplete_matrix' if fail else 'observed_graph_component'})
 notes.append('A graph edge requires both directions blocked; use the directed matrix to assess competition, including missing and asymmetric pairs.')
 if asymmetric:notes.append('Asymmetric blocking observed; inspect order, avidity, detection and incomplete saturation before interpretation.')
 if any(not c['all_pairs_reciprocally_blocked'] for c in communities):notes.append('A graph community is not a mutually blocking clique; transitive edges do not establish identical epitopes.')
 clustering=cluster_profiles(matrix,ids,cfg['clustering']) if np.isfinite(matrix).all() and len(ids)>=3 else {'status':'withheld','reason':'Incomplete controlled matrix or fewer than three antibodies; no imputation'}
 if len(exps)<2:notes.append('One independent experiment: no experiment-level interval; feature-bootstrap stability does not replace replication.')
 if fail:notes.append('Missing or failed controls/pairs retained; community membership is provisional and incomplete clustering is withheld.')
 r={'analysis_type':'epitope_binning','primary':{'status':'limited' if fail else 'estimated','pairs':pairs,'communities':communities,'asymmetric_pairs':asymmetric,'assay_context':{k:cfg[k] for k in ('format','valency','readout','assay_rationale','controls_share_response_scale','independent_experiments','independence_source','blocked_min','nonblocked_max','threshold_source','minimum_control_span','control_span_source','self_block_min','asymmetry_difference','graph_rule','background_handling','background_source','clustering')}},'matrix':{'antibodies':ids,'blocking':[[None if not np.isfinite(v) else float(v) for v in row] for row in matrix]},'clustering':clustering,'must_mention':list(dict.fromkeys(notes)),'failing_items':fail,'limitations':['Thresholds and cluster cut declared before analysis.','Equal weights across declared experiments; technical wells are not independent experiments.','Conditional panel bootstrap BP, not pvclust multiscale AU support.','Canonical normalized-response table only; no unverified vendor parser.']}
 return r
