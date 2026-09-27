"""Categorical-visit marginal repeated-measures models: Python REML/Satterthwaite.

Kenward-Roger uses optional R mmrm through the shared fixed-script bridge.
No subject random intercept is added to the marginal residual covariance.
"""
from copy import deepcopy
import numpy as np
from scipy.optimize import minimize
from . import nonparametric
from .groups import read_unit_table
from .mixed_inference import RandomInterceptModel, fit_reml, satterthwaite_unconstrained, _gradient
DEFAULTS={"analysis_type":"mmrm","schema_version":1,"input":None,"source":"User-supplied repeated measurements",
    "comparison":{"groups":[],"arms":[],"control_arm":None,"outcome":None,"unit":None,"rationale":None,
        "independent_units":None,"covariance_appropriate":None,"covariance_rationale":None,
        "covariates":[],"covariance":"unstructured","inference":"satterthwaite","missing_policy":"require_complete",
        "missingness_rationale":None,"confidence_level":.95,"post_hoc":"arm_vs_control_each_condition"},
    "report":{"plot_style":"prism_like"}}


def resolve_config(raw):
    if not isinstance(raw,dict) or set(raw)-set(DEFAULTS):raise ValueError('Unsupported MMRM config')
    c=deepcopy(DEFAULTS)
    for k,v in raw.items():
        if isinstance(c[k],dict):
            if not isinstance(v,dict) or set(v)-set(c[k]):raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:c[k]=v
    q=c['comparison']
    if c['analysis_type']!='mmrm' or type(c['schema_version']) is not int or c['schema_version']!=1:raise ValueError('Unsupported MMRM schema')
    for v in [c['input'],c['source']]+[q[k] for k in ('outcome','unit','rationale','covariance_rationale')]:
        if not isinstance(v,str) or not v.strip():raise ValueError('Declare input, source, outcome, unit and design/covariance rationale')
    for k in ('groups','arms'):
        v=q[k]
        if not isinstance(v,list) or len(v)<2 or any(not isinstance(s,str) or not s.strip() for s in v) or len(set(v))!=len(v):raise ValueError('Declare distinct ordered visits and arms (at least two each)')
    if any(q[k] is not True for k in ('independent_units','covariance_appropriate')):raise ValueError('MMRM applicability flags must be literal true')
    if not isinstance(q['covariates'],list) or any(not isinstance(v,str) or not v.strip() for v in q['covariates']) or len(set(q['covariates']))!=len(q['covariates']):raise ValueError('covariates must be distinct numeric baseline column names')
    if q['covariance'] not in ('unstructured','ar1') or q['inference'] not in ('satterthwaite','kenward_roger'):raise ValueError('Unsupported covariance/inference')
    if q['missing_policy'] not in ('require_complete','available_case'):raise ValueError('Unsupported missing policy')
    if q['missing_policy']=='available_case' and (not isinstance(q['missingness_rationale'],str) or not q['missingness_rationale'].strip()):raise ValueError('Available-case inference requires a MAR rationale')
    if q['post_hoc'] not in ('none','arm_vs_control_each_condition'):raise ValueError('Unsupported family')
    if (q['post_hoc']=='none' and q['control_arm'] is not None) or (q['post_hoc']!='none' and q['control_arm'] not in q['arms']):raise ValueError('Declare control_arm only for a vs-control family')
    if type(q['confidence_level']) not in (int,float) or not .5<q['confidence_level']<1:raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard','prism_like'):raise ValueError('Unsupported style')
    return c


def load_data(path,cfg):
    d=read_unit_table(path);q=cfg['comparison'];u=d[~d.exclude]
    if 'arm' not in d or set(d.arm)!=set(q['arms']) or set(d.group)!=set(q['groups']) or set(d.outcome)!={q['outcome']} or set(d.unit)!={q['unit']}:raise ValueError('MMRM columns differ from declaration')
    if d.groupby('independent_unit_id').arm.nunique().gt(1).any():raise ValueError('Each unit belongs to one arm')
    if u.duplicated(['independent_unit_id','group']).any():raise ValueError('Duplicate unit/visit')
    if q['missing_policy']=='require_complete' and not u.groupby('independent_unit_id').size().eq(len(q['groups'])).all():raise ValueError('Incomplete units; no silent deletion')
    counts=u.groupby(['arm','group']).size().reindex(__import__('pandas').MultiIndex.from_product([q['arms'],q['groups']])).fillna(0)
    if counts.min()<6:raise ValueError('At least six observed units per arm/visit required')
    presence=u.assign(present=1).pivot(index='independent_unit_id',columns='group',values='present').fillna(0).to_numpy()
    if np.min(presence.T@presence)<3:raise ValueError('At least three observed units for every covariance pair required')
    for col in q['covariates']:
        if col not in d or col in ('value','arm','group','observation_id','independent_unit_id','outcome','unit','exclude','exclusion_reason'):raise ValueError('Invalid baseline covariate column')
        d[col]=__import__('pandas').to_numeric(d[col],errors='raise')
        if not np.isfinite(d[col]).all() or d.groupby('independent_unit_id')[col].nunique().gt(1).any():raise ValueError('Covariates must be finite and constant within a unit')
    y,x,*_=matrices(d,q)
    if np.linalg.matrix_rank(x)<x.shape[1]:raise ValueError('MMRM design matrix is rank deficient')
    return d


class MarginalModel(RandomInterceptModel):
    def __init__(self,y,x,ids,visits,k,structure):
        super().__init__(y,x,ids)
        self.visits,self.k,self.structure=np.asarray(visits,int),k,structure
        self.tril=np.tril_indices(k)
    def unpack(self,u):
        if self.structure=='ar1':
            rho=np.tanh(u[1]);return np.exp(u[0])*rho**abs(np.arange(self.k)[:,None]-np.arange(self.k))
        lower=np.zeros((self.k,self.k));lower[self.tril]=u
        lower[np.diag_indices(self.k)]=np.exp(np.diag(lower))
        return lower@lower.T
    def start(self):
        scale=max(float(np.var(self.y)),1e-6)
        if self.structure=='ar1':return np.array([np.log(scale),np.arctanh(.3)])
        lower=np.zeros((self.k,self.k));lower[np.diag_indices(self.k)]=np.log(np.sqrt(scale));return lower[self.tril]
    def _inverse_blocks(self,u):
        sigma=self.unpack(u)
        for ix in self.blocks:
            v=sigma[np.ix_(self.visits[ix],self.visits[ix])]
            chol=np.linalg.cholesky(v)
            yield ix,np.linalg.inv(v),2*np.log(np.diag(chol)).sum()
    def deviance(self,u):
        try:
            _,_,ld,info,quad=self.gls(u)
            result=ld+np.linalg.slogdet(info)[1]+quad
            return float(result) if np.isfinite(result) else np.inf
        except (np.linalg.LinAlgError,ValueError,OverflowError):return np.inf


def matrices(d,q):
    u=d[~d.exclude];k=len(q['groups']);a=len(q['arms'])
    vi=u.group.map({g:i for i,g in enumerate(q['groups'])}).to_numpy(int)
    ai=u.arm.map({g:i for i,g in enumerate(q['arms'])}).to_numpy(int)
    x=np.eye(a*k)[ai*k+vi]
    eye=np.eye(a*k);joint=np.array([eye[j*k+t]-eye[j*k]-eye[t]+eye[0] for j in range(1,a) for t in range(1,k)])
    pairs=[];contrasts=[]
    if q['post_hoc']!='none':
        control=q['arms'].index(q['control_arm'])
        for j in range(a):
            if j==control:continue
            for t in range(k):
                contrasts.append(eye[j*k+t]-eye[control*k+t]);pairs.append({'arm':q['arms'][j],'control_arm':q['control_arm'],'condition':q['groups'][t]})
    contrasts=np.array(contrasts).reshape(-1,a*k)
    if q['covariates']:
        x=np.column_stack([x,u[q['covariates']].to_numpy(float)])
        joint=np.pad(joint,((0,0),(0,len(q['covariates']))));contrasts=np.pad(contrasts,((0,0),(0,len(q['covariates']))))
    return u.value.to_numpy(float),x,u.independent_unit_id.to_numpy(str),vi,joint,contrasts,pairs


def fit_python(y,x,ids,vi,k,structure,contrasts,joint,level=.95):
    model=MarginalModel(y,x,ids,vi,k,structure)
    fit=fit_reml(model,[model.start()]);beta,cov,*_=model.gls(fit.x)
    if np.max(abs(_gradient(model.deviance,fit.x)))>1e-3:raise ArithmeticError('MMRM REML did not converge to an interior optimum')
    rows,test,_=satterthwaite_unconstrained(model.deviance,lambda u:model.gls(u)[1],fit.x,beta,contrasts,joint,level)
    return {'beta':beta.tolist(),'beta_covariance':cov.tolist(),'covariance_matrix':model.unpack(fit.x).tolist(),
        'tests':test,'contrasts':rows,'deviance_without_constant':float(fit.fun)}


def compare(d,cfg):
    q=cfg['comparison'];y,x,ids,vi,joint,contrasts,pairs=matrices(d,q)
    diagnostics=['available_case_requires_MAR_not_testable_from_observed_data'] if q['missing_policy']=='available_case' else []
    diagnostics+=['covariance_and_contrast_family_were_predeclared','AR1_uses_ordered_visit_lag_not_elapsed_time'] if q['covariance']=='ar1' else ['unstructured_covariance_requires_sufficient_information_for_every_visit_pair']
    environment=None
    try:
        if q['inference']=='kenward_roger':
            from .r_bridge import call
            answer=call('mmrm',{'y':y.tolist(),'x':x.tolist(),'ids':ids.tolist(),'visits':vi.tolist(),'k':len(q['groups']),
                'covariance':q['covariance'],'contrasts':contrasts.tolist(),'joint':joint.tolist(),'level':q['confidence_level']})
            environment=answer.pop('r_environment')
        else:answer=fit_python(y,x,ids,vi,len(q['groups']),q['covariance'],contrasts,joint,q['confidence_level'])
        if not np.all(np.isfinite(answer['beta'])) or not np.isfinite(answer['tests']['p_value']):raise ArithmeticError('Nonfinite fit')
        if any(not np.isfinite([r['standard_error'],r['df'],r['p_unadjusted'],r['ci_low'],r['ci_high']]).all() or r['df']<=0 or r['standard_error']<=0 for r in answer['contrasts']):raise ArithmeticError('Invalid contrast inference')
        diagnostics += answer.get('warnings',[])
        fit={'status':'estimated','reportable':True,'diagnostics':diagnostics,**answer['tests']}
        rows=answer['contrasts'];adjusted=nonparametric.adjust([r['p_unadjusted'] for r in rows],'holm')
        for row,pair,p in zip(rows,pairs,adjusted):row.update(**pair,p_adjusted=float(p),reportable=True)
    except (ArithmeticError,np.linalg.LinAlgError) as exc:
        environment=getattr(exc,'r_environment',environment)
        answer={};rows=[];fit={'status':'withheld','reportable':False,'diagnostics':diagnostics+[str(exc)],'p_value':None}
    fit.update(n_total=len(y),n_units=len(set(ids)),n_groups=len(q['groups']),n_arms=len(q['arms']),
               covariance=q['covariance'],inference=q['inference'])
    return {'schema_version':1,'analysis_type':'mmrm','fits':[fit],'contrasts':rows,
        'estimates':{k:v for k,v in answer.items() if k not in ('tests','contrasts')},'r_environment':environment}
