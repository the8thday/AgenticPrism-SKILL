"""Independent-unit factorial inference, categorical outcomes and paired methods.

All inference runs in Python. Error-ratio convention is Var(error Y)/Var(error X).
"""
from itertools import combinations, permutations
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.stats.contingency import odds_ratio
from statsmodels.stats.proportion import proportion_confint, confint_proportions_2indep
from statsmodels.stats.multitest import multipletests
from .cmc_common import merge, common, text, number, table

TYPES=('factorial_anova','contingency','correlation_regression','method_comparison')
DEFAULTS={'schema_version':1,'analysis_type':None,'input':None,'source':None,
 'assay':{'independent_units':False,'unit_definition':None,'rationale':None,
          'predeclared':False,'same_scale':False,'linear_relation_supported':False},
 'method':None,'confidence_level':.95,
 'factorial':{'ss_type':None,'ss_reason':None,'interaction_alpha':.05,'factor':None,
              'scope':None,'adjustment':None,'control':None},
 'categorical':{'correction':False,'seed':20260930,'resamples':9999,'scores':None,'null_proportion':.5},
 'regression':{'weighting':'unweighted','weight_basis':None},
 'comparison':{'error_variance_ratio':None,'ratio_source':None,'acceptance_limits':None,'acceptance_source':None},
 'report':{'plot_style':'prism_like'}}
METHODS={'factorial_anova':('two_way',),'contingency':('fisher','chi_square','trend','mcnemar','mcnemar_exact','wilson','clopper_pearson'),
 'correlation_regression':('pearson','spearman','linear'), 'method_comparison':('deming','passing_bablok','bland_altman')}

def resolve_config(raw):
 c=merge(raw,DEFAULTS);k=c['analysis_type']
 if k not in TYPES: raise ValueError('Unsupported routine analysis_type')
 common(c,k)
 for f in ('independent_units','predeclared'):
  if c['assay'][f] is not True: raise ValueError(f+' must be literally true; wells are not subjects')
 for f in ('unit_definition','rationale'): text(c['assay'][f],f)
 if c['method'] not in METHODS[k]: raise ValueError('Declare a supported method')
 number(c['confidence_level'],'confidence_level');
 if not .5<c['confidence_level']<1: raise ValueError('Invalid confidence level')
 if k=='factorial_anova':
  f=c['factorial']
  if f['ss_type'] not in ('II','III'): raise ValueError('Declare SS type II or III')
  text(f['ss_reason'],'ss_reason')
  if f['factor'] not in ('factor_a','factor_b') or f['scope'] not in ('simple','marginal'): raise ValueError('Declare contrast factor and simple/marginal scope')
  if f['adjustment'] not in ('tukey','dunnett','sidak','holm'): raise ValueError('Declare contrast family')
  number(f['interaction_alpha'],'interaction_alpha')
  if not 0<f['interaction_alpha']<1: raise ValueError('Invalid interaction threshold')
  if f['adjustment']=='dunnett': text(f['control'],'control')
 if k=='contingency':
  q=c['categorical']
  if type(q['correction']) is not bool or type(q['seed']) is not int or type(q['resamples']) is not int or q['resamples']<999: raise ValueError('Declare valid correction/seed/resamples')
  number(q['null_proportion'],'null_proportion')
  if not 0<q['null_proportion']<1: raise ValueError('Invalid null proportion')
  if c['method']=='trend' and (not isinstance(q['scores'],dict) or not q['scores']): raise ValueError('Declare ordered group scores')
 if k=='correlation_regression':
  if c['regression']['weighting'] not in ('unweighted','inverse_variance'): raise ValueError('Invalid weighting')
  if c['regression']['weighting']!='unweighted':
   if c['method']!='linear': raise ValueError('Only linear regression supports weights')
   text(c['regression']['weight_basis'],'weight_basis')
 if k=='method_comparison':
  if c['assay']['same_scale'] is not True: raise ValueError('Method comparison requires same measurement scale')
  if c['method'] in ('deming','passing_bablok') and c['assay']['linear_relation_supported'] is not True: raise ValueError('Declare plausible linear relation before fitting')
  q=c['comparison']
  if c['method']=='deming':
   number(q['error_variance_ratio'],'Var(error Y)/Var(error X)',True);text(q['ratio_source'],'ratio_source')
  if q['acceptance_limits'] is not None:
   v=q['acceptance_limits']
   if not isinstance(v,list) or len(v)!=2: raise ValueError('Two acceptance limits required')
   for x in v:number(x,'acceptance limit')
   if v[0]>=v[1]:raise ValueError('Ordered limits required')
   text(q['acceptance_source'],'acceptance_source')
 return c

def load_data(path,c):
 k=c['analysis_type'];m=c['method']; required=['observation_id','independent_unit_id'];numeric=[]
 if k=='factorial_anova':required+=['factor_a','factor_b','value'];numeric=['value']
 elif k=='contingency':
  required+=['outcome']
  if m in ('mcnemar','mcnemar_exact'):required+=['before']
  elif m not in ('wilson','clopper_pearson'):required+=['group']
 else:required+=['x','y'];numeric=['x','y']
 if c['regression']['weighting']=='inverse_variance':required+=['weight'];numeric+=['weight']
 d=table(path,required,numeric)
 if d.independent_unit_id.duplicated().any():raise ValueError('Repeated independent_unit_id; use repeated-measures or aggregate with provenance')
 if k=='factorial_anova':
  counts=d.groupby(['factor_a','factor_b']).size()
  if min(d.factor_a.nunique(),d.factor_b.nunique())<2 or len(counts)!=d.factor_a.nunique()*d.factor_b.nunique() or min(counts)<2:raise ValueError('Every crossed cell needs at least two independent units')
 elif k=='contingency':
  if m!='fisher' and m!='chi_square' and not set(d.outcome)<= {'0','1'}:raise ValueError('Binary outcomes must be 0/1')
  if m.startswith('mcnemar') and not set(d.before)<= {'0','1'}:raise ValueError('Binary before must be 0/1')
 else:
  if len(d)<5 or np.ptp(d.x)==0 or np.ptp(d.y)==0:raise ValueError('At least five nonconstant independent pairs required')
  if 'weight' in d and (d.weight<=0).any():raise ValueError('Positive inverse variance weights required')
 return d

def proportion(k,n,method,level=.95):
 return [float(v) for v in proportion_confint(k,n,alpha=1-level,method={'wilson':'wilson','clopper_pearson':'beta'}[method])]

def deming(x,y,ratio=1.,level=.95):
 x=np.asarray(x);y=np.asarray(y);n=len(x)
 def fit(a,b):
  xx=np.var(a,ddof=1);yy=np.var(b,ddof=1);xy=np.cov(a,b,ddof=1)[0,1]
  if abs(xy)<1e-15:raise ValueError('Deming slope unidentified')
  slope=(yy-ratio*xx+np.sqrt((yy-ratio*xx)**2+4*ratio*xy*xy))/(2*xy)
  return np.array([np.mean(b)-slope*np.mean(a),slope])
 theta=fit(x,y);jk=np.array([fit(np.delete(x,i),np.delete(y,i)) for i in range(n)])
 se=np.sqrt((n-1)/n*np.sum((jk-jk.mean(axis=0))**2,axis=0));q=stats.t.ppf((1+level)/2,n-2)
 return {'coefficients':theta.tolist(),'intervals':np.stack([theta-q*se,theta+q*se],axis=1).tolist(),'standard_errors':se.tolist(),'interval_method':'leave-one-subject-out jackknife t'}

def passing_bablok(x,y,level=.95):
 x=np.asarray(x);y=np.asarray(y);n=len(x)
 if min(x.min(),y.min())<0:raise ValueError('Passing-Bablok requires nonnegative measurements')
 if np.corrcoef(x,y)[0,1]<=0:raise ValueError('Passing-Bablok method comparison requires positive association')
 i,j=np.triu_indices(n,1);dx=x[j]-x[i];dy=y[j]-y[i]
 dx[np.abs(dx)<1e-12*(np.abs(x[j])+np.abs(x[i]))/2]=0
 dy[np.abs(dy)<1e-12*(np.abs(y[j])+np.abs(y[i]))/2]=0
 keep=(dx!=0)|(dy!=0);dx=dx[keep];dy=dy[keep]
 slopes=np.sort(np.divide(dy,dx,out=np.full(len(dx),np.inf),where=dx!=0))
 N=len(slopes);offset=int(np.sum(slopes< -1)+np.sum(slopes<= -1))
 if N<3:raise ValueError('Too few informative pairs')
 def rankmedian(numerator):
  center=(numerator-1)/2
  low,high=int(np.floor(center)),int(np.ceil(center))
  if low<0 or high>=N:return np.nan
  return (slopes[low]+slopes[high])/2
 slope=rankmedian(N+offset)
 C=int(np.floor(stats.norm.ppf((1+level)/2)*np.sqrt(n*(n-1)*(2*n+5)/18)+.5))
 lo=rankmedian(N+offset-C);hi=rankmedian(N+offset+C)
 intercept=float(np.median(y-slope*x));a_lo=float(np.median(y-hi*x));a_hi=float(np.median(y-lo*x))
 if not np.isfinite([slope,lo,hi,a_lo,a_hi]).all():raise ValueError('Passing-Bablok interval unbounded; design insufficient')
 residual=y-intercept-slope*x;pos=residual>0;neg=residual<0
 np_,nn=int(sum(pos)),int(sum(neg));score=np.where(pos,np.sqrt(nn/max(np_,1)),np.where(neg,-np.sqrt(np_/max(nn,1)),0.))
 order=np.argsort(x+slope*y,kind='stable');cusum=float(max(abs(np.cumsum(score[order])))/np.sqrt(np_+nn)) if np_ and nn else 0.
 return {'coefficients':[intercept,float(slope)],'intervals':[[a_lo,a_hi],[float(lo),float(hi)]],
  'interval_method':'Passing-Bablok rank interval; tangent median; mcr rank and relative-tie convention',
  'cusum_statistic':cusum,'cusum_5pct_critical':1.36,'linearity_rejected':cusum>1.36,'informative_pairs':N}

def loa(d,level=.95):
 """Exact normal-theory individual quantile CIs by noncentral-t inversion."""
 d=np.asarray(d);n=len(d);mean=float(np.mean(d));s=float(np.std(d,ddof=1));alpha=1-level;z=stats.norm.ppf(.975)
 limits=[]
 for sign in (-1,1):
  # (mean-q_p)/(s/sqrt(n)) ~ nct(n-1,-z_p*sqrt(n)).
  nc=-sign*z*np.sqrt(n)
  ci=[mean-stats.nct.ppf(1-alpha/2,n-1,nc)*s/np.sqrt(n),mean-stats.nct.ppf(alpha/2,n-1,nc)*s/np.sqrt(n)]
  limits.append({'estimate':mean+sign*z*s,'interval':ci})
 return {'bias':mean,'bias_interval':list(stats.t.interval(level,n-1,loc=mean,scale=s/np.sqrt(n))),
  'limits':limits,'interval_method':'exact normal quantile noncentral-t; individual, not simultaneous'}

def factorial(d,c):
 import statsmodels.formula.api as smf
 from statsmodels.stats.anova import anova_lm
 from patsy import build_design_matrices
 fit=smf.ols('value ~ C(factor_a, Sum)*C(factor_b, Sum)',d).fit()
 if np.linalg.matrix_rank(fit.model.exog)!=fit.model.exog.shape[1]:raise ValueError('Rank deficient design')
 an=anova_lm(fit,typ={'II':2,'III':3}[c['factorial']['ss_type']]);names=['C(factor_a, Sum):C(factor_b, Sum)','C(factor_a, Sum)','C(factor_b, Sum)']
 tests=[{'term':n,'F':float(an.loc[n,'F']),'p_value':float(an.loc[n,'PR(>F)']),'df':float(an.loc[n,'df'])} for n in names]
 q=c['factorial'];fac=q['factor'];other='factor_b' if fac=='factor_a' else 'factor_a';levels=sorted(d[fac].unique());others=sorted(d[other].unique())
 material=tests[0]['p_value']<q['interaction_alpha'];allow=not material or q['scope']=='simple';contrasts=[]
 scopes=others if q['scope']=='simple' else [None]
 for scope in scopes:
  X=[]
  for level in levels:
   grid=pd.DataFrame([{fac:level,other:v} for v in (others if scope is None else [scope])])
   X.append(np.asarray(build_design_matrices([fit.model.data.model_spec],grid)[0]).mean(axis=0))
  pairs=list(combinations(range(len(levels)),2))
  if q['adjustment']=='dunnett':
   if q['control'] not in levels:raise ValueError('Control absent')
   ci=levels.index(q['control']);pairs=[(ci,j) for j in range(len(levels)) if j!=ci]
  L=np.array([X[b]-X[a] for a,b in pairs]);est=L@fit.params;cov=L@fit.cov_params()@L.T;se=np.sqrt(np.diag(cov));t=est/se;df=fit.df_resid;p=2*stats.t.sf(abs(t),df);m=len(p);alpha=1-c['confidence_level']
  adj=q['adjustment'];crit=None
  if adj=='tukey':p_adj=stats.studentized_range.sf(abs(t)*np.sqrt(2),len(levels),df);crit=stats.studentized_range.ppf(1-alpha,len(levels),df)/np.sqrt(2)
  elif adj=='dunnett':
   corr=np.ascontiguousarray(cov/np.outer(se,se),dtype=float)
   def prob(v):return stats.multivariate_t.cdf(np.full(m,v),shape=corr,df=df,lower_limit=np.full(m,-v),maxpts=200000,random_state=np.random.default_rng(731))
   p_adj=np.array([1-prob(float(abs(v))) for v in t]);crit=optimize.brentq(lambda v:prob(v)-(1-alpha),.01,30)
  else:
   p_adj=multipletests(p,method='sidak' if adj=='sidak' else 'holm')[1]
   crit=stats.t.ppf(1-(1-(1-alpha)**(1/m))/2,df) if adj=='sidak' else stats.t.ppf(1-alpha/(2*m),df)
  for idx,(a,b) in enumerate(pairs):
   contrasts.append({'a':levels[a],'b':levels[b],'stratum':scope,'estimate':float(est[idx]) if allow else None,'audit_estimate':float(est[idx]),'se':float(se[idx]),'p_adjusted':float(p_adj[idx]) if allow else None,'interval':[float(est[idx]-crit*se[idx]),float(est[idx]+crit*se[idx])] if allow else None,'reportable':allow})
 return {'tests_interaction_first':tests,'contrasts':contrasts,'df_residual':float(fit.df_resid),'interaction_material':material,'contrast_reportable':allow,
 'family_scope':'separate family per declared stratum','interval_method':'Bonferroni simultaneous for Holm; matching single-step intervals otherwise',
 'residual_sd':float(np.sqrt(fit.mse_resid)),'variance_diagnostic_p':float(stats.levene(*[g.value for _,g in d.groupby(['factor_a','factor_b'])]).pvalue)}

def categorical(d,c):
 m=c['method'];q=c['categorical'];level=c['confidence_level'];notes=[]
 if m in ('wilson','clopper_pearson'):
  k=int((d.outcome=='1').sum());n=len(d)
  return {'n':n,'successes':k,'proportion':k/n,'interval':proportion(k,n,m,level),'p_value':float(stats.binomtest(k,n,q['null_proportion']).pvalue)},notes
 if m.startswith('mcnemar'):
  tab=pd.crosstab(d.before,d.outcome).reindex(index=['0','1'],columns=['0','1'],fill_value=0).to_numpy();b,cc=tab[0,1],tab[1,0];nd=b+cc
  p=stats.binomtest(int(b),int(nd),.5).pvalue if nd and m=='mcnemar_exact' else stats.chi2.sf((max(0,abs(b-cc)-(1 if q['correction'] else 0)))**2/nd,1) if nd else 1.
  return {'table':tab.tolist(),'p_value':float(p),'discordant':int(nd)},notes
 groups=sorted(d.group.unique());outcomes=sorted(d.outcome.unique());tab=pd.crosstab(d.group,d.outcome).reindex(index=groups,columns=outcomes,fill_value=0).to_numpy()
 if min(tab.shape)<2:raise ValueError('At least two occupied rows and columns required')
 chi,p,df,expected=stats.chi2_contingency(tab,correction=q['correction']);r={'table':tab.tolist(),'groups':groups,'outcomes':outcomes,'expected':expected.tolist()}
 if np.min(expected)<5:notes.append('Expected count below 5; chi-square approximation may be unreliable; report exact/Monte Carlo method when predeclared.')
 if m=='chi_square':r.update(statistic=float(chi),p_value=float(p),df=int(df))
 elif m=='fisher':
  if tab.shape==(2,2):
   ft=stats.fisher_exact(tab);oo=odds_ratio(tab,kind='conditional');ci=oo.confidence_interval(confidence_level=level)
   r.update(p_value=float(ft.pvalue),conditional_odds_ratio=float(oo.statistic) if np.isfinite(oo.statistic) else None,odds_ratio_interval=[float(v) if np.isfinite(v) else None for v in ci],odds_ratio_upper_open=not np.isfinite(ci.high))
  else:
   ft=stats.fisher_exact(tab,method=stats.MonteCarloMethod(n_resamples=q['resamples'],rng=np.random.default_rng(q['seed'])))
   r.update(p_value=float(ft.pvalue),method='conditional fixed-margins Monte Carlo Fisher',seed=q['seed'],resamples=q['resamples'],monte_carlo_se=float(np.sqrt(ft.pvalue*(1-ft.pvalue)/(q['resamples']+1))))
 elif m=='trend':
  if set(outcomes)!={'0','1'} or set(q['scores'])!=set(groups):raise ValueError('Trend needs binary outcomes and scores for each group')
  scores=np.array([q['scores'][g] for g in groups],float)
  if not np.isfinite(scores).all() or len(set(scores))!=len(scores):raise ValueError('Distinct finite ordered scores required')
  n=tab.sum(axis=1);k=tab[:,outcomes.index('1')];pp=k.sum()/n.sum();center=np.average(scores,weights=n);z=np.sum((scores-center)*k)/np.sqrt(pp*(1-pp)*np.sum(n*(scores-center)**2));r.update(statistic=float(z*z),p_value=float(stats.chi2.sf(z*z,1)))
 if tab.shape==(2,2) and set(outcomes)=={'0','1'}:
  n=tab.sum(axis=1);k=tab[:,outcomes.index('1')];p0,p1=k/n
  rd=confint_proportions_2indep(int(k[1]),int(n[1]),int(k[0]),int(n[0]),method='newcomb',compare='diff',alpha=1-level)
  # Katz log-RR; zeros are withheld rather than silently corrected.
  rr=None;rrci=None
  if min(k)>0:
   rr=p1/p0;se=np.sqrt(1/k[1]-1/n[1]+1/k[0]-1/n[0]);z=stats.norm.ppf((1+level)/2);rrci=[float(np.exp(np.log(rr)-z*se)),float(np.exp(np.log(rr)+z*se))]
  else:notes.append('Risk ratio interval withheld at zero events; no automatic continuity correction.')
  r.update(risk_difference_group1_minus_group0=float(p1-p0),risk_difference_interval=list(map(float,rd)),risk_ratio=float(rr) if rr is not None else None,risk_ratio_interval=rrci,risk_ratio_method='Katz log normal')
 return r,notes

def regression(d,c):
 x=d.x.to_numpy();y=d.y.to_numpy();n=len(x);level=c['confidence_level'];m=c['method']
 if m=='pearson':
  rr=stats.pearsonr(x,y);z=stats.norm.ppf((1+level)/2)/np.sqrt(n-3);ci=np.tanh(np.arctanh(np.clip(rr.statistic,-1+1e-15,1-1e-15))+np.array([-z,z]))
  return {'correlation':float(rr.statistic),'p_value':float(rr.pvalue),'interval':ci.tolist(),'interval_method':'Fisher z','agreement_claim':False}
 if m=='spearman':
  rr=stats.spearmanr(x,y);ties=len(set(x))<n or len(set(y))<n
  if n<=9 and not ties:
   a=stats.rankdata(x);b=stats.rankdata(y);obs=np.sum((a-b)**2);den=0;ext=0
   for perm in permutations(range(1,n+1)):
    v=1-6*np.sum((a-np.array(perm))**2)/(n*(n*n-1));ext+=abs(v)>=abs(rr.statistic)-1e-12;den+=1
   p=ext/den;method='exact two-sided permutation'
  else:p=float(rr.pvalue);method='R cor.test exact=FALSE t approximation; average ranks for ties'
  return {'correlation':float(rr.statistic),'p_value':p,'p_method':method,'ties':ties,'interval':None,'agreement_claim':False}
 import statsmodels.api as sm
 X=sm.add_constant(x);fit=sm.WLS(y,X,weights=d.weight.to_numpy() if 'weight' in d else np.ones(n)).fit();pred=fit.get_prediction(X,weights=d.weight.to_numpy() if 'weight' in d else None).summary_frame(alpha=1-level)
 return {'coefficients':fit.params.tolist(),'intervals':fit.conf_int(alpha=1-level).tolist(),'r_squared':float(fit.rsquared),'agreement_claim':False,
 'bands':pred.to_dict('records'),'x':x.tolist(),'residuals':fit.resid.tolist(),'normality_p':float(stats.shapiro(fit.resid).pvalue),'residual_df':float(fit.df_resid),'weighting':c['regression']['weighting']}

def compute(d,c):
 k=c['analysis_type'];m=c['method'];notes=['The independent unit is '+c['assay']['unit_definition']+'. No automatic outlier removal.']
 if k=='factorial_anova':
  result=factorial(d,c);notes+=['Lead with interaction; main effects average over the other factor. Homoscedastic Gaussian errors are assumed.']
  if not result['contrast_reportable']:notes+=['Material interaction: marginal contrasts withheld; predeclare simple effects.']
 elif k=='contingency':result,n=categorical(d,c);notes+=n
 elif k=='correlation_regression':result=regression(d,c);notes+=['Correlation and R squared do not establish method agreement.']
 elif m=='deming':result=deming(d.x,d.y,c['comparison']['error_variance_ratio'],c['confidence_level'])
 elif m=='passing_bablok':
  result=passing_bablok(d.x,d.y,c['confidence_level'])
  if result['linearity_rejected']:notes+=['CUSUM rejects linearity: Passing-Bablok coefficients are audit-only.']
 else:
  result=loa(d.y-d.x,c['confidence_level']);trend=stats.linregress((d.x+d.y)/2,d.y-d.x);result['proportional_bias_diagnostic']={'slope':float(trend.slope),'p_value':float(trend.pvalue)};notes+=['Exact LoA intervals assume independent normal differences with constant variance; trend is a diagnostic, not a correction.']
 if k=='method_comparison':
  notes+=['Constant bias is the intercept/bias; proportional bias is slope departure from one. Agreement requires predeclared acceptance limits.']
  result['acceptance_limits']=c['comparison']['acceptance_limits']
  if m!='bland_altman':result['constant_bias_interval_excludes_zero']=not result['intervals'][0][0]<=0<=result['intervals'][0][1];result['proportional_bias_interval_excludes_one']=not result['intervals'][1][0]<=1<=result['intervals'][1][1]
 reportable=not result.get('linearity_rejected',False)
 return {'schema_version':1,'analysis_type':k,'method':m,'primary':{'status':'estimated' if reportable else 'withheld','reportable':reportable,'result':result if reportable else None},'detail':result,'must_mention':notes,'failing_items':[] if reportable else ['nonlinear_Passing_Bablok'],'limitations':['Inference is conditional on the declared design and measurement model.'],'fits':[{'status':'estimated' if reportable else 'limited','reportable':reportable}]}
