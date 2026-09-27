"""Lower confidence limits for upper ADA percentiles, not upper FPR bounds.

Parametric mixed-panel limit: noncentral-t moment approximation with crossed
subject/run ANOVA, effective n = marginal variance / variance of grand mean,
and Satterthwaite df for the positive mean-square total. Exact for iid normal
samples only; mixed-panel confidence is assessed by calibration.
Nonparametric: binomial order-statistic bound on explicitly independent pairs, or
(0.9.3) a two-way "pigeonhole" bootstrap (Owen 2007, Ann. Appl. Stat. 1:2) that
resamples subjects and runs independently and uses every panel cell.
"""
import numpy as np
from scipy import stats
from .variance_components import crossed_anova


def normal_lower(mean, variance, mean_variance, df, probability, confidence):
    if not np.isfinite([mean,variance,mean_variance,df,probability,confidence]).all() or min(variance,mean_variance,df)<=0 or not 0<probability<1 or not .5<confidence<1:raise ValueError('Invalid lower-percentile bound parameters')
    effective_n=variance/mean_variance
    factor=float(stats.nct.ppf(1-confidence,df,stats.norm.ppf(probability)*np.sqrt(effective_n)))
    value=float(mean+np.sqrt(mean_variance)*factor)
    if not np.isfinite(value):raise ValueError('Noncentral-t bound outside numerical range')
    return {'lower_bound':value,'effective_n':float(effective_n),'effective_df':float(df),'noncentral_t_factor':factor,'confidence_level':confidence,'percentile_probability':probability}


def parametric_panel_lower(values, probability, confidence):
    y=np.asarray(values,float);a,b=y.shape;vc=crossed_anova(y);ms=np.array(vc['mean_squares']);dfs=np.array(vc['df'])
    w=np.array([1/b,1/a,1-1/a-1/b]);total=float(w@ms)
    df=total**2/np.sum((w*ms)**2/dfs)
    variances=vc['nonnegative_variances'];mean_var=variances['factor_a']/a+variances['factor_b']/b+variances['residual']/(a*b)
    r=normal_lower(y.mean(),total,mean_var,df,probability,confidence)
    r.update(method='normal mixed-panel noncentral-t moment approximation',mean=float(y.mean()),marginal_variance=total,mean_variance=mean_var,
             point_percentile=float(y.mean()+stats.norm.ppf(probability)*np.sqrt(total)),mean_squares=ms.tolist(),component_variances=variances,
             confidence_direction='P(conditional FPR >= target) >= confidence; approximate for crossed panels',
             limitations=['Marginal over a new independent subject AND random run; no per-realized-run FPR guarantee.',
                          'Effective df and effective n are fitted; mixed-panel noncentral-t is approximate, not an exact iid normal result.',
                          'Negative ANOVA components are clipped only for grand-mean variance; total uses positive mean-square coefficients.'])
    return r


def order_statistic_lower(values, probability, confidence):
    y=np.sort(np.asarray(values,float));n=len(y)
    if y.ndim!=1 or n<3 or not np.isfinite(y).all() or not 0<probability<1 or not .5<confidence<1:raise ValueError('At least three independent finite values and valid probabilities required')
    ranks=np.arange(1,n+1);coverage=stats.binom.sf(ranks-1,n,probability);eligible=ranks[coverage>=confidence]
    if not len(eligible):raise ValueError('No finite order statistic attains declared confidence; collect more independent pairs')
    k=int(eligible[-1]);return {'lower_bound':float(y[k-1]),'order':k,'n_independent_pairs':n,'attained_confidence_continuous':float(coverage[k-1]),
        'confidence_level':confidence,'percentile_probability':probability,'method':'binomial lower order statistic',
        'confidence_direction':'P(conditional FPR >= target) >= confidence',
        'limitations':['Requires prespecified mutually independent subject/run pairs, not pooled repeated cells.',
                      'A few independent runs can give a very low cut point and high mean FPR, especially for the titer tail.',
                      'Attained confidence is exact for continuous iid responses; ties are conservative for the >= positive rule.']}


MIN_EXPECTED_TAIL_SUBJECTS = 3  # fixed before the 0.9.3 calibration was run


def two_way_bootstrap_lower(values, probability, confidence, reps, seed):
    """Percentile-bootstrap lower bound of the pooled empirical percentile.

    Subjects and runs are resampled independently with replacement (pigeonhole
    bootstrap for crossed random factors); each resample's type-7 percentile of
    all cells forms the bootstrap distribution, whose (1-confidence) quantile is
    the bound. Owen (2007) shows the scheme is mildly conservative for crossed
    designs; attained confidence is assessed by calibration, not claimed exact.
    """
    y=np.asarray(values,float)
    if y.ndim!=2 or min(y.shape)<3 or not np.isfinite(y).all() or not 0<probability<1 or not .5<confidence<1:raise ValueError('Complete finite subject x run panel with >=3 levels each and valid probabilities required')
    if type(reps) is not int or not 999<=reps<=100000 or type(seed) is not int or seed<0:raise ValueError('Declare integer bootstrap_reps (999-100000) and a nonnegative integer bootstrap_seed')
    a,b=y.shape;expected=a*(1-probability)
    rng=np.random.default_rng(seed)
    si=rng.integers(0,a,(reps,a));ri=rng.integers(0,b,(reps,b))
    q=np.quantile(y[si[:,:,None],ri[:,None,:]].reshape(reps,a*b),probability,axis=1,method='linear')
    bound=float(np.quantile(q,1-confidence,method='linear'))
    return {'lower_bound':bound,'point_percentile':float(np.quantile(y.ravel(),probability,method='linear')),
        'bootstrap_reps':reps,'bootstrap_seed':seed,'bootstrap_percentile_sd':float(q.std(ddof=1)),
        'n_subjects':a,'n_runs':b,'expected_subjects_beyond_percentile':float(expected),
        'tail_support_sufficient':bool(expected>=MIN_EXPECTED_TAIL_SUBJECTS),
        'confidence_level':confidence,'percentile_probability':probability,
        'method':'two-way (subject x run) pigeonhole percentile bootstrap of the pooled type-7 percentile',
        'confidence_direction':'P(conditional FPR >= target) >= confidence; approximate, conservative in the 0.9.3 calibration',
        'limitations':['Targets a new independent subject in a random run (marginal FPR), not every realized run.',
                       'Bootstrap confidence is approximate; the 0.9.3 calibration covered Gaussian crossed panels only.',
                       f'Requires at least {MIN_EXPECTED_TAIL_SUBJECTS} expected subjects beyond the percentile (n_subjects x FPR); extreme titer tails need far larger panels.',
                       'Monte Carlo variation of the bound is controlled by the declared seed and replicate count; report both.']}
