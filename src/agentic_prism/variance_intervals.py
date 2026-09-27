"""Mean-square MLS and correlated quadratic-form MOVER precision intervals.

Graybill & Wang (1980); Burdick & Graybill (1984); Zou et al. (2009).
The unbalanced extension uses REML plug-in Gaussian quadratic-form moments,
not a claim that unbalanced sums of squares are independent chi-squares.
"""
import numpy as np
from scipy import stats, linalg, optimize


def mls_limits(estimates, dfs, weights, confidence=.95, correlation=None, side='two-sided', shortest=False):
    s, df, w = map(lambda a: np.asarray(a, float), (estimates, dfs, weights))
    if s.ndim != 1 or df.shape != s.shape or w.shape != s.shape or not np.isfinite([s, df, w]).all() or np.any(s < 0) or np.any(df <= 0):
        raise ValueError('Finite mean squares, positive df and matching weights required')
    if not .5 < confidence < 1 or side not in ('two-sided','upper'):raise ValueError('Invalid interval request')
    alpha=1-confidence
    tail=np.full(len(s),alpha/2 if side=='two-sided' else alpha)
    if shortest:
        if side!='two-sided':raise ValueError('Shortest allocation applies to two-sided intervals')
        for i,d in enumerate(df):
            tail[i]=optimize.minimize_scalar(lambda t: 1/stats.chi2.ppf(t,d)-1/stats.chi2.ppf(1-alpha+t,d),bounds=(1e-12,alpha-1e-12),method='bounded',options={'xatol':1e-14}).x
    lo=s*df/stats.chi2.ppf(1-alpha+tail if side=='two-sided' else confidence,df)
    hi=s*df/stats.chi2.ppf(tail,df)
    # For a negative coefficient, exchange the component's lower/upper limit.
    left=w*np.where(w>=0,s-lo,hi-s);right=w*np.where(w>=0,hi-s,s-lo)
    corr=np.eye(len(s)) if correlation is None else np.asarray(correlation,float)
    if corr.shape!=(len(s),len(s)) or not np.isfinite(corr).all():raise ValueError('Invalid mean-square correlation')
    center=float(w@s)
    lower=center-np.sqrt(max(0.,float(left@corr@left)))
    upper=center+np.sqrt(max(0.,float(right@corr@right)))
    return {'center':center,'variance_interval':[0. if side=='upper' else max(0.,float(lower)),max(0.,float(upper))],
            'confidence_level':confidence,'side':side,'mean_square_weights':w.tolist()}


def quadratic_design(x,kernels):
    """Ordered orthogonal ANOVA strata; no response-dependent stratum selection.

    Random terms are taken in declared order. E(MS)=M v. In a balanced
    orthogonal design A_i V A_j=0 and each stratum has scalar covariance;
    otherwise moments and correlations are evaluated at the REML estimate.
    """
    x=np.asarray(x,float);n=len(x);q=linalg.orth(x);projectors=[];dfs=[]
    for k in kernels:
        # K is a membership kernel; its column span equals the random design.
        z=linalg.orth(np.asarray(k,float));res=z-q@(q.T@z)
        u,sv,_=linalg.svd(res,full_matrices=False)
        add=u[:,sv>1e-8]
        if not add.shape[1]:raise ValueError('Ordered mean-square strata are confounded; use Satterthwaite')
        projectors.append(add@add.T);dfs.append(add.shape[1]);q=np.c_[q,add]
    p=np.eye(n)-q@q.T;df=n-q.shape[1]
    if df<=0:raise ValueError('No residual mean-square degrees of freedom')
    projectors.append(p);dfs.append(df)
    a=np.asarray(projectors)/np.asarray(dfs)[:,None,None]
    ks=np.asarray([*kernels,np.eye(n)])
    ems=np.einsum('ijk,lkj->il',a,ks)
    if np.linalg.cond(ems)>1e10:raise ValueError('Mean-square component mapping is not identifiable')
    # Balance/orthogonality is a property of every kernel, not fitted variances.
    exact=True
    for p,d in zip(projectors,dfs):
        for k in ks:
            coeff=np.trace(p@k)/d
            if not np.allclose(p@k,coeff*p,rtol=1e-8,atol=1e-8):exact=False
    return {'a':a,'dfs':np.asarray(dfs,float),'ems':ems,'inverse':np.linalg.inv(ems),'exact':exact,'kernels':ks}


def precision_intervals(y,design,variances,confidence=.95,sums=None):
    a=design['a'];v=np.asarray(variances);y=np.asarray(y,float)
    ms=np.einsum('i,kij,j->k',y,a,y)
    # Center first: A annihilates X, but subtraction avoids catastrophic cancellation.
    ms=np.maximum(ms,0.)
    cov=np.einsum('i,ijk->jk',v,design['kernels']);av=a@cov
    moments=2*np.einsum('kij,lji->kl',av,av)
    means=design['ems']@v
    dfs=design['dfs'] if design['exact'] else 2*means**2/np.diag(moments)
    corr=np.eye(len(v)) if design['exact'] else moments/np.sqrt(np.outer(np.diag(moments),np.diag(moments)))
    inv=design['inverse'];out={}
    for label,w in (sums or {'total':np.ones(len(v))}).items():
        coefficients=np.asarray(w)@inv
        out[label]=mls_limits(ms,dfs,coefficients,confidence,corr)
        out[label]['method']=('MLS independent mean squares' if np.all(coefficients>=-1e-12) else 'MOVER signed independent mean squares') if design['exact'] else 'MOVER correlated moment-matched mean squares (REML plug-in)'
    upper=[mls_limits(ms,dfs,inv[i],confidence,corr,side='upper') for i in range(len(v))]
    return {'sums':out,'component_upper_bounds':upper,'mean_squares':ms.tolist(),'df':dfs.tolist(),
            'mean_square_correlation':corr.tolist(),'expected_mean_square_coefficients':design['ems'].tolist(),
            'orthogonal_exact_mean_squares':design['exact'],
            'limitation':'MLS/MOVER coverage is approximate. Unbalanced moment matching and correlation use fitted variances; stratum order is declared, not selected by fit.'}
