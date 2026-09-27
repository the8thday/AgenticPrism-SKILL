"""Balanced two-factor crossed random-effects ANOVA, one result per cell.

Y_ij = mu + a_i + b_j + e_ij. Interaction and measurement error are inseparable.
Raw ANOVA estimates can be negative; clipping is labelled, not called REML.
Conservative intervals use Bonferroni simultaneous chi-square MS intervals.
"""
import numpy as np
from scipy import stats


def crossed_anova(values, confidence=.95):
    y = np.asarray(values, dtype=float)
    if y.ndim != 2 or min(y.shape) < 3 or not np.isfinite(y).all():
        raise ValueError("Complete balanced matrix with at least three levels per factor required")
    if not .5 < confidence < 1:
        raise ValueError("confidence must lie between .5 and 1")
    a, b = y.shape
    mean = float(y.mean())
    residual = y-y.mean(1)[:, None]-y.mean(0)[None, :]+mean
    dfs = np.array([a-1, b-1, (a-1)*(b-1)])
    ms = np.array([b*np.sum((y.mean(1)-mean)**2)/(a-1),
                   a*np.sum((y.mean(0)-mean)**2)/(b-1), np.sum(residual**2)/dfs[2]])
    if ms[2] <= np.finfo(float).eps * max(1., float(np.var(y))):
        raise ValueError("Residual variation is zero or numerically degenerate")
    raw = np.array([(ms[0]-ms[2])/b, (ms[1]-ms[2])/a, ms[2]])
    # Three simultaneous expectation-of-MS intervals give conservative component CIs.
    tail = (1-confidence)/6
    lo = dfs*ms/stats.chi2.ppf(1-tail, dfs)
    hi = dfs*ms/stats.chi2.ppf(tail, dfs)
    intervals = [[max(0., (lo[0]-hi[2])/b), max(0., (hi[0]-lo[2])/b)],
                 [max(0., (lo[1]-hi[2])/a), max(0., (hi[1]-lo[2])/a)], [lo[2], hi[2]]]
    return {"method": "balanced_crossed_random_anova", "n_factor_a": a, "n_factor_b": b,
            "mean": mean, "df": dfs.tolist(), "mean_squares": ms.tolist(),
            "raw_variances": dict(zip(("factor_a", "factor_b", "residual"), raw.tolist())),
            "nonnegative_variances": dict(zip(("factor_a", "factor_b", "residual"), np.maximum(raw, 0).tolist())),
            "confidence_level": confidence, "interval_method": "Bonferroni chi-square mean-square bounds; conservative under additive Gaussian random effects",
            "variance_intervals": dict(zip(("factor_a", "factor_b", "residual"), np.asarray(intervals).tolist())),
            "boundary_components": [k for k, v in zip(("factor_a", "factor_b", "residual"), raw) if v < 0],
            "factor_b_f": float(ms[1]/ms[2]), "factor_b_p": float(stats.f.sf(ms[1]/ms[2], dfs[1], dfs[2])),
            "limitations": ["Balanced crossed additive model only; no separable interaction variance.",
                            "Negative raw ANOVA components retained; clipped values are not constrained REML.",
                            "Intervals require Gaussian independent random effects and common residual variance."]}


def _satterthwaite_limits(variance, df, confidence):
    """None means unavailable; only an upper endpoint may be unbounded."""
    with np.errstate(divide='ignore', over='ignore'):
        endpoints=df*variance/stats.chi2.ppf([(1+confidence)/2,(1-confidence)/2],df)
    if not np.isfinite(endpoints[0]):
        return None
    return [float(endpoints[0]), float(endpoints[1]) if np.isfinite(endpoints[1]) else None]


def fit_components(y, x, kernels, names, confidence=.95, cv_reference=None, *, _centered_retry=False):
    """Independent scalar random-effect REML, common independent residual.

    K_j = Z_j Z_j'; residual I is added here. Expected REML information
    I_jk = trace(P K_j P K_k)/2 (Giesbrecht--Burns, as VCA::getGB).
    Satterthwaite df = 2*v^2 / Var(v), including covariance for sums.
    Boundary component CIs are unavailable, never a degenerate [0,0] claim.
    """
    from scipy.linalg import cho_factor, cho_solve
    from scipy.optimize import minimize
    y, x = np.asarray(y, float), np.asarray(x, float)
    n=len(y); p=x.shape[1]
    ks=np.asarray([*kernels,np.eye(n)],float); labels=[*names,'residual']
    if x.shape[0]!=n or n<=p or np.linalg.matrix_rank(x)!=p:
        raise ValueError('Fixed-effect design must have full rank and residual degrees of freedom')
    if not np.isfinite(y).all() or not np.isfinite(x).all() or not np.isfinite(ks).all():
        raise ValueError('Non-finite model data')
    if not .5 < confidence < 1:raise ValueError('confidence must lie between .5 and 1')
    if len(set(labels))!=len(labels):raise ValueError('Random component names must be unique and not residual')
    scale=float(np.var(y-x@np.linalg.lstsq(x,y,rcond=None)[0]))
    if scale<=np.finfo(float).eps*max(1.,float(np.mean(y*y))):raise ValueError('Residual variation is degenerate')
    # Scale response and variances for invariant optimizer tolerances.
    ys=y/np.sqrt(scale)
    def state(v):
        cov=np.einsum('i,ijk->jk',v,ks)
        cf=cho_factor(cov,lower=True,check_finite=False)
        vi=cho_solve(cf,np.eye(n),check_finite=False)
        vx=vi@x; info=x.T@vx; ci=np.linalg.inv(info)
        proj=vi-vx@ci@vx.T; py=proj@ys
        dev=2*np.log(np.diag(cf[0])).sum()+np.linalg.slogdet(info)[1]+ys@py
        grad=np.einsum('ij,kji->k',proj,ks)-np.einsum('i,kij,j->k',py,ks,py)
        return float(dev),grad,proj,ci,ci@vx.T@ys
    def objective(z):
        try:return state(z)[:2]
        except np.linalg.LinAlgError:
            # Reject a numerically non-positive-definite line-search trial;
            # never perturb covariance or accept an invalid factorization.
            return np.inf,np.zeros(len(z))
    m=len(ks);bounds=[(0,None)]*(m-1)+[(1e-10,None)]
    candidates=[]
    for fraction in (.2,.8):
        v=np.full(m,(1-fraction)/max(1,m-1));v[-1]=fraction
        opt=minimize(objective,v,jac=True,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-8,'maxiter':1000,'maxls':50})
        candidates.append(opt)
    opt=min(candidates,key=lambda z:z.fun);v=opt.x
    # Polish interior scores: relative objective convergence can stop before a
    # printed worked example's last decimal is resolved. Reuse the existing
    # Richardson Hessian primitive for this final small parameter-space step.
    from .mixed_inference import _hessian
    if np.all(v>1e-7):
        for _ in range(3):
            dev,g,_,_,_=state(v)
            if np.max(np.abs(g))<1e-10:break
            h=_hessian(lambda z:state(z)[0],v)
            if np.linalg.eigvalsh(h)[0]<=0:break
            candidate=v-np.linalg.solve(h,g)
            if np.any(candidate<=0):break
            if state(candidate)[0]>dev+1e-10:break
            v=candidate
    dev,g,P,bc,beta=state(v)
    kkt=g.copy();kkt[(v<1e-8)&(g>0)]=0
    if np.max(np.abs(kkt))>1e-4:
        # Independent constrained optimizer fallback for the rare active-boundary
        # solution where L-BFGS-B stopped on relative objective change.
        alt=minimize(objective,v,jac=True,bounds=bounds,method='SLSQP',
                     options={'ftol':1e-13,'maxiter':1000})
        if alt.fun<=dev+1e-9:v=alt.x
        active_ix=np.flatnonzero(v>1e-7)
        for _ in range(5):
            dev,g,_,_,_=state(v)
            if np.max(np.abs(g[active_ix]))<1e-8:break
            def restricted(z):
                full=v.copy();full[active_ix]=z;return state(full)[0]
            h=_hessian(restricted,v[active_ix])
            if np.linalg.eigvalsh(h)[0]<=0:break
            candidate=v.copy();candidate[active_ix]-=np.linalg.solve(h,g[active_ix])
            if np.any(candidate[active_ix]<=0) or state(candidate)[0]>dev+1e-10:break
            v=candidate
        dev,g,P,bc,beta=state(v);kkt=g.copy();kkt[(v<1e-8)&(g>0)]=0
    if np.max(np.abs(kkt))>1e-4:
        if not _centered_retry:
            # Translation by a vector in col(X) leaves the restricted likelihood
            # unchanged. Center only a failed fit to avoid cancellation of a
            # large intercept in y' P y; preserve successful legacy numerics.
            offset=np.linalg.lstsq(x,y,rcond=None)[0]
            retry=fit_components(y-x@offset,x,kernels,names,confidence,cv_reference,_centered_retry=True)
            retry['fixed_coefficients']=(np.asarray(retry['fixed_coefficients'])+offset).tolist()
            return retry
        raise ArithmeticError('REML optimization did not satisfy the constrained score tolerance')
    boundary=v<1e-8
    pk=np.einsum('ij,kjl->kil',P,ks)
    info=.5*np.einsum('kij,lji->kl',pk,pk)
    # Check full identifiability even at a constrained boundary.
    if np.linalg.eigvalsh(info)[0]<=1e-10*np.linalg.eigvalsh(info)[-1]:raise ValueError('Variance components are not separately identifiable in this design')
    active=~boundary;vcov=np.zeros((m,m));vcov[np.ix_(active,active)]=np.linalg.inv(info[np.ix_(active,active)])*scale**2
    v=v*scale; beta=beta*np.sqrt(scale);bc=bc*scale
    def summary(weights):
        variance=float(weights@v); var=float(weights@vcov@weights)
        df=2*variance**2/var if var>0 and variance>0 else None
        # Sums containing boundary estimates also need caution, but are estimable.
        ci=None
        if df is not None:
            ci=_satterthwaite_limits(variance,df,confidence)
        sd=float(np.sqrt(variance))
        return {'variance':variance,'sd':sd,'cv_percent':100*sd/cv_reference if cv_reference is not None else None,
                'df':df,'interval_status':'unavailable_boundary' if df is None else 'unavailable_numeric_range' if ci is None else 'unbounded_upper' if ci[1] is None else 'estimated','variance_interval':ci,'sd_interval':[float(np.sqrt(z)) if z is not None else None for z in ci] if ci else None,
                'cv_interval':[100*np.sqrt(z)/cv_reference if z is not None else None for z in ci] if ci and cv_reference is not None else None}
    components={name:summary(np.eye(m)[i]) for i,name in enumerate(labels)}
    return {'method':'REML expected-information Satterthwaite','n':n,'fixed_coefficients':beta.tolist(),'fixed_covariance':bc.tolist(),
            'components':components,'repeatability':components['residual'],'intermediate_precision':summary(np.ones(m)),
            'component_covariance':vcov.tolist(),'component_order':labels,'confidence_level':confidence,'cv_reference':cv_reference,
            'boundary_components':[labels[i] for i in range(m) if boundary[i]],'max_constrained_score':float(np.max(np.abs(kkt))),
            'limitations':['Independent Gaussian random intercepts and common independent residual variance; no random slopes or residual correlation.',
                'Satterthwaite intervals are approximate and marginal, not simultaneous; covariance of components is included in the total.',
                'Intervals with an unrepresentable lower endpoint are unavailable, not [null,null] or a zero lower bound.',
                'Boundary component intervals are unavailable; total intervals condition the information calculation on active components.',
                'CV uses the declared positive reference as fixed; its uncertainty is not propagated.']}
