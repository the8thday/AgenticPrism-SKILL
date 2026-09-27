"""Time-to-event estimation: Kaplan-Meier, log-rank, Cox proportional hazards and the PH test.

Conventions follow R survival 3.x so results can be checked against survfit,
survdiff, coxph (ties = "efron") and cox.zph (transform = "km"):
- Kaplan-Meier with Greenwood standard error of the cumulative hazard,
  confidence limits on the log or log-log scale; medians and their limits from
  the first time the (limit) curve reaches 0.5, with R's midpoint rule when the
  curve equals 0.5 exactly over an interval.
- Log-rank (rho = 0): chi-square on k-1 df from the hypergeometric variance.
  Optional permutation p value (arm labels permuted, Monte Carlo, fixed seed).
- Cox: Efron partial likelihood by Newton-Raphson from beta = 0; Wald, likelihood
  ratio and score tests; hazard ratios with Wald limits.
- PH test: score test for X * g(t) added to the model at (beta_hat, 0), with
  g = 1 - KM(t-) centred on the event times (cox.zph, transform "km").
"""
import numpy as np
from scipy.stats import chi2, norm


def kaplan_meier(time, event, level=.95, conf_type="log-log"):
    time, event = np.asarray(time, float), np.asarray(event, int)
    grid = np.unique(time)
    n_risk = np.array([(time >= t).sum() for t in grid])
    n_event = np.array([((time == t) & (event == 1)).sum() for t in grid])
    n_censor = np.array([((time == t) & (event == 0)).sum() for t in grid])
    surv = np.cumprod(1 - n_event / n_risk)
    with np.errstate(divide="ignore", invalid="ignore"):
        increments = np.where(n_event > 0, n_event / (n_risk * (n_risk - n_event)), 0.)
    se = np.sqrt(np.cumsum(increments))  # standard error of the cumulative hazard (R std.err)
    z = norm.ppf((1 + level) / 2)
    lower, upper = np.full(len(grid), np.nan), np.full(len(grid), np.nan)
    ok = (surv > 0) & np.isfinite(se)
    if conf_type == "log":
        lower[ok] = np.exp(np.log(surv[ok]) - z * se[ok])
        upper[ok] = np.minimum(np.exp(np.log(surv[ok]) + z * se[ok]), 1)
    elif conf_type == "log-log":
        inner = ok & (surv < 1)
        theta = np.log(-np.log(surv[inner]))
        s_theta = se[inner] / np.log(surv[inner])
        lower[inner] = np.exp(-np.exp(theta - z * s_theta))
        upper[inner] = np.exp(-np.exp(theta + z * s_theta))
        lower[ok & (surv == 1)] = upper[ok & (surv == 1)] = 1
    else:
        raise ValueError("conf_type must be log or log-log")
    return {"time": grid, "n_risk": n_risk, "n_event": n_event, "n_censor": n_censor, "surv": surv,
            "std_err": se, "lower": lower, "upper": upper}


def _first_below(times, values, p=.5, tol=np.sqrt(np.finfo(float).eps)):
    """R quantile.survfit rule: first time the curve is <= p; midpoint if it sits exactly at p."""
    values = np.asarray(values, float)
    hit = np.flatnonzero(np.nan_to_num(values, nan=np.inf) <= p + tol)
    if not len(hit):
        return None
    i = hit[0]
    if abs(values[i] - p) < tol:
        later = np.flatnonzero(np.nan_to_num(values, nan=np.inf) < p - tol)
        if len(later):
            return float((times[i] + times[later[0]]) / 2)
        return float(times[i])
    return float(times[i])


def median_survival(km):
    """Median and its limits from the KM curve and its pointwise limits; None means not reached."""
    events = km["n_event"] > 0
    t = km["time"][events]
    return (_first_below(t, km["surv"][events]), _first_below(t, km["lower"][events]),
            _first_below(t, km["upper"][events]))


def survival_at(km, landmark):
    """Step-function value (and limits) at a landmark time."""
    idx = np.searchsorted(km["time"], landmark, side="right") - 1
    if idx < 0:
        return 1., 1., 1.
    return float(km["surv"][idx]), float(km["lower"][idx]), float(km["upper"][idx])


def logrank(time, event, group, n_groups):
    time, event, group = np.asarray(time, float), np.asarray(event, int), np.asarray(group, int)
    observed, expected = np.zeros(n_groups), np.zeros(n_groups)
    variance = np.zeros((n_groups, n_groups))
    for t in np.unique(time[event == 1]):
        at_risk = time >= t
        n = at_risk.sum()
        d = ((time == t) & (event == 1)).sum()
        n_k = np.bincount(group[at_risk], minlength=n_groups).astype(float)
        d_k = np.bincount(group[(time == t) & (event == 1)], minlength=n_groups).astype(float)
        observed += d_k
        expected += d * n_k / n
        if n > 1:
            variance += d * (n - d) / (n - 1) * (np.diag(n_k / n) - np.outer(n_k, n_k) / n ** 2)
    diff = (observed - expected)[:-1]
    reduced = variance[:-1, :-1]
    # Generalized inverse, as needed when some arm never shares a risk set with an event (R survdiff
    # then also loses degrees of freedom). With no information at all the statistic is undefined.
    rank = int(np.linalg.matrix_rank(reduced, tol=1e-10 * max(1., float(np.abs(reduced).max(initial=0.)))))
    if rank == 0:
        return {"chisq": 0., "df": 0, "p_value": 1., "observed": observed.tolist(), "expected": expected.tolist(),
                "computable": False}
    stat = float(diff @ np.linalg.pinv(reduced) @ diff)
    return {"chisq": stat, "df": rank, "p_value": float(chi2.sf(stat, rank)),
            "observed": observed.tolist(), "expected": expected.tolist(), "computable": rank == n_groups - 1}


def permutation_logrank(time, event, group, n_groups, draws, seed):
    """Monte Carlo permutation p value of the log-rank chi-square (arm labels exchangeable under H0)."""
    observed = logrank(time, event, group, n_groups)["chisq"]
    rng = np.random.default_rng(seed)
    exceed = sum(logrank(time, event, rng.permutation(group), n_groups)["chisq"] >= observed - 1e-12 for _ in range(draws))
    return float((1 + exceed) / (draws + 1))


def _efron_terms(beta, time, event, x):
    """Log partial likelihood, score and information (Efron ties) at beta."""
    eta = x @ beta
    risk = np.exp(eta - eta.max())
    loglik, score, info = 0., np.zeros(x.shape[1]), np.zeros((x.shape[1], x.shape[1]))
    for t in np.unique(time[event == 1]):
        at_risk = time >= t
        dead = (time == t) & (event == 1)
        d = dead.sum()
        s0, s1 = risk[at_risk].sum(), risk[at_risk] @ x[at_risk]
        s2 = (x[at_risk] * risk[at_risk, None]).T @ x[at_risk]
        e0, e1 = risk[dead].sum(), risk[dead] @ x[dead]
        e2 = (x[dead] * risk[dead, None]).T @ x[dead]
        loglik += (eta[dead] - eta.max()).sum()
        score += x[dead].sum(0)
        for r in range(d):
            f = r / d
            a0, a1, a2 = s0 - f * e0, s1 - f * e1, s2 - f * e2
            loglik -= np.log(a0)
            mean = a1 / a0
            score -= mean
            info += a2 / a0 - np.outer(mean, mean)
    return loglik, score, info


def cox(time, event, x, max_iter=30, tol=1e-9):
    time, event, x = np.asarray(time, float), np.asarray(event, int), np.asarray(x, float)
    x = x - x.mean(0)  # centring, as in coxph; coefficients are unchanged
    beta = np.zeros(x.shape[1])
    ll0, score0, info0 = _efron_terms(beta, time, event, x)
    ll = ll0
    converged = False
    for _ in range(max_iter):
        _, score, info = _efron_terms(beta, time, event, x)
        step = np.linalg.solve(info, score)
        candidate = beta + step
        ll_new = _efron_terms(candidate, time, event, x)[0]
        halving = 0
        while ll_new < ll - 1e-12 and halving < 30:
            step /= 2
            candidate = beta + step
            ll_new = _efron_terms(candidate, time, event, x)[0]
            halving += 1
        beta = candidate
        if abs(ll_new - ll) <= tol * max(abs(ll), 1):
            ll = ll_new
            converged = True
            break
        ll = ll_new
    ll, score, info = _efron_terms(beta, time, event, x)
    cov = np.linalg.inv(info)
    return {"beta": beta, "covariance": cov, "loglik": [float(ll0), float(ll)], "converged": converged,
            "lr": float(2 * (ll - ll0)), "wald": float(beta @ info @ beta),
            "score": float(score0 @ np.linalg.solve(info0, score0)), "df": x.shape[1],
            "linear_predictor": x @ beta, "centred_x": x}


def cox_zph(time, event, fit, terms=None):
    """cox.zph(fit, transform = "km"): per-term and global score tests of beta(t) = beta + gamma g(t)."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    x, eta = fit["centred_x"], fit["linear_predictor"]
    km = kaplan_meier(time, event, conf_type="log")
    idx = np.searchsorted(km["time"], time, side="left")  # KM times strictly before t
    ttimes = 1 - np.concatenate([[1.], km["surv"]])[idx]
    g = ttimes - ttimes[event == 1].mean()
    p = x.shape[1]
    risk = np.exp(eta - eta.max())
    u = np.zeros(p)
    imat = np.zeros((2 * p, 2 * p))
    for t in np.unique(time[event == 1]):
        at_risk = time >= t
        dead = (time == t) & (event == 1)
        d = dead.sum()
        gt = g[dead][0]
        s0, s1 = risk[at_risk].sum(), risk[at_risk] @ x[at_risk]
        s2 = (x[at_risk] * risk[at_risk, None]).T @ x[at_risk]
        e0, e1 = risk[dead].sum(), risk[dead] @ x[dead]
        e2 = (x[dead] * risk[dead, None]).T @ x[dead]
        u += gt * x[dead].sum(0)
        for r in range(d):
            f = r / d
            a0, a1, a2 = s0 - f * e0, s1 - f * e1, s2 - f * e2
            mean = a1 / a0
            v = a2 / a0 - np.outer(mean, mean)
            u -= gt * mean
            imat[:p, :p] += v
            imat[:p, p:] += gt * v
            imat[p:, :p] += gt * v
            imat[p:, p:] += gt * gt * v
    terms = terms or [[j] for j in range(p)]
    out = []
    for jj in terms + [list(range(p))]:
        kk = list(range(p)) + [p + j for j in jj]
        uu = np.concatenate([np.zeros(p), u[jj]])
        stat = float(uu @ np.linalg.solve(imat[np.ix_(kk, kk)], uu))
        out.append({"chisq": stat, "df": len(jj), "p_value": float(chi2.sf(stat, len(jj)))})
    return out
