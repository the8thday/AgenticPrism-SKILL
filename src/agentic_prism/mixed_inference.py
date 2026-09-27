"""Satterthwaite small-sample inference for Gaussian linear mixed models (lmerTest-style).

Variance parameters theta enter V_i(theta) = Z_i G(theta) Z_i' + sigma^2 I. For a
contrast l, df = 2 v^2 / (g' A g), with v = l' C(theta) l, C = (X' V^-1 X)^-1,
g = dv/dtheta and A = 2 H^-1, H the Hessian of the REML deviance (-2 log L_R).
Multi-df tests follow Fai and Cornelius (1996) as implemented in lmerTest:
eigen-decompose L C L', compute one df per component and combine them.
Derivatives are central finite differences on log(theta); at an interior REML
optimum the resulting df do not depend on the parameterization.
"""
import numpy as np
from scipy.stats import f as f_dist, t as t_dist


class RandomInterceptModel:
    """V_i = tau 11' + sigma^2 I for each unit block; theta = (tau, sigma^2)."""

    def __init__(self, y, x, groups):
        self.y, self.x = np.asarray(y, float), np.asarray(x, float)
        self.blocks = [np.flatnonzero(groups == g) for g in dict.fromkeys(groups)]

    def _inverse_blocks(self, theta):
        tau, sigma = theta
        for ix in self.blocks:
            m = len(ix)
            inv = (np.eye(m) - tau / (sigma + m * tau) * np.ones((m, m))) / sigma
            logdet = (m - 1) * np.log(sigma) + np.log(sigma + m * tau)
            yield ix, inv, logdet

    def gls(self, theta):
        p = self.x.shape[1]
        info, rhs, logdet = np.zeros((p, p)), np.zeros(p), 0.
        cache = []
        for ix, inv, ld in self._inverse_blocks(theta):
            xi = self.x[ix]
            info += xi.T @ inv @ xi
            rhs += xi.T @ inv @ self.y[ix]
            logdet += ld
            cache.append((ix, inv))
        cov = np.linalg.inv(info)
        beta = cov @ rhs
        quad = sum(float(r @ inv @ r) for ix, inv in cache for r in [self.y[ix] - self.x[ix] @ beta])
        return beta, cov, logdet, info, quad

    def deviance(self, theta):
        """-2 REML log-likelihood up to a constant."""
        if np.any(np.asarray(theta) <= 0):
            return np.inf
        _, _, logdet, info, quad = self.gls(theta)
        return logdet + np.linalg.slogdet(info)[1] + quad


def _base_steps(z, d=.1, eps=1e-4, zero_tol=np.sqrt(np.finfo(float).eps / 7e-7)):
    """numDeriv-style initial steps: d*|z| (eps where |z| is essentially zero)."""
    z = np.asarray(z, float)
    return np.where(np.abs(z) < zero_tol, eps, d * np.abs(z))


def _richardson(values, v=2, order=2):
    """Richardson extrapolation of estimates computed with steps h, h/v, h/v^2, ... (error in even powers of h)."""
    a = list(values)
    for m in range(1, len(a)):
        factor = v ** (order * m)
        a = [(a[i + 1] * factor - a[i]) / (factor - 1) for i in range(len(a) - 1)]
    return a[0]


def _gradient(fun, z, r=4, v=2):
    z = np.asarray(z, float)
    h0 = _base_steps(z)
    grad = np.empty(len(z))
    for j in range(len(z)):
        estimates, h = [], h0[j]
        for _ in range(r):
            e = np.zeros(len(z)); e[j] = h
            estimates.append((fun(z + e) - fun(z - e)) / (2 * h))
            h /= v
        grad[j] = _richardson(estimates, v)
    return grad


def _hessian(fun, z, r=4, v=2):
    """Richardson-extrapolated central-difference Hessian (as numDeriv::hessian, used by lmerTest)."""
    z = np.asarray(z, float)
    k, h0, f0 = len(z), _base_steps(z), fun(z)
    hess = np.empty((k, k))
    for i in range(k):
        for j in range(i, k):
            estimates, hi, hj = [], h0[i], h0[j]
            for _ in range(r):
                ei, ej = np.zeros(k), np.zeros(k)
                ei[i], ej[j] = hi, hj
                if i == j:
                    estimates.append((fun(z + ei) - 2 * f0 + fun(z - ei)) / (hi * hi))
                else:
                    estimates.append((fun(z + ei + ej) - fun(z + ei - ej) - fun(z - ei + ej) + fun(z - ei - ej)) / (4 * hi * hj))
                hi, hj = hi / v, hj / v
            hess[i, j] = hess[j, i] = _richardson(estimates, v)
    return hess


class RandomSlopeModel:
    """V_i = Z_i G Z_i' + sigma^2 I with an unstructured q x q G (q = number of random-effect columns).

    Unconstrained parameters u: the lower-triangular Cholesky factor of G with log
    diagonal, followed by log sigma^2. At an interior REML optimum Satterthwaite df
    do not depend on this parameterization.
    """

    def __init__(self, y, x, groups, z):
        self.y, self.x, self.z = np.asarray(y, float), np.asarray(x, float), np.asarray(z, float)
        # Centre y when the constant lies in the column space of X: the REML deviance is unchanged,
        # and the Woodbury sums below avoid cancellation (beta is shifted back in gls()).
        coef, *_ = np.linalg.lstsq(self.x, np.ones(len(self.y)), rcond=None)
        self.offset, self.shift = 0., np.zeros(self.x.shape[1])
        if np.allclose(self.x @ coef, 1, atol=1e-10):
            self.offset, self.shift = float(self.y.mean()), coef
            self.y = self.y - self.offset
        self.blocks = [np.flatnonzero(groups == g) for g in dict.fromkeys(groups)]
        self.q = self.z.shape[1]
        self.tril = np.tril_indices(self.q)
        # Per-unit sufficient statistics; with Woodbury the deviance costs O(units x q^3).
        self.ztz = np.stack([self.z[ix].T @ self.z[ix] for ix in self.blocks])
        self.ztx = np.stack([self.z[ix].T @ self.x[ix] for ix in self.blocks])
        self.zty = np.stack([self.z[ix].T @ self.y[ix] for ix in self.blocks])
        self.xtx, self.xty, self.yty = self.x.T @ self.x, self.x.T @ self.y, float(self.y @ self.y)
        self.sizes = np.array([len(ix) for ix in self.blocks], float)

    def unpack(self, u):
        lower = np.zeros((self.q, self.q))
        lower[self.tril] = u[:len(self.tril[0])]
        lower[np.diag_indices(self.q)] = np.exp(np.diag(lower))
        return lower @ lower.T, float(np.exp(u[-1]))

    def pack(self, g, sigma):
        lower = np.linalg.cholesky(g)
        lower[np.diag_indices(self.q)] = np.log(np.diag(lower))
        return np.concatenate([lower[self.tril], [np.log(sigma)]])

    def gls(self, u):
        """GLS beta, C = (X'V^-1X)^-1, log|V|, X'V^-1X and r'V^-1r via Woodbury on per-unit statistics.

        V_i^-1 = (I - Z_i M_i^-1 Z_i') / sigma^2 with M_i = sigma^2 G^-1 + Z_i'Z_i, and
        log|V_i| = n_i log sigma^2 + log|M_i| + log|G| - q log sigma^2.
        """
        g, sigma = self.unpack(u)
        m = sigma * np.linalg.inv(g)[None] + self.ztz
        sign, logdet_m = np.linalg.slogdet(m)
        sign_g, logdet_g = np.linalg.slogdet(g)
        if np.any(sign <= 0) or sign_g <= 0:
            raise np.linalg.LinAlgError("V not positive definite")
        mx = np.linalg.solve(m, self.ztx)
        my = np.linalg.solve(m, self.zty[..., None])[..., 0]
        info = (self.xtx - np.einsum("kqp,kqr->pr", self.ztx, mx)) / sigma
        rhs = (self.xty - np.einsum("kqp,kq->p", self.ztx, my)) / sigma
        yvy = (self.yty - np.einsum("kq,kq->", self.zty, my)) / sigma
        logdet = float(np.sum(self.sizes * np.log(sigma) + logdet_m + logdet_g - self.q * np.log(sigma)))
        cov = np.linalg.inv(info)
        beta = cov @ rhs
        quad = float(yvy - beta @ rhs)
        return beta + self.offset * self.shift, cov, logdet, info, quad

    def deviance(self, u):
        try:
            _, _, logdet, info, quad = self.gls(u)
        except np.linalg.LinAlgError:
            return np.inf
        return logdet + np.linalg.slogdet(info)[1] + quad

    def blups(self, u, beta):
        g, sigma = self.unpack(u)
        out = np.zeros((len(self.blocks), self.q))
        for k, ix in enumerate(self.blocks):
            zi = self.z[ix]
            v = zi @ g @ zi.T + sigma * np.eye(len(ix))
            out[k] = g @ zi.T @ np.linalg.solve(v, self.y[ix] + self.offset - self.x[ix] @ beta)
        return out


def fit_reml(model, starts):
    """Minimize the REML deviance over unconstrained parameters from several starts, then polish."""
    from scipy.optimize import minimize
    import warnings
    best = None
    for u0 in starts:
        with warnings.catch_warnings(), np.errstate(invalid="ignore"):
            warnings.simplefilter("ignore", RuntimeWarning)  # deviance is +inf outside the PD region
            res = minimize(model.deviance, u0, method="L-BFGS-B", options={"maxiter": 2000})
        if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
            best = res
    if best is None:
        raise ArithmeticError("REML optimization failed")
    # Newton refinement with Richardson derivatives: quadratic convergence to the interior optimum.
    u, f = best.x.copy(), best.fun
    for _ in range(25):
        grad, hess = _gradient(model.deviance, u), _hessian(model.deviance, u)
        if not (np.all(np.isfinite(grad)) and np.all(np.isfinite(hess))) or np.linalg.eigvalsh(hess).min() <= 0:
            break
        step = np.linalg.solve(hess, grad)
        for _ in range(30):
            trial = u - step
            f_trial = model.deviance(trial)
            if f_trial <= f:
                break
            step /= 2
        else:
            break
        converged = f - f_trial < 1e-12 and np.max(np.abs(step)) < 1e-9
        u, f = trial, f_trial
        if converged:
            break
    else:
        pass
    if f < best.fun or np.allclose(u, best.x):
        best.x, best.fun = u, f
        return best
    polished = minimize(model.deviance, best.x, method="Nelder-Mead",
                        options={"xatol": 1e-11, "fatol": 1e-13, "maxiter": 20000, "maxfev": 40000})
    return polished if polished.fun <= best.fun else best


def satterthwaite_unconstrained(deviance, covariance, u0, beta, contrasts, joint=None, level=.95, n_family=None):
    """Satterthwaite t/F given the REML deviance and C(u) as functions of unconstrained parameters u."""
    u0 = np.asarray(u0, float)
    hess = _hessian(deviance, u0)
    eig = np.linalg.eigvalsh(hess)
    if not np.all(np.isfinite(hess)) or eig.min() <= 0:
        raise ArithmeticError("REML deviance Hessian is not positive definite")
    a = 2 * np.linalg.inv(hess)
    cov = covariance(u0)

    def df_for(l):
        v = float(l @ cov @ l)
        g = _gradient(lambda u: float(l @ covariance(u) @ l), u0)
        return v, 2 * v * v / float(g @ a @ g)

    rows = []
    m = n_family or len(contrasts)
    for l in contrasts:
        v, df = df_for(np.asarray(l, float))
        se, estimate = np.sqrt(v), float(np.asarray(l) @ beta)
        stat = estimate / se
        margin = float(t_dist.ppf(1 - (1 - level) / (2 * m), df) * se)
        rows.append({"estimate": estimate, "standard_error": float(se), "df": float(df), "statistic": float(stat),
                     "p_unadjusted": float(2 * t_dist.sf(abs(stat), df)), "ci_low": estimate - margin,
                     "ci_high": estimate + margin})
    omnibus = None
    if joint is not None:
        L = np.atleast_2d(np.asarray(joint, float))
        vl = L @ cov @ L.T
        values, vectors = np.linalg.eigh(vl)
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
        q = L.shape[0]
        lb = vectors.T @ (L @ beta)
        fstat = float(np.sum(lb ** 2 / values) / q)
        nus = np.array([df_for(row)[1] for row in vectors.T @ L])
        if q == 1:
            ddf = float(nus[0])
        elif np.all(np.abs(np.diff(nus)) < 1e-8):
            ddf = float(nus.mean())
        elif np.any(nus <= 2):
            ddf = 2.
        else:
            e = float(np.sum(nus / (nus - 2)))
            ddf = 2 * e / (e - q)
        omnibus = {"f_statistic": fstat, "df_numerator": q, "df_denominator": ddf,
                   "p_value": float(f_dist.sf(fstat, q, ddf)), "component_df": nus.tolist()}
    return rows, omnibus, {"variance_parameter_covariance_log_scale": a.tolist()}


def satterthwaite(model, theta, beta, contrasts, joint=None, level=.95, n_family=None):
    """Random-intercept convenience wrapper: theta = (tau, sigma^2) on the natural scale, derivatives in log theta."""
    return satterthwaite_unconstrained(lambda z: model.deviance(np.exp(z)), lambda z: model.gls(np.exp(z))[1],
                                       np.log(np.asarray(theta, float)), beta, contrasts, joint, level, n_family)
