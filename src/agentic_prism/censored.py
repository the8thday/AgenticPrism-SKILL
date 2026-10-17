"""Group comparison of values censored at detection or quantification limits (0.14.0).

Rows below a lower limit (or above an upper limit) are kept as censored at that
limit; nothing is substituted. Methods:
- censored_lognormal: maximum likelihood normal regression on natural-log values
  with left/right censoring (Tobit / survreg-type), one mean per group and a
  common SD. Wald intervals for group geometric means and ratios to the reference,
  likelihood-ratio test of equal means.
- peto_peto: Peto-Peto (Peto-Prentice) k-sample rank test for left-censored data,
  computed on flipped values as R survdiff(rho = 1); no effect size.
Descriptive per group: censored fraction, limits, Kaplan-Meier (flipped) quartiles.

Paired designs (0.14.1, opt-in `paired`), following Helsel (2011) and the NADA2
reference code: censored_paired_difference (interval-censored normal ML of the
paired differences, original or log scale; R cen_paired / survreg interval2) and
censored_sign_test (R cen_signtest, with Fong's tie correction).
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats
from .survival import kaplan_meier

TYPES = {'censored_comparison'}
METHODS = ('censored_lognormal', 'peto_peto')
PAIRED_METHODS = ('censored_paired_difference', 'censored_sign_test')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'censored_comparison', 'input': None, 'source': None, 'method': None,
            'design': {'unit': None, 'rationale': None, 'outcome': None, 'outcome_unit': None,
                       'limit_source': None},
            'comparison': {'groups': None, 'reference': None, 'confidence_level': .95},
            'paired': None,
            'report': {'plot_style': 'prism_like'}}
MIN_DETECTED = 3


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported censored-comparison configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    pr = c['paired']
    if c['schema_version'] != 1 or c['analysis_type'] != 'censored_comparison' or c['method'] not in METHODS + PAIRED_METHODS:
        raise ValueError('Declare method censored_lognormal or peto_peto (independent units), or censored_paired_difference or censored_sign_test (paired)')
    if (pr is not None) != (c['method'] in PAIRED_METHODS):
        raise ValueError('Paired methods need the paired block and the paired block needs a paired method (censored_paired_difference or censored_sign_test)')
    if pr is not None:
        if not isinstance(pr, dict) or not {'scale', 'scale_rationale', 'pairing_rationale'} <= set(pr) <= {'scale', 'scale_rationale', 'pairing_rationale', 'inference'}:
            raise ValueError('paired must be {scale, scale_rationale, pairing_rationale[, inference]}')
        pr.setdefault('inference', 'small_sample_t')
        if pr['inference'] not in ('small_sample_t', 'wald_z'): raise ValueError('paired inference must be small_sample_t or wald_z (NADA2 cen_paired)')
        if pr['scale'] not in ('original', 'log'): raise ValueError('paired scale must be original or log')
        if not _text(pr['scale_rationale']) or not _text(pr['pairing_rationale']): raise ValueError('Declare scale_rationale and pairing_rationale')
        if len(c['comparison']['groups'] or []) != 2: raise ValueError('A paired comparison has exactly two groups')
    else:
        c.pop('paired')  # unused 0.14.1 key: resolved configurations of independent designs stay as in 0.14.0
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    q = c['comparison']; g = q['groups']
    if not isinstance(g, list) or len(g) < 2 or len(set(g)) != len(g) or not all(_text(x) for x in g):
        raise ValueError('Declare at least two distinct groups')
    if q['reference'] not in g: raise ValueError('Declare the reference group (one of the groups)')
    lv = q['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'unit_id': str, 'group': str, 'censored': str})
    need = {'unit_id', 'group', 'value', 'censored'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded unit needs an exclusion_reason')
    else: d['exclude'] = False
    d['censored'] = d.censored.str.strip().str.lower()
    if not d.censored.isin(['none', 'left', 'right']).all():
        raise ValueError('censored must be none, left (true value below the stated limit) or right (above the limit)')
    v = pd.to_numeric(d.value, errors='coerce')
    if v.isna().any() or not np.isfinite(v).all():
        raise ValueError('value must be numeric: for censored rows give the limit itself, never text such as "<LLOQ" or a substituted number')
    d['value'] = v
    paired = c.get('paired') is not None
    if d.unit_id.isna().any() or (d.duplicated(['unit_id', 'group']).any() if paired else d.unit_id.duplicated().any()):
        raise ValueError('One row per unit and group in a paired design' if paired else 'One row per independent unit; unit_id must be unique')
    if paired:
        k = d.groupby('unit_id').agg(n=('group', 'size'), ex=('exclude', 'sum'))
        if (k.n != 2).any(): raise ValueError('Every unit must appear once in each of the two groups (incomplete pairs: ' + ', '.join(k.index[k.n != 2][:5]) + '); no silent pair deletion')
        if ((k.ex != 0) & (k.ex != 2)).any(): raise ValueError('Excluding one value of a pair requires excluding both rows of that unit, each with a reason')
        if (d.censored == 'right').any(): raise ValueError('Paired methods handle left-censored values only (as the NADA2 reference code)')
        u = d[~d.exclude]
        if set(u.group) != set(c['comparison']['groups']): raise ValueError('Groups in the data do not match the declared groups')
        if c['paired']['scale'] == 'log' and (u.value <= 0).any(): raise ValueError('The log scale needs positive values and limits; no pseudocounts')
        if (u.value < 0).any(): raise ValueError('Paired censored methods need non-negative values (a censored value spans 0 to its limit)')
        return d
    u = d[~d.exclude]
    if set(u.group) != set(c['comparison']['groups']): raise ValueError('Groups in the data do not match the declared groups')
    if c['method'] == 'censored_lognormal' and (u.value <= 0).any(): raise ValueError('censored_lognormal needs positive values and limits; no pseudocounts')
    if c['method'] == 'peto_peto' and (u.censored == 'right').any():
        raise ValueError('peto_peto handles left-censored values only; data with upper-limit censoring need censored_lognormal')
    return d


def tobit(y, kind, X):
    """MLE of y ~ N(X b, s^2) with left/right censoring at the stated values; returns estimates and covariance.

    Analytic gradient; convergence is judged by the gradient at the solution, not by
    optimizer status flags (0.14.0 post-registration fix: strict numerical-gradient
    tolerances had reported converged optima as failures).
    """
    det, left, right = kind == 'none', kind == 'left', kind == 'right'

    def parts(theta):
        b, ls = theta[:-1], theta[-1]; s = np.exp(ls); z = (y - X @ b)/s
        return b, ls, s, z

    def nll(theta):
        b, ls, s, z = parts(theta)
        return -(stats.norm.logpdf(z[det]).sum() - det.sum()*ls + stats.norm.logcdf(z[left]).sum() + stats.norm.logsf(z[right]).sum())

    def grad(theta):
        b, ls, s, z = parts(theta); dmu = np.zeros_like(z); dls = np.zeros_like(z)
        dmu[det] = z[det]/s; dls[det] = z[det]**2 - 1
        hl = np.exp(stats.norm.logpdf(z[left]) - stats.norm.logcdf(z[left])); dmu[left] = -hl/s; dls[left] = -z[left]*hl
        hr = np.exp(stats.norm.logpdf(z[right]) - stats.norm.logsf(z[right])); dmu[right] = hr/s; dls[right] = z[right]*hr
        return -np.r_[X.T @ dmu, dls.sum()]

    start_b = np.linalg.lstsq(X, y, rcond=None)[0]; resid = y - X @ start_b
    start = np.r_[start_b, np.log(max(resid[det].std(ddof=1) if det.sum() > 1 else 1., 1e-3))]
    fit = optimize.minimize(nll, start, jac=grad, method='BFGS', options={'gtol': 1e-9, 'maxiter': 5000})
    theta = fit.x
    for _ in range(20):  # Newton polish on the numerical Hessian of the analytic gradient
        g = grad(theta)
        if np.max(np.abs(g)) < 1e-9: break
        try: step = np.linalg.solve(_jacobian(grad, theta), g)
        except np.linalg.LinAlgError: break
        trial = theta - step
        if not np.isfinite(nll(trial)) or nll(trial) > nll(theta) + 1e-12: break
        theta = trial
    h = _jacobian(grad, theta); h = (h + h.T)/2
    try: cov = np.linalg.inv(h)
    except np.linalg.LinAlgError: cov = np.full_like(h, np.nan)
    gmax = float(np.max(np.abs(grad(theta))))
    return {'theta': theta, 'cov': cov, 'loglik': -float(nll(theta)), 'converged': bool(np.isfinite(gmax) and gmax < 1e-6), 'max_abs_gradient': gmax}


def _jacobian(g, x):
    n = len(x); J = np.zeros((n, n)); e = np.maximum(np.abs(x), 1.)*1e-6
    for i in range(n):
        up, dn = x.copy(), x.copy(); up[i] += e[i]; dn[i] -= e[i]
        J[:, i] = (g(up) - g(dn))/(2*e[i])
    return J


def _hessian(f, x):
    n = len(x); h = np.zeros((n, n)); e = np.maximum(np.abs(x), 1.)*1e-4
    for i in range(n):
        for j in range(i, n):
            def at(si, sj):
                z = x.copy(); z[i] += si*e[i]; z[j] += sj*e[j]; return f(z)
            v = (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1))/(4*e[i]*e[j])
            h[i, j] = h[j, i] = v
    return h


def peto_peto(value, kind, group, groups):
    """Peto-Peto k-sample test for left-censored data via flipped values (R survdiff, rho = 1)."""
    flip = value.max() + 1. - value; event = (kind == 'none').astype(int)
    order = np.argsort(flip, kind='mergesort'); t, e, g = flip[order], event[order], group[order]
    times = np.unique(t[e == 1]); k = len(groups)
    O = np.zeros(k); E = np.zeros(k); V = np.zeros((k, k)); S = 1.
    for tt in times:
        at_risk = t >= tt; n = at_risk.sum(); dmask = (t == tt) & (e == 1); d = dmask.sum()
        nj = np.array([(at_risk & (g == x)).sum() for x in groups]); dj = np.array([(dmask & (g == x)).sum() for x in groups])
        w = S
        O += w*dj; E += w*d*nj/n
        if n > 1: V += w*w*d*(n-d)/(n*n*(n-1))*(np.diag(nj*n) - np.outer(nj, nj))
        S *= 1 - d/n
    diff = (O - E)[:-1]; Vr = V[:-1, :-1]
    chi = float(diff @ np.linalg.solve(Vr, diff)); df = k - 1
    return {'chi_square': chi, 'df': df, 'p_value': float(stats.chi2.sf(chi, df)), 'weighted_observed': O.tolist(), 'weighted_expected': E.tolist()}


def interval_normal(lo, hi):
    """Intercept-only normal MLE for interval-censored values [lo, hi] (lo == hi exact; +-inf open), as survreg interval2."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    use = np.isfinite(lo) | np.isfinite(hi); lo, hi = lo[use], hi[use]
    ex = np.isfinite(lo) & np.isfinite(hi) & (lo == hi); lft = ~np.isfinite(lo); rgt = ~np.isfinite(hi); itv = ~(ex | lft | rgt)

    def log_p(zl, zh):
        # log(Phi(zh) - Phi(zl)) without underflow in either tail
        up = zl > 0
        a = np.where(up, stats.norm.logsf(zl), stats.norm.logcdf(zh)); b = np.where(up, stats.norm.logsf(zh), stats.norm.logcdf(zl))
        return a + np.log1p(-np.exp(np.minimum(b - a, 0.)))

    def nll(theta):
        mu, ls = theta; s = np.exp(ls)
        v = stats.norm.logpdf((lo[ex]-mu)/s).sum() - ex.sum()*ls + stats.norm.logcdf((hi[lft]-mu)/s).sum() + stats.norm.logsf((lo[rgt]-mu)/s).sum()
        return -float(v + log_p((lo[itv]-mu)/s, (hi[itv]-mu)/s).sum())

    def grad(theta):
        mu, ls = theta; s = np.exp(ls); gm = gs = 0.
        z = (lo[ex]-mu)/s; gm += (z/s).sum(); gs += (z**2 - 1).sum()
        z = (hi[lft]-mu)/s; h = np.exp(stats.norm.logpdf(z) - stats.norm.logcdf(z)); gm += (-h/s).sum(); gs += (-z*h).sum()
        z = (lo[rgt]-mu)/s; h = np.exp(stats.norm.logpdf(z) - stats.norm.logsf(z)); gm += (h/s).sum(); gs += (z*h).sum()
        zl, zh = (lo[itv]-mu)/s, (hi[itv]-mu)/s; lp = log_p(zl, zh)
        fl, fh = np.exp(stats.norm.logpdf(zl) - lp), np.exp(stats.norm.logpdf(zh) - lp)
        gm += ((fl - fh)/s).sum(); gs += (zl*fl - zh*fh).sum()
        return -np.array([gm, gs])

    # Start from the midpoints of all bounded differences (finite end for one-sided ones), not only the exact ones.
    e = np.where(np.isfinite(lo) & np.isfinite(hi), (lo + hi)/2, np.where(np.isfinite(lo), lo, hi))
    start = np.array([e.mean(), np.log(max(e.std(ddof=1), 1e-3))])
    fit = optimize.minimize(nll, start, jac=grad, method='BFGS', options={'gtol': 1e-9, 'maxiter': 5000})
    theta = fit.x
    for _ in range(50):
        g = grad(theta)
        if np.max(np.abs(g)) < 1e-9: break
        try: step = np.linalg.solve(_jacobian(grad, theta), g)
        except np.linalg.LinAlgError: break
        trial = theta - step
        if not np.isfinite(nll(trial)) or nll(trial) > nll(theta) + 1e-12: break
        theta = trial
    h = _jacobian(grad, theta); h = (h + h.T)/2
    try: cov = np.linalg.inv(h)
    except np.linalg.LinAlgError: cov = np.full_like(h, np.nan)
    gmax = float(np.max(np.abs(grad(theta))))
    return {'theta': theta, 'cov': cov, 'loglik': -nll(theta), 'converged': bool(np.isfinite(gmax) and gmax < 1e-6), 'max_abs_gradient': gmax,
            'n_used': int(use.sum()), 'n_exact': int(ex.sum()), 'n_left': int(lft.sum()), 'n_right': int(rgt.sum()), 'n_interval': int(itv.sum())}


def paired_intervals(xv, xc, yv, yc, scale):
    """Interval for x - y per pair as NADA2 cen_paired (censored values span [0, limit]); log scale uses (-inf, ln limit]."""
    xv, yv = np.asarray(xv, float), np.asarray(yv, float); xc, yc = np.asarray(xc, bool), np.asarray(yc, bool)
    if scale == 'log':
        with np.errstate(divide='ignore'):
            xl, yl = np.log(xv), np.log(yv)
        xlo = np.where(xc, -np.inf, xl); ylo = np.where(yc, -np.inf, yl)
        lo = xlo - yl; hi = xl - ylo
        lo = np.where(np.isnan(lo), -np.inf, lo); hi = np.where(np.isnan(hi), np.inf, hi)
        return lo, hi
    d1 = xv - np.where(yc, 0., yv); d2 = np.where(xc, 0., xv) - yv
    return np.minimum(d1, d2), np.maximum(d1, d2)


def sign_test(xv, xc, yv, yc):
    """NADA2 cen_signtest (two-sided): undeterminable directions are ties; exact binomial and Fong-corrected p."""
    x, y = np.asarray(xv, float).copy(), np.asarray(yv, float).copy(); xc, yc = np.asarray(xc, bool).copy(), np.asarray(yc, bool).copy()
    for i in range(len(x)):
        if xc[i] and yc[i]: x[i] = y[i] = max(x[i], y[i])
        elif xc[i] and x[i] > y[i]: y[i] = x[i]; yc[i] = True
        elif yc[i] and y[i] > x[i]: x[i] = y[i]; xc[i] = True
    x = x - xc*.001*x; y = y - yc*.001*y
    npl, nmi = int((x > y).sum()), int((y > x).sum()); N = len(x); ties = N - npl - nmi
    out = {'n_pairs': N, 'n_plus': npl, 'n_minus': nmi, 'ties': ties}
    if npl + nmi == 0: return out
    out['p_value'] = float(stats.binomtest(npl, npl + nmi, .5).pvalue)
    if ties:
        num = 1 - stats.binom.cdf(max(npl, nmi) - 1, N, .5); den = 1 - stats.binom.cdf(np.floor((N - ties + 1)/2) - 1, N, .5)
        out['p_value_fong'] = float(num/den)
    return out


def compute_paired(d, c):
    q = c['comparison']; ref = q['reference']; other = [g for g in q['groups'] if g != ref][0]; lv = q['confidence_level']; m = c['method']
    sc = c['paired']['scale']; z = stats.norm.ppf((1+lv)/2)
    u = d[~d.exclude]
    w = u.pivot(index='unit_id', columns='group', values=['value', 'censored'])
    xv, yv = w['value'][other].to_numpy(float), w['value'][ref].to_numpy(float)
    xc, yc = (w['censored'][other] == 'left').to_numpy(), (w['censored'][ref] == 'left').to_numpy()
    n = len(xv); desc = {g: describe(u[u.group == g], lv) for g in q['groups']}
    flags, failing = [], []
    for g, s in desc.items():
        if s['censored_fraction'] > .5: flags.append(f'heavy_censoring:{g}')
    if n < 10: flags.append('fewer_than_10_pairs')
    both_det = int((~xc & ~yc).sum()); both_cen = int((xc & yc).sum())
    res = {'method': m, 'design': 'paired', 'groups': q['groups'], 'reference': ref, 'difference': f'{other} - {ref}', 'confidence_level': lv,
           'n_pairs': n, 'pairs_both_detected': both_det, 'pairs_both_censored': both_cen, 'pairs_one_censored': n - both_det - both_cen,
           'group_summaries': desc, 'diagnostics': flags}
    if m == 'censored_sign_test':
        st = sign_test(xv, xc, yv, yc)
        if 'p_value' not in st:
            res.update(status='withheld', reportable=False, reason='no_pair_with_a_determinable_direction', test=st)
        else:
            res.update(status='estimated', reportable=True, test=st)
    else:
        res['scale'] = sc
        lo, hi = paired_intervals(xv, xc, yv, yc, sc)
        if sc == 'log': res['pairs_uninformative'] = int((~np.isfinite(lo) & ~np.isfinite(hi)).sum())
        if res['pairs_uninformative'] if sc == 'log' else False: flags.append('uninformative_pairs_both_censored')
        if both_det < MIN_DETECTED:
            res.update(status='withheld', reportable=False, reason='fewer_than_3_pairs_with_both_values_detected')
        else:
            f = interval_normal(lo, hi); cov = f['cov']
            if not f['converged'] or not np.isfinite(cov).all() or (np.diag(cov) <= 0).any():
                res.update(status='withheld', reportable=False, reason='likelihood_not_converged_or_singular_information')
            else:
                mu = float(f['theta'][0]); se_ml = float(np.sqrt(cov[0, 0])); s = float(np.exp(f['theta'][1])); k = f['n_used']
                if c['paired']['inference'] == 'wald_z':
                    se, df = se_ml, None; q_ = z; pv = float(2*stats.norm.sf(abs(mu/se)))
                else:
                    # 0.14.1 post-registration: ML SE scaled by sqrt(k/(k-1)) with t(k-1), k informative pairs; equals the paired t test without censoring.
                    se, df = se_ml*np.sqrt(k/(k - 1)), k - 1; q_ = float(stats.t.ppf((1+lv)/2, df)); pv = float(2*stats.t.sf(abs(mu/se), df))
                zz = mu/se
                est = {'mean_difference': mu, 'se': se, 'se_ml': se_ml, 'df': df, 'inference': c['paired']['inference'], 'ci': [mu - q_*se, mu + q_*se],
                       'statistic': zz, 'p_two_sided': pv,
                       'sd_difference': s, 'sd_difference_se': float(s*np.sqrt(cov[1, 1])), 'loglik': f['loglik'],
                       'intervals_used': {k: f[k] for k in ('n_used', 'n_exact', 'n_left', 'n_right', 'n_interval')},
                       'interval_method': ('Wald z' if df is None else f'small-sample t ({df} df; ML SE x sqrt(k/(k-1)))') + ' on the ' + sc + ' scale from the observed information (survreg-type, gaussian, interval2)'}
                if sc == 'log': est.update(ratio=float(np.exp(mu)), ratio_ci=[float(np.exp(mu - q_*se)), float(np.exp(mu + q_*se))])
                res.update(status='estimated', reportable=True, estimate=est)
    if not res['reportable']: failing.append({'reason': res['reason']})
    notes = ['Independent unit: ' + c['design']['unit'] + f'; each unit was measured in both {ref} and {other} (' + c['paired']['pairing_rationale'] + '). Values below the stated limit stayed censored (' + c['design']['limit_source'] + '); nothing was substituted.',
             'Censored rows carry information only that the true value lies below the stated limit; limits must be the sample-specific reporting limits after dilution.']
    if m == 'censored_sign_test':
        notes.append('Censored sign test (Helsel; NADA2 cen_signtest): pairs whose direction cannot be determined are ties and are left out of the exact binomial test; Fong\'s tie-corrected p value is shown in addition. It tests the direction of paired differences and gives no effect size; a non-significant result does not show equivalence.')
    elif sc == 'log':
        notes.append('Paired differences of natural-log values were fitted by normal maximum likelihood with censored differences as intervals (Helsel cen_paired on the log scale); the ratio is of paired geometric means. Pairs with both values censored carry no information about the log difference and were not used (' + str(res.get('pairs_uninformative', 0)) + ').')
    else:
        notes.append('Paired differences on the original scale were fitted by normal maximum likelihood with each censored value spanning 0 to its limit (Helsel cen_paired / NADA2), assuming the differences are normal.')
    if m == 'censored_paired_difference':
        notes.append('Wald z (NADA2 cen_paired) was anticonservative with 20 pairs in the registered calibration (coverage 0.93); it is kept only to match the reference software.' if c['paired']['inference'] == 'wald_z' else
                     'Intervals and p values use a small-sample t correction (ML SE scaled by sqrt(k/(k-1)), k - 1 df over k informative pairs), which equals the paired t test without censoring; Helsel\'s Wald z (NADA2) is available as inference: wald_z.')
    for f in flags:
        if f.startswith('heavy_censoring'): notes.append(f"More than half the values in group {f.split(':')[1]} are censored; results lean heavily on the assumed distribution (parametric) or have little information (sign test).")
        if f == 'fewer_than_10_pairs': notes.append('Fewer than 10 pairs.')
    if not res['reportable']: notes.append('Inference withheld: ' + res['reason'])
    return {'analysis_type': 'censored_comparison', 'primary': res, 'must_mention': notes, 'failing_items': failing,
            'preprocessing': {'excluded_units': d[d.exclude][['unit_id', 'group'] + (['exclusion_reason'] if 'exclusion_reason' in d else [])].to_dict('records')},
            'limitations': ['Two paired groups with left censoring only; no repeated measures beyond pairs, covariates or right/interval-censored inputs.',
                            'No imputation, substitution or regression-on-order-statistics.',
                            'No Prism, NADA or software-equivalence claim beyond the recorded benchmarks.']}


def describe(part, level):
    v = part.value.to_numpy(float); kind = part.censored.to_numpy()
    out = {'n': int(len(v)), 'n_detected': int((kind == 'none').sum()), 'n_below_limit': int((kind == 'left').sum()),
           'n_above_limit': int((kind == 'right').sum()), 'censored_fraction': float((kind != 'none').mean()),
           'lower_limits': sorted(set(v[kind == 'left'].tolist())), 'upper_limits': sorted(set(v[kind == 'right'].tolist())),
           'detected_min': float(v[kind == 'none'].min()) if (kind == 'none').any() else None,
           'detected_max': float(v[kind == 'none'].max()) if (kind == 'none').any() else None}
    if (kind == 'right').any() or not (kind == 'none').any():
        out['km_quartiles'] = None; return out
    M = v.max() + 1.; km = kaplan_meier(M - v, (kind == 'none').astype(int), level)
    q = {}
    for p, name in ((.75, 'q25'), (.5, 'median'), (.25, 'q75')):
        # Flipped survival S(M - x) = P(X < x); the original p-quantile is where flipped survival first drops to <= p.
        hit = [M - t for t, s in zip(km['time'], km['surv']) if s <= 1 - p + 1e-12]
        q[name] = float(hit[0]) if hit else None
    out['km_quartiles'] = q
    return out


def compute(d, c):
    if c.get('paired') is not None: return compute_paired(d, c)
    q = c['comparison']; groups = q['groups']; ref = q['reference']; lv = q['confidence_level']; m = c['method']
    u = d[~d.exclude].copy(); z = stats.norm.ppf((1+lv)/2)
    desc = {g: describe(u[u.group == g], lv) for g in groups}
    flags, failing = [], []
    for g, s in desc.items():
        if s['censored_fraction'] > .5: flags.append(f'heavy_censoring:{g}')
        if s['n'] < 5: flags.append(f'fewer_than_5_units:{g}')
    res = {'method': m, 'groups': groups, 'reference': ref, 'confidence_level': lv, 'group_summaries': desc, 'diagnostics': flags}
    if m == 'peto_peto':
        if sum(s['n_detected'] for s in desc.values()) < 2:
            res.update(status='withheld', reportable=False, reason='fewer_than_two_detected_values'); failing.append({'reason': res['reason']})
        else:
            res.update(status='estimated', reportable=True, test=peto_peto(u.value.to_numpy(float), u.censored.to_numpy(), u.group.to_numpy(), groups))
    else:
        few = [g for g, s in desc.items() if s['n_detected'] < MIN_DETECTED]
        if few:
            res.update(status='withheld', reportable=False, reason='fewer_than_3_detected_values_in:' + ','.join(few))
            failing.append({'reason': res['reason']})
        else:
            order = [ref] + [g for g in groups if g != ref]
            X = np.column_stack([(u.group == g).to_numpy(float) for g in order]); y = np.log(u.value.to_numpy(float)); kind = u.censored.to_numpy()
            full = tobit(y, kind, X); null = tobit(y, kind, np.ones((len(y), 1)))
            cov = full['cov']; k = len(order)
            if not full['converged'] or not np.isfinite(cov).all() or (np.diag(cov) <= 0).any():
                res.update(status='withheld', reportable=False, reason='likelihood_not_converged_or_singular_information'); failing.append({'reason': res['reason']})
            else:
                b = full['theta'][:k]; sigma = float(np.exp(full['theta'][-1]))
                gm = []
                for i, g in enumerate(order):
                    se = float(np.sqrt(cov[i, i])); gm.append({'group': g, 'geometric_mean': float(np.exp(b[i])), 'log_mean': float(b[i]), 'log_se': se,
                                                              'ci': [float(np.exp(b[i]-z*se)), float(np.exp(b[i]+z*se))]})
                ratios = []
                for i, g in enumerate(order[1:], 1):
                    est = b[i]-b[0]; se = float(np.sqrt(cov[i, i]+cov[0, 0]-2*cov[i, 0])); w = est/se
                    ratios.append({'group': g, 'versus': ref, 'ratio_of_geometric_means': float(np.exp(est)), 'log_ratio': float(est), 'log_se': se,
                                   'ci': [float(np.exp(est-z*se)), float(np.exp(est+z*se))], 'wald_z': float(w), 'p_two_sided': float(2*stats.norm.sf(abs(w)))})
                lr = 2*(full['loglik'] - null['loglik'])
                res.update(status='estimated', reportable=True, sigma_log=sigma, sigma_log_se=float(sigma*np.sqrt(cov[-1, -1])),
                           loglik=full['loglik'], geometric_means=gm, ratios=ratios,
                           likelihood_ratio_test={'statistic': float(lr), 'df': k-1, 'p_value': float(stats.chi2.sf(max(lr, 0.), k-1))},
                           interval_method='Wald on the log scale from the observed information (likelihood-based, common SD)')
    notes = ['Independent unit: ' + c['design']['unit'] + '. Values below (or above) the stated limit were kept as censored at that limit (' + c['design']['limit_source'] + '); nothing was replaced by zero, the limit, half the limit or any other number.',
             'Censored rows carry information only that the true value lies beyond the stated limit; limits must be the sample-specific reporting limits after dilution.']
    if m == 'censored_lognormal':
        notes.append('Censored lognormal maximum likelihood assumes log values are normal with a common SD across groups, censoring depends only on the limit, and units are independent; ratios are of geometric means.')
    else:
        notes.append('Peto-Peto is a rank test on flipped values (left-censoring only); it tests for a difference in distributions and gives no effect size. A non-significant result does not show equivalence.')
    for f in flags:
        if f.startswith('heavy_censoring'): notes.append(f"More than half the values in group {f.split(':')[1]} are censored; results lean heavily on the assumed distribution (parametric) or have little information (rank test).")
        if f.startswith('fewer_than_5'): notes.append(f"Group {f.split(':')[1]} has fewer than 5 units.")
    if not res['reportable']: notes.append('Inference withheld: ' + res['reason'])
    return {'analysis_type': 'censored_comparison', 'primary': res, 'must_mention': notes, 'failing_items': failing,
            'preprocessing': {'excluded_units': d[d.exclude][['unit_id'] + (['exclusion_reason'] if 'exclusion_reason' in d else [])].to_dict('records')},
            'limitations': ['Independent units with one value each; no repeated measures, paired designs or covariates.',
                            'No imputation, substitution or regression-on-order-statistics; interval-censored values (between two limits) are not supported.',
                            'No Prism, NADA or software-equivalence claim beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    d = pd.read_csv(Path(run)/'input.csv', dtype={'unit_id': str, 'group': str, 'censored': str})
    if 'exclude' in d: d = d[~d.exclude.astype(str).str.lower().eq('true')]
    if p.get('design') == 'paired': return _render_paired(d, p, token, out, pstyle, font_setup, plt)
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(1.8+1.1*len(p['groups']), 3.6), layout='constrained'); shown = set()
        for i, g in enumerate(p['groups']):
            part = d[d.group == g]; jit = np.linspace(-.15, .15, len(part)) if len(part) > 1 else np.zeros(1)
            for kind, marker, label in (('none', 'o', 'detected'), ('left', 'v', 'below limit (plotted at limit)'), ('right', '^', 'above limit (plotted at limit)')):
                sel = (part.censored.str.lower() == kind).to_numpy()
                if sel.any():
                    ax.scatter(i + jit[sel], part.value[sel], s=20, marker=marker, color=pstyle.color(token, 0 if kind == 'none' else 1),
                               facecolors=None if kind == 'none' else 'none', label=None if label in shown else label); shown.add(label)
        if p['method'] == 'censored_lognormal' or (d.value > 0).all(): ax.set_yscale('log')
        ax.set(xticks=range(len(p['groups'])), xticklabels=p['groups'], ylabel='Value', title='Censored values kept at their limits')
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'censored_values', token, 300)
    return ['censored_values']


def _render_paired(d, p, token, out, pstyle, font_setup, plt):
    ref = p['reference']; other = [g for g in p['groups'] if g != ref][0]
    w = d.pivot(index='unit_id', columns='group', values=['value', 'censored'])
    x, y = w['value'][ref].astype(float), w['value'][other].astype(float)
    cx, cy = w['censored'][ref].str.lower() == 'left', w['censored'][other].str.lower() == 'left'
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(4.2, 4), layout='constrained')
        for sel, marker, label, col in (((~cx) & (~cy), 'o', 'both detected', 0), (cx & ~cy, '<', f'{ref} below limit', 1),
                                        (~cx & cy, 'v', f'{other} below limit', 1), (cx & cy, 'x', 'both below limit', 2)):
            if sel.any(): ax.scatter(x[sel], y[sel], s=22, marker=marker, color=pstyle.color(token, col), label=label)
        lim = [min(x.min(), y.min()), max(x.max(), y.max())]; ax.plot(lim, lim, ls='--', color=pstyle.MUTED, lw=1)
        if (x > 0).all() and (y > 0).all(): ax.set_xscale('log'); ax.set_yscale('log')
        ax.set(xlabel=f'{ref} (censored at limit)', ylabel=f'{other} (censored at limit)', title='Paired values; censored kept at limits')
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'censored_paired', token, 300)
    return ['censored_paired']
