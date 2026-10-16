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
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats
from .survival import kaplan_meier

TYPES = {'censored_comparison'}
METHODS = ('censored_lognormal', 'peto_peto')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'censored_comparison', 'input': None, 'source': None, 'method': None,
            'design': {'unit': None, 'rationale': None, 'outcome': None, 'outcome_unit': None,
                       'limit_source': None},
            'comparison': {'groups': None, 'reference': None, 'confidence_level': .95},
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
    if c['schema_version'] != 1 or c['analysis_type'] != 'censored_comparison' or c['method'] not in METHODS:
        raise ValueError('Declare method censored_lognormal or peto_peto')
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
    if d.unit_id.isna().any() or d.unit_id.duplicated().any(): raise ValueError('One row per independent unit; unit_id must be unique')
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
