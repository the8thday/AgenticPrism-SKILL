"""Ordinal scores (histopathology grades, clinical scores) by cumulative-logit proportional odds (0.15.3).

logit P(Y <= j) = zeta_j - eta, eta = sum beta x (MASS::polr parameterisation): exp(beta) is the
odds ratio for a higher score. One row per independent subject; a declared group factor with a
reference level and optional declared numeric pre-treatment covariates. Maximum likelihood with
an analytic gradient; Wald and profile-likelihood intervals; likelihood-ratio test for the group
factor; proportional-odds check by a nominal (category-specific group effect) likelihood-ratio
test. Repeated scores per subject are refused (use a per-subject summary or a mixed model).
"""
from copy import deepcopy
import math
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'ordinal_regression'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'ordinal_regression', 'input': None, 'source': None,
            'design': {'outcome': None, 'categories': None, 'subject_unit': None, 'independence_rationale': None},
            'group': None, 'covariates': None, 'proportional_odds_check': None,
            'confidence_level': .95, 'report': {'plot_style': 'prism_like'}}
SEPARATION_LIMIT = 10.


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported ordinal configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict) and k in ('design', 'report'):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'ordinal_regression': raise ValueError('Unsupported ordinal schema')
    for name, v in [('input', c['input']), ('source', c['source'])] + [(k, v) for k, v in c['design'].items() if k != 'categories']:
        if not _text(v): raise ValueError('Declare ' + (name if name in ('input', 'source') else 'design.' + name))
    cat = c['design']['categories']
    if not isinstance(cat, list) or len(cat) < 3 or len(set(map(str, cat))) != len(cat): raise ValueError('design.categories must list at least 3 distinct ordered categories, lowest first')
    c['design']['categories'] = [str(x) for x in cat]
    g = c['group']
    if not isinstance(g, dict) or set(g) != {'column', 'reference'} or not _text(g['column']) or not _text(g['reference']): raise ValueError('Declare group {column, reference}')
    cv = c['covariates']
    if cv is None: c['covariates'] = []
    elif not isinstance(cv, list) or not all(isinstance(x, dict) and set(x) == {'name', 'rationale'} and _text(x['name']) and _text(x['rationale']) for x in cv):
        raise ValueError('covariates must be [{name, rationale}] for declared numeric pre-treatment covariates')
    po = c['proportional_odds_check']
    if not isinstance(po, dict) or set(po) != {'method', 'alpha'} or po['method'] != 'nominal_lr' or not 0 < po['alpha'] < .5:
        raise ValueError('Declare proportional_odds_check {method: nominal_lr, alpha}')
    if type(c['confidence_level']) not in (int, float) or not .5 < c['confidence_level'] < 1: raise ValueError('Invalid confidence_level')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype=str)
    gc = c['group']['column']; covs = [x['name'] for x in c['covariates']]
    need = {'subject', 'score', gc, *covs}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.str.strip().str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()): raise ValueError('Every excluded row needs an exclusion_reason')
    else: d['exclude'] = False
    if d.subject.isna().any() or d.subject.duplicated().any(): raise ValueError('One row per independent subject: repeated scores need a per-subject summary or a mixed model (not available here)')
    cats = c['design']['categories']
    bad = sorted(set(d.score.dropna()) - set(cats))
    if bad or d.score.isna().any(): raise ValueError('Scores must be one of the declared categories; unexpected: ' + ', '.join(bad[:5]))
    for col in covs:
        v = pd.to_numeric(d[col], errors='coerce')
        if v.isna().any(): raise ValueError(f'Covariate {col} must be numeric on every row')
        d[col] = v
    if c['group']['reference'] not in set(d[gc]): raise ValueError('group.reference is not a level in the data')
    return d


def _design(d, c):
    gc = c['group']['column']; ref = c['group']['reference']; levels = [ref] + sorted(set(d[gc]) - {ref})
    G = np.column_stack([(d[gc] == l).to_numpy(float) for l in levels[1:]])
    Z = np.column_stack([d[x['name']].to_numpy(float) - d[x['name']].mean() for x in c['covariates']]) if c['covariates'] else np.zeros((len(d), 0))
    y = d.score.map({k: i for i, k in enumerate(c['design']['categories'])}).to_numpy(int)
    return G, Z, y, levels


def _F(z): return 1/(1 + np.exp(-z))
def _f(z): e = np.exp(-np.abs(z)); return e/(1 + e)**2


def negll(theta, X, y, J, nominal=None):
    """Negative log-likelihood and gradient. theta = zeta (J-1) + beta (p) [+ nominal block (J-2) x q]."""
    p = X.shape[1]; zeta = theta[:J - 1]; beta = theta[J - 1:J - 1 + p]; eta = X @ beta
    # Cut-point-specific extra effects for the nominal columns (first cut point carries none).
    n = len(y); add = np.zeros((n, J - 1))
    if nominal is not None:
        q = nominal.shape[1]; B = theta[J - 1 + p:].reshape(J - 2, q)
        add[:, 1:] = nominal @ B.T
    lo = np.full(n, -np.inf); hi = np.full(n, np.inf)
    upper = zeta[None, :] - eta[:, None] - add  # (n, J-1)
    hi_idx = y < J - 1; lo_idx = y > 0
    hi[hi_idx] = upper[hi_idx, y[hi_idx]]; lo[lo_idx] = upper[lo_idx, y[lo_idx] - 1]
    Fh = np.where(np.isfinite(hi), _F(hi), 1.); Fl = np.where(np.isfinite(lo), _F(lo), 0.)
    pr = np.clip(Fh - Fl, 1e-300, None); ll = np.sum(np.log(pr))
    fh = np.where(np.isfinite(hi), _f(hi), 0.); fl = np.where(np.isfinite(lo), _f(lo), 0.)
    gz = np.zeros(J - 1); gu = np.zeros((n, J - 1))
    gu[hi_idx, y[hi_idx]] += fh[hi_idx]/pr[hi_idx]; gu[lo_idx, y[lo_idx] - 1] -= fl[lo_idx]/pr[lo_idx]
    gz = gu.sum(0); gb = -(gu.sum(1)) @ X
    grad = np.r_[gz, gb]
    if nominal is not None: grad = np.r_[grad, -(gu[:, 1:].T @ nominal).ravel()]
    return -ll, -grad, upper


def fit(X, y, J, nominal=None, fixed=None):
    """ML fit; `fixed` = (index, value) holds one coefficient for profile likelihood."""
    counts = np.bincount(y, minlength=J); cum = np.cumsum(counts)[:-1]/len(y)
    z0 = np.log(np.clip(cum, 1e-3, 1 - 1e-3)/(1 - np.clip(cum, 1e-3, 1 - 1e-3)))
    k = J - 1 + X.shape[1] + ((J - 2)*nominal.shape[1] if nominal is not None else 0)
    th0 = np.r_[z0, np.zeros(k - (J - 1))]
    free = np.ones(k, bool)
    if fixed is not None: free[fixed[0]] = False; th0[fixed[0]] = fixed[1]
    def obj(t):
        full = th0.copy(); full[free] = t; v, g, _ = negll(full, X, y, J, nominal); return v, g[free]
    r = optimize.minimize(obj, th0[free], jac=True, method='BFGS', options={'gtol': 1e-10, 'maxiter': 5000})
    full = th0.copy(); full[free] = r.x
    return full, -r.fun, r


def hessian(theta, X, y, J, nominal=None, h=1e-5):
    k = len(theta); H = np.zeros((k, k))
    for i in range(k):
        e = np.zeros(k); e[i] = h
        H[i] = (negll(theta + e, X, y, J, nominal)[1] - negll(theta - e, X, y, J, nominal)[1])/(2*h)
    return (H + H.T)/2


def profile_ci(theta, ll, X, y, J, idx, level):
    crit = stats.chi2.ppf(level, 1)/2; est = theta[idx]
    def dev(b): return ll - fit(X, y, J, fixed=(idx, b))[1] - crit
    out = []
    for sgn in (-1, 1):
        step = .5; b = est
        for _ in range(60):
            b2 = est + sgn*step
            if dev(b2) > 0: out.append(optimize.brentq(dev, min(b, b2), max(b, b2), xtol=1e-10)); break
            b, step = b2, step*2
        else: out.append(None)
    return out


def compute(d, c):
    use = d[~d.exclude].reset_index(drop=True); cats = c['design']['categories']; J = len(cats)
    G, Z, y, levels = _design(use, c); X = np.column_stack([G, Z])
    counts = np.bincount(y, minlength=J)
    if (counts == 0).any(): raise ValueError('Declared categories with no observations: ' + ', '.join(cats[i] for i in np.where(counts == 0)[0]) + ' (merge categories only by a rule declared before seeing the data)')
    lv = c['confidence_level']; z = stats.norm.ppf(.5 + lv/2)
    theta, ll, res = fit(X, y, J); cov = np.linalg.inv(hessian(theta, X, y, J)); se = np.sqrt(np.diag(cov))
    names = [f'{c["group"]["column"]}={l}' for l in levels[1:]] + [x['name'] for x in c['covariates']]
    coefs = []; failing = []
    for i, nm in enumerate(names):
        k = J - 1 + i; b = theta[k]
        rec = {'term': nm, 'log_odds_ratio': float(b), 'se': float(se[k]), 'odds_ratio': float(math.exp(b)),
               'wald_ci': [float(math.exp(b - z*se[k])), float(math.exp(b + z*se[k]))], 'reasons': []}
        if abs(b) > SEPARATION_LIMIT or not math.isfinite(se[k]): rec['reasons'].append('separation or non-convergence (|log OR| > 10)')
        else:
            pci = profile_ci(theta, ll, X, y, J, k, lv)
            rec['profile_ci'] = [float(math.exp(v)) if v is not None else None for v in pci]
            if None in pci: rec['reasons'].append('profile-likelihood limit not found (flat likelihood)')
        rec['status'] = 'withheld' if rec['reasons'] else 'reportable'
        coefs.append(rec)
        if rec['status'] == 'withheld': failing.append({'term': nm, 'reason': '; '.join(rec['reasons'])})
    # Group LR test against the model without the group factor.
    _, ll0, _ = fit(Z if Z.shape[1] else np.zeros((len(y), 0)), y, J)
    lr = 2*(ll - ll0); gdf = len(levels) - 1
    test = {'statistic': float(lr), 'df': gdf, 'p_value': float(stats.chi2.sf(lr, gdf))}
    # Proportional-odds check: category-specific group effects.
    po = {'method': 'nominal_lr', 'alpha': c['proportional_odds_check']['alpha']}
    try:
        thn, lln, rn = fit(X, y, J, nominal=G)
        _, _, upper = negll(thn, X, y, J, nominal=G)
        mono = bool(np.all(np.diff(upper, axis=1) > 0))
        stat = 2*(lln - ll); df = (J - 2)*G.shape[1]
        po.update(statistic=float(stat), df=df, p_value=float(stats.chi2.sf(max(stat, 0.), df)), monotone=mono)
        po['status'] = 'unavailable' if not mono else 'violated' if po['p_value'] < po['alpha'] else 'not_rejected'
    except (np.linalg.LinAlgError, ValueError) as exc:
        po.update(status='unavailable', reason=str(exc)[:200])
    if po['status'] == 'violated': failing.append({'term': 'proportional_odds', 'reason': f"nominal LR p = {po['p_value']:.3g} < {po['alpha']}: a single odds ratio does not summarise the group effect"})
    sparse = [cats[i] for i in np.where(counts < 5)[0]]
    table = pd.crosstab(use[c['group']['column']], use.score).reindex(columns=cats, fill_value=0)
    pred = []
    for i, l in enumerate(levels):
        xg = np.zeros(X.shape[1]);
        if i: xg[i - 1] = 1
        cum = _F(theta[:J - 1] - xg @ theta[J - 1:]); probs = np.diff(np.r_[0, cum, 1])
        pred.append({'group': l, 'category_probabilities': dict(zip(cats, map(float, probs)))})
    primary = {'categories': cats, 'groups': levels, 'reference': levels[0], 'coefficients': coefs,
               'thresholds': dict(zip([f'{cats[j]}|{cats[j + 1]}' for j in range(J - 1)], map(float, theta[:J - 1]))),
               'group_lr_test': test, 'proportional_odds_check': po, 'loglik': float(ll), 'n_subjects': int(len(y)),
               'counts': {str(k): {str(kk): int(vv) for kk, vv in v.items()} for k, v in table.to_dict('index').items()},
               'sparse_categories': sparse, 'predicted_probabilities_at_covariate_means': pred,
               'excluded': d[d.exclude][['subject', 'exclusion_reason']].to_dict('records') if d.exclude.any() else []}
    return {'analysis_type': 'ordinal_regression', 'primary': primary, 'must_mention': _notes(c, primary), 'failing_items': failing,
            'limitations': ['Independent subjects only: repeated scores over time, cage or litter clustering are not modelled.',
                            'Proportional odds is an assumption; the nominal LR test has low power in small studies, so not rejecting it does not establish it.',
                            'Scores are treated as ordered categories, not as interval numbers; means of scores are not estimated.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks (MASS::polr, ordinal::clm).']}


def _notes(c, p):
    d = c['design']; n = [f"Ordinal outcome {d['outcome']} with categories {' < '.join(p['categories'])}; one score per {d['subject_unit']} ({d['independence_rationale']}); cumulative-logit proportional-odds model, {p['n_subjects']} subjects; odds ratios above 1 mean higher scores than {p['reference']}."]
    lv = int(100*c['confidence_level'])
    for x in p['coefficients']:
        n.append(f"{x['term']}: OR {x['odds_ratio']:.3g}" + (f" ({lv}% profile CI {x['profile_ci'][0]:.3g}-{x['profile_ci'][1]:.3g}; Wald {x['wald_ci'][0]:.3g}-{x['wald_ci'][1]:.3g})" if x['status'] == 'reportable' else f" withheld: {'; '.join(x['reasons'])}") + '.')
    t = p['group_lr_test']; n.append(f"Group LR test {t['statistic']:.3g} on {t['df']} df, p {t['p_value']:.3g}.")
    po = p['proportional_odds_check']
    n.append({'not_rejected': f"Proportional odds not rejected (nominal LR p {po.get('p_value', float('nan')):.3g}); with few subjects this test has little power.",
              'violated': f"Proportional odds rejected (nominal LR p {po.get('p_value', float('nan')):.3g}): report category-wise comparisons; a single odds ratio does not summarise the effect.",
              'unavailable': 'Proportional-odds check unavailable (category-specific model not estimable or non-monotone).'}[po['status']])
    if p['sparse_categories']: n.append('Sparse categories (< 5 subjects): ' + ', '.join(p['sparse_categories']) + '; estimates involving them are imprecise.')
    if c['covariates']: n.append('Adjusted for declared pre-treatment covariates (centred): ' + ', '.join(f"{x['name']} ({x['rationale']})" for x in c['covariates']) + '.')
    return n


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    cats = p['categories']; groups = p['groups']
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.4, .5 + .45*len(groups)), layout='constrained')
        for gi, g in enumerate(groups):
            row = p['counts'].get(g, {}); tot = sum(row.values()) or 1; left = 0.
            for k, cat in enumerate(cats):
                w = row.get(cat, 0)/tot
                ax.barh(gi, w, left=left, color=pstyle.color(token, k), label=cat if gi == 0 else None); left += w
        ax.set_yticks(range(len(groups))); ax.set_yticklabels(groups); ax.set_xlim(0, 1); ax.invert_yaxis()
        ax.set(xlabel='Proportion of subjects', title='Observed score distribution'); ax.legend(frameon=False, fontsize=7, ncol=min(len(cats), 5), loc='lower center', bbox_to_anchor=(.5, 1.08))
        pstyle.save(fig, out/'ordinal_distribution', token, 300)
    return ['ordinal_distribution']
