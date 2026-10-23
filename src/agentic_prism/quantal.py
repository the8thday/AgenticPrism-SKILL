"""Quantal (responder / mortality) dose response: probit, logit or cloglog with Fieller limits (0.15.3).

P(dose) = c + (1 - c) F(a + b x), x = log10(dose) or dose, F the declared link's distribution
function; c is the natural (background) response: 0, a declared known rate (Abbott) or estimated
from control groups. Maximum likelihood by Fisher scoring with expected information (as R glm).
EDp solves F(a + b x) = p (additional response over background) with Fieller limits for
x_p = (F^-1(p) - a)/b; g = t^2 V(b)/b^2 >= 1 gives unbounded limits and the EDp is withheld.
Heterogeneity: Pearson chi-square on k - parameters degrees of freedom; Finney's policy scales
the covariance by chi2/df and uses t on those df. Relative potency (log10 dose only): parallel
lines with common slope, rho = D_reference / D_test at equal effect = 10^((a_T - a_R)/b), Fieller
limits; parallelism by a likelihood-ratio test or declared slope-ratio equivalence margins.
"""
from copy import deepcopy
import math
import numpy as np
import pandas as pd
from scipy import special, stats

TYPES = {'quantal_dose_response'}
LINKS = {
    'probit': (stats.norm.cdf, stats.norm.pdf, stats.norm.ppf),
    'logit': (lambda z: 1/(1 + np.exp(-z)), lambda z: np.exp(-z)/(1 + np.exp(-z))**2, lambda p: np.log(p/(1 - p))),
    'cloglog': (lambda z: 1 - np.exp(-np.exp(z)), lambda z: np.exp(z - np.exp(z)), lambda p: np.log(-np.log(1 - p))),
}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'quantal_dose_response', 'input': None, 'source': None,
            'design': {'endpoint': None, 'response_definition': None, 'subject_unit': None, 'dose_unit': None, 'independence_rationale': None},
            'data_format': None,
            'model': {'link': None, 'metameter': None, 'rationale': None},
            'natural_response': None,
            'heterogeneity': None,
            'effective_doses': {'percent': [50], 'confidence_level': .95},
            'relative_potency': None,
            'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported quantal configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict) and k in ('design', 'model', 'effective_doses', 'report'):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'quantal_dose_response': raise ValueError('Unsupported quantal schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + (name if name in ('input', 'source') else 'design.' + name))
    if c['data_format'] not in ('grouped', 'individual'): raise ValueError('data_format must be grouped (dose, n, responders) or individual (dose, response 0/1)')
    m = c['model']
    if m['link'] not in LINKS: raise ValueError('model.link must be probit, logit or cloglog')
    if m['metameter'] not in ('log10_dose', 'dose'): raise ValueError('model.metameter must be log10_dose or dose')
    if not _text(m['rationale']): raise ValueError('Declare model.rationale (link and metameter chosen before fitting)')
    nr = c['natural_response']
    if not isinstance(nr, dict) or nr.get('method') not in ('none', 'abbott_known', 'estimate'): raise ValueError('Declare natural_response.method: none, abbott_known (with rate and source) or estimate (needs control groups)')
    if nr['method'] == 'abbott_known':
        if set(nr) != {'method', 'rate', 'source'} or type(nr['rate']) not in (int, float) or not 0 < nr['rate'] < 1 or not _text(nr['source']): raise ValueError('abbott_known needs rate in (0, 1) and source')
    elif set(nr) - {'method', 'rationale'} or not _text(nr.get('rationale')): raise ValueError(f"natural_response {nr['method']} needs a rationale")
    h = c['heterogeneity']
    if not isinstance(h, dict) or h.get('policy') not in ('finney_if_significant', 'always', 'never') or not _text(h.get('rationale')): raise ValueError('Declare heterogeneity {policy: finney_if_significant, always or never, rationale}')
    if h['policy'] == 'finney_if_significant':
        if set(h) != {'policy', 'alpha', 'rationale'} or type(h['alpha']) not in (int, float) or not 0 < h['alpha'] < .5: raise ValueError('finney_if_significant needs alpha in (0, 0.5)')
    elif set(h) != {'policy', 'rationale'}: raise ValueError('heterogeneity always/never take only policy and rationale')
    e = c['effective_doses']
    if not isinstance(e['percent'], list) or not e['percent'] or not all(type(p) in (int, float) and 0 < p < 100 for p in e['percent']) or len(set(e['percent'])) != len(e['percent']):
        raise ValueError('effective_doses.percent must be distinct values in (0, 100)')
    if type(e['confidence_level']) not in (int, float) or not .5 < e['confidence_level'] < 1: raise ValueError('Invalid effective_doses.confidence_level')
    rp = c['relative_potency']
    if rp is None: c.pop('relative_potency')
    else:
        if not isinstance(rp, dict) or not {'reference', 'parallelism', 'rationale'} <= set(rp) <= {'reference', 'parallelism', 'rationale', 'confidence_level'} or not _text(rp['reference']) or not _text(rp['rationale']):
            raise ValueError('relative_potency needs reference, parallelism and rationale (optional confidence_level)')
        rp.setdefault('confidence_level', .95)
        pa = rp['parallelism']
        if not isinstance(pa, dict) or pa.get('method') not in ('lr_test', 'slope_ratio_equivalence'): raise ValueError('parallelism.method must be lr_test or slope_ratio_equivalence')
        if pa['method'] == 'lr_test' and (set(pa) != {'method', 'alpha'} or not 0 < pa['alpha'] < .5): raise ValueError('lr_test needs alpha')
        if pa['method'] == 'slope_ratio_equivalence':
            if set(pa) != {'method', 'margin', 'confidence_level', 'source'} or not isinstance(pa['margin'], list) or len(pa['margin']) != 2 or not 0 < pa['margin'][0] < 1 < pa['margin'][1] \
                    or not .5 < pa['confidence_level'] < 1 or not _text(pa['source']):
                raise ValueError('slope_ratio_equivalence needs margin [low < 1 < high], confidence_level and source')
        if m['metameter'] != 'log10_dose': raise ValueError('relative potency needs metameter log10_dose (parallel lines on log dose)')
        if nr['method'] == 'estimate': raise ValueError('relative potency with an estimated natural response is not available')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'preparation': str})
    need = {'dose', 'n', 'responders'} if c['data_format'] == 'grouped' else {'dose', 'response'}
    if c.get('relative_potency'): need.add('preparation')
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'preparation', 'subject', 'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()): raise ValueError('Every excluded row needs an exclusion_reason')
    else: d['exclude'] = False
    for col in need - {'preparation'}:
        v = pd.to_numeric(d[col], errors='coerce')
        if v.isna().any() or not np.isfinite(v).all(): raise ValueError(f'{col} must be numeric')
        d[col] = v
    if (d.dose < 0).any(): raise ValueError('dose cannot be negative')
    if c['data_format'] == 'grouped':
        if (d.n <= 0).any() or (d.n != d.n.round()).any() or (d.responders < 0).any() or (d.responders > d.n).any() or (d.responders != d.responders.round()).any():
            raise ValueError('n must be positive integers and 0 <= responders <= n')
    else:
        if not d.response.isin([0, 1]).all(): raise ValueError('response must be 0 or 1')
        if 'subject' in d and d.subject.duplicated().any(): raise ValueError('Each subject appears once (one quantal response per independent subject)')
    if 'preparation' not in d: d['preparation'] = 'single'
    return d


def _groups(d, c):
    use = d[~d.exclude]
    if c['data_format'] == 'individual':
        g = use.groupby(['preparation', 'dose']).response.agg(['size', 'sum']).reset_index().rename(columns={'size': 'n', 'sum': 'responders'})
    else:
        g = use.groupby(['preparation', 'dose'])[['n', 'responders']].sum().reset_index()
    return g


def fisher(x_design, n, r, link, c_fixed=None, estimate_c=False, start=None, tol=1e-12, maxit=200):
    """ML by Fisher scoring. x_design: (k, p) linear predictor design. Returns params, covariance, loglik."""
    F, f, Finv = LINKS[link]
    p_lin = x_design.shape[1]
    if start is None:
        pe = np.clip((r + .5)/(n + 1), .02, .98); y = Finv(pe)
        dosed = np.all(np.abs(x_design) < 1e5, axis=1)  # control rows (predictor -1e6) inform only the background
        beta, *_ = np.linalg.lstsq(x_design[dosed], y[dosed], rcond=None)
        theta = np.r_[beta, [math.log(.05/.95)] if estimate_c else []]
    else: theta = np.asarray(start, float)
    c0 = c_fixed or 0.
    def parts(th):
        eta = x_design @ th[:p_lin]
        cc = float(special.expit(th[p_lin])) if estimate_c else c0
        Fz, fz = F(eta), f(eta); P = np.clip(cc + (1 - cc)*Fz, 1e-15, 1 - 1e-15)
        G = (1 - cc)*fz[:, None]*x_design
        if estimate_c: G = np.column_stack([G, (1 - Fz)*cc*(1 - cc)])  # derivative w.r.t. logit c
        return P, G
    ll_old = -np.inf
    for it in range(maxit):
        P, G = parts(theta)
        w = n/(P*(1 - P)); score = G.T @ ((r - n*P)/(P*(1 - P))); info = G.T @ (G*w[:, None])
        step = np.linalg.solve(info, score)
        lam = 1.
        ll = float(np.sum(r*np.log(P) + (n - r)*np.log(1 - P)))
        while lam > 1e-6:
            P2, _ = parts(theta + lam*step); ll2 = float(np.sum(r*np.log(P2) + (n - r)*np.log(1 - P2)))
            if ll2 >= ll - 1e-12: break
            lam /= 2
        theta = theta + lam*step
        if np.max(np.abs(lam*step)) < tol and abs(ll2 - ll_old) < 1e-12: break
        ll_old = ll2
    P, G = parts(theta); w = n/(P*(1 - P)); info = G.T @ (G*w[:, None])
    cov = np.linalg.inv(info); ll = float(np.sum(r*np.log(P) + (n - r)*np.log(1 - P)))
    return theta, cov, ll, P, it + 1


def fieller(num, den, v_nn, v_dd, v_nd, t):
    """Fieller limits for num/den; None when unbounded (g >= 1)."""
    g = t*t*v_dd/den**2
    if g >= 1: return None, g
    est = num/den; a = den*den - t*t*v_dd; b = -2*(num*den - t*t*v_nd); cq = num*num - t*t*v_nn
    disc = b*b - 4*a*cq
    if disc < 0: return None, g
    r1, r2 = (-b - math.sqrt(disc))/(2*a), (-b + math.sqrt(disc))/(2*a)
    return [min(r1, r2), max(r1, r2)], g


def _xform(c, dose):
    return np.log10(dose) if c['model']['metameter'] == 'log10_dose' else np.asarray(dose, float)


def _back(c, x):
    return 10**x if c['model']['metameter'] == 'log10_dose' else x


def _fit_single(g, c):
    nr = c['natural_response']; link = c['model']['link']
    pos = g[g.dose > 0] if c['model']['metameter'] == 'log10_dose' else g
    ctrl = g[g.dose == 0]
    use = pos if nr['method'] != 'estimate' else pd.concat([pos, ctrl]) if c['model']['metameter'] == 'log10_dose' else g
    if nr['method'] == 'estimate' and ctrl.empty: raise ValueError('natural_response estimate needs dose-0 control groups')
    if nr['method'] == 'none' and not ctrl.empty and ctrl.responders.sum() > 0:
        raise ValueError(f'Controls show {int(ctrl.responders.sum())} responder(s) but natural_response is none; declare abbott_known or estimate')
    x = _xform(c, use.dose.where(use.dose > 0, 1.)) if c['model']['metameter'] == 'log10_dose' else use.dose.to_numpy(float)
    X = np.column_stack([np.ones(len(use)), x])
    if c['model']['metameter'] == 'log10_dose' and nr['method'] == 'estimate':
        # Controls inform only the background: their predictor is set far below the doses tested.
        X[(use.dose == 0).to_numpy(), 1] = -1e6
    n, r = use.n.to_numpy(float), use.responders.to_numpy(float)
    diag = []
    try:
        th, cov, ll, P, iters = fisher(X, n, r, link, c_fixed=nr.get('rate') if nr['method'] == 'abbott_known' else None, estimate_c=nr['method'] == 'estimate')
        boundary = nr['method'] == 'estimate' and float(special.expit(th[2])) < 1e-4
    except (np.linalg.LinAlgError, OverflowError, FloatingPointError):
        if nr['method'] != 'estimate': raise
        boundary = True
    if boundary:
        # The natural-response estimate is on the boundary (0): its information vanishes, so the
        # model is refitted with c = 0 and the boundary is reported.
        keep = (use.dose > 0).to_numpy() if c['model']['metameter'] == 'log10_dose' else np.ones(len(use), bool)
        th, cov, ll, P0, iters = fisher(X[keep], n[keep], r[keep], link)
        P = np.full(len(use), 1e-12); P[keep] = P0; th = np.r_[th, -np.inf]; cov = np.pad(cov, ((0, 1), (0, 1))); diag.append('natural_response_at_boundary_zero')
    k = len(use); npar = len(th); P = np.clip(P, 1e-12, 1 - 1e-12)
    chi2 = float(np.sum((r - n*P)**2/(n*P*(1 - P)))); df = k - npar
    return {'theta': th, 'cov': cov, 'loglik': ll, 'P': P, 'iterations': iters, 'k': k, 'df': df, 'chi2': chi2, 'groups': use, 'X': X, 'diagnostics': diag}


def _hetero(c, chi2, df):
    h = c['heterogeneity']; p = float(stats.chi2.sf(chi2, df)) if df > 0 else None
    scale = (h['policy'] == 'always' or (h['policy'] == 'finney_if_significant' and p is not None and p < h['alpha'])) and df > 0
    return {'pearson_chi2': chi2, 'df': df, 'p_value': p, 'scaled': bool(scale), 'heterogeneity_factor': chi2/df if scale else 1.}


def _eds(c, fit, het, prep_label, tested):
    F, f, Finv = LINKS[c['model']['link']]; a, b = fit['theta'][0], fit['theta'][1]; V = fit['cov'][:2, :2]*het['heterogeneity_factor']
    lv = c['effective_doses']['confidence_level']
    t = stats.t.ppf(.5 + lv/2, het['df']) if het['scaled'] else stats.norm.ppf(.5 + lv/2)
    out = []
    partial = int(((fit['groups'].responders > 0) & (fit['groups'].responders < fit['groups'].n) & (fit['groups'].dose > 0)).sum())
    pd_ok = bool(np.all(np.isfinite(V)) and np.all(np.linalg.eigvalsh(V) > 0))
    for pct in c['effective_doses']['percent']:
        q = float(Finv(pct/100)); x = (q - a)/b
        lim, g = fieller(q - a, b, V[0, 0], V[1, 1], -V[0, 1], t)
        var_d = V[0, 0] + 2*x*V[0, 1] + x*x*V[1, 1]
        se_d = math.sqrt(var_d)/abs(b) if math.isfinite(var_d) and var_d > 0 else None
        row = {'preparation': prep_label, 'percent': pct, 'x': float(x), 'dose': float(_back(c, x)), 'g': float(g), 't': float(t),
               'fieller_x': [float(v) for v in lim] if lim else None, 'fieller_dose': [float(_back(c, v)) for v in lim] if lim else None,
               'delta_se_x': se_d, 'delta_dose_audit': [float(_back(c, x - t*se_d)), float(_back(c, x + t*se_d))] if se_d else None, 'reasons': []}
        if not pd_ok: row['reasons'].append('covariance of (a, b) not positive definite (near separation)')
        if partial == 0: row['reasons'].append('no dose group with a partial response (separation): maximum likelihood estimate not reliable')
        if b <= 0: row['reasons'].append('non-positive slope: response does not increase with dose')
        if lim is None: row['reasons'].append(f'Fieller limits unbounded (g = {g:.3g} >= 1)')
        if not tested[0] <= row['dose'] <= tested[1]: row['reasons'].append(f"ED{pct:g} outside the tested dose range {tested[0]:g}-{tested[1]:g}")
        row['status'] = 'withheld' if row['reasons'] else 'reportable'
        if row['status'] == 'withheld': row['fieller_dose_audit'] = row.pop('fieller_dose'); row['fieller_dose'] = None
        out.append(row)
    return out, partial


def _relative_potency(g, c):
    rp = c['relative_potency']; link = c['model']['link']; preps = sorted(g.preparation.unique())
    if rp['reference'] not in preps: raise ValueError('relative_potency.reference is not a preparation in the data')
    g = g[g.dose > 0]; others = [p for p in preps if p != rp['reference']]; order = [rp['reference']] + others
    x = np.log10(g.dose.to_numpy(float)); n, r = g.n.to_numpy(float), g.responders.to_numpy(float)
    D = np.column_stack([(g.preparation == p).to_numpy(float) for p in order])
    par, cov, ll_par, P, _ = fisher(np.column_stack([D, x]), n, r, link)
    sep, cov_sep, ll_sep, _, _ = fisher(np.column_stack([D, D*x[:, None]]), n, r, link)
    k = len(g); df = k - len(par); chi2 = float(np.sum((r - n*P)**2/(n*P*(1 - P)))); het = _hetero(c, chi2, df)
    lr = 2*(ll_sep - ll_par); lr_df = len(order) - 1; lr_p = float(stats.chi2.sf(lr, lr_df))
    lv = rp['confidence_level']; t = stats.t.ppf(.5 + lv/2, df) if het['scaled'] else stats.norm.ppf(.5 + lv/2)
    V = cov*het['heterogeneity_factor']; b = par[-1]; nb = len(order); res = []
    pa = rp['parallelism']
    for j, p in enumerate(others, start=1):
        num = par[j] - par[0]; v_nn = V[j, j] + V[0, 0] - 2*V[j, 0]; v_dd = V[-1, -1]; v_nd = V[j, -1] - V[0, -1]
        lim, gq = fieller(num, b, v_nn, v_dd, v_nd, t)
        ratio = sep[nb + j]/sep[nb]; se_lr = math.sqrt(cov_sep[nb + j, nb + j]/sep[nb + j]**2 + cov_sep[nb, nb]/sep[nb]**2 - 2*cov_sep[nb + j, nb]/(sep[nb + j]*sep[nb])) if sep[nb + j] > 0 and sep[nb] > 0 else float('nan')
        row = {'test': p, 'reference': rp['reference'], 'log10_potency': float(num/b), 'potency_ratio': float(10**(num/b)),
               'fieller_potency_ratio': [float(10**v) for v in lim] if lim else None, 'g': float(gq), 'slope_ratio_test_over_reference': float(ratio), 'reasons': []}
        if pa['method'] == 'slope_ratio_equivalence':
            z = stats.norm.ppf(.5 + pa['confidence_level']/2); ci = [float(ratio*math.exp(-z*se_lr)), float(ratio*math.exp(z*se_lr))] if math.isfinite(se_lr) else None
            row['slope_ratio_ci'] = ci
            if ci is None or not (pa['margin'][0] <= ci[0] and ci[1] <= pa['margin'][1]): row['reasons'].append(f"slope-ratio interval not within the declared margin {pa['margin']}")
        elif lr_p < pa['alpha']: row['reasons'].append(f"parallelism rejected (LR p = {lr_p:.3g} < {pa['alpha']})")
        if lim is None: row['reasons'].append(f'Fieller limits unbounded (g = {gq:.3g})')
        row['status'] = 'withheld' if row['reasons'] else 'reportable'
        if row['status'] == 'withheld': row['fieller_potency_ratio_audit'] = row.pop('fieller_potency_ratio'); row['fieller_potency_ratio'] = None
        res.append(row)
    return {'common_slope': float(b), 'parallelism_lr': {'statistic': float(lr), 'df': lr_df, 'p_value': lr_p, 'method': pa['method']},
            'heterogeneity': het, 'comparisons': res}


def compute(d, c):
    g = _groups(d, c); results = []; failing = []
    for prep, part in g.groupby('preparation', sort=True):
        try:
            fit = _fit_single(part, c)
        except np.linalg.LinAlgError:
            # Singular information: all-or-none responses or separation, no finite ML estimate.
            reason = 'information matrix singular: no finite maximum likelihood estimate (separation or all-or-none responses)'
            results.append({'preparation': prep, 'status': 'not_estimable', 'reason': reason,
                            'groups': part[['dose', 'n', 'responders']].to_dict('records'),
                            'effective_doses': [{'preparation': prep, 'percent': p, 'status': 'withheld', 'reasons': [reason]} for p in c['effective_doses']['percent']]})
            failing += [{'preparation': prep, 'item': f'ED{p:g}', 'reason': reason} for p in c['effective_doses']['percent']]
            continue
        het = _hetero(c, fit['chi2'], fit['df'])
        pos = part[part.dose > 0].dose; tested = [float(pos.min()), float(pos.max())]
        eds, partial = _eds(c, fit, het, prep, tested)
        groups = fit['groups'].assign(expected=fit['n_exp'] if 'n_exp' in fit else fit['groups'].n.to_numpy()*fit['P'])
        rec = {'preparation': prep, 'intercept': float(fit['theta'][0]), 'slope': float(fit['theta'][1]),
               'covariance_ab': fit['cov'][:2, :2].tolist(), 'loglik': fit['loglik'], 'iterations': fit['iterations'],
               'natural_response': (float(special.expit(fit['theta'][2])) if c['natural_response']['method'] == 'estimate'
                                    else c['natural_response'].get('rate', 0.)), 'heterogeneity': het,
               'partial_response_groups': partial, 'tested_dose_range': tested, 'effective_doses': eds,
               'groups': groups[['dose', 'n', 'responders', 'expected']].to_dict('records')}
        if c['natural_response']['method'] == 'estimate' and math.isfinite(fit['theta'][2]):
            se = math.sqrt(fit['cov'][2, 2]); lc = fit['theta'][2]; z = stats.norm.ppf(.975)
            rec['natural_response_ci95'] = [float(special.expit(lc - z*se)), float(special.expit(lc + z*se))]
        rec['diagnostics'] = fit['diagnostics'] + (['only_one_partial_response_group'] if partial == 1 else [])
        results.append(rec)
        failing += [{'preparation': prep, 'item': f"ED{e['percent']:g}", 'reason': '; '.join(e['reasons'])} for e in eds if e['status'] == 'withheld']
    primary = {'link': c['model']['link'], 'metameter': c['model']['metameter'], 'fits': results}
    if c.get('relative_potency'):
        try:
            primary['relative_potency'] = _relative_potency(g, c)
        except np.linalg.LinAlgError:
            primary['relative_potency'] = {'status': 'not_estimable', 'comparisons': [{'test': p, 'reference': c['relative_potency']['reference'], 'status': 'withheld',
                                           'reasons': ['information matrix singular in the parallel-line fit (separation)']} for p in sorted(set(g.preparation) - {c['relative_potency']['reference']})]}
        failing += [{'item': f"relative potency {x['test']} vs {x['reference']}", 'reason': '; '.join(x['reasons'])} for x in primary['relative_potency']['comparisons'] if x['status'] == 'withheld']
    return {'analysis_type': 'quantal_dose_response', 'primary': primary, 'must_mention': _notes(c, primary), 'failing_items': failing,
            'limitations': ['Independent subjects per dose group are assumed (litter, cage or tank effects are not modelled beyond the heterogeneity factor).',
                            'Fieller limits are approximate (asymptotic normality of a and b); in small groups they can undercover (see the registered calibration).',
                            'No model is chosen by fit: the link and metameter are declared; multistage, Weibull or benchmark-dose (BMD) models are not part of this analysis type.',
                            'No Prism, SAS PROBIT or software-equivalence claim beyond the recorded benchmarks.']}


def _notes(c, p):
    m = c['model']; nr = c['natural_response']; h = c['heterogeneity']; d = c['design']
    notes = [f"Quantal dose response for {d['endpoint']} ({d['response_definition']}); unit {d['subject_unit']} ({d['independence_rationale']}); {m['link']} on {m['metameter']} declared before fitting ({m['rationale']}); maximum likelihood with expected information.",
             {'none': 'No natural response (controls declared without responders).',
              'abbott_known': f"Natural response fixed at {nr.get('rate', 0):g} (Abbott, {nr.get('source')}); EDp is the dose giving p% additional response beyond background.",
              'estimate': 'Natural response estimated jointly from control groups; EDp is the dose giving p% additional response beyond background.'}[nr['method']],
             f"Heterogeneity policy {h['policy']}" + (f" (alpha {h['alpha']})" if 'alpha' in h else '') + f" ({h['rationale']}); Fieller limits use t on the residual df when the covariance is scaled, otherwise normal quantiles."]
    for f in p['fits']:
        if f.get('status') == 'not_estimable':
            notes.append(f"{f['preparation']}: not estimable ({f['reason']}); no EDp is reported."); continue
        he = f['heterogeneity']
        notes.append(f"{f['preparation']}: slope {f['slope']:.4g}, Pearson chi2 {he['pearson_chi2']:.3g} on {he['df']} df" + (f" (p {he['p_value']:.3g})" if he['p_value'] is not None else '') +
                     ('; covariance scaled by ' + f"{he['heterogeneity_factor']:.3g}" if he['scaled'] else '; not scaled') + '. ' +
                     '; '.join(f"ED{e['percent']:g} {e['dose']:.4g} {d['dose_unit']}" + (f" ({int(100*c['effective_doses']['confidence_level'])}% Fieller {e['fieller_dose'][0]:.4g}-{e['fieller_dose'][1]:.4g})" if e['status'] == 'reportable' else f" withheld: {'; '.join(e['reasons'])}") for e in f['effective_doses']) + '.')
    rp = p.get('relative_potency')
    if rp and rp.get('status') == 'not_estimable':
        notes.append('Relative potency not estimable: ' + rp['comparisons'][0]['reasons'][0] + '.')
    elif rp:
        notes.append(f"Relative potency by parallel {m['link']} lines (common slope {rp['common_slope']:.4g}); parallelism {rp['parallelism_lr']['method']}: LR {rp['parallelism_lr']['statistic']:.3g} on {rp['parallelism_lr']['df']} df, p {rp['parallelism_lr']['p_value']:.3g}" +
                     (' (a non-significant test does not prove parallelism)' if rp['parallelism_lr']['method'] == 'lr_test' else '') + '. ' +
                     '; '.join(f"{x['test']} vs {x['reference']}: {x['potency_ratio']:.4g}" + (f" (Fieller {x['fieller_potency_ratio'][0]:.4g}-{x['fieller_potency_ratio'][1]:.4g})" if x['status'] == 'reportable' else f" withheld: {'; '.join(x['reasons'])}") for x in rp['comparisons']) + '. A ratio above 1 means the test preparation is more potent.')
    return notes


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    import json
    c = json.loads((Path(run)/'config.resolved.json').read_text()); p = result['primary']; F = LINKS[p['link']][0]
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5, 3.4), layout='constrained')
        for k, f in enumerate(p['fits']):
            if f.get('status') == 'not_estimable': continue
            gr = pd.DataFrame(f['groups']); gp = gr[gr.dose > 0]
            ax.scatter(gp.dose, gp.responders/gp.n, s=16, color=pstyle.color(token, k), label=f['preparation'])
            xs = np.geomspace(gp.dose.min(), gp.dose.max(), 200) if p['metameter'] == 'log10_dose' else np.linspace(gp.dose.min(), gp.dose.max(), 200)
            xx = np.log10(xs) if p['metameter'] == 'log10_dose' else xs
            cc = f['natural_response']; ax.plot(xs, cc + (1 - cc)*F(f['intercept'] + f['slope']*xx), color=pstyle.color(token, k), lw=1)
        if p['metameter'] == 'log10_dose': ax.set_xscale('log')
        ax.set(xlabel=f"Dose ({c['design']['dose_unit']})", ylabel='Proportion responding', title=f"{p['link']} fit (tested range)", ylim=(-.02, 1.02))
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'quantal_fit', token, 300)
    return ['quantal_fit']
