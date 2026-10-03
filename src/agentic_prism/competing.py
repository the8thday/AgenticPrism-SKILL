"""Competing risks: cumulative incidence, Gray's test and Fine-Gray regression (0.13.1).

Numerics follow cmprsk 2.2-12 (Robert Gray) so they can be checked value by value:
`cinc` (Aalen-Johansen cumulative incidence and its variance), `crst` (Gray's K-sample
test, rho = 0, one stratum) and `crr` (Fine-Gray subdistribution proportional hazards:
the same Newton iterations with backtracking, gtol 1e-6, maxiter 10, and the `crrvv`
sandwich variance that includes the estimated censoring distribution; one censoring
group, fixed covariates). Cause-specific Cox hazard ratios use the package Efron Cox with
other causes treated as censored. Pointwise CIF intervals use the log(-log) transform.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from . import survival as sv

TYPES = {'competing_risks'}
TIME_UNITS = ('hour', 'day', 'week', 'month', 'year')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'competing_risks', 'input': None, 'source': None,
            'study': {'endpoint': None, 'time_origin': None, 'time_unit': None, 'censoring_rationale': None,
                      'competing_events_rationale': None, 'rationale': None},
            'causes': None, 'cause_of_interest': None,
            'comparison': {'arms': None, 'control_arm': None, 'covariates': [], 'landmarks': [], 'confidence_level': .95},
            'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported competing-risks configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'competing_risks' or c['schema_version'] != 1:
        raise ValueError('Unsupported competing-risks schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    s = c['study']
    for k in ('endpoint', 'time_origin', 'censoring_rationale', 'competing_events_rationale', 'rationale'):
        _text(s[k], 'study.' + k + ' (competing_events_rationale: why the other events preclude the event of interest rather than censor it)')
    if s['time_unit'] not in TIME_UNITS:
        raise ValueError(f'study.time_unit must be one of {TIME_UNITS}')
    causes = c['causes']
    if not isinstance(causes, dict) or len(causes) < 2 or any(not k.isdigit() or k == '0' for k in causes) or any(not isinstance(v, str) or not v.strip() for v in causes.values()):
        raise ValueError('causes must map at least two positive integer codes ("1", "2", ...) to labels; 0 is censored')
    if c['cause_of_interest'] not in causes:
        raise ValueError('cause_of_interest must be one of the declared cause codes')
    q = c['comparison']
    arms = q['arms']
    if not isinstance(arms, list) or len(arms) < 2 or len(set(arms)) != len(arms) or any(not isinstance(a, str) or not a.strip() for a in arms):
        raise ValueError('Declare at least two distinct arms')
    if q['control_arm'] not in arms:
        raise ValueError('control_arm must be one of the arms (reference for hazard ratios)')
    if not isinstance(q['covariates'], list) or len(set(q['covariates'])) != len(q['covariates']) or set(q['covariates']) & {'subject_id', 'arm', 'time', 'status'}:
        raise ValueError('covariates must be distinct numeric column names')
    if not isinstance(q['landmarks'], list) or any(isinstance(t, bool) or not isinstance(t, (int, float)) or t <= 0 for t in q['landmarks']):
        raise ValueError('landmarks must be positive times')
    lv = q['confidence_level']
    if isinstance(lv, bool) or not isinstance(lv, (int, float)) or not .5 < lv < 1:
        raise ValueError('confidence_level must lie in (0.5, 1)')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    q = cfg['comparison']
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'subject_id', 'arm', 'time', 'status'} | set(q['covariates'])
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    if d.subject_id.duplicated().any():
        raise ValueError('One row per subject (first event only); recurrent events are not supported')
    d['time'] = pd.to_numeric(d.time, errors='raise')
    if not np.isfinite(d.time).all() or (d.time < 0).any():
        raise ValueError('Times must be finite and nonnegative')
    allowed = {'0'} | set(cfg['causes'])
    if set(d.status) - allowed:
        raise ValueError(f'status codes {sorted(set(d.status) - allowed)} are not declared; 0 = censored')
    d['status'] = d.status.astype(int)
    for k in q['covariates']:
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'Covariate {k} must be finite for every subject')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    if set(d.arm) != set(q['arms']):
        raise ValueError('Input arms do not match declared arms')
    return d


# ---------- cmprsk ports ----------

def cinc(time, failed, cause):
    """Port of cmprsk cinc: step-function corners (x, f, v) for one group and cause.
    time sorted ascending; failed = any-cause failure indicator; cause = cause-of-interest indicator."""
    y, ic, icc = np.asarray(time, float), np.asarray(failed, int), np.asarray(cause, int)
    n = len(y); fk = 1.; v1 = v2 = v3 = 0.; x, f, v = [0.], [0.], [0.]
    rs = float(n); ll = 0
    while ll < n:
        lu = ll
        while lu + 1 < n and y[lu + 1] == y[ll]:
            lu += 1
        nd1 = int(icc[ll:lu + 1].sum()); nd2 = int((ic[ll:lu + 1] - icc[ll:lu + 1]).sum()); nd = nd1 + nd2
        if nd:
            fkn = fk * (rs - nd) / rs
            if nd1 > 0:
                f.append(f[-1]); f.append(f[-1] + fk * nd1 / rs)
            if nd2 > 0 and fkn > 0:
                t5 = 1. - (nd2 - 1.) / (rs - 1.) if nd2 > 1 else 1.
                t6 = fk * fk * t5 * nd2 / (rs * rs); t3 = 1. / fkn; t4 = f[-1] / fkn
                v1 += t4 * t4 * t6; v2 += t3 * t4 * t6; v3 += t3 * t3 * t6
            if nd1 > 0:
                t5 = 1. - (nd1 - 1.) / (rs - 1.) if nd1 > 1 else 1.
                t6 = fk * fk * t5 * nd1 / (rs * rs); t3 = 1. / fkn if fkn > 0 else 0.; t4 = 1. + t3 * f[-1]
                v1 += t4 * t4 * t6; v2 += t3 * t4 * t6; v3 += t3 * t3 * t6
                t2 = f[-1]
                x += [y[lu], y[lu]]; v += [v[-1], v1 + t2 * t2 * v3 - 2 * t2 * v2]
            fk = fkn
        rs = n - (lu + 1); ll = lu + 1
    x.append(y[-1]); f.append(f[-1]); v.append(v[-1])
    return np.array(x), np.array(f), np.array(v)


def gray_test(time, code, group, ng, rho=0.):
    """Port of cmprsk crst (one stratum). code: 0 censored, 1 cause of interest, 2 other cause; group 0..ng-1."""
    order = np.argsort(time, kind='stable'); y = np.asarray(time, float)[order]; m = np.asarray(code, int)[order]; ig = np.asarray(group, int)[order]
    n = len(y); ng1 = ng - 1
    rs = np.bincount(ig, minlength=ng).astype(float)
    s = np.zeros(ng1); v = np.zeros((ng1, ng1)); f1m = np.zeros(ng); f1 = np.zeros(ng); skmm = np.ones(ng); skm = np.ones(ng)
    v3 = np.zeros(ng); v2 = np.zeros((ng1, ng)); c = np.zeros((ng, ng)); fm = 0.; f = 0.
    ll = 0
    while ll < n:
        lu = ll
        while lu + 1 < n and y[lu + 1] == y[ll]:
            lu += 1
        d = np.zeros((3, ng))
        for i in range(ll, lu + 1):
            d[m[i], ig[i]] += 1
        nd1, nd2 = d[1].sum(), d[2].sum()
        if nd1 or nd2:
            tr = tq = 0.
            for i in range(ng):
                if rs[i] <= 0:
                    continue
                td = d[1, i] + d[2, i]
                skm[i] = skmm[i] * (rs[i] - td) / rs[i]
                f1[i] = f1m[i] + skmm[i] * d[1, i] / rs[i]
                tr += rs[i] / skmm[i]; tq += rs[i] * (1 - f1m[i]) / skmm[i]
            f = fm + nd1 / tr; fb = (1 - fm) ** rho
            a = np.zeros((ng, ng))
            for i in range(ng):
                if rs[i] <= 0:
                    continue
                t1 = rs[i] / skmm[i]
                a[i, i] = fb * t1 * (1 - t1 / tr)
                if a[i, i] != 0:
                    c[i, i] += a[i, i] * nd1 / (tr * (1 - fm))
                for j in range(i + 1, ng):
                    if rs[j] <= 0:
                        continue
                    a[i, j] = -fb * t1 * rs[j] / (skmm[j] * tr)
                    if a[i, j] != 0:
                        c[i, j] += a[i, j] * nd1 / (tr * (1 - fm))
            for i in range(1, ng):
                for j in range(i):
                    a[i, j] = a[j, i]; c[i, j] = c[j, i]
            for i in range(ng1):
                if rs[i] > 0:
                    s[i] += fb * (d[1, i] - nd1 * rs[i] * (1 - f1m[i]) / (skmm[i] * tq))
            if nd1 > 0:
                for k in range(ng):
                    if rs[k] <= 0:
                        continue
                    t4 = 1 - (1 - f) / skm[k] if skm[k] > 0 else 1.
                    t5 = 1 - (nd1 - 1) / (tr * skmm[k] - 1) if nd1 > 1 else 1.
                    t3 = t5 * skmm[k] * nd1 / (tr * rs[k])
                    v3[k] += t4 * t4 * t3
                    for i in range(ng1):
                        t1 = a[i, k] - t4 * c[i, k]; v2[i, k] += t1 * t4 * t3
                        for j in range(i + 1):
                            v[i, j] += t1 * (a[j, k] - t4 * c[j, k]) * t3
            if nd2 > 0:
                for k in range(ng):
                    if skm[k] <= 0 or d[2, k] <= 0:
                        continue
                    t4 = (1 - f) / skm[k]
                    t5 = 1 - (d[2, k] - 1.) / (rs[k] - 1.) if d[2, k] > 1 else 1.
                    t3 = t5 * (skmm[k] ** 2 * d[2, k]) / rs[k] ** 2
                    v3[k] += t4 * t4 * t3
                    for i in range(ng1):
                        t1 = t4 * c[i, k]; v2[i, k] -= t1 * t4 * t3
                        for j in range(i + 1):
                            v[i, j] += t1 * (t4 * c[j, k]) * t3
        if lu >= n - 1:
            break
        for i in range(ll, lu + 1):
            rs[ig[i]] -= 1
        fm = f; f1m = f1.copy(); skmm = skm.copy(); ll = lu + 1
    for i in range(ng1):
        for j in range(i + 1):
            v[i, j] += np.sum(c[i] * c[j] * v3) + np.sum(c[i] * v2[j]) + np.sum(c[j] * v2[i])
    v = np.tril(v) + np.tril(v, -1).T
    stat = float(s @ np.linalg.solve(v, s)) if np.linalg.matrix_rank(v) == ng1 else None
    return {'statistic': stat, 'df': ng1, 'p_value': None if stat is None else float(chi2.sf(stat, ng1)), 'score': s.tolist()}


def _censoring_weights(t, cens):
    """Kaplan-Meier of the censoring distribution evaluated at each t-, as crr's uuu."""
    out = np.empty(len(t)); g = 1.
    times = np.unique(t)
    gmap = {}
    for u in times:
        gmap[u] = g
        at = t >= u; c = int(((t == u) & (cens == 1)).sum())
        g *= 1 - c / at.sum() if c else 1.
    for i, u in enumerate(t):
        out[i] = gmap[u]
    return out


def _risk(t, ici, x, b, w, i):
    """Weights and covariates of crr's risk set at the time of subject i."""
    eta = np.exp(x @ b)
    before = t < t[i]
    weight = np.where(before, np.where(ici > 1, eta * w[i] / w, 0.), eta)
    return weight


def _crrf(t, ici, x, b, w):
    lik = 0.; s = np.zeros(x.shape[1]); v = np.zeros((x.shape[1],) * 2)
    for u in np.unique(t[ici == 1]):
        events = np.flatnonzero((t == u) & (ici == 1)); first = int(np.flatnonzero(t == u)[0])
        wt = _risk(t, ici, x, b, w, first); sw = wt.sum(); xbar = wt @ x / sw
        lik += -float((x[events] @ b).sum()) + len(events) * np.log(sw)
        s += -x[events].sum(0) + len(events) * xbar
        dx = x - xbar; v += len(events) * (dx.T * wt) @ dx / sw
    return lik, s, v


def fine_gray(time, code, x, gtol=1e-6, maxiter=10):
    """Port of cmprsk crr (one censoring group, fixed covariates): Newton with backtracking and crrvv variance.
    code: 0 censored, 1 cause of interest, 2 competing."""
    order = np.argsort(time, kind='stable'); t = np.asarray(time, float)[order]; ici = np.asarray(code, int)[order]
    x = np.asarray(x, float)[order]; n, p = x.shape
    w = _censoring_weights(t, (ici == 0).astype(int))
    b = np.zeros(p); converged = False
    for ll in range(maxiter + 1):
        lik, s, h = _crrf(t, ici, x, b, w)
        if np.max(np.abs(s) * np.maximum(np.abs(b), 1)) < max(abs(lik), 1) * gtol:
            converged = True; break
        if ll == maxiter:
            break
        sc = -np.linalg.solve(h, s); bn = b + sc; fbn = _crrf(t, ici, x, bn, w)[0]; i = 0
        while not np.isfinite(fbn) or fbn > lik + 1e-4 * float(sc @ s):
            i += 1; sc = sc * .5; bn = b + sc; fbn = _crrf(t, ici, x, bn, w)[0]
            if i > 20:
                break
        if i > 20:
            break
        b = bn
    # crrvv: information v and the sandwich middle v2 (with the censoring-distribution term)
    eta = np.exp(x @ b)
    xb0 = np.zeros(n); xb = np.zeros((n, p))
    for i in np.flatnonzero(ici == 1):
        wt = np.where(t < t[i], np.where(ici > 1, eta * w[i] / w, 0.), eta)
        xb0[i] = wt.sum(); xb[i] = wt @ x
    icrsk = float(n); ss2 = np.zeros(p); qu = np.zeros(p); v = np.zeros((p, p)); v2 = np.zeros((p, p)); lc = 0
    fail = np.flatnonzero(ici == 1)
    for i in range(n):
        # d Lambda hat portion of eta_i over failures j at which i is at risk
        at = fail[(t[fail] <= t[i])]
        st1 = -np.sum((x[i] - xb[at] / xb0[at, None]) * (eta[i] / xb0[at])[:, None], axis=0)
        if ici[i] > 1:
            later = fail[t[fail] > t[i]]
            st1 -= np.sum((x[i] - xb[later] / xb0[later, None]) * (eta[i] * w[later] / w[i] / xb0[later])[:, None], axis=0)
        if ici[i] == 1:
            st1 += x[i] - xb[i] / xb0[i]
            wt = np.where(t < t[i], np.where(ici > 1, eta * w[i] / w, 0.), eta)
            dx = x - xb[i] / xb0[i]; v += (dx.T * wt) @ dx / xb0[i]
        if i == 0 or t[i] > t[i - 1]:
            block = np.flatnonzero(t == t[i])
            if (ici[block] == 0).any():
                # q(u): failures j1 at or after t_i; type-2 subjects j2 before t_i contribute w(j1)/w(j2) (separable)
                prior = (t < t[i]) & (ici > 1)
                a_ = np.sum(eta[prior] / w[prior]); bvec = (eta[prior] / w[prior]) @ x[prior]
                j1 = fail[t[fail] >= t[i]]
                qu = np.sum((w[j1, None] * bvec - xb[j1] * (w[j1] * a_ / xb0[j1])[:, None]) / xb0[j1, None], axis=0)
                ss2 = ss2 - qu / icrsk ** 2 * int((ici[block] == 0).sum())
        st2 = ss2.copy()
        if ici[i] == 0:
            st2 += qu / icrsk
        st1 = st1 + st2; v2 += np.outer(st1, st1)
        if i < n - 1 and t[i + 1] > t[i]:
            icrsk -= (i - lc + 1); lc = i + 1
    hinv = np.linalg.inv(v); cov = hinv @ v2 @ hinv.T
    return {'beta': b, 'covariance': cov, 'information': v, 'converged': converged, 'neg_log_partial_likelihood': float(_crrf(t, ici, x, b, w)[0])}


# ---------- analysis ----------

def _cif_ci(f, var, z):
    if f <= 0 or f >= 1 or var <= 0:
        return [None, None]
    se = np.sqrt(var) / (f * abs(np.log(f)))
    return [float(f ** np.exp(z * se)), float(f ** np.exp(-z * se))]


def _at(x, f, v, landmark):
    k = int(np.searchsorted(x, landmark, side='right') - 1)
    return float(f[k]), float(v[k])


def compute(d, cfg):
    q = cfg['comparison']; arms = q['arms']; control = q['control_arm']; lv = q['confidence_level']; z = norm.ppf((1 + lv) / 2)
    used = d[~d.exclude].sort_values('time', kind='stable').reset_index(drop=True)
    t, status = used.time.to_numpy(float), used.status.to_numpy(int)
    ci_code = int(cfg['cause_of_interest']); failed = (status > 0).astype(int)
    curves, landmarks = [], []
    for arm in arms:
        g = used[used.arm == arm]
        for code, label in cfg['causes'].items():
            x, f, v = cinc(g.time.to_numpy(float), (g.status > 0).astype(int).to_numpy(), (g.status == int(code)).astype(int).to_numpy())
            curves.append({'arm': arm, 'cause': code, 'label': label, 'time': x.tolist(), 'estimate': f.tolist(), 'variance': v.tolist()})
            for lm in q['landmarks']:
                fe, fv = _at(x, f, v, lm)
                beyond = lm > g.time.max()
                landmarks.append({'arm': arm, 'cause': code, 'label': label, 'time': lm, 'cif': None if beyond else fe,
                                  'se': None if beyond else float(np.sqrt(fv)), 'ci': [None, None] if beyond else _cif_ci(fe, fv, z),
                                  'status': 'beyond_follow_up' if beyond else 'estimated', 'at_risk': int((g.time >= lm).sum())})
    group = used.arm.map({a: i for i, a in enumerate(arms)}).to_numpy()
    gray = []
    for code, label in cfg['causes'].items():
        m = np.where(status == int(code), 1, np.where(status > 0, 2, 0))
        gray.append({'cause': code, 'label': label, **gray_test(t, m, group, len(arms)),
                     'method': "Gray's K-sample test (rho = 0), unstratified; cmprsk cuminc"})
    others = [a for a in arms if a != control]
    design = np.column_stack([(used.arm == a).astype(float) for a in others] + [used[c].to_numpy(float) for c in q['covariates']])
    names = [f'arm[{a}] vs {control}' for a in others] + list(q['covariates'])
    code = np.where(status == ci_code, 1, np.where(status > 0, 2, 0))
    fg = fine_gray(t, code, design); se = np.sqrt(np.diag(fg['covariance']))
    subdistribution = [{'term': nm, 'log_shr': float(b), 'se': float(s_), 'shr': float(np.exp(b)), 'ci': [float(np.exp(b - z * s_)), float(np.exp(b + z * s_))],
                        'p_value': float(2 * norm.sf(abs(b / s_)))} for nm, b, s_ in zip(names, fg['beta'], se)]
    cs = sv.cox(t, (status == ci_code).astype(int), design); cse = np.sqrt(np.diag(cs['covariance']))
    cause_specific = [{'term': nm, 'log_hr': float(b), 'se': float(s_), 'hr': float(np.exp(b)), 'ci': [float(np.exp(b - z * s_)), float(np.exp(b + z * s_))],
                       'p_value': float(2 * norm.sf(abs(b / s_)))} for nm, b, s_ in zip(names, cs['beta'], cse)]
    counts = used.groupby(['arm', 'status']).size().unstack(fill_value=0).reindex(index=arms, fill_value=0)
    events = {a: {('censored' if int(k) == 0 else cfg['causes'][str(k)]): int(v_) for k, v_ in counts.loc[a].items()} for a in arms}
    must = ['Cumulative incidence accounts for competing events; 1 - Kaplan-Meier treating them as censored would overstate the incidence of the event of interest.',
            'The subdistribution hazard ratio (Fine-Gray) describes the cumulative incidence; the cause-specific hazard ratio describes the event rate among those still event-free. Report both and do not read the Fine-Gray ratio as a rate ratio.',
            'Censoring must be independent of both event types; administrative censoring is assumed, informative dropout is not handled.']
    failing = []
    if not fg['converged']:
        failing.append({'item': 'fine_gray', 'reason': 'Newton iterations did not converge (cmprsk rule)'}); must.append('Fine-Gray did not converge; its estimates are not reportable.')
    if not cs['converged']:
        failing.append({'item': 'cause_specific_cox', 'reason': 'Cox iterations did not converge'})
    n_interest = int((status == ci_code).sum())
    if n_interest < 10 * len(names):
        must.append(f'Only {n_interest} events of interest for {len(names)} regression terms; hazard ratios are imprecise and may be biased.')
    for lm in landmarks:
        if lm['status'] == 'beyond_follow_up':
            must.append(f"Landmark {lm['time']} is beyond follow-up in arm {lm['arm']}; no estimate.")
    return {'analysis_type': 'competing_risks', 'primary': {'cause_of_interest': {'code': cfg['cause_of_interest'], 'label': cfg['causes'][cfg['cause_of_interest']]},
            'events_by_arm': events, 'gray_tests': gray, 'fine_gray': {'converged': fg['converged'], 'terms': subdistribution, 'variance': 'cmprsk crr sandwich including the censoring-distribution term'},
            'cause_specific_cox': {'converged': cs['converged'], 'terms': cause_specific, 'ties': 'Efron', 'other_causes': 'treated as censored'},
            'landmarks': landmarks, 'confidence_level': lv, 'study': cfg['study']},
            'curves': curves, 'must_mention': must, 'failing_items': failing,
            'limitations': ['One censoring group and no strata; time-varying covariates, frailty and recurrent events are not supported.',
                            'Pointwise log(-log) CIF intervals; no simultaneous bands. Wald intervals for hazard ratios.']}


def render(run, result, style):
    """Cumulative incidence step curves from saved values only."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    causes = list(dict.fromkeys((c['cause'], c['label']) for c in result['curves'])); arms = list(dict.fromkeys(c['arm'] for c in result['curves']))
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, axes = plt.subplots(1, len(causes), figsize=(4.2 * len(causes), 3.6), sharey=True, layout='constrained')
        for ax, (code, label) in zip(np.atleast_1d(axes), causes):
            for j, arm in enumerate(arms):
                c = next(c for c in result['curves'] if c['arm'] == arm and c['cause'] == code)
                ax.plot(c['time'], c['estimate'], '-', drawstyle='steps-post', color=pstyle.color(token, j), lw=1.5, label=arm)
            star = ' (event of interest)' if code == result['primary']['cause_of_interest']['code'] else ''
            ax.set(title=f'{label}{star}', xlabel=f"Time ({result['primary']['study']['time_unit']})", ylim=(0, 1))
        np.atleast_1d(axes)[0].set_ylabel('Cumulative incidence'); np.atleast_1d(axes)[0].legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'cumulative_incidence', token, 300)
    return ['cumulative_incidence']
