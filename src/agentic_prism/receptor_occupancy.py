"""Receptor occupancy (RO) from flow-cytometry MFI, with an optional exposure-occupancy fit (0.14.2).

Formats (declared before the data are seen), background MFI subtracted from every signal:
- free_baseline:   RO = 100 (1 - F_t / F_0); F = free-receptor signal (detection reagent competes
                   with the drug), F_0 the same subject's predose sample. Assumes total receptor
                   per cell does not change (no internalization or shedding).
- free_over_total: RO = 100 (1 - (F/T)_t / (F/T)_0); T = total-receptor signal (non-competing
                   reagent); corrects for receptor-level changes.
- bound_saturation: RO = 100 B / B_sat; B = drug-bound signal (anti-drug reagent) and B_sat the
                   same sample spiked ex vivo with saturating drug.
Values are never clipped to 0-100; out-of-range values are kept and counted. A sample whose
denominator is not above background is not computable. Summaries per group and timepoint use
the subject as the unit (one sample per subject, t interval). Optional exposure-occupancy fit:
RO = 100 C^h / (C^h + EC50^h) (Emax 100, baseline 0; h fixed at 1 or estimated), unweighted
least squares over post-dose samples; EC50 interval from a subject-cluster percentile bootstrap,
because repeated samples from one subject are not independent. 0.14.3 adds the opt-in
`interval: cluster_robust_t`: Bell-McCaffrey (CR2) cluster-robust variance of log EC50 (and log h)
at the least-squares estimate, with Satterthwaite degrees of freedom (Pustejovsky-Tipton
working-independence form), because the percentile cluster bootstrap undercovered with 6-8 subjects.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'receptor_occupancy'}
METHODS = {'free_baseline': ['free_mfi', 'free_background_mfi'],
           'free_over_total': ['free_mfi', 'free_background_mfi', 'total_mfi', 'total_background_mfi'],
           'bound_saturation': ['bound_mfi', 'bound_saturated_mfi', 'bound_background_mfi']}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'receptor_occupancy', 'input': None, 'source': None,
            'design': {'receptor': None, 'drug': None, 'cell_population': None, 'subject_unit': None, 'matrix': None},
            'method': {'name': None, 'rationale': None, 'reagent_evidence': None},
            'summary': {'confidence_level': .95},
            'exposure_response': None,
            'report': {'plot_style': 'prism_like'}}
ER_KEYS = {'concentration_unit', 'hill', 'bootstrap_replicates', 'seed', 'rationale', 'interval'}
INTERVALS = ('cluster_bootstrap', 'cluster_robust_t')


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported receptor-occupancy configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'receptor_occupancy': raise ValueError('Unsupported receptor-occupancy schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + name)
    m = c['method']
    if m['name'] not in METHODS: raise ValueError('Declare method.name: free_baseline, free_over_total or bound_saturation')
    if not _text(m['rationale']) or not _text(m['reagent_evidence']):
        raise ValueError('Declare method.rationale and reagent_evidence (competition of the free reagent with the drug, non-competition of the total reagent, or the saturating spike), each with its source')
    lv = c['summary']['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    er = c['exposure_response']
    if er is None:
        c.pop('exposure_response')
    else:
        if not isinstance(er, dict) or not {'concentration_unit', 'rationale'} <= set(er) <= ER_KEYS: raise ValueError('exposure_response needs concentration_unit and rationale (optional hill, interval, bootstrap_replicates, seed)')
        if er.get('interval', 'cluster_bootstrap') not in INTERVALS: raise ValueError('interval must be cluster_bootstrap or cluster_robust_t')
        er.setdefault('hill', 1.)
        if er.get('interval') == 'cluster_robust_t':
            if {'bootstrap_replicates', 'seed'} & set(er): raise ValueError('bootstrap_replicates and seed apply only to interval cluster_bootstrap')
        else:
            er.setdefault('bootstrap_replicates', 1000); er.setdefault('seed', 14220)
        if er['hill'] != 'estimate' and (type(er['hill']) not in (int, float) or er['hill'] <= 0): raise ValueError('hill must be a positive number (fixed) or "estimate"')
        if er.get('interval') != 'cluster_robust_t' and (type(er['bootstrap_replicates']) is not int or er['bootstrap_replicates'] < 200): raise ValueError('bootstrap_replicates must be an integer of at least 200')
        if not _text(er['concentration_unit']) or not _text(er['rationale']): raise ValueError('Declare concentration_unit and rationale')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'subject': str, 'timepoint': str, 'group': str})
    m = c['method']['name']; need = {'subject', 'timepoint', 'group', *METHODS[m]}
    if m != 'bound_saturation': need.add('baseline')
    if c.get('exposure_response'): need.add('drug_concentration')
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason', 'drug_concentration', 'baseline'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded sample needs an exclusion_reason')
    else: d['exclude'] = False
    for col in METHODS[m] + (['drug_concentration'] if 'drug_concentration' in d else []):
        v = pd.to_numeric(d[col], errors='coerce')
        if col != 'drug_concentration' and (v.isna().any() or not np.isfinite(v).all()): raise ValueError(f'{col} must be numeric (MFI)')
        d[col] = v
    if 'baseline' in d:
        b = d.baseline.astype(str).str.lower()
        if not b.isin(['true', 'false']).all(): raise ValueError('baseline must be true or false')
        d['baseline'] = b.eq('true')
    if d[['subject', 'timepoint', 'group']].isna().any().any(): raise ValueError('subject, timepoint and group are required on every row')
    if d.duplicated(['subject', 'timepoint']).any(): raise ValueError('One sample per subject and timepoint; average technical replicate tubes first by a declared rule')
    if d.groupby('subject').group.nunique().max() > 1: raise ValueError('Each subject belongs to one group')
    if m != 'bound_saturation':
        nb = d[~d.exclude].groupby('subject').baseline.sum()
        if (nb != 1).any(): raise ValueError('Each subject needs exactly one baseline (predose) sample for this method: ' + ', '.join(nb.index[nb != 1][:5]))
    if c.get('exposure_response') and (d.loc[~d.exclude & ~d.get('baseline', pd.Series(False, index=d.index)), 'drug_concentration'] < 0).any():
        raise ValueError('drug_concentration cannot be negative')
    return d


def occupancy(d, m):
    """Per-sample RO (percent) and status; never clipped."""
    out = d.copy(); out['ro_percent'] = np.nan; out['ro_status'] = 'computed'
    if m == 'bound_saturation':
        num = d.bound_mfi - d.bound_background_mfi; den = d.bound_saturated_mfi - d.bound_background_mfi
        ok = den > 0; out.loc[ok, 'ro_percent'] = 100*num[ok]/den[ok]; out.loc[~ok, 'ro_status'] = 'saturated_signal_not_above_background'
        return out
    f = d.free_mfi - d.free_background_mfi
    if m == 'free_over_total':
        t = d.total_mfi - d.total_background_mfi; sig = f/t.where(t > 0); out.loc[t <= 0, 'ro_status'] = 'total_signal_not_above_background'
    else:
        sig = f
    for s, part in d.groupby('subject'):
        base = part[part.baseline]
        b = float(sig[base.index[0]]) if len(base) else np.nan
        idx = part.index
        if not np.isfinite(b) or b <= 0:
            out.loc[idx, 'ro_status'] = 'baseline_signal_not_above_background'; continue
        ok = out.loc[idx, 'ro_status'] == 'computed'
        out.loc[idx[ok], 'ro_percent'] = 100*(1 - sig[idx[ok]]/b)
    out.loc[out.baseline & (out.ro_status == 'computed'), 'ro_status'] = 'baseline_reference'
    return out


def hill(c, ec50, h): return 100*c**h/(c**h + ec50**h)


def fit_exposure(cc, ro, h):
    """Least squares on log EC50 (and log h when estimated)."""
    pos = cc > 0; lc = np.log(cc[pos]) if pos.any() else np.array([0.])
    starts = [np.median(lc), np.quantile(lc, .25), np.quantile(lc, .75)]
    best = None
    for s0 in starts:
        z0 = [s0] + ([0.] if h == 'estimate' else [])
        def res(z): return hill(cc, np.exp(z[0]), np.exp(z[1]) if h == 'estimate' else h) - ro
        try: r = optimize.least_squares(res, z0, bounds=([lc.min() - 10] + ([-3.] if h == 'estimate' else []), [lc.max() + 10] + ([3.] if h == 'estimate' else [])), method='trf', xtol=1e-12, ftol=1e-12)
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best


def cluster_robust(cc, ro, subj, h, z, level):
    """CR2 sandwich for the log parameters at the least-squares estimate z; Satterthwaite df per parameter.

    J is the Jacobian of the fitted RO with respect to (log EC50[, log h]); bread (J'J)^-1; each subject's
    residuals are adjusted by (I - H_gg)^(-1/2) (Bell and McCaffrey 2002). The df for a parameter is
    (tr W)^2 / tr(W^2), W = U'U, with u_g = (I - H)_(.,g) A_g J_g (J'J)^-1 c (Pustejovsky and Tipton 2018, working
    independence). Returns None when a subject's I - H_gg is singular (that subject determines the fit).
    """
    ec = np.exp(z[0]); hh = np.exp(z[1]) if h == 'estimate' else float(h)
    pred = hill(cc, ec, hh); e = ro - pred; p = pred/100; cols = [-hh*100*p*(1 - p)]
    if h == 'estimate':
        lr = np.log(np.where(cc > 0, cc, 1.)/ec); cols.append(np.where(cc > 0, hh*lr*100*p*(1 - p), 0.))
    J = np.column_stack(cols); B = np.linalg.inv(J.T@J); H = J@B@J.T; M = np.eye(len(e)) - H; k = J.shape[1]
    meat = np.zeros((k, k)); U = [[] for _ in range(k)]; names = list(dict.fromkeys(subj))
    for g in names:
        ix = np.flatnonzero(subj == g); w, vec = np.linalg.eigh(np.eye(len(ix)) - H[np.ix_(ix, ix)])
        if w.min() <= 1e-10: return None
        A = (vec/np.sqrt(w))@vec.T; s = J[ix].T@A@e[ix]; meat += np.outer(s, s)
        for j in range(k): U[j].append(M[ix, :].T@A@J[ix]@B[:, j])
    V = B@meat@B; a = (1 - level)/2; out = []
    for j in range(k):
        W = np.column_stack(U[j]); W = W.T@W; df = float(np.trace(W)**2/np.trace(W@W)); se = float(np.sqrt(V[j, j])); q = float(stats.t.ppf(1 - a, df))
        out.append({'se_log': se, 'df': df, 'ci': [float(np.exp(z[j] - q*se)), float(np.exp(z[j] + q*se))]})
    return out


def compute(d, c):
    m = c['method']['name']; lv = c['summary']['confidence_level']; u = d[~d.exclude]
    s = occupancy(u, m); post = s[s.ro_status == 'computed'] if m != 'bound_saturation' else s[s.ro_status == 'computed']
    samples = [{'subject': r.subject, 'group': r.group, 'timepoint': r.timepoint, 'ro_percent': None if pd.isna(r.ro_percent) else float(r.ro_percent),
                'status': r.ro_status, **({'drug_concentration': None if pd.isna(r.drug_concentration) else float(r.drug_concentration)} if 'drug_concentration' in s else {})}
               for r in s.itertuples()]
    out_of_range = int(((post.ro_percent < 0) | (post.ro_percent > 100)).sum()); failing = []
    summaries = []
    for (g, tp), part in post.groupby(['group', 'timepoint'], sort=False):
        v = part.ro_percent.to_numpy(float); n = len(v); r = {'group': g, 'timepoint': tp, 'n_subjects': n, 'mean_ro_percent': float(v.mean()), 'median_ro_percent': float(np.median(v)),
                                                              'n_below_0': int((v < 0).sum()), 'n_above_100': int((v > 100).sum())}
        if n >= 3 and v.std(ddof=1) > 0:
            h = float(stats.t.ppf((1+lv)/2, n-1)*v.std(ddof=1)/np.sqrt(n)); r.update(status='estimated', ci=[r['mean_ro_percent']-h, r['mean_ro_percent']+h])
        else: r.update(status='withheld', reason='fewer_than_3_subjects' if n < 3 else 'zero_between_subject_variance')
        summaries.append(r)
    not_comp = s[~s.ro_status.isin(['computed', 'baseline_reference'])]
    for r in not_comp.itertuples(): failing.append({'subject': r.subject, 'timepoint': r.timepoint, 'reason': r.ro_status})
    primary = {'status': 'estimated' if len(post) else 'withheld', 'reportable': bool(len(post)), 'method': m, 'confidence_level': lv,
               'samples': samples, 'summaries': summaries, 'n_samples_out_of_0_100': out_of_range}
    er = c.get('exposure_response')
    if er:
        e = post[post.drug_concentration.notna()] if 'drug_concentration' in post else post.iloc[:0]
        subj = sorted(e.subject.unique()); res = {'concentration_unit': er['concentration_unit'], 'hill': er['hill'], 'n_samples': int(len(e)), 'n_subjects': len(subj)}
        if len(subj) < 6 or len(e) < 8:
            res.update(status='withheld', reason='fewer_than_6_subjects_or_8_samples'); failing.append({'reason': 'exposure_response_' + res['reason']})
        else:
            cc, ro = e.drug_concentration.to_numpy(float), e.ro_percent.to_numpy(float); f = fit_exposure(cc, ro, er['hill'])
            ec50 = float(np.exp(f.x[0])); hh = float(np.exp(f.x[1])) if er['hill'] == 'estimate' else float(er['hill'])
            if er.get('interval') == 'cluster_robust_t':
                _robust(res, e, cc, ro, f, ec50, hh, er, lv, failing)
            else: _bootstrap(res, e, subj, cc, ro, f, ec50, hh, er, lv, failing)
        primary['exposure_response'] = res
    notes = _notes(m, c, out_of_range, er, failing)
    return {'analysis_type': 'receptor_occupancy', 'primary': primary, 'must_mention': notes, 'failing_items': failing,
            'limitations': ['MFI-based formats only, with declared background; no receptor-number (molecules per cell) calibration, gating or compensation.',
                            'Exposure-occupancy fit fixes Emax at 100% and baseline at 0%; no PK model, hysteresis or time-course PK/PD.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks.']}


def _range_diag(cc, ro, ec50):
    cmin, cmax = float(cc[cc > 0].min()) if (cc > 0).any() else 0., float(cc.max()); diag = []
    if not cmin <= ec50 <= cmax: diag.append('EC50_outside_observed_concentrations')
    if ro.max() < 50 or ro.min() > 50: diag.append('half_occupancy_not_bracketed_by_data')
    return cmin, cmax, diag


def _robust(res, e, cc, ro, f, ec50, hh, er, lv, failing):
    cmin, cmax, diag = _range_diag(cc, ro, ec50); cr = cluster_robust(cc, ro, e.subject.to_numpy(), er['hill'], f.x, lv)
    if cr is None: diag.append('subject_leverage_singular')
    ok = not diag
    res.update(status='estimated' if ok else 'withheld', reportable=ok, interval='cluster_robust_t', ec50=ec50,
               ec50_ci_cluster_robust_t=cr[0]['ci'] if cr else None, ec50_ci_df=cr[0]['df'] if cr else None, ec50_se_log=cr[0]['se_log'] if cr else None,
               hill_estimate=hh if er['hill'] == 'estimate' else None, residual_ss=float(2*f.cost), diagnostics=diag, observed_concentration_range=[cmin, cmax])
    if er['hill'] == 'estimate' and cr: res.update(hill_ci_cluster_robust_t=cr[1]['ci'], hill_ci_df=cr[1]['df'])
    if not ok: failing.append({'reason': 'exposure_response_withheld:' + ','.join(diag)})


def _bootstrap(res, e, subj, cc, ro, f, ec50, hh, er, lv, failing):
    rng = np.random.default_rng(er['seed']); groups = {k: v.index.to_numpy() for k, v in e.groupby('subject')}; boot = []; failed = 0
    for _ in range(er['bootstrap_replicates']):
        pick = rng.choice(subj, len(subj), replace=True); idx = np.concatenate([groups[k] for k in pick]); sub = e.loc[idx]
        if sub.drug_concentration.nunique() < 3: failed += 1; continue
        b = fit_exposure(sub.drug_concentration.to_numpy(float), sub.ro_percent.to_numpy(float), er['hill'])
        if b is None: failed += 1; continue
        boot.append(b.x)
    boot = np.array(boot); a = (1 - lv)/2
    lo, hi = np.quantile(np.exp(boot[:, 0]), [a, 1-a]) if len(boot) else (None, None)
    cmin, cmax, diag = _range_diag(cc, ro, ec50)
    if failed > .05*er['bootstrap_replicates']: diag.append('more_than_5pct_bootstrap_failures')
    ok = not diag
    res.update(status='estimated' if ok else 'withheld', reportable=ok, ec50=ec50, ec50_ci_cluster_bootstrap=[float(lo), float(hi)] if len(boot) else None,
               hill_estimate=hh if er['hill'] == 'estimate' else None, residual_ss=float(2*f.cost), bootstrap_failed=failed, diagnostics=diag,
               observed_concentration_range=[cmin, cmax])
    if er['hill'] == 'estimate' and len(boot): res['hill_ci_cluster_bootstrap'] = [float(x) for x in np.quantile(np.exp(boot[:, 1]), [a, 1-a])]
    if not ok: failing.append({'reason': 'exposure_response_withheld:' + ','.join(diag)})


def _notes(m, c, out_of_range, er, failing):
    notes = ['Receptor occupancy format declared before the data: ' + m + ' (' + c['method']['rationale'] + '; reagent evidence: ' + c['method']['reagent_evidence'] + '). Background MFI was subtracted from every signal.',
             {'free_baseline': 'RO = 100 (1 - free/free at the subject\'s predose sample); this assumes receptor levels per cell did not change (internalization, shedding or cell-population shifts lower the free signal and are counted as occupancy).',
              'free_over_total': 'RO = 100 (1 - (free/total)/(free/total) predose); normalizing to total receptor corrects for receptor-level changes but requires a total-receptor reagent that does not compete with the drug.',
              'bound_saturation': 'RO = 100 bound/bound after ex vivo saturation of the same sample; this requires a truly saturating spike and an anti-drug reagent that detects bound drug equally in both aliquots.'}[m],
             'RO values were not clipped to 0-100%; ' + str(out_of_range) + ' post-dose sample(s) fell outside that range and are kept (assay noise around 0% or 100%).',
             'Subject is the unit for summaries (one sample per subject and timepoint); technical tubes are not subjects. RO is a pharmacodynamic marker of target engagement at the sampled site and cell population, not of tissue occupancy.']
    if er and er.get('interval') == 'cluster_robust_t':
        notes.append('Exposure-occupancy model RO = 100 C^h/(C^h + EC50^h) with Emax 100 and baseline 0 (' + er['rationale'] + '); EC50 interval from the cluster-robust (CR2, Bell-McCaffrey) variance of log EC50 with subject as the cluster and Satterthwaite degrees of freedom, so the interval reflects the number of subjects rather than the number of samples. The EC50 is an in vivo apparent value for this matrix and population, not an intrinsic KD.')
    elif er:
        notes.append('Exposure-occupancy model RO = 100 C^h/(C^h + EC50^h) with Emax 100 and baseline 0 (' + er['rationale'] + '); EC50 interval from a subject-cluster percentile bootstrap (' + str(er['bootstrap_replicates']) + ' replicates, seed ' + str(er['seed']) + '). The EC50 is an in vivo apparent value for this matrix and population, not an intrinsic KD.')
    if failing: notes.append('Not computable or withheld: ' + '; '.join(f"{x.get('subject', '')} {x.get('timepoint', '')} {x['reason']}".strip() for x in failing[:10]) + ('; ...' if len(failing) > 10 else '') + '.')
    return notes


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    smp = pd.DataFrame([x for x in p['samples'] if x['status'] == 'computed'])
    if smp.empty: return []
    er = p.get('exposure_response')
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        n = 2 if er and smp.get('drug_concentration') is not None and smp.drug_concentration.notna().any() else 1
        fig, axes = plt.subplots(1, n, figsize=(4.4*n, 3.3), layout='constrained', squeeze=False); ax = axes[0, 0]
        for k, (g, part) in enumerate(smp.groupby('group', sort=False)):
            tps = list(dict.fromkeys(smp.timepoint)); xs = [tps.index(t) for t in part.timepoint]
            ax.scatter(np.array(xs) + (k - .5)*.12, part.ro_percent, s=14, color=pstyle.color(token, k), label=g)
        ax.axhline(0, ls=':', color=pstyle.MUTED); ax.axhline(100, ls=':', color=pstyle.MUTED)
        tps = list(dict.fromkeys(smp.timepoint)); ax.set_xticks(range(len(tps))); ax.set_xticklabels(tps, rotation=30)
        ax.set(ylabel='Receptor occupancy (%)', title=f"RO ({p['method']}), not clipped"); ax.legend(frameon=False, fontsize=7)
        if n == 2:
            ax = axes[0, 1]; e = smp[smp.drug_concentration.notna() & (smp.drug_concentration > 0)]
            ax.scatter(e.drug_concentration, e.ro_percent, s=14, color=pstyle.color(token, 0)); ax.set_xscale('log')
            if er.get('ec50'):
                cc = np.geomspace(e.drug_concentration.min(), e.drug_concentration.max(), 200)
                ax.plot(cc, hill(cc, er['ec50'], er['hill_estimate'] or er['hill']), color=pstyle.color(token, 1), lw=1.2)
            ax.set(xlabel=f"Drug concentration ({er['concentration_unit']})", ylabel='RO (%)', title='Exposure-occupancy (observed range)')
        pstyle.save(fig, out/'receptor_occupancy', token, 300)
    return ['receptor_occupancy']
