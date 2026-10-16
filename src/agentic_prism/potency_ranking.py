"""Many candidates against one benchmark tested in the same experiments (0.14.0).

Each candidate's log10 EC50 is compared with the reference's log10 EC50 from the
same independent experiment; the per-experiment differences are summarized with a
Student interval (n experiments - 1 df). Relative potency RP = 10^-(mean difference)
(> 1 means more potent than the reference for an EC50/IC50-type quantity). The
family of m candidate-versus-reference comparisons is controlled together: Holm
adjusted p values and Bonferroni simultaneous intervals at the declared level. An
optional declared fold threshold classifies candidates only when the simultaneous
interval lies wholly beyond it. Candidates are not compared with each other.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'potency_ranking'}
STATUSES = ('estimated', 'above_range', 'below_range', 'failed')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'potency_ranking', 'input': None, 'source': None,
            'design': {'unit': None, 'rationale': None, 'quantity': None, 'quantity_unit': None, 'estimate_source': None},
            'comparison': {'reference': None, 'confidence_level': .95, 'fold_threshold': None, 'fold_threshold_source': None,
                           'minimum_experiments': 3},
            'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported potency-ranking configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'potency_ranking': raise ValueError('Unsupported potency-ranking schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items(), ('reference', c['comparison']['reference'])]:
        if not _text(v): raise ValueError('Declare ' + name)
    q = c['comparison']; lv = q['confidence_level']
    if type(lv) not in (int, float) or not .5 < lv < 1: raise ValueError('Invalid confidence level')
    if q['fold_threshold'] is not None:
        if type(q['fold_threshold']) not in (int, float) or q['fold_threshold'] <= 1: raise ValueError('fold_threshold must exceed 1')
        if not _text(q['fold_threshold_source']): raise ValueError('Declare fold_threshold_source (set before results)')
    elif q['fold_threshold_source'] is not None: raise ValueError('fold_threshold_source needs a fold_threshold')
    if type(q['minimum_experiments']) is not int or q['minimum_experiments'] < 3: raise ValueError('minimum_experiments must be an integer of at least 3')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'candidate': str, 'experiment_id': str, 'status': str, 'unit': str})
    need = {'candidate', 'experiment_id', 'ec50', 'unit', 'status'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'source_run', 'note'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if d[['candidate', 'experiment_id', 'unit', 'status']].isna().any().any(): raise ValueError('candidate, experiment_id, unit and status are required on every row')
    d['status'] = d.status.str.strip().str.lower()
    if not d.status.isin(STATUSES).all(): raise ValueError('status must be one of ' + ', '.join(STATUSES) + ' (keep non-estimable runs; never drop them)')
    if d.duplicated(['candidate', 'experiment_id']).any(): raise ValueError('One estimate per candidate and experiment; summarize technical replicates within the experiment first')
    if d.unit.nunique() != 1 or d.unit.iloc[0] != c['design']['quantity_unit']: raise ValueError('All estimates must share the declared quantity_unit')
    est = d.status == 'estimated'; v = pd.to_numeric(d.ec50, errors='coerce')
    if (est & (v.isna() | ~(v > 0))).any(): raise ValueError('Estimated EC50 values must be positive numbers')
    d['ec50'] = v
    if c['comparison']['reference'] not in set(d.candidate): raise ValueError('The reference candidate is not in the data')
    if d.candidate.nunique() < 2: raise ValueError('At least one candidate besides the reference is required')
    return d


def holm(p):
    p = np.asarray(p, float); m = len(p); order = np.argsort(p); adj = np.empty(m); run = 0.
    for rank, i in enumerate(order):
        run = max(run, min(1., (m - rank)*p[i])); adj[i] = run
    return adj


def compute(d, c):
    q = c['comparison']; ref = q['reference']; lv = q['confidence_level']; mmin = q['minimum_experiments']
    refrows = d[d.candidate == ref].set_index('experiment_id')
    ref_ok = set(refrows.index[refrows.status == 'estimated']); ref_failed = sorted(set(refrows.index) - ref_ok)
    cands = [x for x in dict.fromkeys(d.candidate) if x != ref]; rows = []
    for cand in cands:
        part = d[d.candidate == cand].set_index('experiment_id'); shared = [e for e in part.index if e in ref_ok]
        status = part.loc[shared, 'status'] if shared else pd.Series(dtype=str)
        r = {'candidate': cand, 'n_experiments_with_reference': len(shared), 'experiments': shared,
             'non_estimable_runs': {k: int(v) for k, v in status[status != 'estimated'].value_counts().items()},
             'runs_without_valid_reference': sorted(set(part.index) - set(shared))}
        if (status != 'estimated').any():
            r.update(status='withheld', reportable=False, reason='non_estimable_runs_present')
        elif len(shared) < mmin:
            r.update(status='withheld', reportable=False, reason=f'fewer_than_{mmin}_experiments_with_reference')
        else:
            diff = np.log10(part.loc[shared, 'ec50'].to_numpy(float)) - np.log10(refrows.loc[shared, 'ec50'].to_numpy(float))
            n = len(diff); mean = float(diff.mean()); se = float(diff.std(ddof=1)/np.sqrt(n))
            if se <= 0:
                r.update(status='withheld', reportable=False, reason='zero_between_experiment_variance')
            else:
                t = mean/se
                r.update(status='estimated', reportable=True, log10_ec50_difference=mean, se=se, df=n-1, t=float(t),
                         p_unadjusted=float(2*stats.t.sf(abs(t), n-1)), per_experiment_differences=diff.tolist(),
                         relative_potency=float(10**(-mean)))
        rows.append(r)
    est = [r for r in rows if r['reportable']]; m = len(est)
    if m:
        adj = holm([r['p_unadjusted'] for r in est]); simult = 1 - (1-lv)/m
        for r, p in zip(est, adj):
            h = stats.t.ppf((1+simult)/2, r['df'])*r['se']; hp = stats.t.ppf((1+lv)/2, r['df'])*r['se']
            r.update(p_holm=float(p), simultaneous_level=simult,
                     relative_potency_simultaneous_ci=[float(10**(-(r['log10_ec50_difference']+h))), float(10**(-(r['log10_ec50_difference']-h)))],
                     relative_potency_unadjusted_ci=[float(10**(-(r['log10_ec50_difference']+hp))), float(10**(-(r['log10_ec50_difference']-hp)))])
            f = q['fold_threshold']
            if f:
                lo, hi = r['relative_potency_simultaneous_ci']
                r['tier'] = (f'at_least_{f:g}x_more_potent' if lo >= f else f'at_least_{f:g}x_less_potent' if hi <= 1/f
                             else 'not_resolved_against_threshold')
    ranked = sorted(est, key=lambda r: -r['relative_potency'])
    for i, r in enumerate(ranked, 1): r['point_rank'] = i
    notes = ['Independent experiment is the unit; each candidate is compared with the reference ' + ref + ' from the same experiments (paired log10 differences). '
             'Relative potency > 1 means a lower ' + c['design']['quantity'] + ' than the reference.',
             f'The {m} candidate-versus-reference comparisons form one family: Holm adjusted p values and Bonferroni simultaneous {lv:.0%} intervals. Unadjusted intervals are shown for context only; do not pick the candidate with the smallest unadjusted p value.',
             'Point ranks are descriptive. Candidates were not compared with each other; overlapping intervals do not establish an order between candidates.',
             'Estimates come from the declared source runs (' + c['design']['estimate_source'] + '); their own fit quality and withholding apply first.']
    if q['fold_threshold']: notes.append(f"Fold threshold {q['fold_threshold']:g} ({q['fold_threshold_source']}) is applied to simultaneous intervals only.")
    if ref_failed: notes.append('Experiments where the reference was not estimated were excluded for every candidate: ' + ', '.join(ref_failed) + '.')
    wh = [r for r in rows if not r['reportable']]
    if wh: notes.append('Withheld candidates (non-estimable runs are informative, e.g. EC50 above the top dose, so they are not dropped): ' +
                        '; '.join(f"{r['candidate']} ({r['reason']})" for r in wh) + '.')
    primary = {'status': 'estimated' if m else 'withheld', 'reportable': m > 0, 'reference': ref, 'confidence_level': lv,
               'n_candidates': len(cands), 'n_compared': m, 'reference_failed_experiments': ref_failed, 'candidates': rows,
               'quantity': c['design']['quantity'], 'quantity_unit': c['design']['quantity_unit']}
    return {'analysis_type': 'potency_ranking', 'primary': primary, 'must_mention': notes,
            'failing_items': [{'candidate': r['candidate'], 'reason': r['reason']} for r in wh] + ([] if m else [{'reason': 'no_candidate_compared'}]),
            'limitations': ['Paired Student intervals assume log differences are roughly normal and experiments independent; no mixed model, Dunnett or between-candidate comparisons.',
                            'EC50 estimates are taken as given; their within-experiment uncertainty is reflected only through between-experiment variation.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    est = sorted([r for r in p['candidates'] if r['reportable']], key=lambda r: r['relative_potency'])
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.5, 1.2 + .32*max(len(est), 1)), layout='constrained')
        for i, r in enumerate(est):
            lo, hi = r['relative_potency_simultaneous_ci']
            ax.plot([lo, hi], [i, i], color=pstyle.color(token, 0), lw=1.6)
            ax.plot(r['relative_potency'], i, 'o', color=pstyle.color(token, 0), ms=4)
        ax.axvline(1, ls='--', color=pstyle.MUTED, lw=1)
        ax.set_xscale('log'); ax.set_yticks(range(len(est))); ax.set_yticklabels([r['candidate'] for r in est])
        ax.set(xlabel=f"Relative potency vs {p['reference']} (simultaneous {p['confidence_level']:.0%} intervals)", title='Candidates vs reference')
        pstyle.save(fig, out/'potency_ranking', token, 300)
    return ['potency_ranking']
