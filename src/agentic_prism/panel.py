"""Two-group comparisons across a declared panel of outcomes with family-wise or FDR control (0.14.0).

Every measured outcome must belong to one declared family; nothing is dropped
after seeing results. Per outcome: Welch t (independent groups) or paired t (same
units), on the original or declared log scale. Per family: Holm adjusted p values
with Bonferroni simultaneous intervals (family-wise error), or Benjamini-Hochberg
adjusted p values (false discovery rate; intervals stay unadjusted and say so).

0.14.1 opt-in `censoring` block: outcomes with values below (or above) a quantification
limit are kept censored and analysed on the log scale by censored normal ML with a
small-sample t correction (independent: two-group Tobit, common SD, SE x sqrt(N/(N-2)),
N - 2 df; paired: interval-censored paired log differences as censored.compute_paired,
SE x sqrt(k/(k-1)), k - 1 df). Without censoring these equal Student's pooled and paired
t tests. Outcomes without censored values keep the Welch or paired t test.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'panel_comparison'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'panel_comparison', 'input': None, 'source': None,
            'design': {'unit': None, 'rationale': None, 'pairing': None, 'scale': 'original', 'scale_rationale': None},
            'comparison': {'group_a': None, 'group_b': None, 'confidence_level': .95, 'alpha': .05},
            'families': None, 'multiplicity': {'procedure': None, 'rationale': None}, 'censoring': None,
            'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported panel configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'panel_comparison': raise ValueError('Unsupported panel schema')
    d = c['design']; q = c['comparison']
    for name, v in [('input', c['input']), ('source', c['source']), ('unit', d['unit']), ('rationale', d['rationale']),
                    ('group_a', q['group_a']), ('group_b', q['group_b']), ('multiplicity rationale', c['multiplicity']['rationale'])]:
        if not _text(v): raise ValueError('Declare ' + name)
    if q['group_a'] == q['group_b']: raise ValueError('Declare two different groups')
    if d['pairing'] not in ('independent', 'paired'): raise ValueError('Declare pairing: independent or paired (same units in both groups)')
    if d['scale'] not in ('original', 'log'): raise ValueError('scale must be original or log')
    if d['scale'] == 'log' and not _text(d['scale_rationale']): raise ValueError('Declare why the log scale is used (before results)')
    if c['multiplicity']['procedure'] not in ('holm', 'benjamini_hochberg'): raise ValueError('Declare the multiplicity procedure: holm or benjamini_hochberg')
    f = c['families']
    if not isinstance(f, dict) or not f or not all(_text(k) and isinstance(v, list) and v and all(_text(x) for x in v) for k, v in f.items()):
        raise ValueError('Declare families as {name: [outcomes]}; every measured outcome must be in exactly one family')
    allout = [x for v in f.values() for x in v]
    if len(allout) != len(set(allout)): raise ValueError('An outcome appears in more than one family')
    for k in ('confidence_level', 'alpha'):
        v = q[k]
        if type(v) not in (int, float) or not 0 < v < 1: raise ValueError('Invalid ' + k)
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    cz = c['censoring']
    if cz is None:
        c.pop('censoring')  # unused 0.14.1 key: resolved configurations stay as in 0.14.0
    else:
        if not isinstance(cz, dict) or set(cz) != {'method', 'limit_source'} or cz['method'] != 'censored_lognormal' or not _text(cz['limit_source']):
            raise ValueError('censoring must be {method: censored_lognormal, limit_source}')
        if d['scale'] != 'log': raise ValueError('Censored outcomes are analysed on the log scale; declare scale: log with a rationale')
    return c


def load_data(path, c):
    d = pd.read_csv(path, dtype={'unit_id': str, 'group': str, 'outcome': str, 'unit': str, 'censored': str})
    need = {'unit_id', 'group', 'outcome', 'value', 'unit'}
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    cens = c.get('censoring') is not None
    if cens and 'censored' not in d: raise ValueError('The censoring block needs a censored column (none, left, right)')
    if 'censored' in d and not cens: raise ValueError('A censored column needs the censoring block (method censored_lognormal, limit_source); values below a limit are never substituted')
    extra = set(d.columns) - need - {'exclude', 'exclusion_reason'} - ({'censored'} if cens else set())
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if 'exclude' in d:
        ex = d.exclude.astype(str).str.lower()
        if not ex.isin(['true', 'false']).all(): raise ValueError('exclude must be true or false')
        d['exclude'] = ex.eq('true')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded row needs an exclusion_reason')
    else: d['exclude'] = False
    v = pd.to_numeric(d.value, errors='coerce')
    if v.isna().any() or not np.isfinite(v).all():
        raise ValueError('value must be numeric; values below a quantification limit belong in censored-data (no substitution)')
    d['value'] = v
    if cens:
        d['censored'] = d.censored.str.strip().str.lower()
        if not d.censored.isin(['none', 'left', 'right']).all(): raise ValueError('censored must be none, left (below the stated limit) or right (above it); give the limit as value')
        if c['design']['pairing'] == 'paired' and (d.censored == 'right').any(): raise ValueError('Paired censored outcomes handle left censoring only')
    declared = {x for vals in c['families'].values() for x in vals}; measured = set(d.outcome)
    if measured - declared: raise ValueError('Measured outcomes missing from the declared families: ' + ', '.join(sorted(measured - declared)) + ' (no outcome may be dropped)')
    if declared - measured: raise ValueError('Declared outcomes not in the data: ' + ', '.join(sorted(declared - measured)))
    q = c['comparison']
    if set(d.group) != {q['group_a'], q['group_b']}: raise ValueError('Groups in the data do not match the declared groups')
    if d.duplicated(['unit_id', 'group', 'outcome']).any(): raise ValueError('One value per unit, group and outcome; average technical replicates first')
    if d.groupby('outcome').unit.nunique().max() > 1: raise ValueError('Each outcome must use a single unit')
    if c['design']['pairing'] == 'independent' and d.groupby('unit_id').group.nunique().max() > 1:
        raise ValueError('Independent groups share unit IDs; declare pairing: paired if the same units were measured in both groups')
    if c['design']['scale'] == 'log' and (d.loc[~d.exclude, 'value'] <= 0).any(): raise ValueError('Log scale requires positive values; no pseudocounts')
    return d


def holm(p):
    p = np.asarray(p, float); m = len(p); order = np.argsort(p); adj = np.empty(m); run = 0.
    for r, i in enumerate(order):
        run = max(run, min(1., (m - r)*p[i])); adj[i] = run
    return adj


def benjamini_hochberg(p):
    p = np.asarray(p, float); m = len(p); order = np.argsort(p)[::-1]; adj = np.empty(m); run = 1.
    for k, i in enumerate(order):
        rank = m - k; run = min(run, p[i]*m/rank); adj[i] = run
    return adj


def censored_outcome(part, c):
    from . import censored as cz
    q = c['comparison']; lv = q['confidence_level']; a_, b_ = q['group_a'], q['group_b']
    r = {'n_a': int((part.group == a_).sum()), 'n_b': int((part.group == b_).sum()), 'method': 'censored_lognormal_small_sample_t',
         'n_censored_a': int(((part.group == a_) & (part.censored != 'none')).sum()), 'n_censored_b': int(((part.group == b_) & (part.censored != 'none')).sum())}
    if c['design']['pairing'] == 'paired':
        a = part[part.group == a_].set_index('unit_id'); b = part[part.group == b_].set_index('unit_id')
        shared = sorted(set(a.index) & set(b.index)); r['n_pairs'] = len(shared); r['unpaired_units_ignored'] = sorted(set(a.index) ^ set(b.index))
        if r['unpaired_units_ignored']: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'incomplete_pairs_no_silent_deletion'}
        xc, yc = (b.loc[shared, 'censored'] == 'left').to_numpy(), (a.loc[shared, 'censored'] == 'left').to_numpy()
        if (~xc & ~yc).sum() < 3: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'fewer_than_3_pairs_with_both_values_detected'}
        lo, hi = cz.paired_intervals(b.loc[shared, 'value'], xc, a.loc[shared, 'value'], yc, 'log'); f = cz.interval_normal(lo, hi)
        k = f['n_used']; r['pairs_uninformative'] = int(len(shared) - k)
        if not f['converged'] or not np.isfinite(f['cov']).all() or f['cov'][0, 0] <= 0: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'likelihood_not_converged_or_singular_information'}
        est = float(f['theta'][0]); se_ml = float(np.sqrt(f['cov'][0, 0])); se = se_ml*np.sqrt(k/(k - 1)); df = k - 1
    else:
        det = part[part.censored == 'none']
        if min((det.group == a_).sum(), (det.group == b_).sum()) < 3: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'fewer_than_3_detected_values_in_a_group'}
        X = np.column_stack([np.ones(len(part)), (part.group == b_).to_numpy(float)])
        f = cz.tobit(np.log(part.value.to_numpy(float)), part.censored.to_numpy(), X); cov = f['cov']; N = len(part)
        if not f['converged'] or not np.isfinite(cov).all() or (np.diag(cov) <= 0).any(): return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'likelihood_not_converged_or_singular_information'}
        est = float(f['theta'][1]); se_ml = float(np.sqrt(cov[1, 1])); se = se_ml*np.sqrt(N/(N - 2)); df = N - 2
    t = est/se; h = stats.t.ppf((1+lv)/2, df)*se
    r.update(status='estimated', reportable=True, estimate=est, se=float(se), se_ml=se_ml, df=float(df), t=float(t),
             p_unadjusted=float(2*stats.t.sf(abs(t), df)), ci_unadjusted=[est-h, est+h])
    return r


def one_outcome(part, c):
    if c.get('censoring') is not None and 'censored' in part and (part.censored != 'none').any(): return censored_outcome(part, c)
    q = c['comparison']; lv = q['confidence_level']; log = c['design']['scale'] == 'log'
    a = part[part.group == q['group_a']].set_index('unit_id').value; b = part[part.group == q['group_b']].set_index('unit_id').value
    if log: a, b = np.log(a), np.log(b)
    r = {'n_a': int(len(a)), 'n_b': int(len(b))}
    if c['design']['pairing'] == 'paired':
        shared = sorted(set(a.index) & set(b.index)); r['n_pairs'] = len(shared); r['unpaired_units_ignored'] = sorted(set(a.index) ^ set(b.index))
        if r['unpaired_units_ignored']: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'incomplete_pairs_no_silent_deletion'}
        if len(shared) < 3: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'fewer_than_three_pairs'}
        dlt = (b[shared] - a[shared]).to_numpy(float); n = len(dlt); est = float(dlt.mean()); se = float(dlt.std(ddof=1)/np.sqrt(n)); df = n - 1
    else:
        if min(len(a), len(b)) < 3: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'fewer_than_three_units_per_group'}
        va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b); se = float(np.sqrt(va + vb)); est = float(b.mean() - a.mean())
        df = float((va + vb)**2/(va**2/(len(a)-1) + vb**2/(len(b)-1))) if se > 0 else 1.
    if se <= 0: return {**r, 'status': 'withheld', 'reportable': False, 'reason': 'zero_variance'}
    t = est/se; h = stats.t.ppf((1+lv)/2, df)*se
    r.update(status='estimated', reportable=True, estimate=est, se=se, df=float(df), t=float(t), p_unadjusted=float(2*stats.t.sf(abs(t), df)),
             ci_unadjusted=[est-h, est+h])
    return r


def compute(d, c):
    u = d[~d.exclude]; q = c['comparison']; lv = q['confidence_level']; alpha = q['alpha']; proc = c['multiplicity']['procedure']
    log = c['design']['scale'] == 'log'; rows = []
    for fam, outs in c['families'].items():
        res = []
        for o in outs:
            part = u[u.outcome == o]; r = one_outcome(part, c)
            r.update(family=fam, outcome=o, unit=str(part.unit.iloc[0]) if len(part) else None); res.append(r)
        ok = [r for r in res if r['reportable']]; m = len(ok)
        if m:
            p = [r['p_unadjusted'] for r in ok]; adj = holm(p) if proc == 'holm' else benjamini_hochberg(p)
            for r, pa in zip(ok, adj):
                r['p_adjusted'] = float(pa); r['significant_after_adjustment'] = bool(pa < alpha); r['family_size'] = m
                if proc == 'holm':
                    simult = 1 - (1-lv)/m; h = stats.t.ppf((1+simult)/2, r['df'])*r['se']; r['ci_simultaneous'] = [r['estimate']-h, r['estimate']+h]
                    r['simultaneous_level'] = simult
        for r in res:
            if log and r['reportable']:
                r['ratio_b_over_a'] = float(np.exp(r['estimate'])); r['ratio_ci_unadjusted'] = [float(np.exp(x)) for x in r['ci_unadjusted']]
                if 'ci_simultaneous' in r: r['ratio_ci_simultaneous'] = [float(np.exp(x)) for x in r['ci_simultaneous']]
        rows += res
    wh = [r for r in rows if not r['reportable']]
    notes = ['Independent unit: ' + c['design']['unit'] + '. Every measured outcome was assigned to a declared family before analysis; none was dropped.',
             ('Holm adjusted p values control the family-wise error within each family; Bonferroni simultaneous intervals match that control.' if proc == 'holm' else
              'Benjamini-Hochberg adjusted p values control the expected false discovery rate within each family; some declared discoveries may be false, and the intervals are unadjusted.') +
             ' Rationale: ' + c['multiplicity']['rationale'] + '.',
             f"Comparison {q['group_b']} versus {q['group_a']} ({c['design']['pairing']}); " + ("ratios of geometric means on the natural-log scale." if log else "differences of means on the original scale.") +
             ' Welch/paired t intervals assume roughly normal values (or paired differences) on that scale.',
             'Outcomes measured on the same units are correlated; Holm and Benjamini-Hochberg remain valid for the dependence studied in calibration, which is not every possible dependence.',
             'A non-significant adjusted p value does not show that an outcome is unaffected.']
    if c.get('censoring') is not None:
        cen = [r['outcome'] for r in rows if r.get('method') == 'censored_lognormal_small_sample_t']
        notes.append(('Outcomes with values beyond a quantification limit (' + ', '.join(cen) + ') were kept censored (' + c['censoring']['limit_source'] + ') and analysed by censored lognormal maximum likelihood with a small-sample t correction'
                      + (' (common SD in both groups)' if c['design']['pairing'] == 'independent' else ' on paired log differences; pairs with both values censored carry no information') + '; other outcomes used the ' + ('Welch' if c['design']['pairing'] == 'independent' else 'paired') + ' t test. Nothing was substituted.') if cen else
                     'No outcome had censored values; every outcome used the ' + ('Welch' if c['design']['pairing'] == 'independent' else 'paired') + ' t test.')
    if wh: notes.append('Withheld outcomes (not tested, still counted as measured): ' + '; '.join(f"{r['outcome']} ({r['reason']})" for r in wh) + '.')
    primary = {'status': 'estimated' if len(wh) < len(rows) else 'withheld', 'reportable': len(wh) < len(rows), 'procedure': proc, 'alpha': alpha,
               'confidence_level': lv, 'group_a': q['group_a'], 'group_b': q['group_b'], 'pairing': c['design']['pairing'], 'scale': c['design']['scale'],
               'outcomes': rows, 'n_significant_after_adjustment': sum(r.get('significant_after_adjustment', False) for r in rows)}
    return {'analysis_type': 'panel_comparison', 'primary': primary, 'must_mention': notes,
            'failing_items': [{'outcome': r['outcome'], 'reason': r['reason']} for r in wh],
            'limitations': [('Two groups per panel; no multi-arm, covariate-adjusted, censored or nonparametric per-outcome tests; no gatekeeping or hierarchical procedures.' if c.get('censoring') is None else
                             'Two groups per panel; censored outcomes only by censored lognormal ML (left censoring for paired designs); no multi-arm, covariate-adjusted or nonparametric per-outcome tests; no gatekeeping or hierarchical procedures.'),
                            'Families and procedure are user declarations; the module cannot verify they were fixed before the data were seen.',
                            'No Prism or software-equivalence claim beyond the recorded benchmarks.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    rows = [r for r in p['outcomes'] if r['reportable']]
    log = p['scale'] == 'log'
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.8, 1.2 + .3*max(len(rows), 1)), layout='constrained')
        for i, r in enumerate(rows):
            ci = r.get('ci_simultaneous', r['ci_unadjusted']); est = r['estimate']
            if log: ci, est = [np.exp(x) for x in ci], np.exp(est)
            col = pstyle.color(token, 1 if r['significant_after_adjustment'] else 0)
            ax.plot(ci, [i, i], color=col, lw=1.5); ax.plot(est, i, 'o', color=col, ms=4)
        ax.axvline(1 if log else 0, ls='--', color=pstyle.MUTED, lw=1)
        if log: ax.set_xscale('log')
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{r['family']}: {r['outcome']}" for r in rows])
        kind = 'simultaneous (Holm/Bonferroni)' if p['procedure'] == 'holm' else 'unadjusted (BH controls FDR, not intervals)'
        ax.set(xlabel=f"{p['group_b']} vs {p['group_a']} ({'ratio' if log else 'difference'}; {kind} intervals)", title='Panel comparison')
        pstyle.save(fig, out/'panel_comparison', token, 300)
    return ['panel_comparison']
