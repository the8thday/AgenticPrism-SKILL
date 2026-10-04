"""qPCR relative quantification: delta-delta Cq with declared efficiencies (Pfaffl / geNorm form) (0.13.2).

Per biological replicate: technical replicate Cq values are averaged per gene; with amplification
factor E (2 = 100 %), log2 NRQ = -Cq_target * log2 E_target + mean_r(Cq_ref,r * log2 E_ref,r), the log of
the target quantity relative to the geometric mean of the reference genes (Hellemans et al. 2007;
with E = 2 and one reference this is -delta Cq of Livak & Schmittgen 2001). Each condition is compared
with the control condition on log2 NRQ: Welch t for independent biological replicates or paired t when
the same biological source is measured in every condition; fold change = 2^difference with its interval.
Holm within each target gene. Undetermined Cq values are never imputed.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'qpcr_relative'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'qpcr_relative', 'input': None, 'source': None,
            'assay': {'chemistry': None, 'instrument': None, 'cq_method': None, 'rationale': None},
            'genes': {'targets': None, 'references': None, 'reference_rationale': None},
            'efficiency': {'mode': None, 'values': None, 'source': None},
            'design': {'conditions': None, 'control_condition': None, 'pairing': None, 'pairing_rationale': None},
            'qc': {'cq_max': None, 'technical_sd_flag': .3, 'reference_shift_flag_cycles': .5, 'source': None},
            'uncertainty': {'level': .95}, 'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def _list(v, name, minimum=1):
    if not isinstance(v, list) or len(v) < minimum or len(set(v)) != len(v) or any(not isinstance(x, str) or not x.strip() for x in v):
        raise ValueError(f'Declare {name} as a list of at least {minimum} distinct names')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported qPCR configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'qpcr_relative' or c['schema_version'] != 1:
        raise ValueError('Unsupported qPCR schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    for k in ('chemistry', 'instrument', 'cq_method', 'rationale'):
        _text(c['assay'][k], 'assay.' + k)
    g = c['genes']; _list(g['targets'], 'genes.targets'); _list(g['references'], 'genes.references')
    if set(g['targets']) & set(g['references']):
        raise ValueError('A gene cannot be both target and reference')
    _text(g['reference_rationale'], 'genes.reference_rationale (evidence the references are stable across these conditions)')
    e = c['efficiency']
    if e['mode'] not in ('assume_100_percent', 'declared'):
        raise ValueError('efficiency.mode must be assume_100_percent or declared')
    _text(e['source'], 'efficiency.source')
    genes = g['targets'] + g['references']
    if e['mode'] == 'declared':
        if not isinstance(e['values'], dict) or set(e['values']) != set(genes) or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not 1.6 <= v <= 2.2 for v in e['values'].values()):
            raise ValueError('efficiency.values must give every gene an amplification factor between 1.6 and 2.2 (2 = 100 %)')
    elif e['values'] is not None:
        raise ValueError('efficiency.values apply only to mode declared')
    d = c['design']; _list(d['conditions'], 'design.conditions', 2)
    if d['control_condition'] not in d['conditions']:
        raise ValueError('control_condition must be one of the conditions')
    if d['pairing'] not in ('independent', 'paired'):
        raise ValueError('design.pairing must be independent or paired')
    _text(d['pairing_rationale'], 'design.pairing_rationale')
    q = c['qc']
    for k, lo, hi in (('cq_max', 20, 45), ('technical_sd_flag', 0, 2), ('reference_shift_flag_cycles', 0, 5)):
        v = q[k]
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not lo < v <= hi:
            raise ValueError(f'qc.{k} out of range')
    _text(q['source'], 'qc.source')
    lv = c['uncertainty']['level']
    if isinstance(lv, bool) or not isinstance(lv, (int, float)) or not .5 < lv < 1:
        raise ValueError('uncertainty.level must lie in (0.5, 1)')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'sample_id', 'condition', 'gene', 'cq'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    d['role'] = d.get('role', pd.Series('sample', index=d.index)).replace('', 'sample')
    if set(d.role) - {'sample', 'ntc'}:
        raise ValueError('role must be sample or ntc')
    undetermined = d.cq.str.strip().str.lower().isin(('', 'undetermined', 'nan', 'n/a'))
    d['cq'] = pd.to_numeric(d.cq.where(~undetermined, None), errors='raise')
    d['undetermined'] = undetermined
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    genes = set(cfg['genes']['targets'] + cfg['genes']['references']); s = d[d.role == 'sample']
    if set(s.gene) - genes or genes - set(s.gene):
        raise ValueError('Genes in the file must match the declared targets and references')
    if set(s.condition) != set(cfg['design']['conditions']):
        raise ValueError('Conditions in the file must match the declared conditions')
    if cfg['design']['pairing'] == 'independent' and (s.groupby('sample_id').condition.nunique() > 1).any():
        raise ValueError('A sample_id appears in several conditions; declare pairing: paired if the same biological source was measured in each')
    if cfg['design']['pairing'] == 'paired':
        if 'pair_id' not in d:
            raise ValueError('Paired designs need a pair_id column naming the biological source')
        if (s.groupby(['pair_id', 'condition']).sample_id.nunique() > 1).any():
            raise ValueError('Each pair_id must have one biological sample per condition')
    return d


def _log2e(cfg, gene):
    return 1. if cfg['efficiency']['mode'] == 'assume_100_percent' else float(np.log2(cfg['efficiency']['values'][gene]))


def compute(d, cfg):
    g, des, qc = cfg['genes'], cfg['design'], cfg['qc']; lv = cfg['uncertainty']['level']; must, diag, failing = [], [], []
    s = d[(d.role == 'sample') & ~d.exclude]
    ntc = d[(d.role == 'ntc') & ~d.exclude & ~d.undetermined & (d.cq < qc['cq_max'])]
    if len(ntc):
        diag.append('ntc_amplification'); must.append(f"No-template controls amplified for {sorted(set(ntc.gene))} (Cq < {qc['cq_max']}); contamination or primer dimers must be ruled out before interpretation.")
    keys = ['sample_id', 'condition'] + (['pair_id'] if des['pairing'] == 'paired' else [])
    wells = []
    for key, w in s.groupby(keys + ['gene'], sort=False):
        info = dict(zip(keys + ['gene'], key)); v = w.cq.dropna().to_numpy(float)
        row = {**info, 'technical_n': int(len(w)), 'undetermined_n': int(w.undetermined.sum()), 'mean_cq': float(v.mean()) if len(v) else None,
               'technical_sd': float(v.std(ddof=1)) if len(v) > 1 else None}
        row['status'] = 'no_determined_cq' if not len(v) else 'above_cq_max' if row['mean_cq'] >= qc['cq_max'] else 'ok'
        if row['technical_sd'] is not None and row['technical_sd'] > qc['technical_sd_flag']:
            row['flag'] = 'technical_sd_above_flag'
        wells.append(row)
    wt = pd.DataFrame(wells)
    if (wt.status != 'ok').any():
        bad = wt[wt.status != 'ok']; diag.append('cq_not_quantifiable')
        must.append(f'{len(bad)} sample-gene combination(s) have no determined Cq or a mean Cq at or above {qc["cq_max"]}; they are not imputed and those samples drop out of that gene.')
    partial = wt[(wt.undetermined_n > 0) & (wt.status == 'ok')]
    if len(partial):
        must.append(f'{len(partial)} sample-gene combination(s) mix determined and undetermined technical replicates; the mean uses the determined wells only, which can bias Cq low near the detection limit.')
    if 'flag' in wt and wt.flag.notna().any():
        must.append(f"{int(wt.flag.notna().sum())} sample-gene combination(s) have technical SD above {qc['technical_sd_flag']} cycles.")
    pivot = wt[wt.status == 'ok'].pivot_table(index=keys, columns='gene', values='mean_cq').reset_index()
    refs = g['references']
    have_refs = pivot[refs].notna().all(axis=1) if set(refs) <= set(pivot) else pd.Series(False, index=pivot.index)
    pivot['ref_log2'] = np.where(have_refs, np.mean([pivot[r] * _log2e(cfg, r) for r in refs], axis=0) if set(refs) <= set(pivot) else np.nan, np.nan)
    # reference stability: reference log2 quantity by condition
    ref_shift = []
    for cond in des['conditions']:
        if cond == des['control_condition']:
            continue
        a = pivot.loc[pivot.condition == des['control_condition'], 'ref_log2'].dropna(); b = pivot.loc[pivot.condition == cond, 'ref_log2'].dropna()
        if len(a) > 1 and len(b) > 1:
            shift = float(b.mean() - a.mean()); p = float(stats.ttest_ind(b, a, equal_var=False).pvalue)
            ref_shift.append({'condition': cond, 'reference_cycle_shift': shift, 'p_value': p})
            if abs(shift) > qc['reference_shift_flag_cycles']:
                diag.append('reference_shift'); must.append(f"Reference genes shift by {shift:+.2f} cycles in {cond} versus {des['control_condition']}; a treatment effect on the references would bias every fold change.")
    results = []
    for tgt in g['targets']:
        sub = pivot[[c for c in keys] + ([tgt] if tgt in pivot else []) + ['ref_log2']].copy()
        if tgt not in sub:
            continue
        sub['log2_nrq'] = -sub[tgt] * _log2e(cfg, tgt) + sub['ref_log2']; sub = sub.dropna(subset=['log2_nrq'])
        ctrl = sub[sub.condition == des['control_condition']]; rows = []
        for cond in des['conditions']:
            if cond == des['control_condition']:
                continue
            tr = sub[sub.condition == cond]; row = {'target': tgt, 'condition': cond, 'control': des['control_condition'], 'n_condition': int(len(tr)), 'n_control': int(len(ctrl)),
                                                   'log2_fold_change': None, 'fold_change': None, 'fold_change_ci': [None, None], 'p_value': None}
            if des['pairing'] == 'paired':
                m = tr.merge(ctrl, on='pair_id', suffixes=('', '_c')); diff = (m.log2_nrq - m.log2_nrq_c).to_numpy(float); row['n_pairs'] = int(len(diff))
                if len(diff) < 2:
                    row['status'] = 'withheld_fewer_than_two_pairs'
                else:
                    est, se, df = float(diff.mean()), float(diff.std(ddof=1) / np.sqrt(len(diff))), len(diff) - 1; method = 'paired t on log2 NRQ'
            else:
                x, y = ctrl.log2_nrq.to_numpy(float), tr.log2_nrq.to_numpy(float)
                if len(x) < 2 or len(y) < 2:
                    row['status'] = 'withheld_fewer_than_two_biological_replicates'
                else:
                    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y); se = float(np.sqrt(vx + vy)); est = float(y.mean() - x.mean())
                    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)); method = 'Welch t on log2 NRQ'
            if 'status' not in row:
                q = stats.t.ppf((1 + lv) / 2, df)
                row.update(status='estimated', log2_fold_change=est, standard_error=se, df=float(df), log2_ci=[est - q * se, est + q * se],
                           fold_change=float(2 ** est), fold_change_ci=[float(2 ** (est - q * se)), float(2 ** (est + q * se))],
                           p_value=float(2 * stats.t.sf(abs(est / se), df)), method=method)
            rows.append(row)
        est = [r for r in rows if r['status'] == 'estimated']; order = np.argsort([r['p_value'] for r in est]); running = 0.
        for rank, i in enumerate(order):
            running = max(running, min(1., (len(est) - rank) * est[i]['p_value'])); est[i]['p_holm'] = running
        results += rows
    failing = [r for r in results if r['status'] != 'estimated']
    must.insert(0, f"Fold changes are relative to {des['control_condition']} and normalised to {', '.join(refs)} ({'geometric mean' if len(refs) > 1 else 'single reference'}); "
                f"statistics use biological replicates ({des['pairing']}), never technical wells.")
    if cfg['efficiency']['mode'] == 'assume_100_percent':
        must.append('Amplification efficiency was assumed to be 100 % for every gene; unequal real efficiencies bias fold changes, increasingly with the Cq difference.')
    else:
        must.append('Declared efficiencies are treated as known; their uncertainty is not propagated into the intervals.')
    if len(refs) == 1:
        must.append('A single reference gene was used; MIQE recommends at least two validated references.')
    return {'analysis_type': 'qpcr_relative', 'primary': {'comparisons': results, 'reference_stability': ref_shift, 'replicate_values': pivot.to_dict('records'),
            'wells': wells, 'genes': g, 'efficiency': cfg['efficiency'], 'design': des, 'level': lv},
            'must_mention': must, 'failing_items': failing,
            'limitations': ['Relative quantification only; absolute copy numbers need a standard curve (standard_curve).',
                            'Holm within each target gene; multiplicity across genes is not adjusted.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    rows = [r for r in result['primary']['comparisons'] if r['status'] == 'estimated']
    if not rows:
        return []
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.6, .5 + .45 * len(rows)), layout='constrained')
        for i, r in enumerate(rows):
            ax.plot(r['log2_ci'], [i, i], '-', color=pstyle.color(token, 0), lw=1.6); ax.plot(r['log2_fold_change'], i, 'o', color=pstyle.color(token, 0), ms=5)
        ax.axvline(0, color=pstyle.MUTED, lw=.8, ls=(0, (3, 2)))
        ax.set(yticks=range(len(rows)), yticklabels=[f"{r['target']}: {r['condition']}" for r in rows], xlabel=f"log2 fold change vs {rows[0]['control']} (CI)")
        ax.invert_yaxis(); pstyle.save(fig, out/'qpcr_fold_change', token, 300)
    return ['qpcr_fold_change']
