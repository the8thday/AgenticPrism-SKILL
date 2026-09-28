"""Analytical comparability / biosimilarity of one quality attribute across lots.

The lot is the unit: replicate measurements are averaged within a lot first.
The attribute's tier, the method and every margin or k are user declarations
with a source; nothing is inferred from the data or from a tier number.

- equivalence_tost: two one-sided t tests on the difference of lot means (test -
  reference) against +/-margin, i.e. the 100(1-2 alpha)% interval inside the
  margin. Welch (default) or pooled variance. The margin is either an absolute
  declared value or a declared multiple of the reference-lot SD (as in the
  withdrawn FDA 2017 draft, 1.5 sigma_R); an estimated margin carries sampling
  uncertainty that the nominal test size does not account for.
- quality_range: reference mean +/- k * reference SD; the fraction of test lots
  inside, optionally against a declared required fraction.
- descriptive: summaries only, no decision (for example for a tier-3 attribute).
"""
import numpy as np
from scipy import stats
from .cmc_common import merge, text, number, common, table

DEFAULTS = {'analysis_type': 'comparability', 'schema_version': 1, 'input': None, 'source': None,
    'attribute': {'name': None, 'unit': None, 'tier': None, 'tier_source': None, 'rationale': None},
    'design': {'lot_is_unit': None, 'independent_lots': None, 'lot_selection': None},
    'method': None,
    'equivalence': {'margin_type': None, 'margin': None, 'margin_source': None, 'alpha': .05, 'variance': 'welch'},
    'quality_range': {'k': None, 'k_source': None, 'required_fraction': None},
    'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    c = merge(raw, DEFAULTS); common(c, 'comparability')
    a, d, m = c['attribute'], c['design'], c['method']
    for k in ('name', 'unit', 'tier_source', 'rationale'):
        text(a[k], 'attribute.' + k)
    if a['tier'] not in (1, 2, 3) or isinstance(a['tier'], bool):
        raise ValueError('attribute.tier must be declared as 1, 2 or 3 from the risk assessment; it is never inferred')
    if d['lot_is_unit'] is not True or d['independent_lots'] is not True:
        raise ValueError('Declare lot_is_unit and independent_lots as literal true')
    text(d['lot_selection'], 'design.lot_selection')
    if m not in ('equivalence_tost', 'quality_range', 'descriptive'):
        raise ValueError('Declare method: equivalence_tost, quality_range or descriptive (not inferred from the tier)')
    e, q = c['equivalence'], c['quality_range']
    if m == 'equivalence_tost':
        if e['margin_type'] not in ('absolute', 'reference_sd_multiple'):
            raise ValueError('equivalence.margin_type must be absolute or reference_sd_multiple')
        number(e['margin'], 'equivalence.margin', positive=True); text(e['margin_source'], 'equivalence.margin_source')
        number(e['alpha'], 'equivalence.alpha', positive=True)
        if not e['alpha'] < .5:
            raise ValueError('alpha must be below .5')
        if e['variance'] not in ('welch', 'pooled'):
            raise ValueError('equivalence.variance must be welch or pooled')
    elif any(e[k] is not None for k in ('margin_type', 'margin', 'margin_source')):
        raise ValueError('equivalence settings apply only to equivalence_tost')
    if m == 'quality_range':
        number(q['k'], 'quality_range.k', positive=True); text(q['k_source'], 'quality_range.k_source')
        if q['required_fraction'] is not None:
            number(q['required_fraction'], 'quality_range.required_fraction', positive=True)
            if q['required_fraction'] > 1:
                raise ValueError('required_fraction cannot exceed 1')
    elif any(q[k] is not None for k in q):
        raise ValueError('quality_range settings apply only to quality_range')
    return c


def load_data(path, cfg):
    d = table(path, ['observation_id', 'lot_id', 'product', 'value', 'exclude'], ['value'])
    if set(d['product']) - {'reference', 'test'} or set(d['product']) != {'reference', 'test'}:
        raise ValueError('product must be reference or test, and both must be present')
    if set(d.exclude) - {'true', 'false'}:
        raise ValueError('exclude must be true/false')
    d['exclude'] = d.exclude == 'true'
    if 'exclusion_reason' not in d or (d.exclude & (d.exclusion_reason.str.strip() == '')).any():
        raise ValueError('Excluded rows need an exclusion_reason column entry')
    if (d.groupby('lot_id')['product'].nunique() > 1).any():
        raise ValueError('A lot cannot be both reference and test')
    return d


def compute(d, cfg):
    used = d[~d.exclude]
    lots = used.groupby(['product', 'lot_id']).value.agg(['mean', 'count']).reset_index()
    ref, tst = lots[lots['product'] == 'reference']['mean'].to_numpy(), lots[lots['product'] == 'test']['mean'].to_numpy()
    nr, nt = len(ref), len(tst)
    must, failing = [], []
    summary = {'reference': {'lots': nr, 'mean': float(ref.mean()), 'sd': float(ref.std(ddof=1)) if nr > 1 else None,
                             'min': float(ref.min()), 'max': float(ref.max())},
               'test': {'lots': nt, 'mean': float(tst.mean()), 'sd': float(tst.std(ddof=1)) if nt > 1 else None,
                        'min': float(tst.min()), 'max': float(tst.max())}}
    if (lots['count'] > 1).any():
        must.append('Replicate measurements were averaged within each lot before analysis; the lot is the unit.')
    if nr < 10 or nt < 6:
        must.append(f'Few lots ({nr} reference, {nt} test): the reference SD and any SD-based margin or range are imprecise.')
    m = cfg['method']; primary = {'method': m, 'tier': cfg['attribute']['tier'], 'tier_source': cfg['attribute']['tier_source']}
    detail = None
    if m == 'equivalence_tost':
        e = cfg['equivalence']
        if nr < 2 or nt < 2:
            raise ValueError('Equivalence testing needs at least two lots per product')
        sr, st_ = ref.std(ddof=1), tst.std(ddof=1)
        margin = e['margin'] if e['margin_type'] == 'absolute' else e['margin']*sr
        diff = float(tst.mean()-ref.mean())
        if e['variance'] == 'welch':
            se = float(np.sqrt(st_**2/nt+sr**2/nr)); df = float((st_**2/nt+sr**2/nr)**2/((st_**2/nt)**2/(nt-1)+(sr**2/nr)**2/(nr-1)))
        else:
            sp2 = ((nt-1)*st_**2+(nr-1)*sr**2)/(nt+nr-2); se = float(np.sqrt(sp2*(1/nt+1/nr))); df = float(nt+nr-2)
        t = stats.t.ppf(1-e['alpha'], df); ci = [diff-t*se, diff+t*se]
        p_low, p_up = float(stats.t.sf((diff+margin)/se, df)), float(stats.t.cdf((diff-margin)/se, df))
        equivalent = bool(-margin < ci[0] and ci[1] < margin)
        detail = {'difference_test_minus_reference': diff, 'se': se, 'df': df, 'variance': e['variance'],
                  'interval': ci, 'interval_level': 1-2*e['alpha'], 'margin': float(margin), 'margin_type': e['margin_type'],
                  'margin_declared': e['margin'], 'margin_source': e['margin_source'],
                  'p_lower': p_low, 'p_upper': p_up, 'tost_p': max(p_low, p_up), 'equivalent': equivalent}
        primary.update(equivalent=equivalent, difference=diff, interval=ci, margin=float(margin))
        if e['margin_type'] == 'reference_sd_multiple':
            must.append(f'The margin {margin:.4g} is {e["margin"]} x the reference-lot SD estimated from {nr} lots; its sampling uncertainty is not reflected in the nominal test size.')
        if not equivalent:
            failing.append({'reason': 'interval_outside_margin', 'interval': ci, 'margin': float(margin)})
    elif m == 'quality_range':
        q = cfg['quality_range']
        if nr < 3:
            raise ValueError('A quality range needs at least three reference lots')
        sr = ref.std(ddof=1); lo, hi = ref.mean()-q['k']*sr, ref.mean()+q['k']*sr
        inside = (tst >= lo) & (tst <= hi); frac = float(inside.mean())
        outside = lots[(lots['product'] == 'test')][~inside].lot_id.tolist()
        detail = {'range': [float(lo), float(hi)], 'k': q['k'], 'k_source': q['k_source'], 'test_lots_inside': int(inside.sum()),
                  'test_lots': nt, 'fraction_inside': frac, 'test_lots_outside': outside,
                  'reference_lots_inside': int(((ref >= lo) & (ref <= hi)).sum()), 'required_fraction': q['required_fraction']}
        primary.update(range=[float(lo), float(hi)], fraction_inside=frac)
        if q['required_fraction'] is not None:
            detail['meets_required_fraction'] = bool(frac >= q['required_fraction']-1e-12)
            primary['meets_required_fraction'] = detail['meets_required_fraction']
            if not detail['meets_required_fraction']:
                failing.append({'reason': 'fraction_inside_below_required', 'fraction_inside': frac, 'required': q['required_fraction']})
        if outside:
            must.append('Test lots outside the quality range: ' + ', '.join(map(str, outside)) + '.')
        must.append('A quality range compares individual lots with the reference spread; it is not an equivalence test of means and has no nominal error rate.')
    else:
        must.append('Descriptive comparison only; no decision rule was declared.')
    if cfg['attribute']['tier'] == 3 and m != 'descriptive':
        must.append('The attribute was declared tier 3 but a statistical decision method was requested; confirm this against the risk assessment.')
    return {'analysis_type': 'comparability', 'schema_version': 1, 'attribute': cfg['attribute'], 'summary': summary,
            'lots': lots.rename(columns={'mean': 'lot_mean', 'count': 'n_measurements'}).to_dict('records'),
            'primary': primary, 'detail': detail, 'must_mention': must, 'failing_items': failing,
            'limitations': ['Tier, method, margin and k are user declarations; the analysis does not establish them.',
                            'Lots are assumed independent and representative of each manufacturing process.',
                            'One attribute per analysis; totality-of-evidence judgements across attributes are outside the tool.']}
