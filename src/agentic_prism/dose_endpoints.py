"""Relative ECx/ICx from symmetric 4PL; profile the requested endpoint, not C50's interval."""
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import f as f_dist
from .dose_fit import _bounds, _evaluate, _noise_floor, summarize_dose
from .dose_schema import DOSE_UNITS


def effect_endpoints(group, fit, cfg):
    rows = []
    levels = cfg.get('effect_levels', {}).get('percentages', [])
    if not levels:
        return rows
    used = group.loc[~group.exclude]
    x, y = used.concentration_canonical.to_numpy(float), used.response.to_numpy(float)
    factor = DOSE_UNITS[fit['input_unit']]
    prefix = cfg['assay']['endpoint'][:2]
    for percent in levels:
        name = prefix + format(percent, 'g')
        row = {'curve_id': fit['curve_id'], 'sample_id': fit['sample_id'], 'experiment_id': fit['experiment_id'],
               'endpoint': name, 'percent': percent, 'basis': 'relative_low_to_high_dose_effect',
               'canonical_unit': fit['canonical_unit'], 'input_unit': fit['input_unit'],
               'status': 'withheld', 'reportable': False, 'estimate_canonical': None, 'estimate_input_unit': None,
               'ci_canonical': [None, None], 'ci_input_unit': [None, None],
               'ci_method': cfg['uncertainty']['method'], 'ci_level': cfg['uncertainty']['level'],
               'ci_status': 'not_computed', 'range_status': 'unknown', 'diagnostics': []}
        if fit['half_response_canonical'] is None:
            row['diagnostics'].append('parent_fit_failed'); rows.append(row); continue
        p = percent / 100.
        shift = float(np.log10(p/(1-p)))
        hill = abs(fit['hill_slope_signed'])
        centre = fit['log10_half_response_canonical'] + shift/hill
        estimate = float(10.**centre)
        low, high = fit['lowest_tested_positive_canonical'], fit['highest_tested_canonical']
        row.update(estimate_canonical=estimate, estimate_input_unit=estimate/factor,
                   log10_estimate_canonical=centre,
                   response_at_endpoint=fit['bottom']+(fit['top']-fit['bottom'])*(p if cfg['assay']['direction']=='increasing' else 1-p),
                   range_status='within_range' if low <= estimate <= high else 'below_range' if estimate < low else 'above_range')
        if not fit['reportable']:
            row['diagnostics'].append('parent_fit_not_reportable')
        if row['range_status'] != 'within_range':
            row['diagnostics'].append('effect_endpoint_outside_tested_range')
        if percent == 50:
            row.update(ci_canonical=list(fit['ci_canonical']), ci_status=fit['ci_status'])
        elif cfg['uncertainty']['method'] == 'none':
            row['ci_status'] = 'not_requested'
        elif fit['objective_sse'] <= _noise_floor(y, cfg):
            row['ci_status'] = 'noise_scale_not_estimable'
        elif not row['diagnostics']:
            mid_bounds, hb, _, _ = _bounds(cfg, [x[x > 0]])
            threshold = fit['objective_sse']*(1+float(f_dist.ppf(cfg['uncertainty']['level'], 1, fit['residual_df']))/fit['residual_df'])

            @lru_cache(maxsize=512)
            def profile(z):
                # Restrict the nuisance interval to the original C50 search domain.
                cuts = list(hb)
                for bound in mid_bounds:
                    if z != bound and shift/(z-bound) > 0:
                        h = np.log10(shift/(z-bound))
                        if hb[0] < h < hb[1]: cuts.append(float(h))
                cuts = sorted(cuts)
                def objective(lh):
                    l50 = z-shift/10.**lh
                    r = _evaluate([(x, y)], [l50], [lh], [[0]], cfg)[2]
                    return float(r@r)
                candidates = []
                for left, right in zip(cuts[:-1], cuts[1:]):
                    l50 = z-shift/10.**((left+right)/2)
                    if mid_bounds[0]-1e-12 <= l50 <= mid_bounds[1]+1e-12:
                        opt = minimize_scalar(objective, bounds=(left,right), method='bounded', options={'xatol':1e-10})
                        candidates += [(objective(left),left), (objective(right),right), (float(opt.fun),float(opt.x))]
                if not candidates: return float('inf'), True
                ss, lh = min(candidates)
                l50 = z-shift/10.**lh
                at_bound = min(abs(lh-hb[0]), abs(lh-hb[1]), abs(l50-mid_bounds[0]), abs(l50-mid_bounds[1])) < 1e-4
                return ss, at_bound

            limits = [mid_bounds[0]+min(shift/10.**hb[0],shift/10.**hb[1]),
                      mid_bounds[1]+max(shift/10.**hb[0],shift/10.**hb[1])]
            def root(bound):
                previous = centre
                # Local steps resolve short confidence intervals even when Hill's
                # lower search bound gives a very wide endpoint search domain.
                distance = abs(bound-centre)
                for step in np.geomspace(min(.01,distance), distance, 55):
                    z = centre + np.sign(bound-centre)*step
                    if profile(float(z))[0] > threshold:
                        return float(brentq(lambda v: profile(float(v))[0]-threshold, min(previous,z),max(previous,z),xtol=1e-9))
                    previous = z
                return None
            roots = [root(b) for b in limits]
            if None in roots:
                row['ci_status'] = 'profile_interval_open'
                row['diagnostics'].append('effect_endpoint_profile_interval_open')
            else:
                row.update(ci_canonical=[float(10.**v) for v in roots], ci_status='two_sided_profile_f')
                if any(profile(v)[1] for v in roots):
                    row['diagnostics'].append('effect_endpoint_profile_nuisance_at_boundary')
        else:
            row['ci_status'] = 'withheld_before_profiling'
        row['ci_input_unit'] = [None if v is None else v/factor for v in row['ci_canonical']]
        row['reportable'] = not row['diagnostics'] and row['ci_status'] in ('two_sided_profile_f','not_requested')
        row['status'] = 'estimated' if row['reportable'] else 'withheld'
        # Numeric candidates remain auditable but are removed by safe_copy in facts.
        rows.append(row)
    return rows


def summarize_endpoints(rows, cfg):
    summaries = []
    for name in dict.fromkeys(r['endpoint'] for r in rows):
        fs = [{**r, 'half_response_canonical': r['estimate_canonical'],
               'log10_half_response_canonical': r.get('log10_estimate_canonical')} for r in rows if r['endpoint']==name]
        for s in summarize_dose(fs, cfg):
            s['endpoint'] = name
            summaries.append(s)
    return summaries
