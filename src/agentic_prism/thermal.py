"""Thermal unfolding: apparent Tm from nanoDSF / DSF / CD melting curves (0.13.1).

Model with k declared transitions (two-state van 't Hoff form, T in kelvin):
    f_i(T) = 1 / (1 + exp(-(dH_i / R) * (1/Tm_i - 1/T)))
    y(T)   = a + b*(T - T0) + sum_i A_i * f_i(T) + c*(T - T0)*f_k(T)
a + b(T - T0) is the native baseline; the last term lets the unfolded baseline have its
own slope (exactly the sloping-baseline two-state model when k = 1). Only Tm_1, the
log gaps Tm_{i+1} - Tm_i and log10 dH_i are optimized; the k + 3 linear coefficients are
solved by least squares for each proposal. Each Tm gets a profile-F interval with all
other parameters re-optimized. Thermal unfolding of antibodies is usually irreversible
and scan-rate dependent, so Tm and dH are apparent quantities.

A Savitzky-Golay first-derivative extremum is reported alongside as a descriptive
comparison with instrument software; it has no interval.
"""
from copy import deepcopy
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.signal import savgol_filter
from scipy.special import expit
from scipy.stats import f as f_dist, t as t_dist

TYPES = {'thermal_unfolding'}
R_KJ = 8.314462618e-3  # kJ / (mol K)
TECHNIQUES = ('nanodsf_ratio', 'intrinsic_fluorescence', 'extrinsic_dye', 'circular_dichroism', 'static_light_scattering')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'thermal_unfolding', 'input': None, 'source': None,
            'assay': {'technique': None, 'readout': None, 'scan_rate_c_per_min': None, 'buffer': None,
                      'reversibility': None, 'rationale': None},
            'model': {'transitions': None, 'transitions_source': None, 'fit_range_c': None, 'dh_bounds_kj_mol': [50., 5000.]},
            'derivative': {'window_c': 3.0, 'polyorder': 2},
            'replicates': {'independent_unit': 'none', 'independence_source': None},
            'comparisons': [], 'uncertainty': {'level': .95}, 'report': {'plot_style': 'prism_like'}}
MIN_BASELINE_POINTS = 5
STRUCTURE_ALPHA = .01  # k transitions must beat k - 1 at this level, or every Tm is withheld


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported thermal-unfolding configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'thermal_unfolding' or c['schema_version'] != 1:
        raise ValueError('Unsupported thermal-unfolding schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    a, m = c['assay'], c['model']
    if a['technique'] not in TECHNIQUES:
        raise ValueError('assay.technique must be one of ' + ', '.join(TECHNIQUES) + ' (DSC heat capacity needs a different model)')
    if a['technique'] == 'static_light_scattering':
        raise ValueError('Light scattering reports aggregation onset (Tagg), not unfolding; a sigmoid Tm here would be mislabelled')
    for k in ('readout', 'buffer', 'rationale'):
        _text(a[k], 'assay.' + k)
    if isinstance(a['scan_rate_c_per_min'], bool) or not isinstance(a['scan_rate_c_per_min'], (int, float)) or not a['scan_rate_c_per_min'] > 0:
        raise ValueError('Declare assay.scan_rate_c_per_min: apparent Tm depends on the heating rate')
    if a['reversibility'] not in ('irreversible_or_unknown', 'reversible_demonstrated'):
        raise ValueError('Declare assay.reversibility: irreversible_or_unknown or reversible_demonstrated (cooling/reheating evidence)')
    if m['transitions'] not in (1, 2, 3) or type(m['transitions']) is not int:
        raise ValueError('Declare model.transitions (1-3) from the molecule and prior data, before fitting')
    _text(m['transitions_source'], 'model.transitions_source')
    r = m['fit_range_c']
    if r is not None and (not isinstance(r, list) or len(r) != 2 or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in r) or not r[0] < r[1]):
        raise ValueError('model.fit_range_c must be null or [low, high] in degrees C, declared before fitting')
    b = m['dh_bounds_kj_mol']
    if not isinstance(b, list) or len(b) != 2 or not 0 < b[0] < b[1]:
        raise ValueError('dh_bounds_kj_mol must be [low, high] > 0')
    dv = c['derivative']
    if isinstance(dv['window_c'], bool) or not isinstance(dv['window_c'], (int, float)) or not dv['window_c'] > 0 or dv['polyorder'] not in (2, 3):
        raise ValueError('derivative.window_c > 0 and polyorder 2 or 3')
    rep = c['replicates']
    if rep['independent_unit'] not in ('none', 'replicate_id'):
        raise ValueError('replicates.independent_unit must be none or replicate_id')
    if rep['independent_unit'] == 'replicate_id':
        _text(rep['independence_source'], 'replicates.independence_source (separate preparations or runs, not capillaries of one mix)')
    if not isinstance(c['comparisons'], list):
        raise ValueError('comparisons must be a list')
    ids = set()
    for item in c['comparisons']:
        if not isinstance(item, dict) or set(item) != {'id', 'reference_sample', 'test_sample'} or item['id'] in ids \
                or item['reference_sample'] == item['test_sample']:
            raise ValueError('Each comparison needs a unique id, reference_sample and a different test_sample')
        ids.add(item['id'])
    if c['comparisons'] and rep['independent_unit'] != 'replicate_id':
        raise ValueError('Tm comparisons need declared independent replicates')
    lv = c['uncertainty']['level']
    if isinstance(lv, bool) or not isinstance(lv, (int, float)) or not .5 < lv < 1:
        raise ValueError('uncertainty.level must lie in (0.5, 1)')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'sample_id', 'replicate_id', 'curve_id', 'temperature_c', 'signal'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    for k in ('sample_id', 'replicate_id', 'curve_id'):
        if d[k].str.strip().eq('').any():
            raise ValueError(f'Missing {k}')
    for k in ('temperature_c', 'signal'):
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'Nonfinite {k}')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    for cid, g in d.groupby('curve_id'):
        if g.sample_id.nunique() != 1 or g.replicate_id.nunique() != 1:
            raise ValueError(f'Curve {cid} mixes samples or replicates')
        if g.temperature_c.duplicated().any():
            raise ValueError(f'Curve {cid} repeats a temperature; do not merge capillaries into one curve')
    if d.groupby(['sample_id', 'replicate_id']).curve_id.nunique().max() > 1:
        raise ValueError('One curve per sample replicate; technical capillaries of one preparation are not separate replicates')
    samples = set(d.sample_id)
    for item in cfg['comparisons']:
        if {item['reference_sample'], item['test_sample']} - samples:
            raise ValueError(f"Comparison {item['id']} names an unknown sample")
    return d


# ---------- model ----------

def fractions(t_c, tms, dhs):
    tk = np.asarray(t_c, float)[:, None] + 273.15
    return expit((np.asarray(dhs)[None, :] / R_KJ) * (1 / tk - 1 / (np.asarray(tms)[None, :] + 273.15)) * -1)


def basis(t_c, tms, dhs, t0):
    f = fractions(t_c, tms, dhs); x = np.asarray(t_c, float) - t0
    return np.column_stack([np.ones_like(x), x, f, x * f[:, -1]])


def _unpack(z, k):
    tms = [z[0]]
    for g in z[1:k]:
        tms.append(tms[-1] + 10 ** g)
    return np.array(tms), 10 ** np.asarray(z[k:2 * k])


def _width_c(tm, dh):
    """10-90 % transition width in degrees C: 2 ln 9 R Tm^2 / dH."""
    return float(2 * np.log(9) * R_KJ * (tm + 273.15) ** 2 / dh)


def derivative_extrema(t, y, cfg, k):
    step = float(np.median(np.diff(t)))
    grid = np.arange(t[0], t[-1] + step / 2, step)
    yi = np.interp(grid, t, y)
    w = max(cfg['derivative']['polyorder'] + 2, int(round(cfg['derivative']['window_c'] / step)) | 1)
    if w >= len(grid):
        return {'status': 'window_too_wide', 'inflections_c': []}
    dy = savgol_filter(yi, w, cfg['derivative']['polyorder'], deriv=1, delta=step)
    sign = 1. if np.median(dy) >= 0 else -1.
    s = sign * dy
    peaks = [i for i in range(1, len(s) - 1) if s[i] >= s[i - 1] and s[i] > s[i + 1] and s[i] > .2 * s.max()]
    peaks = sorted(sorted(peaks, key=lambda i: -s[i])[:k])
    out = []
    for i in peaks:  # parabolic refinement of the extremum
        a, b, c = s[i - 1], s[i], s[i + 1]; den = a - 2 * b + c
        out.append(float(grid[i] + (.5 * (a - c) / den * step if den != 0 else 0.)))
    return {'status': 'computed', 'inflections_c': out, 'window_points': int(w), 'grid_step_c': step,
            'method': 'Savitzky-Golay first derivative on a uniform grid; local extrema above 20% of the largest, top k by height'}


def fit_curve(g, cfg, _profile=True):
    k = cfg['model']['transitions']; lv = cfg['uncertainty']['level']
    used = g[~g.exclude].sort_values('temperature_c')
    r = cfg['model']['fit_range_c']
    if r is not None:
        used = used[(used.temperature_c >= r[0]) & (used.temperature_c <= r[1])]
    t, y = used.temperature_c.to_numpy(float), used.signal.to_numpy(float)
    res = {'curve_id': str(g.curve_id.iloc[0]), 'sample_id': str(g.sample_id.iloc[0]), 'replicate_id': str(g.replicate_id.iloc[0]),
           'n_points': int(len(t)), 'status': 'failed', 'reportable': False, 'transitions': [], 'diagnostics': []}
    n_par = 2 * k + k + 3
    if len(t) - n_par < 10:
        res['diagnostics'].append('too_few_points_for_declared_transitions'); return res, [], []
    t0 = float(t.mean()); scale = max(float(np.ptp(y)), 1e-12)
    dhb = np.log10(cfg['model']['dh_bounds_kj_mol']); lo_t, hi_t = float(t.min()), float(t.max()); span = hi_t - lo_t
    lower = np.r_[lo_t, [np.log10(.5)] * (k - 1), [dhb[0]] * k]
    upper = np.r_[hi_t, [np.log10(span)] * (k - 1), [dhb[1]] * k]

    def evaluate(z):
        tms, dhs = _unpack(z, k); X = basis(t, tms, dhs, t0)
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        return X @ coef, coef, (X @ coef - y) / scale

    def as_start(tm_list, dh):
        tm_list = sorted(tm_list)
        gaps = [np.log10(max(.6, b_ - a_)) for a_, b_ in zip(tm_list, tm_list[1:])]
        return np.r_[tm_list[0], gaps, [np.log10(dh)] * k]
    candidates = []
    seeds = derivative_extrema(t, y, cfg, k).get('inflections_c', [])
    if seeds:  # derivative extrema seed the Tm values; missing transitions are spread over the remaining range
        seeds = sorted(seeds)
        while len(seeds) < k:
            seeds = sorted(seeds + [min(hi_t - 1., max(seeds) + (hi_t - max(seeds)) / 2)])
        candidates.append(seeds[:k])
    for frac in np.linspace(.15, .85, 8 if k == 1 else 6):
        first = lo_t + frac * span * (1. if k == 1 else .75)
        for spread in ((.5, 1.) if k > 1 else (1.,)):
            step = spread * (hi_t - first) / (k + .5)
            candidates.append([first + j * step for j in range(k)])
    starts = [as_start(c_, dh) for c_ in candidates for dh in (250., 600., 1200.)]
    best = None
    for s in starts:
        try:
            sol = least_squares(lambda z: evaluate(z)[2], np.clip(s, lower + 1e-6, upper - 1e-6), bounds=(lower, upper), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
        except ValueError:
            continue
        if np.isfinite(sol.fun).all() and (best is None or sol.cost < best.cost):
            best = sol
    if best is None:
        res['diagnostics'].append('optimizer_failed'); return res, [], []
    if k > 1:  # local restarts: shift each Tm of the best solution by +/-2 and +/-4 C
        tms0, dhs0 = _unpack(best.x, k)
        for i in range(k):
            for shift in (-4., -2., 2., 4.):
                moved = list(tms0); moved[i] += shift
                if np.all(np.diff(moved) > .6) and lo_t < moved[0] and moved[-1] < hi_t:
                    s0 = np.r_[moved[0], [np.log10(b_ - a_) for a_, b_ in zip(moved, moved[1:])], np.log10(dhs0)]
                    try:
                        sol = least_squares(lambda z: evaluate(z)[2], np.clip(s0, lower + 1e-6, upper - 1e-6), bounds=(lower, upper), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)
                    except ValueError:
                        continue
                    if np.isfinite(sol.fun).all() and sol.cost < best.cost:
                        best = sol
    pred, coef, resid = evaluate(best.x)
    tms, dhs = _unpack(best.x, k)
    sse = float(np.sum((pred - y) ** 2)); df = len(t) - n_par; rmse = float(np.sqrt(sse / df))
    rank = int(np.linalg.matrix_rank(best.jac))
    boundary = bool(np.any(np.minimum(best.x - lower, upper - best.x) < 1e-4))
    d = res['diagnostics']
    if boundary:
        d.append('numerical_boundary_hit')
    if rank < 2 * k:
        d.append('transition_parameters_not_locally_identifiable')
    rows = []
    for i, (tm, dh) in enumerate(zip(tms, dhs)):
        width = _width_c(tm, dh); amp = float(coef[2 + i])
        flags = []
        if not lo_t <= tm <= hi_t:
            flags.append('tm_outside_measured_range')
        if abs(amp) <= 3 * rmse:
            flags.append('amplitude_not_above_noise')
        if i == 0 and (t < tm - width).sum() < MIN_BASELINE_POINTS:
            flags.append('native_baseline_not_covered')
        if i == k - 1 and (t > tm + width).sum() < MIN_BASELINE_POINTS:
            flags.append('unfolded_baseline_not_covered')
        if i > 0 and tm - tms[i - 1] < (width + _width_c(tms[i - 1], dhs[i - 1])) / 2:
            flags.append('transitions_overlap')
            if 'transitions_overlap' not in rows[i - 1]['diagnostics']:
                rows[i - 1]['diagnostics'].append('transitions_overlap')
        rows.append({'transition': i + 1, 'tm_c': float(tm), 'apparent_dh_kj_mol': float(dh), 'width_10_90_c': width,
                     'amplitude': amp, 'ci_c': [None, None], 'ci_status': 'not_computed', 'diagnostics': flags})
    if not _profile:  # reduced-model fit for the structure test only
        res.update(sse=sse, residual_df=df, status='estimated')
        return res, [], []
    structure = None
    if k > 1:  # do the data support k transitions rather than k - 1? (3 parameters: Tm, dH, amplitude)
        reduced = deepcopy(cfg); reduced['model']['transitions'] = k - 1
        r_fit, _, _ = fit_curve(g, reduced, _profile=False)
        if r_fit.get('sse') is not None and sse > 0:
            f_stat = max(0., (r_fit['sse'] - sse) / 3) / (sse / df)
            structure = {'test': f'extra_sum_of_squares_F_{k}_vs_{k - 1}_transitions', 'f_statistic': float(f_stat), 'df': [3, df],
                         'p_value': float(f_dist.sf(f_stat, 3, df)), 'alpha': STRUCTURE_ALPHA}
            if structure['p_value'] >= STRUCTURE_ALPHA:
                d.append('declared_transitions_not_supported_by_data')
    # profile-F interval for each Tm (other nonlinear parameters re-optimized)
    threshold = sse * (1 + float(f_dist.ppf(lv, 1, df)) / df)
    for i in range(k):
        def profile(value, i=i):
            fixed = list(tms)

            @lru_cache(maxsize=256)
            def at(v):
                def fun(w):
                    tm_list = list(w[:k - 1]); dh_list = 10 ** np.asarray(w[k - 1:])
                    full = tm_list[:i] + [v] + tm_list[i:]
                    if np.any(np.diff(full) <= .05):
                        return np.full(len(t), 1e3)
                    X = basis(t, full, dh_list, t0); cf = np.linalg.lstsq(X, y, rcond=None)[0]
                    return (X @ cf - y) / scale
                others = [x for j, x in enumerate(fixed) if j != i]
                lb = np.r_[[lo_t] * (k - 1), [dhb[0]] * k]; ub = np.r_[[hi_t] * (k - 1), [dhb[1]] * k]
                # Several nuisance starts: a single start under-optimized some noisy profiles
                # (0.13.1 post-hoc supplement), which narrows the interval.
                starts = [np.r_[others, np.log10(dhs)]] + [np.r_[others, [np.log10(dh)] * k] for dh in (300., 900.)]
                starts += [np.r_[[x + shift for x in others], np.log10(dhs)] for shift in (-2., 2.)] if k > 1 else []
                best = np.inf
                for w0 in starts:
                    sol = least_squares(fun, np.clip(w0, lb + 1e-6, ub - 1e-6), bounds=(lb, ub), xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=3000)
                    best = min(best, float(np.sum((sol.fun * scale) ** 2)))
                return best
            return at(float(value))

        centre = float(tms[i]); ends = []
        for bound in (lo_t, hi_t):
            prev, root = centre, None
            for cur in np.linspace(centre, bound, 41)[1:]:
                if profile(cur) > threshold:
                    root = float(brentq(lambda v: profile(v) - threshold, min(prev, cur), max(prev, cur), xtol=1e-7)); break
                prev = cur
            ends.append(root)
        rows[i]['ci_c'] = ends
        rows[i]['ci_status'] = 'two_sided_profile_f' if None not in ends else 'open'
        if None in ends:
            rows[i]['diagnostics'].append('profile_interval_open')
    blocking = {'tm_outside_measured_range', 'amplitude_not_above_noise', 'native_baseline_not_covered',
                'unfolded_baseline_not_covered', 'transitions_overlap', 'profile_interval_open'}
    curve_ok = not boundary and rank >= 2 * k and 'declared_transitions_not_supported_by_data' not in d
    # A transition next to an unresolved one is not reportable either: a spurious low-amplitude
    # component beside a merged one passed every per-transition check (0.13.1 overlap stress).
    unresolved = {'tm_outside_measured_range', 'amplitude_not_above_noise', 'transitions_overlap'}
    for i, row in enumerate(rows):
        if any(unresolved & set(rows[j]['diagnostics']) for j in (i - 1, i + 1) if 0 <= j < k) and not unresolved & set(row['diagnostics']):
            row['diagnostics'].append('adjacent_transition_unresolved')
    blocking.add('adjacent_transition_unresolved')
    for row in rows:
        row['reportable'] = bool(curve_ok and not blocking & set(row['diagnostics']))
    res.update(status='estimated' if curve_ok else 'limited', reportable=bool(curve_ok and all(r['reportable'] for r in rows)),
               transitions=rows, structure_test=structure, sse=sse, residual_df=df, rmse=rmse, n_parameters=n_par, temperature_range_c=[lo_t, hi_t],
               centre_c=t0, linear_coefficients=dict(zip(['native_intercept', 'native_slope'] + [f'amplitude_{i+1}' for i in range(k)] + ['unfolded_slope_change'], map(float, coef))),
               derivative=derivative_extrema(t, y, cfg, k))
    grid_t = np.linspace(lo_t, hi_t, 300); X = basis(grid_t, tms, dhs, t0)
    grid = [{'curve_id': res['curve_id'], 'temperature_c': float(a), 'predicted': float(b)} for a, b in zip(grid_t, X @ coef)]
    obs = [{'curve_id': res['curve_id'], 'temperature_c': float(a), 'signal': float(b), 'predicted': float(p), 'residual': float(b - p)}
           for a, b, p in zip(t, y, pred)]
    return res, grid, obs


def compute(d, cfg):
    fits, grids, obs = [], [], []
    for _, g in d.groupby('curve_id', sort=False):
        f, gr, ob = fit_curve(g, cfg); fits.append(f); grids += gr; obs += ob
    k, lv = cfg['model']['transitions'], cfg['uncertainty']['level']
    independent = cfg['replicates']['independent_unit'] == 'replicate_id'
    summaries, comps, failing = [], [], []
    for f in fits:
        if not f['reportable']:
            failing.append({'curve_id': f['curve_id'], 'diagnostics': f['diagnostics'] + [x for r in f['transitions'] for x in r['diagnostics']]})
    for sample in dict.fromkeys(f['sample_id'] for f in fits):
        fs = [f for f in fits if f['sample_id'] == sample]
        for i in range(k):
            row = {'sample_id': sample, 'transition': i + 1, 'n_replicates': len(fs), 'mean_tm_c': None, 'sd_c': None, 'ci_c': [None, None]}
            if not independent:
                row['status'] = 'independent_replicates_not_declared'
            elif any(not f['transitions'] or not f['transitions'][i]['reportable'] for f in fs):
                row['status'] = 'withheld_incomplete_estimates'
            else:
                v = np.array([f['transitions'][i]['tm_c'] for f in fs]); row['mean_tm_c'] = float(v.mean())
                if len(v) > 1:
                    sd = float(v.std(ddof=1)); half = float(t_dist.ppf((1 + lv) / 2, len(v) - 1) * sd / np.sqrt(len(v)))
                    row.update(sd_c=sd, ci_c=[row['mean_tm_c'] - half, row['mean_tm_c'] + half], status='summarized_t')
                else:
                    row['status'] = 'one_replicate_no_interval'
            summaries.append(row)
    pvals = []
    for item in cfg['comparisons']:
        for i in range(k):
            ref = [s for s in summaries if s['sample_id'] == item['reference_sample'] and s['transition'] == i + 1][0]
            tst = [s for s in summaries if s['sample_id'] == item['test_sample'] and s['transition'] == i + 1][0]
            row = {'id': item['id'], 'transition': i + 1, 'reference_sample': item['reference_sample'], 'test_sample': item['test_sample'],
                   'delta_tm_c': None, 'ci_c': [None, None], 'p_value': None, 'holm_p': None}
            if ref['status'] != 'summarized_t' or tst['status'] != 'summarized_t':
                row['status'] = 'withheld_summary_not_available'
            else:
                a = np.array([f['transitions'][i]['tm_c'] for f in fits if f['sample_id'] == item['reference_sample']])
                b = np.array([f['transitions'][i]['tm_c'] for f in fits if f['sample_id'] == item['test_sample']])
                va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b); se = np.sqrt(va + vb)
                dfw = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1)) if va + vb > 0 else np.inf
                delta = float(b.mean() - a.mean()); q = float(t_dist.ppf((1 + lv) / 2, dfw))
                p = float(2 * t_dist.sf(abs(delta / se), dfw)) if se > 0 else 0.
                row.update(delta_tm_c=delta, ci_c=[delta - q * se, delta + q * se], welch_df=float(dfw), p_value=p, status='estimated')
                pvals.append((len(comps), p))
            comps.append(row)
    if pvals:  # Holm across every comparison x transition in the declared family
        order = sorted(pvals, key=lambda z: z[1]); m = len(order); running = 0.
        for rank_, (idx, p) in enumerate(order):
            running = max(running, min(1., (m - rank_) * p)); comps[idx]['holm_p'] = running
    must = ['Tm is apparent: it depends on scan rate, buffer, concentration and irreversible aggregation; compare only curves measured under the same declared conditions.',
            'Each transition is assigned by order of Tm; mapping a transition to a domain (CH2, Fab, CH3) needs separate evidence.']
    if cfg['assay']['reversibility'] == 'irreversible_or_unknown':
        must.append('Apparent dH is a curve-shape parameter, not a thermodynamic enthalpy, because reversibility was not demonstrated.')
    if not independent:
        must.append('Replicates were not declared independent; no between-replicate summary or comparison was made.')
    if any('transitions_overlap' in r['diagnostics'] for f in fits for r in f['transitions']):
        must.append('Some transitions overlap; their individual Tm values depend on the model and are withheld.')
    return {'analysis_type': 'thermal_unfolding', 'primary': {'fits': fits, 'summaries': summaries, 'comparisons': comps,
            'assay_context': {**cfg['assay'], 'transitions': k, 'transitions_source': cfg['model']['transitions_source'], 'fit_range_c': cfg['model']['fit_range_c']}},
            'curve_grid': grids, 'observations': obs, 'must_mention': must, 'failing_items': failing,
            'limitations': ['Two-state van t Hoff transitions with linear native baseline and one unfolded-slope change; kinetic (irreversible) models are not fitted.',
                            'Derivative inflections are descriptive and depend on the smoothing window; instrument software may differ.',
                            'Comparisons use Welch intervals on independent replicate Tm values with Holm-adjusted p values across the declared family.']}


def render(run, result, style):
    """Per-sample melting curves, fits and residuals from saved values only."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); names = []
    obs, grid, fits = pd.DataFrame(result['observations']), pd.DataFrame(result['curve_grid']), result['primary']['fits']
    for n, sample in enumerate(dict.fromkeys(f['sample_id'] for f in fits)):
        fs = [f for f in fits if f['sample_id'] == sample]
        with plt.rc_context(pstyle.rc(token, *font_setup())):
            fig, (ax, res) = plt.subplots(2, 1, figsize=(6.4, 4.8), sharex=True, height_ratios=[2.4, 1], layout='constrained')
            for j, f in enumerate(fs):
                o, g = obs[obs.curve_id == f['curve_id']], grid[grid.curve_id == f['curve_id']]
                colour = pstyle.color(token, j)
                ax.plot(o.temperature_c, o.signal, '.', ms=2.5, color=colour, alpha=.6)
                ax.plot(g.temperature_c, g.predicted, '-', lw=1.3, color=colour, label=f['replicate_id'])
                res.plot(o.temperature_c, o.residual, '.', ms=2, color=colour, alpha=.7)
                for t_ in f['transitions']:
                    ax.axvline(t_['tm_c'], color=colour, lw=.7, ls=(0, (3, 2)) if t_['reportable'] else (0, (1, 2)))
            res.axhline(0, color=pstyle.MUTED, lw=.8)
            ax.set(ylabel='Signal', title=f'{sample}: apparent Tm (dashed; dotted = withheld)'); res.set(xlabel='Temperature (°C)', ylabel='Residual')
            ax.legend(frameon=False, fontsize=7)
            name = f'melting_{n+1}'; pstyle.save(fig, out/name, token, 300); names.append(name)
    return names
