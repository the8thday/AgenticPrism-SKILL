"""Four-parameter logistic fits, replicate summaries and curve comparisons.

Only log10 C50 and log10 |Hill| are optimized numerically. Plateaus are solved
for each proposal: bounded linear least squares (Top - Bottom >= 0), refined to
the Prism-style relative objective sum(((Y - Ycurve) / Ycurve)^2) when relative
weighting is selected. Top is always the higher plateau. Fixed plateaus remove
their parameter from the fit and from the residual degrees of freedom.
"""
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq, least_squares, lsq_linear, minimize_scalar
from scipy.special import expit
from scipy.stats import f as f_dist, spearmanr, t as t_dist
from .dose_schema import DOSE_UNITS

PENALTY = 1e3  # relative residual assigned where the curve is not positive


def response(concentration, bottom, top, log50, hill, direction):
    x = np.asarray(concentration, dtype=float)
    if np.any(x < 0):
        raise ValueError("Negative dose")
    signed = 1. if direction == "increasing" else -1.
    z = np.empty_like(x)
    positive = x > 0
    z[positive] = expit(np.log(10.) * signed * hill * (np.log10(x[positive]) - log50))
    z[~positive] = 0. if signed > 0 else 1.
    return bottom + (top - bottom) * z


def free_plateaus(cfg):
    return (cfg["fit"]["fixed_bottom"] is None) + (cfg["fit"]["fixed_top"] is None)


def _weighted(prediction, y, cfg):
    if cfg["fit"]["weighting"] == "relative":
        safe = np.where(prediction > 0, prediction, 1.)
        return np.where(prediction > 0, (prediction - y) / safe, PENALTY)
    return prediction - y


def _plateaus(shapes, ys, cfg):
    """Best (Bottom, Top) shared by the given curves for fixed curve shapes."""
    fixed_bottom, fixed_top = cfg["fit"]["fixed_bottom"], cfg["fit"]["fixed_top"]
    if fixed_bottom is not None and fixed_top is not None:
        return float(fixed_bottom), float(fixed_top)
    s, y = np.concatenate(shapes), np.concatenate(ys)
    # Free variables are (Bottom, span), span, or span, with span = Top - Bottom >= 0.
    if fixed_bottom is None and fixed_top is None:
        design, offset, lower = np.column_stack((np.ones_like(s), s)), 0., [-np.inf, 0.]
    elif fixed_bottom is not None:
        design, offset, lower = s[:, None], float(fixed_bottom), [0.]
    else:
        design, offset, lower = -(1. - s)[:, None], float(fixed_top), [0.]
    upper = [np.inf] * len(lower)
    relative = cfg["fit"]["weighting"] == "relative"
    w = 1. / np.abs(y) if relative else np.ones_like(y)
    v = lsq_linear(design * w[:, None], (y - offset) * w, bounds=(lower, upper), method="bvls").x
    if relative:
        v = least_squares(lambda vv: _weighted(offset + design @ vv, y, cfg), v, bounds=(lower, upper),
                          xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200).x
    if fixed_bottom is None and fixed_top is None:
        return float(v[0]), float(v[0] + v[1])
    if fixed_bottom is not None:
        return float(fixed_bottom), float(fixed_bottom + v[0])
    return float(fixed_top - v[0]), float(fixed_top)


def _evaluate(curves, log50s, log_hills, groups, cfg):
    """Predictions, per-curve plateaus and weighted residuals; groups share plateaus."""
    direction = cfg["assay"]["direction"]
    shapes = [response(x, 0., 1., l, 10. ** h, direction) for (x, _), l, h in zip(curves, log50s, log_hills)]
    plateaus = [None] * len(curves)
    for group in groups:
        pair = _plateaus([shapes[i] for i in group], [curves[i][1] for i in group], cfg)
        for i in group:
            plateaus[i] = pair
    predictions = [b + (t - b) * s for (b, t), s in zip(plateaus, shapes)]
    residual = np.concatenate([_weighted(p, y, cfg) for p, (_, y) in zip(predictions, curves)])
    return predictions, plateaus, residual


def _bounds(cfg, positives):
    lo = float(min(np.log10(p.min()) for p in positives))
    hi = float(max(np.log10(p.max()) for p in positives))
    mid = cfg["fit"]["log50_bounds"] or [lo - 3., hi + 3.]
    return mid, list(np.log10(cfg["fit"]["hill_bounds"])), lo, hi


def _best(fun, starts, lower, upper, max_nfev):
    best = None
    for start in starts:
        try:
            s = least_squares(fun, np.clip(start, lower + 1e-8, upper - 1e-8), bounds=(lower, upper),
                              max_nfev=max_nfev, ftol=1e-11, xtol=1e-11, gtol=1e-11)
        except ValueError:
            continue
        if s.success and np.isfinite(s.fun).all() and (best is None or s.fun @ s.fun < best.fun @ best.fun):
            best = s
    return best


def _noise_floor(y, cfg):
    # Objective below this is treated as noiseless: response units squared, or relative residuals.
    return 1e-22 * (max(1., float(y @ y)) if cfg["fit"]["weighting"] == "unweighted" else len(y))


def _at_bound(z, lower, upper):
    return bool(np.any(np.minimum(np.asarray(z) - lower, upper - np.asarray(z)) < 1e-4))


def fit_dose_curve(group, cfg, *, shape_diagnostics=True):
    curve_id = str(group.curve_id.iloc[0])
    used = group.loc[~group.exclude]
    f = cfg["fit"]
    n_parameters = 2 + free_plateaus(cfg)
    result = {"curve_id": curve_id, "sample_id": str(group.sample_id.iloc[0]),
              "experiment_id": str(group.experiment_id.iloc[0]), "endpoint": cfg["assay"]["endpoint"],
              "direction": cfg["assay"]["direction"], "status": "failed", "reportable": False,
              "canonical_unit": str(group.canonical_unit.iloc[0]), "input_unit": str(group.concentration_unit.iloc[0]),
              "half_response_canonical": None, "half_response_input_unit": None,
              "log10_half_response_canonical": None, "bottom": None, "top": None,
              "hill_slope_signed": None, "n_input": len(group), "n_fit": len(used),
              "n_distinct_positive_doses": int(used.loc[used.concentration_canonical > 0, "concentration_canonical"].nunique()),
              "weighting": f["weighting"], "fixed_bottom": f["fixed_bottom"], "fixed_top": f["fixed_top"],
              "n_parameters": n_parameters, "objective_sse": None, "residual_df": None, "rmse": None,
              "range_status": "unknown", "ci_method": cfg["uncertainty"]["method"],
              "ci_level": cfg["uncertainty"]["level"], "ci_status": "not_computed",
              "ci_canonical": [None, None], "diagnostics": []}

    def fail(reason):
        result["diagnostics"].append(reason)
        return result, [], []

    x = used.concentration_canonical.to_numpy(float)
    y = used.response.to_numpy(float)
    positive = x[x > 0]
    if len(x) < 10 or len(np.unique(positive)) < 6:
        return fail("insufficient_distinct_doses_or_residual_df")
    y_span = float(np.ptp(y))
    y_scale = max(float(np.max(np.abs(y))), 1e-12)
    if y_span <= y_scale * 1e-9:
        return fail("flat_response")
    if f["weighting"] == "relative" and np.any(y <= 0):
        return fail("nonpositive_response_for_relative_weighting")
    direction = cfg["assay"]["direction"]
    low_dose, high_dose = float(positive.min()), float(positive.max())
    bounds, (hill_min, hill_max), lo, hi = _bounds(cfg, [positive])
    lower = np.array([bounds[0], hill_min], dtype=float)
    upper = np.array([bounds[1], hill_max], dtype=float)
    if not lower[0] < upper[0]:
        return fail("invalid_midpoint_bounds")
    # Unweighted residuals are divided by the response scale for unit-free optimizer tolerances.
    scale = y_scale if f["weighting"] == "unweighted" else 1.
    curves = [(x, y)]

    def evaluate(log50, log_hill):
        predictions, plateaus, residual = _evaluate(curves, [log50], [log_hill], [[0]], cfg)
        return predictions[0], plateaus[0], residual

    m = int(np.ceil(np.sqrt(f["multistart"])))
    starts = [(a, b) for a in np.linspace(lo, hi, m)
              for b in np.linspace(max(hill_min, np.log10(.5)), min(hill_max, np.log10(2.)), m)][:f["multistart"]]
    best = _best(lambda z: evaluate(*z)[2] / scale, starts, lower, upper, f["max_nfev"])
    if best is None:
        return fail("optimizer_failed")
    prediction, (bottom, top), weighted = evaluate(*best.x)
    sse = float(weighted @ weighted)
    df = len(x) - n_parameters
    rmse = float(np.sqrt(np.sum((prediction - y) ** 2) / df))
    midpoint = float(10. ** best.x[0])
    amp = top - bottom
    rank = int(np.linalg.matrix_rank(best.jac))
    boundary = bool(_at_bound(best.x, lower, upper) or amp <= y_scale * 1e-8)
    result.update(half_response_canonical=midpoint,
                  half_response_input_unit=midpoint / DOSE_UNITS[str(group.concentration_unit.iloc[0])],
                  log10_half_response_canonical=float(best.x[0]), bottom=bottom, top=top,
                  hill_slope_signed=float((1 if direction == "increasing" else -1) * 10. ** best.x[1]),
                  objective_sse=sse, residual_df=df, rmse=rmse, rate_jacobian_rank=rank,
                  numerical_boundary_hit=boundary,
                  range_status="within_range" if low_dose <= midpoint <= high_dose else "below_range" if midpoint < low_dose else "above_range",
                  lowest_tested_positive_canonical=low_dose, highest_tested_canonical=high_dose)
    if f["weighting"] == "relative":
        result["weighted_rmse"] = float(np.sqrt(sse / df))
    if boundary:
        result["diagnostics"].append("numerical_boundary_hit")
    if rank < 2:
        result["diagnostics"].append("midpoint_or_hill_not_locally_identifiable")
    if result["range_status"] != "within_range":
        result["diagnostics"].append("half_response_outside_tested_range")
    if amp <= 3 * rmse:
        result["diagnostics"].append("low_response_range_relative_to_noise")
    shape_edges = response(np.array([low_dose, high_dose]), 0., 1., best.x[0], 10. ** best.x[1], direction)
    if min(shape_edges) > .1 or max(shape_edges) < .9:
        result["diagnostics"].append("plateau_not_well_covered")
    intermediate = sum(.1 < response(np.array([dose]), 0., 1., best.x[0], 10. ** best.x[1], direction)[0] < .9
                       for dose in np.unique(positive))
    if intermediate < 3:
        result["diagnostics"].append("few_transition_doses_for_variable_slope")
    if (x == 0).sum() > len(x) / 4:
        result["diagnostics"].append("zero_controls_dominate_unweighted_fit")
    rho, p_value = spearmanr(np.log10(positive), y[x > 0])
    opposite = bool(np.isfinite(rho) and p_value < .05 and rho * (1 if direction == "increasing" else -1) < 0)
    result["observed_trend_spearman"] = float(rho) if np.isfinite(rho) else None
    if opposite:
        result["diagnostics"].append("observed_direction_opposite_to_declared")

    # Only the default 4PL workflow changes. Internal 4PL reference fits used by
    # other declared models retain their original scientific artifacts.
    shape_passed = True
    if shape_diagnostics and cfg.get('model', 'relative_four_parameter_logistic') == 'relative_four_parameter_logistic':
        from .dose_diagnostics import shape_checks
        checks = shape_checks(x, y, prediction, n_parameters, f['weighting'])
        result['shape_diagnostics'] = checks
        result['diagnostics'].extend(checks['diagnostics'])
        shape_passed = checks['passed']

    if cfg["uncertainty"]["method"] == "profile_f":
        if sse <= _noise_floor(y, cfg):
            result["ci_status"] = "noise_scale_not_estimable"
        else:
            threshold = sse * (1. + float(f_dist.ppf(cfg["uncertainty"]["level"], 1, df)) / df)

            @lru_cache(maxsize=512)
            def profile(log_midpoint):
                def objective(log_hill):
                    r = evaluate(log_midpoint, log_hill)[2]
                    return float(r @ r)
                opt = minimize_scalar(objective, bounds=(hill_min, hill_max), method="bounded",
                                      options={"xatol": 1e-9})
                candidates = [(objective(hill_min), hill_min), (objective(hill_max), hill_max),
                              (float(opt.fun), float(opt.x))]
                return min(candidates, key=lambda z: z[0])

            def endpoint(bound):
                grid = np.linspace(float(best.x[0]), bound, 45)
                previous = float(best.x[0])
                for current in grid[1:]:
                    if profile(float(current))[0] > threshold:
                        root = brentq(lambda v: profile(float(v))[0] - threshold, min(previous, current),
                                      max(previous, current), xtol=1e-9)
                        return float(root), profile(float(root))[1]
                    previous = float(current)
                return None, None

            lower_root, lower_hill = endpoint(bounds[0])
            upper_root, upper_hill = endpoint(bounds[1])
            if lower_root is not None and upper_root is not None:
                result["ci_canonical"] = [float(10. ** lower_root), float(10. ** upper_root)]
                result["ci_status"] = "two_sided_profile_f"
                if any(min(abs(h - hill_min), abs(h - hill_max)) < 1e-4 for h in (lower_hill, upper_hill)):
                    result["diagnostics"].append("profile_hill_nuisance_at_boundary")
            else:
                result["ci_status"] = "both_open" if lower_root is None and upper_root is None else "lower_open" if lower_root is None else "upper_open"
                result["diagnostics"].append("profile_interval_open")
    else:
        result["ci_status"] = "not_requested"

    result["reportable"] = (not boundary and rank == 2 and result["range_status"] == "within_range"
                            and amp > 3 * rmse and not opposite and shape_passed
                            and result["ci_status"] in ("two_sided_profile_f", "not_requested")
                            and "profile_hill_nuisance_at_boundary" not in result["diagnostics"])
    result["status"] = "limited" if not result["reportable"] else "estimated_with_diagnostics" if result["diagnostics"] else "estimated"
    result["diagnostics"] = list(dict.fromkeys(result["diagnostics"]))
    dose_grid = np.geomspace(low_dose, high_dose, 200)
    grid = [{"curve_id": curve_id, "concentration_canonical": float(dose), "predicted_response": float(pred)}
            for dose, pred in zip(dose_grid, response(dose_grid, bottom, top, best.x[0], 10. ** best.x[1], direction))]
    all_dose = group.concentration_canonical.to_numpy(float)
    all_pred = response(all_dose, bottom, top, best.x[0], 10. ** best.x[1], direction)
    observations = [{"observation_id": str(row.observation_id), "curve_id": curve_id,
                     "concentration_canonical": float(row.concentration_canonical),
                     "observed_response": float(row.response), "predicted_response": float(pred),
                     "residual": float(row.response - pred), "used_in_fit": not bool(row.exclude)}
                    for (_, row), pred in zip(group.iterrows(), all_pred)]
    return result, grid, observations


def summarize_dose(fits, cfg):
    """Per sample: technical curves averaged on log10 C50 within an experiment,
    then experiments weighted equally. Any failed or limited curve withholds the
    whole sample rather than silently selecting successes."""
    rows = []
    rep, level = cfg["replicates"], cfg["uncertainty"]["level"]
    for sample in dict.fromkeys(f["sample_id"] for f in fits):
        fs = [f for f in fits if f["sample_id"] == sample]
        units = {f["canonical_unit"] for f in fs}
        input_units = {f["input_unit"] for f in fs}
        row = {"sample_id": sample, "endpoint": cfg["assay"]["endpoint"], "n_curves": len(fs),
               "n_failed_or_limited": sum(not f["reportable"] for f in fs),
               "n_experiments": len({f["experiment_id"] for f in fs}),
               "canonical_unit": units.pop() if len(units) == 1 else None,
               "input_unit": input_units.pop() if len(input_units) == 1 else None,
               "geometric_mean_canonical": None, "ci_low_canonical": None, "ci_high_canonical": None,
               "geometric_mean_input_unit": None, "ci_low_input_unit": None, "ci_high_input_unit": None,
               "log10_sd_between_experiments": None}
        if rep["independent_unit"] == "none":
            row["status"] = "independent_units_unconfirmed"
        elif not rep["conditions_comparable"]:
            row["status"] = "comparability_unconfirmed"
        elif row["canonical_unit"] is None:
            row["status"] = "mixed_dose_units_not_summarized"
        elif row["n_failed_or_limited"]:
            row["status"] = "withheld_incomplete_estimates"
        else:
            logs = [np.mean([f["log10_half_response_canonical"] for f in fs if f["experiment_id"] == e])
                    for e in dict.fromkeys(f["experiment_id"] for f in fs)]
            mean = float(np.mean(logs))
            row["geometric_mean_canonical"] = float(10. ** mean)
            row["status"] = "one_experiment_no_between_experiment_ci"
            if len(logs) > 1:
                sd = float(np.std(logs, ddof=1))
                half = float(t_dist.ppf((1 + level) / 2, len(logs) - 1) * sd / np.sqrt(len(logs)))
                row.update(ci_low_canonical=float(10. ** (mean - half)), ci_high_canonical=float(10. ** (mean + half)),
                           log10_sd_between_experiments=sd, status="summarized_log_t")
            if row["input_unit"] is not None:
                factor = DOSE_UNITS[row["input_unit"]]
                for key in ("geometric_mean", "ci_low", "ci_high"):
                    value = row[key + "_canonical"]
                    row[key + "_input_unit"] = None if value is None else value / factor
        rows.append(row)
    return rows


def _f_test(sse_null, df_null, sse_alt, df_alt):
    numerator_df = df_null - df_alt
    statistic = max(sse_null - sse_alt, 0.) / numerator_df / (sse_alt / df_alt)
    return {"F": float(statistic), "df_numerator": int(numerator_df), "df_denominator": int(df_alt),
            "p_value": float(f_dist.sf(statistic, numerator_df, df_alt)),
            "sse_null": float(sse_null), "sse_alternative": float(sse_alt)}


def compare_curves(item, groups, fits, cfg):
    """Relative potency from a parallel 4PL (shared plateaus and Hill, separate C50)
    with a profile-F interval on log10 RP, a parallelism F test against separate
    fits, and a shared-C50 extra-sum-of-squares F test. RP = C50_ref / C50_test."""
    ref_id, test_id = item["reference_curve"], item["test_curve"]
    rf, tf = fits[ref_id], fits[test_id]
    out = {"comparison_id": item["id"], "reference_curve": ref_id, "test_curve": test_id,
           "reference_sample": rf["sample_id"], "test_sample": tf["sample_id"],
           "reference_experiment": rf["experiment_id"], "test_experiment": tf["experiment_id"],
           "endpoint": cfg["assay"]["endpoint"], "canonical_unit": rf["canonical_unit"],
           "status": "failed", "reportable": False, "relative_potency": None, "log10_relative_potency": None,
           "rp_ci": [None, None], "ci_status": "not_computed", "ci_level": cfg["uncertainty"]["level"],
           "individual_c50_ratio": None, "parallel_model": None, "parallelism_f_test": None,
           "shared_c50_f_test": None, "parallelism_method": cfg["comparison_settings"]["parallelism_method"],
           "parallelism_equivalence": None, "rp_acceptance": None, "diagnostics": []}

    def fail(reason):
        out["diagnostics"].append(reason)
        return out, []

    if rf["half_response_canonical"] is None or tf["half_response_canonical"] is None:
        return fail("individual_fit_failed")
    out["individual_c50_ratio"] = rf["half_response_canonical"] / tf["half_response_canonical"]
    frames = [groups[ref_id].loc[~groups[ref_id].exclude], groups[test_id].loc[~groups[test_id].exclude]]
    curves = [(g.concentration_canonical.to_numpy(float), g.response.to_numpy(float)) for g in frames]
    positives = [x[x > 0] for x, _ in curves]
    bounds, hill_bounds, _, _ = _bounds(cfg, positives)
    lower_single, upper_single = bounds[0], bounds[1]
    y_all = np.concatenate([y for _, y in curves])
    scale = max(float(np.max(np.abs(y_all))), 1e-12) if cfg["fit"]["weighting"] == "unweighted" else 1.
    n = len(y_all)
    k = 2 + free_plateaus(cfg)
    df_separate = n - 2 * k
    sse_separate = rf["objective_sse"] + tf["objective_sse"]
    lr, lt = rf["log10_half_response_canonical"], tf["log10_half_response_canonical"]
    hr, ht = np.log10(abs(rf["hill_slope_signed"])), np.log10(abs(tf["hill_slope_signed"]))
    max_nfev = cfg["fit"]["max_nfev"]

    # Shared C50 (plateaus and Hill separate): Prism's "does logEC50 differ" test.
    shared = _best(lambda z: _evaluate(curves, [z[0], z[0]], [z[1], z[2]], [[0], [1]], cfg)[2] / scale,
                   [(l, hr, ht) for l in (lr, lt, (lr + lt) / 2)],
                   np.array([lower_single, hill_bounds[0], hill_bounds[0]]),
                   np.array([upper_single, hill_bounds[1], hill_bounds[1]]), max_nfev)
    if shared is not None:
        r = _evaluate(curves, [shared.x[0]] * 2, list(shared.x[1:]), [[0], [1]], cfg)[2]
        out["shared_c50_f_test"] = _f_test(float(r @ r), df_separate + 1, sse_separate, df_separate)

    # Parallel model: shared plateaus and Hill, separate C50.
    lower = np.array([lower_single, lower_single, hill_bounds[0]])
    upper = np.array([upper_single, upper_single, hill_bounds[1]])

    def parallel(z):
        return _evaluate(curves, [z[0], z[1]], [z[2], z[2]], [[0, 1]], cfg)

    best = _best(lambda z: parallel(z)[2] / scale, [(lr, lt, h) for h in (hr, ht, (hr + ht) / 2)],
                 lower, upper, max_nfev)
    if best is None:
        return fail("parallel_optimizer_failed")
    predictions, plateaus, residual = parallel(best.x)
    sse_parallel = float(residual @ residual)
    k_parallel = 3 + free_plateaus(cfg)
    df_parallel = n - k_parallel
    boundary = _at_bound(best.x, lower, upper) or plateaus[0][1] - plateaus[0][0] <= 1e-8 * max(np.abs(y_all))
    rank = int(np.linalg.matrix_rank(best.jac))
    delta = float(best.x[0] - best.x[1])
    sign = 1 if cfg["assay"]["direction"] == "increasing" else -1
    out.update(relative_potency=float(10. ** delta), log10_relative_potency=delta, status="limited",
               parallel_model={"bottom": plateaus[0][0], "top": plateaus[0][1],
                               "hill_slope_signed": float(sign * 10. ** best.x[2]),
                               "reference_c50_canonical": float(10. ** best.x[0]),
                               "test_c50_canonical": float(10. ** best.x[1]),
                               "objective_sse": sse_parallel, "residual_df": df_parallel,
                               "n_parameters": k_parallel, "numerical_boundary_hit": bool(boundary),
                               "jacobian_rank": rank})
    parallelism = _f_test(sse_parallel, df_parallel, sse_separate, df_separate)
    out["parallelism_f_test"] = parallelism
    if boundary:
        out["diagnostics"].append("numerical_boundary_hit")
    if rank < 3:
        out["diagnostics"].append("parallel_model_not_locally_identifiable")
    if out["parallelism_method"] == "equivalence":
        from .equivalence import parameter_equivalence
        out["parallelism_equivalence"] = parameter_equivalence(rf, tf, curves, cfg)
        if not out["parallelism_equivalence"]["equivalent"]:
            out["diagnostics"].append("parallelism_equivalence_not_demonstrated")
    elif parallelism["p_value"] < cfg["comparison_settings"]["parallelism_alpha"]:
        out["diagnostics"].append("nonparallel_by_f_test")
    if not (rf["reportable"] and tf["reportable"]):
        out["diagnostics"].append("individual_fit_not_reportable")
    if rf["experiment_id"] != tf["experiment_id"]:
        out["diagnostics"].append("curves_from_different_experiments")

    if cfg["uncertainty"]["method"] == "profile_f" and sse_parallel > _noise_floor(y_all, cfg):
        threshold = sse_parallel * (1. + float(f_dist.ppf(cfg["uncertainty"]["level"], 1, df_parallel)) / df_parallel)
        state = {"z": np.array([best.x[0], best.x[2]])}

        def profile(d):
            # log10 C50_test = log10 C50_ref - d; both midpoints stay inside the search bounds.
            lo_ref, hi_ref = max(lower_single, lower_single + d), min(upper_single, upper_single + d)
            if not lo_ref < hi_ref:
                raise ArithmeticError("relative potency outside midpoint search domain")
            lb, ub = np.array([lo_ref, hill_bounds[0]]), np.array([hi_ref, hill_bounds[1]])
            fun = lambda z: parallel([z[0], z[0] - d, z[1]])[2] / scale
            s = _best(fun, [state["z"], np.array([best.x[0], best.x[2]])], lb, ub, max_nfev)
            if s is None:
                raise ArithmeticError("profile nuisance optimization failed")
            state["z"] = s.x
            return float(s.fun @ s.fun) * scale ** 2 - threshold

        def endpoint(direction):
            span = (upper_single - lower_single) / 2
            previous = delta
            state["z"] = np.array([best.x[0], best.x[2]])
            for current in np.linspace(delta, delta + direction * span, 41)[1:]:
                try:
                    value = profile(current)
                except ArithmeticError:
                    return None
                if value > 0:
                    return float(brentq(profile, min(previous, current), max(previous, current), xtol=1e-9))
                previous = current
            return None

        low, high = endpoint(-1), endpoint(1)
        if low is not None and high is not None:
            out["rp_ci"] = [float(10. ** low), float(10. ** high)]
            out["ci_status"] = "two_sided_profile_f"
        else:
            out["ci_status"] = "both_open" if low is None and high is None else "lower_open" if low is None else "upper_open"
            out["diagnostics"].append("profile_interval_open")
    elif cfg["uncertainty"]["method"] == "profile_f":
        out["ci_status"] = "noise_scale_not_estimable"
    else:
        out["ci_status"] = "not_requested"
    blocking = {"numerical_boundary_hit", "parallel_model_not_locally_identifiable", "nonparallel_by_f_test",
                "parallelism_equivalence_not_demonstrated", "individual_fit_not_reportable", "profile_interval_open"}
    out["reportable"] = not blocking & set(out["diagnostics"]) and out["ci_status"] in ("two_sided_profile_f", "not_requested")
    out["status"] = "estimated_with_diagnostics" if out["reportable"] and out["diagnostics"] else "estimated" if out["reportable"] else "limited"
    limits = cfg["comparison_settings"]["rp_acceptance_limits"]
    if limits is not None:
        from .equivalence import rp_acceptance
        out["rp_acceptance"] = rp_acceptance(out["rp_ci"] if out["reportable"] else [None, None], out["ci_level"],
                                             limits, cfg["comparison_settings"]["rp_acceptance_rationale"])
    grid = []
    b, t = plateaus[0]
    for (x, _), curve_id, log50 in zip(curves, (ref_id, test_id), best.x[:2]):
        dose = np.geomspace(x[x > 0].min(), x.max(), 200)
        grid += [{"comparison_id": item["id"], "curve_id": curve_id, "concentration_canonical": float(c),
                  "predicted_response": float(p)}
                 for c, p in zip(dose, response(dose, b, t, log50, 10. ** best.x[2], cfg["assay"]["direction"]))]
    return out, grid


def summarize_potency(comparisons, cfg):
    """Geometric-mean RP per reference/test sample pair across independent
    experiments. Only same-experiment comparisons count; any non-reportable one
    withholds the pair, and one comparison per experiment is required."""
    rows = []
    rep, level = cfg["replicates"], cfg["uncertainty"]["level"]
    pairs = dict.fromkeys((c["reference_sample"], c["test_sample"]) for c in comparisons)
    for reference, test in pairs:
        cs = [c for c in comparisons if (c["reference_sample"], c["test_sample"]) == (reference, test)]
        experiments = [c["reference_experiment"] for c in cs]
        row = {"reference_sample": reference, "test_sample": test, "n_comparisons": len(cs),
               "n_experiments": len(set(experiments)), "n_not_reportable": sum(not c["reportable"] for c in cs),
               "geometric_mean_rp": None, "ci_low": None, "ci_high": None, "log10_sd_between_experiments": None}
        if rep["independent_unit"] == "none":
            row["status"] = "independent_units_unconfirmed"
        elif not rep["conditions_comparable"]:
            row["status"] = "comparability_unconfirmed"
        elif any(c["reference_experiment"] != c["test_experiment"] for c in cs):
            row["status"] = "cross_experiment_comparison_not_summarized"
        elif len(set(experiments)) != len(experiments):
            row["status"] = "multiple_comparisons_per_experiment_not_summarized"
        elif row["n_not_reportable"]:
            row["status"] = "withheld_incomplete_estimates"
        else:
            logs = [c["log10_relative_potency"] for c in cs]
            mean = float(np.mean(logs))
            row["geometric_mean_rp"] = float(10. ** mean)
            row["status"] = "one_experiment_no_between_experiment_ci"
            if len(logs) > 1:
                sd = float(np.std(logs, ddof=1))
                half = float(t_dist.ppf((1 + level) / 2, len(logs) - 1) * sd / np.sqrt(len(logs)))
                row.update(ci_low=float(10. ** (mean - half)), ci_high=float(10. ** (mean + half)),
                           log10_sd_between_experiments=sd, status="summarized_log_t")
        limits = cfg["comparison_settings"]["rp_acceptance_limits"]
        if limits is not None:
            closed = row["ci_low"] is not None
            row["rp_acceptance_limits"] = list(limits)
            row["rp_acceptance_within"] = bool(limits[0] <= row["ci_low"] and row["ci_high"] <= limits[1]) if closed else None
        rows.append(row)
    return rows
