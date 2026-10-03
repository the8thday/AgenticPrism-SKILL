"""Opt-in dose-response shapes beyond the symmetric 4PL (0.13.1).

Asymmetric 5PL. Y = Bottom + (Top - Bottom) * expit(ln10*s*h*(log10 x - log10 C))^g, the
Richards form used by ELISA calibration. The optimizer works on (log10 EC50, log10 h,
log10 g), where EC50 is the relative midpoint between the plateaus and
log10 C = log10 EC50 - logit(0.5^(1/g)) / (ln10*s*h); so the profile-F interval is on
EC50 itself. Plateaus come from the shared bounded linear solve, honoring fixed
plateaus and relative weighting exactly as the 4PL does.

Bell-shaped (Prism's "Bell-shaped dose-response"): two 4PL phases sharing the middle
plateau, Y = Plateau1 + (Dip - Plateau1)*S1 + (Plateau2 - Dip)*S2, where S1, S2 are
rising 4PL shapes with midpoints EC50_1 < EC50_2 and Hill slopes h1, h2; the declared
direction is that of the first phase (increasing: bell; decreasing: U-shape) and sets the
sign of Dip - Plateau. Reported Hill slopes carry the sign of the first phase. Only log10 EC50_1, log10 h1, the log10 gap to EC50_2 and
log10 h2 are optimized; the three plateaus are solved by ordinary least squares.
"Dip" is Prism's name for the middle plateau (the peak of a bell). Each EC50 is a
model parameter: it equals a half-maximal concentration only when the curve reaches
the middle plateau, which is a reportability gate.
"""
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq
from scipy.special import logit
from scipy.stats import f as f_dist, spearmanr
from .calibration import LOG_G_BOUNDS, shape as richards
from .dose_fit import (_at_bound, _best, _bounds, _noise_floor, _plateaus, _weighted, fit_dose_curve,
                       free_plateaus, response)
from .dose_schema import DOSE_UNITS

REACHED_FRACTION = .9  # the fitted curve must get this close to the middle plateau from both sides


def _log_c(log50, hill, g, direction):
    signed = hill if direction == "increasing" else -hill
    return log50 - logit(.5 ** (1. / g)) / (np.log(10.) * signed)


def five_pl_shape(x, log50, hill, g, direction):
    """Richards shape normalized so the value at EC50 (= 10**log50) is exactly 0.5."""
    return richards(x, _log_c(log50, hill, g, direction), hill, g, direction)


def bell_basis(x, log50_1, hill_1, log50_2, hill_2):
    """Columns for (Plateau1, Dip, Plateau2). Phase shapes always rise with dose, so Plateau1 is the
    low-dose plateau for either direction; the declared direction only fixes the sign of Dip - Plateau."""
    s1 = response(x, 0., 1., log50_1, hill_1, "increasing")
    s2 = response(x, 0., 1., log50_2, hill_2, "increasing")
    return np.column_stack((1. - s1, s1 - s2, s2))


def _profile_interval(profile_sse, centre, lower, upper, threshold):
    """Scan from the optimum to each search bound; return (root or None) on both sides."""
    def endpoint(bound):
        previous = centre
        for current in np.linspace(centre, bound, 45)[1:]:
            if profile_sse(float(current)) > threshold:
                return float(brentq(lambda v: profile_sse(float(v)) - threshold, min(previous, current),
                                    max(previous, current), xtol=1e-9))
            previous = float(current)
        return None
    return endpoint(lower), endpoint(upper)


def _base(group, cfg, n_parameters):
    f = cfg["fit"]
    used = group.loc[~group.exclude]
    return {"curve_id": str(group.curve_id.iloc[0]), "sample_id": str(group.sample_id.iloc[0]),
            "experiment_id": str(group.experiment_id.iloc[0]), "endpoint": cfg["assay"]["endpoint"],
            "direction": cfg["assay"]["direction"], "model": cfg["model"], "status": "failed", "reportable": False,
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


def _range(value, low, high):
    return "within_range" if low <= value <= high else "below_range" if value < low else "above_range"


def _outputs(group, result, predict, low_dose, high_dose):
    curve_id = result["curve_id"]
    dose_grid = np.geomspace(low_dose, high_dose, 200)
    grid = [{"curve_id": curve_id, "concentration_canonical": float(c), "predicted_response": float(p)}
            for c, p in zip(dose_grid, predict(dose_grid))]
    all_pred = predict(group.concentration_canonical.to_numpy(float))
    observations = [{"observation_id": str(row.observation_id), "curve_id": curve_id,
                     "concentration_canonical": float(row.concentration_canonical),
                     "observed_response": float(row.response), "predicted_response": float(p),
                     "residual": float(row.response - p), "used_in_fit": not bool(row.exclude)}
                    for (_, row), p in zip(group.iterrows(), all_pred)]
    return grid, observations


def _opposite(x, y, sign):
    positive = x > 0
    rho, p_value = spearmanr(np.log10(x[positive]), y[positive])
    return (float(rho) if np.isfinite(rho) else None), bool(np.isfinite(rho) and p_value < .05 and rho * sign < 0)


def fit_dose_5pl(group, cfg):
    f, direction = cfg["fit"], cfg["assay"]["direction"]
    sign = 1 if direction == "increasing" else -1
    n_parameters = 3 + free_plateaus(cfg)
    result = _base(group, cfg, n_parameters)
    result.update(asymmetry=None, c_parameter_canonical=None, model_comparison=None)
    used = group.loc[~group.exclude]
    x, y = used.concentration_canonical.to_numpy(float), used.response.to_numpy(float)
    positive = x[x > 0]
    if len(x) - n_parameters < 3 or len(np.unique(positive)) < 7:
        result["diagnostics"].append("insufficient_distinct_doses_for_5pl")
        return result, [], []
    y_scale = max(float(np.max(np.abs(y))), 1e-12)
    if float(np.ptp(y)) <= y_scale * 1e-9:
        result["diagnostics"].append("flat_response")
        return result, [], []
    if f["weighting"] == "relative" and np.any(y <= 0):
        result["diagnostics"].append("nonpositive_response_for_relative_weighting")
        return result, [], []
    bounds, (hill_min, hill_max), lo, hi = _bounds(cfg, [positive])
    lower = np.array([bounds[0], hill_min, LOG_G_BOUNDS[0]])
    upper = np.array([bounds[1], hill_max, LOG_G_BOUNDS[1]])
    scale = y_scale if f["weighting"] == "unweighted" else 1.

    def evaluate(z):
        s = five_pl_shape(x, z[0], 10. ** z[1], 10. ** z[2], direction)
        b, t = _plateaus([s], [y], cfg)
        prediction = b + (t - b) * s
        return prediction, (b, t), _weighted(prediction, y, cfg)

    # Steep starts matter: without h = 4 the S.alba benchmark stopped at the g upper bound (0.13.1 initial run).
    starts = [(a, h, g) for a in np.linspace(lo, hi, 5) for h in (np.log10(.6), 0., np.log10(1.8), np.log10(4.))
              for g in (np.log10(.4), 0., np.log10(2.5))]
    best = _best(lambda z: evaluate(z)[2] / scale, starts, lower, upper, f["max_nfev"])
    if best is None:
        result["diagnostics"].append("optimizer_failed")
        return result, [], []
    prediction, (bottom, top), weighted = evaluate(best.x)
    sse = float(weighted @ weighted)
    df = len(x) - n_parameters
    rmse = float(np.sqrt(np.sum((prediction - y) ** 2) / df))
    log50, hill, g = float(best.x[0]), float(10. ** best.x[1]), float(10. ** best.x[2])
    midpoint = float(10. ** log50)
    amp = top - bottom
    rank = int(np.linalg.matrix_rank(best.jac))
    boundary = bool(_at_bound(best.x, lower, upper) or amp <= y_scale * 1e-8)
    low_dose, high_dose = float(positive.min()), float(positive.max())
    result.update(half_response_canonical=midpoint, half_response_input_unit=midpoint / DOSE_UNITS[result["input_unit"]],
                  log10_half_response_canonical=log50, bottom=bottom, top=top, hill_slope_signed=float(sign * hill),
                  asymmetry=g, c_parameter_canonical=float(10. ** _log_c(log50, hill, g, direction)),
                  objective_sse=sse, residual_df=df, rmse=rmse, rate_jacobian_rank=rank, numerical_boundary_hit=boundary,
                  range_status=_range(midpoint, low_dose, high_dose),
                  lowest_tested_positive_canonical=low_dose, highest_tested_canonical=high_dose)
    if f["weighting"] == "relative":
        result["weighted_rmse"] = float(np.sqrt(sse / df))
    d = result["diagnostics"]
    if boundary:
        d.append("numerical_boundary_hit")
    if rank < 3:
        d.append("5pl_parameters_not_locally_identifiable")
    if result["range_status"] != "within_range":
        d.append("half_response_outside_tested_range")
    if amp <= 3 * rmse:
        d.append("low_response_range_relative_to_noise")
    curve = lambda c: five_pl_shape(np.asarray(c, float), log50, hill, g, direction)
    edges = curve([low_dose, high_dose])
    if min(edges) > .1 or max(edges) < .9:
        d.append("plateau_not_well_covered")
    if sum(.1 < v < .9 for v in curve(np.unique(positive))) < 4:
        d.append("few_transition_doses_for_asymmetric_slope")
    rho, opposite = _opposite(x, y, sign)
    result["observed_trend_spearman"] = rho
    if opposite:
        d.append("observed_direction_opposite_to_declared")
    base, _, _ = fit_dose_curve(group, cfg)
    if base["objective_sse"] is not None and sse > 0:
        df4 = base["residual_df"]
        f_stat = max(0., (base["objective_sse"] - sse) / (df4 - df)) / (sse / df)
        n = len(y)
        aicc = lambda ss, k: n * np.log(ss / n) + 2 * (k + 1) + 2 * (k + 1) * (k + 2) / (n - k - 2)
        result["model_comparison"] = {"test": "extra_sum_of_squares_F_nested_g_equals_1", "sse_4pl": base["objective_sse"],
                                      "sse_5pl": sse, "df_4pl": df4, "df_5pl": df, "f_statistic": float(f_stat),
                                      "p_value": float(f_dist.sf(f_stat, df4 - df, df)),
                                      "aicc_4pl": float(aicc(base["objective_sse"], base["n_parameters"])),
                                      "aicc_5pl": float(aicc(sse, n_parameters)),
                                      "interpretation": "Informational; the model is the one declared before analysis."}
        if result["model_comparison"]["p_value"] >= .05:
            d.append("asymmetry_not_supported_over_4pl")

    nuisance_at_bound = False
    if cfg["uncertainty"]["method"] == "profile_f":
        if sse <= _noise_floor(y, cfg):
            result["ci_status"] = "noise_scale_not_estimable"
        else:
            threshold = sse * (1. + float(f_dist.ppf(cfg["uncertainty"]["level"], 1, df)) / df)
            state = {"z": np.array(best.x[1:])}

            @lru_cache(maxsize=512)
            def profile(value):
                fun = lambda z: evaluate([value, z[0], z[1]])[2] / scale
                s = _best(fun, [state["z"], np.array(best.x[1:])], lower[1:], upper[1:], f["max_nfev"])
                if s is None:
                    return np.inf, tuple(state["z"])
                state["z"] = s.x
                return float(s.fun @ s.fun) * scale ** 2, tuple(s.x)

            low, high = _profile_interval(lambda v: profile(v)[0], log50, bounds[0], bounds[1], threshold)
            if low is not None and high is not None:
                result["ci_canonical"] = [float(10. ** low), float(10. ** high)]
                result["ci_status"] = "two_sided_profile_f"
                nuisance_at_bound = any(_at_bound(profile(v)[1], lower[1:], upper[1:]) for v in (low, high))
                if nuisance_at_bound:
                    d.append("profile_nuisance_at_boundary")
            else:
                result["ci_status"] = "both_open" if low is None and high is None else "lower_open" if low is None else "upper_open"
                d.append("profile_interval_open")
    else:
        result["ci_status"] = "not_requested"
    result["reportable"] = (not boundary and rank == 3 and result["range_status"] == "within_range" and amp > 3 * rmse
                            and not opposite and not nuisance_at_bound
                            and result["ci_status"] in ("two_sided_profile_f", "not_requested"))
    result["status"] = "limited" if not result["reportable"] else "estimated_with_diagnostics" if d else "estimated"
    result["diagnostics"] = list(dict.fromkeys(d))
    grid, observations = _outputs(group, result, lambda c: bottom + amp * curve(c), low_dose, high_dose)
    return result, grid, observations


def fit_dose_bell(group, cfg):
    f, direction = cfg["fit"], cfg["assay"]["direction"]
    sign = 1 if direction == "increasing" else -1
    result = _base(group, cfg, 7)
    result.update(plateau_1=None, middle_plateau=None, plateau_2=None, second_phase=None, peak=None)
    used = group.loc[~group.exclude]
    x, y = used.concentration_canonical.to_numpy(float), used.response.to_numpy(float)
    positive = x[x > 0]
    if len(x) - 7 < 3 or len(np.unique(positive)) < 10:
        result["diagnostics"].append("insufficient_distinct_doses_for_bell_shape")
        return result, [], []
    y_scale = max(float(np.max(np.abs(y))), 1e-12)
    if float(np.ptp(y)) <= y_scale * 1e-9:
        result["diagnostics"].append("flat_response")
        return result, [], []
    bounds, (hill_min, hill_max), lo, hi = _bounds(cfg, [positive])
    max_gap = float(bounds[1] - bounds[0])
    # z = (log10 EC50_1, log10 h1, gap = log10 EC50_2 - log10 EC50_1, log10 h2)
    lower = np.array([bounds[0], hill_min, 0., hill_min])
    upper = np.array([bounds[1], hill_max, max_gap, hill_max])

    def evaluate(z):
        basis = bell_basis(x, z[0], 10. ** z[1], z[0] + z[2], 10. ** z[3])
        plateaus = np.linalg.lstsq(basis, y, rcond=None)[0]
        prediction = basis @ plateaus
        return prediction, plateaus, prediction - y

    width = hi - lo
    starts = [(lo + a * width, h, gap * width, h) for a in (.15, .3, .45) for gap in (.35, .6) for h in (0., np.log10(2.))]
    best = _best(lambda z: evaluate(z)[2] / y_scale, starts, lower, upper, f["max_nfev"])
    if best is None:
        result["diagnostics"].append("optimizer_failed")
        return result, [], []
    z = best.x
    prediction, (p1, dip, p2), residual = evaluate(z)
    sse = float(residual @ residual)
    df = len(x) - 7
    rmse = float(np.sqrt(sse / df))
    log50_1, log50_2 = float(z[0]), float(z[0] + z[2])
    h1, h2 = float(10. ** z[1]), float(10. ** z[3])
    low_dose, high_dose = float(positive.min()), float(positive.max())
    rank = int(np.linalg.matrix_rank(best.jac))
    boundary = bool(_at_bound(z, lower, upper))
    factor = DOSE_UNITS[result["input_unit"]]
    predict = lambda c: bell_basis(np.asarray(c, float), log50_1, h1, log50_2, h2) @ np.array([p1, dip, p2])
    fine = np.geomspace(low_dose, high_dose, 2000)
    curve = predict(fine)
    k = int(np.argmax(sign * curve))
    peak_response = float(curve[k])
    reached = [float((peak_response - p) / (dip - p)) if dip != p else 0. for p in (p1, p2)]
    result.update(half_response_canonical=float(10. ** log50_1), half_response_input_unit=float(10. ** log50_1) / factor,
                  log10_half_response_canonical=log50_1, hill_slope_signed=float(sign * h1),
                  plateau_1=float(p1), middle_plateau=float(dip), plateau_2=float(p2),
                  second_phase={"half_response_canonical": float(10. ** log50_2), "half_response_input_unit": float(10. ** log50_2) / factor,
                                "log10_half_response_canonical": log50_2, "hill_slope_signed": float(sign * h2),
                                "range_status": _range(10. ** log50_2, low_dose, high_dose), "ci_canonical": [None, None],
                                "ci_status": "not_computed"},
                  peak={"concentration_canonical": float(fine[k]), "concentration_input_unit": float(fine[k]) / factor,
                        "response": peak_response, "fraction_of_middle_plateau_reached": reached,
                        "note": "Model-based extremum on the tested range; not an observed maximum."},
                  objective_sse=sse, residual_df=df, rmse=rmse, rate_jacobian_rank=rank, numerical_boundary_hit=boundary,
                  range_status=_range(10. ** log50_1, low_dose, high_dose),
                  lowest_tested_positive_canonical=low_dose, highest_tested_canonical=high_dose)
    d = result["diagnostics"]
    if boundary:
        d.append("numerical_boundary_hit")
    if rank < 4:
        d.append("bell_parameters_not_locally_identifiable")
    if result["range_status"] != "within_range" or result["second_phase"]["range_status"] != "within_range":
        d.append("half_response_outside_tested_range")
    consistent = sign * (dip - p1) > 0 and sign * (dip - p2) > 0
    if not consistent:
        d.append("phases_inconsistent_with_declared_bell_direction")
    if min(abs(dip - p1), abs(dip - p2)) <= 3 * rmse:
        d.append("low_response_range_relative_to_noise")
    if consistent and min(reached) < REACHED_FRACTION:
        d.append("middle_plateau_not_reached")
    for log50, hill in ((log50_1, h1), (log50_2, h2)):
        s = response(np.unique(positive), 0., 1., log50, hill, "increasing")
        if sum((.1 < s) & (s < .9)) < 2:
            d.append("few_transition_doses_in_a_phase")
            break

    cis_ok = True
    if cfg["uncertainty"]["method"] == "profile_f":
        if sse <= _noise_floor(y, cfg):
            result["ci_status"] = result["second_phase"]["ci_status"] = "noise_scale_not_estimable"
            cis_ok = False
        else:
            threshold = sse * (1. + float(f_dist.ppf(cfg["uncertainty"]["level"], 1, df)) / df)
            nuisance_lower = lower[1:]
            for phase in (1, 2):
                state = {"z": np.array(z[1:])}

                def full(value, w, phase=phase):
                    # w = (log10 h1, gap, log10 h2); phase 2 fixes log10 EC50_2 = value.
                    return [value, w[0], w[1], w[2]] if phase == 1 else [value - w[1], w[0], w[1], w[2]]

                @lru_cache(maxsize=512)
                def profile(value, phase=phase):
                    gap_max = max_gap if phase == 1 else float(value - bounds[0])
                    ub = np.array([hill_max, gap_max, hill_max])
                    if gap_max <= 0:
                        return np.inf, None
                    start = np.minimum(state["z"], ub - 1e-6)
                    s = _best(lambda w: evaluate(full(value, w))[2] / y_scale, [start, np.minimum(z[1:], ub - 1e-6)],
                              nuisance_lower, ub, f["max_nfev"])
                    if s is None:
                        return np.inf, None
                    state["z"] = s.x
                    return float(s.fun @ s.fun) * y_scale ** 2, tuple(s.x)

                centre = log50_1 if phase == 1 else log50_2
                low, high = _profile_interval(lambda v: profile(v)[0], centre, bounds[0], bounds[1], threshold)
                target = result if phase == 1 else result["second_phase"]
                if low is not None and high is not None:
                    target["ci_canonical"] = [float(10. ** low), float(10. ** high)]
                    target["ci_status"] = "two_sided_profile_f"
                else:
                    target["ci_status"] = "both_open" if low is None and high is None else "lower_open" if low is None else "upper_open"
                    d.append("profile_interval_open")
                    cis_ok = False
    else:
        result["ci_status"] = result["second_phase"]["ci_status"] = "not_requested"
    blocking = {"numerical_boundary_hit", "bell_parameters_not_locally_identifiable", "half_response_outside_tested_range",
                "phases_inconsistent_with_declared_bell_direction", "low_response_range_relative_to_noise",
                "middle_plateau_not_reached", "profile_interval_open"}
    result["reportable"] = bool(cis_ok and not blocking & set(d))
    result["status"] = "limited" if not result["reportable"] else "estimated_with_diagnostics" if d else "estimated"
    result["diagnostics"] = list(dict.fromkeys(d))
    grid, observations = _outputs(group, result, predict, low_dose, high_dose)
    return result, grid, observations


def second_phase_fits(fits):
    """Bell-shaped fits viewed through their second-phase EC50, for the shared summary rules."""
    return [{**f, "log10_half_response_canonical": (f["second_phase"] or {}).get("log10_half_response_canonical")} for f in fits]
