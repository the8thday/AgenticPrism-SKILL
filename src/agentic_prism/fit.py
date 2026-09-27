"""Single-site fitting and KD intervals under explicit residual assumptions.

Profile-F uses SS(theta)/SSmin = 1 + F(1,n-p; level)/(n-p).
For inverse-SD fitting the supplied SDs define relative variance, not a known
absolute noise scale. Baseline controls in zero_control mode are conditioned on.
"""
from copy import deepcopy
import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.stats import f as f_dist, t as t_dist


def response(c, kd, amplitude, baseline):
    return baseline + amplitude * c / (kd + c)


def fit_curve(g, cfg, compute_ci=True):
    h = g.loc[~g.exclude].copy()
    f = cfg["fit"]
    c, y = h.concentration_M.to_numpy(float), h.response.to_numpy(float)
    result = {"curve_id": str(g.curve_id.iloc[0]), "sample_id": str(g.sample_id.iloc[0]),
              "experiment_id": str(g.experiment_id.iloc[0]), "response_unit": str(g.response_unit.iloc[0]),
              "n_total": len(g), "n_included": len(h), "n_fit": 0,
              "kd_M": None, "amplitude": None, "baseline": None,
              "ci_low_M": None, "ci_high_M": None, "ci_status": "not_computed",
              "ci_method": cfg["uncertainty"]["parameter_ci"], "ci_level": cfg["uncertainty"]["level"],
              "optimizer_success": False, "range_status": "unknown", "identifiability": "unknown",
              "reportable": False, "status": "failed", "diagnostics": [], "residual_scale": f["residual_scale"],
              "baseline_mode": f["baseline_mode"], "interpretation": cfg["assay"]["interpretation"]}
    def fail(reason):
        result["diagnostics"].append(reason)
        return result
    if len(y) == 0:
        return fail("no_included_observations")
    if f["residual_scale"] == "log" and np.any(y <= 0):
        return fail("nonpositive_response_for_log")
    if np.ptp(y) <= max(np.max(np.abs(y)), 1e-300) * 1e-10:
        return fail("flat_response")
    positive = c[c > 0]
    if not len(positive):
        return fail("no_positive_concentration")
    result.update(min_positive_M=float(positive.min()), max_concentration_M=float(c.max()))
    fitted_b = f["baseline_mode"] == "fitted"
    if f["baseline_mode"] == "zero_control":
        if not np.any(c == 0):
            return fail("missing_zero_control")
        baseline = float(y[c == 0].mean())
        result["diagnostics"].append("ci_conditional_on_estimated_zero_baseline")
        mask = c > 0
    else:
        baseline = float(f["baseline_value"]) if f["baseline_mode"] == "fixed" else float(y.min())
        mask = np.ones(len(c), dtype=bool)
    x, obs = c[mask], y[mask]
    p = 3 if fitted_b else 2
    result["n_fit"] = len(x)
    if len(x) <= p or len(np.unique(x)) < p + 1:
        return fail("insufficient_distinct_concentrations_or_df")
    if not fitted_b and f["residual_scale"] == "log" and baseline <= 0:
        return fail("log_model_requires_positive_baseline")
    # Scale response to make optimizer tolerances insensitive to response units.
    scale = max(float(np.max(np.abs(obs))), float(np.ptp(obs)), 1e-300)
    amp0 = max(float(np.ptp(obs)) / scale, .01)
    klo, khi = np.log10(f["kd_bounds_M"])
    # For KD >> C the model is B + (A/KD)*C: keeping the observed slope over the
    # whole KD search domain needs A up to ~(KD_max/C_min) response scales.
    amp_hi = max(8., khi - np.log10(positive.min()) + 3.)
    low = [klo, -12.] + ([-8. if f["residual_scale"] == "linear" else -12.] if fitted_b else [])
    high = [khi, amp_hi] + ([8.] if fitted_b else [])
    low, high = np.array(low), np.array(high)
    weights = h.response_sd.to_numpy(float)[mask] / scale if f["weighting"] == "inverse_sd" else np.ones(len(x))

    def decode(z):
        b = (z[2] * scale if f["residual_scale"] == "linear" else 10.**z[2] * scale) if fitted_b else baseline
        return 10.**z[0], 10.**z[1] * scale, b

    def resid(z):
        pred = response(x, *decode(z))
        if f["residual_scale"] == "log":
            return np.log10(pred) - np.log10(obs)
        return (pred / scale - obs / scale) / weights

    base_start = float(np.clip(baseline / scale, -7., 7.)) if f["residual_scale"] == "linear" else np.log10(max(baseline / scale, 1e-10))
    starts = np.linspace(max(klo + .01, np.log10(positive.min()) - 2), min(khi - .01, np.log10(c.max()) + 2), f["multistart"])
    solutions = []
    for k in starts:
        z0 = np.clip(np.array([k, np.log10(amp0)] + ([base_start] if fitted_b else [])), low + 1e-7, high - 1e-7)
        s = least_squares(resid, z0, bounds=(low, high), max_nfev=f["max_nfev"], ftol=1e-11, xtol=1e-11, gtol=1e-11)
        if s.success and np.isfinite(s.fun).all():
            solutions.append(s)
    if not solutions:
        return fail("optimizer_failed")
    best = min(solutions, key=lambda s: float(s.fun @ s.fun))
    kd, amp, b = map(float, decode(best.x))
    ss = float(best.fun @ best.fun)
    df = len(x) - p
    bound = bool(np.any(np.minimum(best.x - low, high - best.x) < 1e-4))
    jac_cond = float(np.linalg.cond(best.jac))
    result.update(kd_M=kd, amplitude=amp, baseline=b, optimizer_success=True,
                  objective_sse=ss, objective_scale=scale if f["residual_scale"] == "linear" and f["weighting"] == "unweighted" else 1.,
                  residual_df=df, numerical_boundary_hit=bound,
                  jacobian_condition=jac_cond if np.isfinite(jac_cond) else None,
                  multistart_successes=len(solutions), multistart_sse_spread=float(np.ptp([s.fun @ s.fun for s in solutions])),
                  range_status="below_range" if kd < positive.min() else "above_range" if kd > c.max() else "within_range",
                  identifiability="local_full_rank" if np.linalg.matrix_rank(best.jac) == p else "rank_deficient")
    if bound:
        result["diagnostics"].append("numerical_boundary_hit")
    if result["range_status"] != "within_range":
        result["diagnostics"].append("out_of_range_no_automatic_statistical_bound")
    if result["identifiability"] == "rank_deficient":
        result["diagnostics"].append("rank_deficient")
    if c.max() / (kd + c.max()) < .8:
        result["diagnostics"].append("limited_upper_plateau_model_based")
    if compute_ci and cfg["uncertainty"]["parameter_ci"] == "profile_f":
        if ss <= 1e-20:
            result["ci_status"] = "noise_scale_not_estimable"
        else:
            threshold = ss * (1 + float(f_dist.ppf(cfg["uncertainty"]["level"], 1, df)) / df)

            def linear_start(k):
                # Exact (weighted) linear solve for A, B at fixed KD. Starting only from
                # the best-fit amplitude stalls on a flat plateau far from the optimum.
                h_k = x / (10.**k + x)
                target = obs / scale - (0. if fitted_b else baseline / scale)
                design = np.c_[h_k, np.ones_like(h_k)] if fitted_b else h_k[:, None]
                coef = np.linalg.lstsq(design / weights[:, None], target / weights, rcond=None)[0]
                amp = np.log10(max(coef[0], 10.**low[1]))
                if not fitted_b:
                    return np.array([amp])
                if f["residual_scale"] == "linear":
                    return np.array([amp, coef[1]])
                return np.array([amp, np.log10(coef[1])]) if coef[1] > 0 else None

            def profile(k):
                def nuisance(z):
                    return resid(np.r_[k, z])
                opt = None
                for z0 in (best.x[1:], linear_start(k)):
                    if z0 is None:
                        continue
                    s = least_squares(nuisance, np.clip(z0, low[1:] + 1e-7, high[1:] - 1e-7), bounds=(low[1:], high[1:]),
                                      max_nfev=f["max_nfev"], ftol=1e-10, xtol=1e-10, gtol=1e-10)
                    if s.success and np.isfinite(s.fun).all() and (opt is None or s.fun @ s.fun < opt.fun @ opt.fun):
                        opt = s
                if opt is None:
                    raise ArithmeticError("profile nuisance optimization failed")
                if high[1] - opt.x[0] < 1e-4:
                    # A bound-limited profile would close an interval the data leave open.
                    raise ArithmeticError("profile constrained by numerical amplitude bound")
                return float(opt.fun @ opt.fun) - threshold

            def endpoint(edge):
                previous = best.x[0]
                for k in np.linspace(best.x[0], edge, 35)[1:]:
                    if profile(k) >= 0:
                        root = brentq(profile, min(previous, k), max(previous, k), xtol=1e-7)
                        return float(10.**root)
                    previous = k
                return None
            try:
                lo, hi = endpoint(klo), endpoint(khi)
                result["ci_low_M"], result["ci_high_M"] = lo, hi
                result["ci_status"] = "two_sided" if lo is not None and hi is not None else "lower_open" if lo is None and hi is not None else "upper_open" if hi is None and lo is not None else "both_open"
                result["profile_sse_threshold"] = threshold
            except (ValueError, ArithmeticError):
                result["ci_status"] = "profile_failed"
            if result["ci_status"] != "two_sided":
                result["diagnostics"].append("interval_" + result["ci_status"])
    precise_ci = result["ci_status"] == "two_sided" or cfg["uncertainty"]["parameter_ci"] == "none" or result["ci_status"] == "noise_scale_not_estimable"
    result["reportable"] = bool(result["range_status"] == "within_range" and not bound and result["identifiability"] == "local_full_rank" and precise_ci)
    result["status"] = "estimated" if result["reportable"] else "limited"
    if not compute_ci:
        result["reportable"] = False
        result["ci_status"] = "not_requested_for_sensitivity"
        result["status"] = "fit_only" if result["range_status"] == "within_range" and not bound and result["identifiability"] == "local_full_rank" else "limited"
    return result


def summarize(fits, cfg):
    """Equal experiment weights on log(KD); technical curves first collapsed.

    Any failed or non-reportable curve withholds its entire sample summary.
    This avoids silently selecting successes or converting limits to points.
    """
    rows = []
    for sid in dict.fromkeys(f["sample_id"] for f in fits):
        fs = [f for f in fits if f["sample_id"] == sid]
        row = {"sample_id": sid, "n_curves": len(fs), "n_failed_or_limited": sum(not f["reportable"] for f in fs),
               "n_experiments": len(set(f["experiment_id"] for f in fs if f["experiment_id"])),
               "geometric_mean_kd_M": None, "ci_low_M": None, "ci_high_M": None}
        rep = cfg["replicates"]
        if rep["independent_unit"] == "none":
            row["status"] = "independent_units_unconfirmed"
        elif not rep["conditions_comparable"]:
            row["status"] = "comparability_unconfirmed"
        elif row["n_failed_or_limited"]:
            row["status"] = "withheld_incomplete_estimates"
        else:
            logs = [np.mean([np.log10(f["kd_M"]) for f in fs if f["experiment_id"] == e]) for e in dict.fromkeys(f["experiment_id"] for f in fs)]
            row["geometric_mean_kd_M"] = float(10.**np.mean(logs))
            row["status"] = "one_experiment_no_between_experiment_ci"
            if len(logs) > 1:
                half = float(t_dist.ppf((1 + cfg["uncertainty"]["level"]) / 2, len(logs) - 1) * np.std(logs, ddof=1) / np.sqrt(len(logs)))
                row.update(ci_low_M=float(10.**(np.mean(logs) - half)), ci_high_M=float(10.**(np.mean(logs) + half)), status="summarized_log_t")
        rows.append(row)
    return rows


def sensitivity(g, cfg, primary):
    rows = []
    for scenario in cfg["sensitivity"]:
        cc, gg = deepcopy(cfg), g.copy()
        if scenario == "alternate_loss":
            if cc["fit"]["weighting"] != "unweighted":
                rows.append({"curve_id": primary["curve_id"], "scenario": scenario, "status": "unsupported_with_weighting"})
                continue
            cc["fit"]["residual_scale"] = "log" if cc["fit"]["residual_scale"] == "linear" else "linear"
        else:
            gg.loc[gg.concentration_M == gg.loc[~gg.exclude, "concentration_M"].max(), "exclude"] = True
        alt = fit_curve(gg, cc, compute_ci=False)
        rows.append({"curve_id": primary["curve_id"], "scenario": scenario,
                     "status": alt["status"], "kd_M": alt["kd_M"],
                     "ratio_to_primary": alt["kd_M"] / primary["kd_M"] if alt["kd_M"] and primary["kd_M"] else None})
    return rows
