"""Calibration-curve helpers for ELISA: 5PL fit, parameter covariance and inverse intervals.

The 5PL (Richards) form is Y = B + (T-B) * expit(ln10*s*h*(log10 x - log10 C))^g, which
is the symmetric 4PL when g = 1. As in the 4PL, only log10 C, log10 h and log10 g are
optimized; plateaus come from the shared bounded linear solve. With g != 1, C is not
the half-response concentration, so the fit also reports the half-response midpoint.

Inverse-prediction intervals use the delta method on log10 concentration, combining
the standard-curve parameter covariance s^2 (J'J)^-1 with the response noise of the
unknown wells implied by the same residual model. They are conditional on that model
and the current plate; they include no plate-to-plate, matrix or dilution error.
"""
import numpy as np
from scipy.special import expit, logit
from scipy.stats import f as f_dist, spearmanr, t as t_dist
from .dose_fit import _best, _at_bound, _plateaus, _weighted, fit_dose_curve

LOG_G_BOUNDS = (np.log10(.05), np.log10(20.))


def shape(x, log50, hill, g, direction):
    x = np.asarray(x, float)
    signed = 1. if direction == "increasing" else -1.
    z = np.empty_like(x)
    positive = x > 0
    z[positive] = expit(np.log(10.) * signed * hill * (np.log10(x[positive]) - log50))
    z[~positive] = 0. if signed > 0 else 1.
    return z ** g


def predict(x, fit, direction):
    """Fitted calibration response at canonical concentrations x."""
    return fit["bottom"] + (fit["top"] - fit["bottom"]) * shape(
        x, fit["log10_c_canonical"], abs(fit["hill_slope_signed"]), fit.get("asymmetry", 1.), direction)


def inverse_parameters(y, bottom, top, log50, hill_signed, g=1.):
    """Canonical concentration at response y for explicit parameters, or None outside the plateaus."""
    fraction = (y - bottom) / (top - bottom)
    if not 0 < fraction < 1:
        return None
    value = 10 ** (log50 + logit(fraction ** (1. / g)) / (np.log(10) * hill_signed))
    return float(value) if np.isfinite(value) and value > 0 else None


def _free(cfg):
    return [name for name, fixed in (("bottom", cfg["fit"]["fixed_bottom"]), ("top", cfg["fit"]["fixed_top"]))
            if fixed is None]


def _vector(fit, cfg, model):
    """Free full-parameter vector (plateaus, log10 C, log10 |h|[, log10 g])."""
    v = [fit[name] for name in _free(cfg)] + [fit["log10_c_canonical"], np.log10(abs(fit["hill_slope_signed"]))]
    return np.array(v + ([np.log10(fit["asymmetry"])] if model == "5pl" else []), float)


def _unpack(v, cfg, model):
    names = _free(cfg)
    values = dict(zip(names, v[:len(names)]))
    bottom = values.get("bottom", cfg["fit"]["fixed_bottom"])
    top = values.get("top", cfg["fit"]["fixed_top"])
    rest = v[len(names):]
    g = 10 ** rest[2] if model == "5pl" else 1.
    return float(bottom), float(top), float(rest[0]), float(10 ** rest[1]), float(g)


def _jacobian(fun, v):
    v = np.asarray(v, float)
    base = np.asarray(fun(v), float)
    jac = np.empty((base.size, v.size))
    for j in range(v.size):
        step = 1e-6 * max(1., abs(v[j]))
        up, down = v.copy(), v.copy()
        up[j] += step
        down[j] -= step
        jac[:, j] = (np.asarray(fun(up)) - np.asarray(fun(down))) / (2 * step)
    return jac


def covariance(fit, x, y, cfg, model, direction):
    """s^2 (J'J)^-1 of the weighted residual objective at the fitted solution, or None."""
    v = _vector(fit, cfg, model)

    def residual(p):
        b, t, log50, hill, g = _unpack(p, cfg, model)
        prediction = b + (t - b) * shape(x, log50, hill, g, direction)
        return _weighted(prediction, y, cfg)

    jac = _jacobian(residual, v)
    df = len(y) - v.size
    if df <= 0 or np.linalg.matrix_rank(jac) < v.size:
        return None
    r = residual(v)
    s2 = float(r @ r) / df
    try:
        cov = s2 * np.linalg.inv(jac.T @ jac)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) < 0):
        return None
    return {"vector": v.tolist(), "covariance": cov.tolist(), "residual_variance": s2, "df": df}


def inverse_interval(responses, dilutions, fit, cfg, model, direction, level):
    """Delta-method interval for mean_i(DF_i * x(y_i)); returns (estimate, low, high, se_log10) or None."""
    cov = fit.get("calibration_covariance")
    if cov is None:
        return None
    v, sigma, s2, df = np.array(cov["vector"]), np.array(cov["covariance"]), cov["residual_variance"], cov["df"]
    responses, dilutions = np.asarray(responses, float), np.asarray(dilutions, float)

    def mean_concentration(p, ys):
        b, t, log50, hill, g = _unpack(p, cfg, model)
        signed = hill if direction == "increasing" else -hill
        xs = [inverse_parameters(yy, b, t, log50, signed, g) for yy in ys]
        if any(value is None for value in xs):
            return np.nan
        return float(np.mean(np.asarray(xs) * dilutions))

    estimate = mean_concentration(v, responses)
    if not np.isfinite(estimate) or estimate <= 0:
        return None
    grad_p = _jacobian(lambda p: [mean_concentration(p, responses)], v)[0]
    grad_y = _jacobian(lambda ys: [mean_concentration(v, ys)], responses)[0]
    if not (np.all(np.isfinite(grad_p)) and np.all(np.isfinite(grad_y))):
        return None
    # Each unknown well's response noise follows the calibration residual model; at the inverse
    # estimate the fitted response equals the observed one, so relative weighting scales by y^2.
    y_var = s2 * (responses ** 2 if cfg["fit"]["weighting"] == "relative" else np.ones_like(responses))
    variance = float(grad_p @ sigma @ grad_p + np.sum(grad_y ** 2 * y_var))
    se_log = np.sqrt(variance) / (estimate * np.log(10))
    half = float(t_dist.ppf((1 + level) / 2, df) * se_log)
    centre = np.log10(estimate)
    return estimate, float(10 ** (centre - half)), float(10 ** (centre + half)), float(se_log)


def fit_5pl(std, cfg, dc, factor):
    """Asymmetric 5PL on the plate's included standards, with the 4PL fitted for comparison."""
    direction, f = cfg["direction"], dc["fit"]
    base, _, _ = fit_dose_curve(std, dc)
    x = std.concentration_canonical.to_numpy(float)
    y = std.response.to_numpy(float)
    positive = x[x > 0]
    result = {**base, "calibration_model": "5pl", "asymmetry": None, "log10_c_canonical": None,
              "reportable": False, "status": "failed", "diagnostics": [d for d in base["diagnostics"]
                                                                        if d == "insufficient_distinct_doses_or_residual_df"]}
    n_parameters = 3 + len(_free(dc))
    result["n_parameters"] = n_parameters
    if len(x) - n_parameters < 3 or len(np.unique(positive)) < 7:
        result["diagnostics"].append("insufficient_distinct_doses_for_5pl")
        return result
    if f["weighting"] == "relative" and np.any(y <= 0):
        result["diagnostics"].append("nonpositive_response_for_relative_weighting")
        return result
    lo, hi = float(np.log10(positive.min())), float(np.log10(positive.max()))
    lower = np.array([lo - 3., np.log10(f["hill_bounds"][0]), LOG_G_BOUNDS[0]])
    upper = np.array([hi + 3., np.log10(f["hill_bounds"][1]), LOG_G_BOUNDS[1]])
    scale = max(float(np.max(np.abs(y))), 1e-12) if f["weighting"] == "unweighted" else 1.

    def evaluate(z):
        s = shape(x, z[0], 10 ** z[1], 10 ** z[2], direction)
        b, t = _plateaus([s], [y], dc)
        prediction = b + (t - b) * s
        return prediction, (b, t), _weighted(prediction, y, dc)

    starts = [(a, h, g) for a in np.linspace(lo, hi, 5) for h in (np.log10(.6), 0., np.log10(1.8))
              for g in (np.log10(.4), 0., np.log10(2.5))]
    best = _best(lambda z: evaluate(z)[2] / scale, starts, lower, upper, f["max_nfev"])
    if best is None:
        result["diagnostics"].append("optimizer_failed")
        return result
    prediction, (bottom, top), weighted = evaluate(best.x)
    sse = float(weighted @ weighted)
    df = len(x) - n_parameters
    rmse = float(np.sqrt(np.sum((prediction - y) ** 2) / df))
    g = float(10 ** best.x[2])
    signed = (1 if direction == "increasing" else -1) * float(10 ** best.x[1])
    midpoint = inverse_parameters(bottom + (top - bottom) / 2, bottom, top, float(best.x[0]), signed, g)
    rank = int(np.linalg.matrix_rank(best.jac))
    boundary = bool(_at_bound(best.x, lower, upper) or top - bottom <= scale * 1e-8)
    diagnostics = []
    if boundary:
        diagnostics.append("numerical_boundary_hit")
    if rank < 3:
        diagnostics.append("5pl_parameters_not_locally_identifiable")
    if top - bottom <= 3 * rmse:
        diagnostics.append("low_response_range_relative_to_noise")
    rho, p_value = spearmanr(np.log10(positive), y[x > 0])
    opposite = bool(np.isfinite(rho) and p_value < .05 and rho * (1 if direction == "increasing" else -1) < 0)
    if opposite:
        diagnostics.append("observed_direction_opposite_to_declared")
    comparison = None
    if base["objective_sse"] is not None and sse > 0:
        df4 = base["residual_df"]
        f_stat = max(0., (base["objective_sse"] - sse) / (df4 - df)) / (sse / df)

        def aicc(ss, k):
            n = len(y)
            return n * np.log(ss / n) + 2 * (k + 1) + 2 * (k + 1) * (k + 2) / (n - k - 2)
        comparison = {"test": "extra_sum_of_squares_F_nested_g_equals_1", "sse_4pl": base["objective_sse"],
                      "sse_5pl": sse, "df_4pl": df4, "df_5pl": df, "f_statistic": float(f_stat),
                      "p_value": float(f_dist.sf(f_stat, df4 - df, df)),
                      "aicc_4pl": float(aicc(base["objective_sse"], base["n_parameters"])),
                      "aicc_5pl": float(aicc(sse, n_parameters)),
                      "interpretation": "Informational; the calibration model is the one declared in the config."}
        if comparison["p_value"] >= .05:
            diagnostics.append("asymmetry_not_supported_over_4pl")
    result.update(calibration_model="5pl", asymmetry=g, log10_c_canonical=float(best.x[0]),
                  c_parameter_canonical=float(10 ** best.x[0]), bottom=bottom, top=top, hill_slope_signed=signed,
                  half_response_canonical=midpoint, log10_half_response_canonical=None if midpoint is None else float(np.log10(midpoint)),
                  half_response_input_unit=None if midpoint is None else midpoint / factor,
                  objective_sse=sse, residual_df=df, rmse=rmse, rate_jacobian_rank=rank, numerical_boundary_hit=boundary,
                  model_comparison=comparison, ci_status="not_computed_for_5pl", ci_canonical=[None, None],
                  diagnostics=diagnostics)
    result["reportable"] = bool(not boundary and rank == 3 and top - bottom > 3 * rmse and not opposite and midpoint is not None)
    result["status"] = "limited" if not result["reportable"] else "estimated_with_diagnostics" if diagnostics else "estimated"
    return result
