"""Parallelism by equivalence for reference-vs-test 4PL curves (USP <1032>/<1034> style).

The unrestricted model fits each curve separately. Hill-slope ratio (test/reference) and
plateau differences (test - reference) get Wald intervals from the full-parameter
covariance with a pooled residual variance, (SSE_ref + SSE_test) / (n - 2k), on the fit
objective's own scale. Parallelism is concluded only when every declared interval lies
inside its predeclared margin (a two one-sided tests decision at (1 - level) / 2 per side).
Margins must come from the laboratory's historical data; the software never proposes them.
"""
import numpy as np
from scipy.stats import t as t_dist
from .calibration import _jacobian
from .dose_fit import _weighted, response


def _curve_jacobian(fit, x, y, cfg):
    f = cfg["fit"]
    free = [name for name, fixed in (("bottom", f["fixed_bottom"]), ("top", f["fixed_top"])) if fixed is None]
    v = np.array([fit[name] for name in free] + [fit["log10_half_response_canonical"],
                                                 np.log10(abs(fit["hill_slope_signed"]))], float)

    def residual(p):
        values = dict(zip(free, p[:len(free)]))
        bottom = values.get("bottom", f["fixed_bottom"])
        top = values.get("top", f["fixed_top"])
        prediction = response(x, bottom, top, p[len(free)], 10 ** p[len(free) + 1], cfg["assay"]["direction"])
        return _weighted(prediction, y, cfg)
    return free, v, _jacobian(residual, v)


def parameter_equivalence(rf, tf, curves, cfg):
    rule = cfg["comparison_settings"]["equivalence"]
    out = {"method": "equivalence_of_unrestricted_parameters_pooled_variance", "confidence_level": rule["confidence_level"],
           "rationale": rule["rationale"], "tests": [], "equivalent": False, "status": "not_estimable"}
    if rf["objective_sse"] is None or tf["objective_sse"] is None:
        return out
    parts = [_curve_jacobian(fit, x, y, cfg) for fit, (x, y) in zip((rf, tf), curves)]
    free = parts[0][0]
    df = sum(len(y) for _, y in curves) - 2 * (len(free) + 2)
    if df <= 0:
        return out
    s2 = (rf["objective_sse"] + tf["objective_sse"]) / df
    covs = []
    for _, v, jac in parts:
        if np.linalg.matrix_rank(jac) < v.size:
            return out
        covs.append(s2 * np.linalg.inv(jac.T @ jac))
    q = float(t_dist.ppf((1 + rule["confidence_level"]) / 2, df))
    index = {name: i for i, name in enumerate(free)}
    index["log10_hill"] = len(free) + 1

    def wald(name):
        i = index[name]
        estimate = parts[1][1][i] - parts[0][1][i]
        se = float(np.sqrt(covs[0][i, i] + covs[1][i, i]))
        return float(estimate), float(estimate - q * se), float(estimate + q * se), se

    rows = []
    d, lo, hi, se = wald("log10_hill")
    rows.append({"parameter": "hill_ratio_test_over_reference", "estimate": float(10 ** d),
                 "ci": [float(10 ** lo), float(10 ** hi)], "limits": rule["hill_ratio_limits"], "standard_error_log10": se})
    for name in ("bottom", "top"):
        limits = rule[f"{name}_difference_limits"]
        if limits is None:
            continue
        if name not in index:
            rows.append({"parameter": f"{name}_difference_test_minus_reference", "estimate": 0., "ci": [None, None],
                         "limits": limits, "within": None, "note": "plateau fixed in both fits; not tested"})
            continue
        d, lo, hi, se = wald(name)
        rows.append({"parameter": f"{name}_difference_test_minus_reference", "estimate": d, "ci": [lo, hi],
                     "limits": limits, "standard_error": se})
    for row in rows:
        if "within" not in row:
            row["within"] = bool(row["limits"][0] <= row["ci"][0] and row["ci"][1] <= row["limits"][1])
    tested = [r["within"] for r in rows if r["within"] is not None]
    out.update(tests=rows, residual_df=df, pooled_residual_variance=float(s2),
               equivalent=bool(tested and all(tested)), status="evaluated")
    return out


def rp_acceptance(interval, level, limits, rationale):
    """Whether a closed RP interval lies entirely inside predeclared acceptance limits."""
    closed = interval[0] is not None and interval[1] is not None
    return {"limits": limits, "rationale": rationale, "interval": list(interval), "interval_level": level,
            "within": bool(limits[0] <= interval[0] and interval[1] <= limits[1]) if closed else None,
            "decision_basis": "entire RP interval inside limits; an open or missing interval is not a pass"}
