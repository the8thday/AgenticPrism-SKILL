"""Plate-specific 4PL/5PL calibration, standard back-calculation QC and bounded inverse.

Defaults screen per-level mean recovery (80-120%, 75-125% at passing edges),
replicate CV and the proportion/count of passing levels. These software rules
are not a complete ICH M10 validation. Legacy LLOQ/ULOQ fields identify this
plate's passing range, not experimentally established assay limits. Failing
levels are reported; exclusions require an explicit reason and a new run.
Unknown intervals (delta method) are conditional on the plate's calibration model.
Dilution linearity compares dilution-corrected results of one sample on one plate.
"""
from copy import deepcopy
import math
import numpy as np
import pandas as pd
from scipy.special import logit
from .dose_schema import DOSE_UNITS, MASS_UNITS, resolve_dose_config
from scipy.stats import linregress
from .dose_fit import fit_dose_curve
from .calibration import covariance, fit_5pl, inverse_interval, inverse_parameters

DEFAULTS = {
    "analysis_type": "elisa_quantification", "schema_version": 1, "input": None,
    "source": "User-supplied ELISA plate", "direction": None,
    "assay": {"analyte": "", "response_definition": "", "matrix": "", "rationale": ""},
    "fit": {"model": "4pl", "weighting": "unweighted", "fixed_bottom": None, "fixed_top": None},
    "uncertainty": {"unknown_interval": "delta_method", "level": .95},
    "dilution_linearity": {"required": False, "summary_policy": "all_wells_required",
                           "recovery_limits_percent": [80., 120.], "max_cv_percent": 20.},
    "qc": {"recovery_limits_percent": [80., 120.], "edge_recovery_limits_percent": [75., 125.],
           "max_cv_percent": 20., "min_passing_fraction": .75, "min_passing_levels": 6},
    "independent_qc": {"required": False, "preparation_independent": None, "rationale": "",
                       "recovery_limits_percent": [80.,120.], "max_cv_percent":20.,
                       "min_levels":3, "min_replicates":2},
    "report": {"plot_style": "prism_like"},
}


def resolve_config(raw):
    if not isinstance(raw, dict):
        raise ValueError("ELISA config must be an object")
    cfg = deepcopy(DEFAULTS)
    if set(raw) - set(cfg):
        raise ValueError(f"Unsupported ELISA setting: {sorted(set(raw)-set(cfg))}")
    for key, value in raw.items():
        if isinstance(cfg[key], dict):
            if not isinstance(value, dict) or set(value) - set(cfg[key]):
                raise ValueError(f"Unsupported ELISA {key} settings")
            cfg[key].update(value)
        else:
            cfg[key] = value
    if cfg["analysis_type"] != "elisa_quantification" or cfg["schema_version"] != 1:
        raise ValueError("Unsupported ELISA analysis type or schema")
    if not isinstance(cfg["input"], str) or not cfg["input"] or cfg["direction"] not in ("increasing", "decreasing"):
        raise ValueError("Declare the input file and observed calibration direction")
    if any(not isinstance(v, str) or not v.strip() for v in cfg["assay"].values()):
        raise ValueError("Document analyte, response, matrix and calibration rationale")
    if cfg["fit"]["weighting"] not in ("unweighted", "relative"):
        raise ValueError("Unsupported weighting")
    if cfg["fit"]["model"] not in ("4pl", "5pl"):
        raise ValueError("fit.model must be 4pl or 5pl")
    u = cfg["uncertainty"]
    if u["unknown_interval"] not in ("delta_method", "none") or isinstance(u["level"], bool) \
            or not isinstance(u["level"], (int, float)) or not .5 < u["level"] < 1:
        raise ValueError("uncertainty.unknown_interval must be delta_method or none, with level in (0.5, 1)")
    dl = cfg["dilution_linearity"]
    if type(dl["required"]) is not bool or dl["summary_policy"] not in ("all_wells_required", "in_range_dilutions_with_linearity"):
        raise ValueError("dilution_linearity.required must be boolean; summary_policy all_wells_required or in_range_dilutions_with_linearity")
    limits = dl["recovery_limits_percent"]
    if not isinstance(limits, list) or len(limits) != 2 or any(type(v) not in (int, float) for v in limits) or not 0 < limits[0] < 100 < limits[1]:
        raise ValueError("dilution_linearity.recovery_limits_percent must be [low, high] bracketing 100")
    if type(dl["max_cv_percent"]) not in (int, float) or not dl["max_cv_percent"] > 0:
        raise ValueError("dilution_linearity.max_cv_percent must be positive")
    if cfg["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    q = cfg["qc"]
    for key in ("recovery_limits_percent", "edge_recovery_limits_percent"):
        v = q[key]
        if not isinstance(v, list) or len(v) != 2 or any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in v) \
                or not 0 < v[0] < 100 < v[1]:
            raise ValueError(f"qc.{key} must be [low, high] percentages bracketing 100")
    inner, edge = q["recovery_limits_percent"], q["edge_recovery_limits_percent"]
    if not (edge[0] <= inner[0] and inner[1] <= edge[1]):
        raise ValueError("qc.edge_recovery_limits_percent must contain recovery_limits_percent")
    if q["max_cv_percent"] is not None and (isinstance(q["max_cv_percent"], bool) or not isinstance(q["max_cv_percent"], (int, float))
                                            or q["max_cv_percent"] <= 0):
        raise ValueError("qc.max_cv_percent must be positive or null")
    if isinstance(q["min_passing_fraction"], bool) or not isinstance(q["min_passing_fraction"], (int, float)) \
            or not 0 < q["min_passing_fraction"] <= 1:
        raise ValueError("qc.min_passing_fraction must be in (0, 1]")
    if type(q["min_passing_levels"]) is not int or q["min_passing_levels"] < 2:
        raise ValueError("qc.min_passing_levels must be an integer >= 2")
    # Delegate fit-option validation to the identical numeric 4PL implementation.
    iq=cfg["independent_qc"]
    if type(iq["required"]) is not bool or iq["preparation_independent"] not in (None, True, False):
        raise ValueError("Invalid independent QC declarations")
    if iq["preparation_independent"] is not None and type(iq["preparation_independent"]) is not bool:
        raise ValueError("QC independence must be a boolean")
    if not isinstance(iq["rationale"],str):
        raise ValueError("QC rationale must be text")
    limits=iq["recovery_limits_percent"]
    if not isinstance(limits,list) or len(limits)!=2 or not all(type(v) in (int,float) and np.isfinite(v) for v in limits) or not 0<limits[0]<100<limits[1]:
        raise ValueError("Invalid independent QC recovery limits")
    if type(iq["max_cv_percent"]) not in (int,float) or not np.isfinite(iq["max_cv_percent"]) or iq["max_cv_percent"]<=0:
        raise ValueError("Invalid independent QC CV limit")
    if any(type(iq[k]) is not int or iq[k]<2 for k in ("min_levels","min_replicates")):
        raise ValueError("Independent QC requires at least two levels and replicates")
    dose_config(cfg)
    return cfg


def dose_config(cfg):
    return resolve_dose_config({
        "analysis_type": "dose_response_4pl", "input": cfg["input"],
        "assay": {"endpoint": "EC50", "direction": cfg["direction"],
                  "tested_agent": cfg["assay"]["analyte"],
                  "response_definition": cfg["assay"]["response_definition"],
                  "relative_half_response_supported": True,
                  "rationale": cfg["assay"]["rationale"]},
        "fit": {k: v for k, v in cfg["fit"].items() if k != "model"}, "uncertainty": {"method": "profile_f", "level": .95},
    })


def load_data(path):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"plate_id", "well_id", "role", "sample_id", "concentration",
                "concentration_unit", "dilution_factor", "response", "response_unit"}
    if d.empty or required - set(d) or d.columns.duplicated().any():
        raise ValueError(f"Invalid ELISA table; missing {sorted(required-set(d))}")
    for key in ("plate_id", "well_id", "role", "sample_id", "concentration_unit", "response_unit"):
        if d[key].str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
    if not d.role.isin(("standard", "unknown", "qc")).all():
        raise ValueError("role must be standard, unknown or qc")
    if d.duplicated(["plate_id", "well_id"]).any():
        raise ValueError("Duplicate plate/well identifier")
    if not d.concentration_unit.isin(DOSE_UNITS).all():
        raise ValueError("Unsupported concentration unit")
    d["response"] = pd.to_numeric(d.response, errors="raise")
    if not np.isfinite(d.response).all():
        raise ValueError("Nonfinite response")
    d["exclude"] = d.get("exclude", pd.Series("false", index=d.index)).str.lower()
    d["exclusion_reason"] = d.get("exclusion_reason", pd.Series("", index=d.index))
    if not d.exclude.isin(("true", "false")).all() or ((d.exclude == "true") & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Exclusions require true/false and a reason")
    d["exclude"] = d.exclude.eq("true")
    standards = d.role == "standard"
    known=d.role.isin(("standard","qc"))
    if d.loc[known, "concentration"].eq("").any() or d.loc[~known, "concentration"].ne("").any():
        raise ValueError("Standard concentration is required; unknown concentration must be blank")
    if d.loc[standards, "dilution_factor"].ne("").any():
        raise ValueError("Standard dilution_factor must be blank")
    d["concentration"] = pd.to_numeric(d.concentration.where(known, "0"), errors="raise")
    d["dilution_factor"] = pd.to_numeric(d.dilution_factor.where(~standards, "1"), errors="raise")
    if (not np.isfinite(d.concentration).all() or (d.concentration < 0).any()
            or not np.isfinite(d.dilution_factor).all() or (d.dilution_factor <= 0).any()):
        raise ValueError("Invalid concentration or dilution factor")
    d["concentration_canonical"] = d.concentration * d.concentration_unit.map(DOSE_UNITS)
    if (d.loc[d.role=="qc","concentration"]<=0).any():
        raise ValueError("QC nominal concentration must be positive (before dilution)")
    d["canonical_unit"] = d.concentration_unit.map(lambda x: "mg/mL" if x in MASS_UNITS else "source_unit" if x == "source_unit" else "M")
    for plate, g in d.groupby("plate_id"):
        if g.concentration_unit.nunique() != 1 or g.response_unit.nunique() != 1:
            raise ValueError(f"Plate {plate} mixes units")
        if not (g.role == "standard").any() or not (g.role == "unknown").any():
            raise ValueError(f"Plate {plate} needs standards and unknowns")
        if not ((g.role == "standard") & ~g.exclude & (g.concentration > 0)).any():
            raise ValueError(f"Plate {plate} needs positive included standards")
    return d


def inverse(y, fit):
    """Concentration (canonical unit) on the fitted curve at response y, or None."""
    if fit["bottom"] is None:
        return None
    return inverse_parameters(y, fit["bottom"], fit["top"], fit.get("log10_c_canonical", fit["log10_half_response_canonical"]),
                              fit["hill_slope_signed"], fit.get("asymmetry") or 1.)


def standards_qc(std, fit, cfg, factor):
    """Back-calculate each positive standard well; judge each level; derive LLOQ/ULOQ."""
    q = cfg["qc"]
    levels = []
    for nominal, g in std[std.concentration_canonical > 0].groupby("concentration_canonical", sort=True):
        back = [inverse(float(y), fit) for y in g.response]
        complete = all(b is not None for b in back)
        mean = float(np.mean(back)) if complete else None
        levels.append({"nominal_input_unit": nominal / factor, "n_wells": len(back),
                       "back_calculated_input_unit": [None if b is None else b / factor for b in back],
                       "mean_back_calculated_input_unit": None if mean is None else mean / factor,
                       "recovery_percent": None if mean is None else 100 * mean / nominal,
                       "cv_percent": float(100 * np.std(back, ddof=1) / mean) if complete and len(back) > 1 else None})

    def passes(level, limits):
        r, cv = level["recovery_percent"], level["cv_percent"]
        return (r is not None and limits[0] <= r <= limits[1]
                and (q["max_cv_percent"] is None or cv is None or cv <= q["max_cv_percent"]))

    edge_passing = [i for i, level in enumerate(levels) if passes(level, q["edge_recovery_limits_percent"])]
    lo, hi = (edge_passing[0], edge_passing[-1]) if edge_passing else (None, None)
    for i, level in enumerate(levels):
        if lo is None or not lo <= i <= hi:
            level.update(role="outside_quantification_range", limits_percent=q["edge_recovery_limits_percent"], passed=False)
        else:
            limits = q["edge_recovery_limits_percent"] if i in (lo, hi) else q["recovery_limits_percent"]
            level.update(role="LLOQ" if i == lo else "ULOQ" if i == hi else "interior",
                         limits_percent=limits, passed=passes(level, limits))
    n_passing = sum(level["passed"] for level in levels)
    required = max(q["min_passing_levels"], math.ceil(q["min_passing_fraction"] * len(levels)))
    accepted = bool(lo is not None and hi > lo and n_passing >= required)
    return {"levels": levels, "n_levels": len(levels), "n_passing": n_passing, "n_required": required,
            "accepted": accepted,
            "lloq_input_unit": levels[lo]["nominal_input_unit"] if accepted else None,
            "uloq_input_unit": levels[hi]["nominal_input_unit"] if accepted else None,
            "criteria": q}


def quantify(d, cfg):
    fits, wells, summaries = [], [], []
    dc = dose_config(cfg)
    for plate, g in d.groupby("plate_id", sort=False):
        std = g[(g.role == "standard") & ~g.exclude].copy()
        std["curve_id"] = plate
        std["sample_id"] = cfg["assay"]["analyte"]
        std["experiment_id"] = plate
        std["observation_id"] = std.well_id
        factor = DOSE_UNITS[str(g.concentration_unit.iloc[0])]
        model = cfg["fit"]["model"]
        if model == "5pl":
            fit = fit_5pl(std, cfg, dc, factor)
        else:
            # ELISA has its own calibration/QC contract; 0.13.4 changes functional
            # dose-response adequacy only, preserving this validated method.
            fit, _, _ = fit_dose_curve(std, dc, shape_diagnostics=False)
            fit.update(calibration_model="4pl", asymmetry=1., log10_c_canonical=fit["log10_half_response_canonical"])
        fit["endpoint"] = "calibration_C50"
        fit["plate_id"] = plate
        fit["calibration_covariance"] = None
        fit["dilution_linearity"] = []
        if fit["bottom"] is not None and cfg["uncertainty"]["unknown_interval"] == "delta_method":
            fit["calibration_covariance"] = covariance(fit, std.concentration_canonical.to_numpy(float),
                                                       std.response.to_numpy(float), dc, model, cfg["direction"])
        qc = standards_qc(std, fit, cfg, factor)
        fit["qc"] = qc
        from .elisa_validation import independent_controls
        fit["independent_qc"]=independent_controls(g,fit,cfg,factor)
        fit["quantification_range_evidence"]="current_plate_standard_back_calculation_only"
        fit["quantification_reportable"]=bool(fit["reportable"] and qc["accepted"] and
            (fit["independent_qc"]["status"]=="passed" or
             (fit["independent_qc"]["status"]=="not_assessed" and not cfg["independent_qc"]["required"])))
        fits.append(fit)
        # Interpolation is limited to the measured standard responses and the QC-derived LLOQ-ULOQ.
        lo_y, hi_y = float(std.response.min()), float(std.response.max())
        lo_c = qc["lloq_input_unit"] * factor if qc["accepted"] else float("nan")
        hi_c = qc["uloq_input_unit"] * factor if qc["accepted"] else float("nan")
        for _, row in g[g.role == "unknown"].iterrows():
            reasons = []
            if row.exclude:
                reasons.append("excluded")
            if not fit["reportable"]:
                reasons.append("standard_fit_not_reportable")
            if not qc["accepted"]:
                reasons.append("calibration_qc_failed")
            if fit["independent_qc"]["status"] not in ("passed","not_assessed"):
                reasons.append("independent_qc_not_passed")
            if cfg["independent_qc"]["required"] and fit["independent_qc"]["status"]!="passed":
                reasons.append("required_independent_qc_missing_or_failed")
            y = float(row.response)
            if not lo_y <= y <= hi_y:
                reasons.append("response_outside_observed_standard_range")
            estimate = None
            diluted = None
            if not reasons:
                diluted = inverse(y, fit)
                if diluted is None:
                    reasons.append("response_outside_fitted_plateaus")
                elif diluted < lo_c:
                    reasons.append("concentration_below_lloq")
                elif diluted > hi_c:
                    reasons.append("concentration_above_uloq")
                else:
                    estimate = diluted / DOSE_UNITS[row.concentration_unit] * float(row.dilution_factor)
            interval = None
            if estimate is not None and cfg["uncertainty"]["unknown_interval"] == "delta_method":
                interval = inverse_interval([y], [float(row.dilution_factor) / DOSE_UNITS[row.concentration_unit]],
                                            fit, dc, model, cfg["direction"], cfg["uncertainty"]["level"])
            wells.append({"plate_id": plate, "well_id": row.well_id, "sample_id": row.sample_id,
                          "response": y, "dilution_factor": float(row.dilution_factor),
                          "concentration_unit": row.concentration_unit,
                          "measured_concentration_in_well": None if diluted is None else diluted / DOSE_UNITS[row.concentration_unit],
                          "estimated_sample_concentration": estimate,
                          "status": "quantified" if not reasons else "withheld",
                          "ci_low": None if interval is None else interval[1],
                          "ci_high": None if interval is None else interval[2],
                          "ci_status": _ci_status(estimate, interval, cfg),
                          "independent_qc_status":fit["independent_qc"]["status"],
                          "interpretation_scope":"exploratory_interpolation_not_assay_validation", "diagnostics": reasons})
        for sample, sub in g[g.role == "unknown"].groupby("sample_id", sort=False):
            found = [w for w in wells if w["plate_id"] == plate and w["sample_id"] == sample]
            linearity, per_dilution = dilution_linearity(found, cfg)
            fit["dilution_linearity"].extend({"sample_id": sample, **row} for row in per_dilution)
            rule = cfg["dilution_linearity"]
            in_range = linearity.pop("_in_range")
            if rule["summary_policy"] == "all_wells_required":
                chosen = found
                complete = bool(found) and all(w["status"] == "quantified" for w in found)
            else:
                chosen = [w for w in found if w["dilution_factor"] in in_range]
                complete = linearity["linearity_status"] == "passed"
            status = "complete" if complete else "incomplete"
            if complete and rule["required"] and linearity["linearity_status"] != "passed":
                status, complete = "withheld_dilution_linearity_not_passed", False
            values = [w["estimated_sample_concentration"] for w in chosen if w["status"] == "quantified"]
            mean = float(np.mean(values)) if complete else None
            sd = float(np.std(values, ddof=1)) if complete and len(values) > 1 else None
            interval = None
            if complete and cfg["uncertainty"]["unknown_interval"] == "delta_method":
                interval = inverse_interval([w["response"] for w in chosen],
                                            [w["dilution_factor"] / DOSE_UNITS[w["concentration_unit"]] for w in chosen],
                                            fit, dc, model, cfg["direction"], cfg["uncertainty"]["level"])
            summaries.append({"plate_id": plate, "sample_id": sample, "n_wells": len(found),
                              "n_quantified": sum(w["status"] == "quantified" for w in found), "unit": sub.concentration_unit.iloc[0],
                              "mean_concentration": mean, "sd_concentration": sd,
                              "cv_percent": None if sd is None else 100 * sd / mean,
                              "status": status, "n_wells_in_estimate": len(values) if complete else 0,
                              "ci_low": None if interval is None else interval[1],
                              "ci_high": None if interval is None else interval[2],
                              "ci_level": cfg["uncertainty"]["level"], "ci_status": _ci_status(mean, interval, cfg),
                              **linearity})
    return fits, wells, summaries


def _ci_status(estimate, interval, cfg):
    if cfg["uncertainty"]["unknown_interval"] == "none":
        return "not_requested"
    if estimate is None:
        return "not_computed_no_reportable_estimate"
    return "calibration_conditional_delta_method" if interval is not None else "not_estimable"


def dilution_linearity(found, cfg):
    """Agreement of dilution-corrected results across a sample's dilutions on one plate.

    A dilution is in range only when every one of its wells is quantified. Each in-range
    dilution's mean is compared with the equal-weight mean across in-range dilutions.
    """
    rule = cfg["dilution_linearity"]
    by_dilution = {}
    for w in found:
        by_dilution.setdefault(w["dilution_factor"], []).append(w)
    rows, in_range = [], []
    for dilution in sorted(by_dilution):
        group = by_dilution[dilution]
        ok = all(w["status"] == "quantified" for w in group)
        mean = float(np.mean([w["estimated_sample_concentration"] for w in group])) if ok else None
        rows.append({"dilution_factor": dilution, "n_wells": len(group), "in_range": ok, "mean_corrected_concentration": mean})
        if ok:
            in_range.append(dilution)
    summary = {"n_dilutions": len(by_dilution), "n_dilutions_in_range": len(in_range), "linearity_status": "not_assessed",
               "dilution_cv_percent": None, "dilution_recovery_min_percent": None, "dilution_recovery_max_percent": None,
               "dilution_trend_log_slope": None, "dilution_trend_p": None, "linearity_diagnostics": "", "_in_range": in_range}
    if len(in_range) < 2:
        if len(by_dilution) >= 2:
            summary["linearity_diagnostics"] = "fewer_than_two_dilutions_in_range"
        return summary, rows
    means = np.array([r["mean_corrected_concentration"] for r in rows if r["in_range"]])
    centre = float(np.mean(means))
    for r in rows:
        r["recovery_vs_dilution_mean_percent"] = None if not r["in_range"] else 100 * r["mean_corrected_concentration"] / centre
    recoveries = 100 * means / centre
    cv = float(100 * np.std(means, ddof=1) / centre)
    low, high = rule["recovery_limits_percent"]
    passed = bool(np.all((recoveries >= low) & (recoveries <= high)) and cv <= rule["max_cv_percent"])
    diagnostics = []
    slope = p_value = None
    increasing = "corrected_concentration_increases_with_dilution_possible_hook_or_matrix_interference"
    decreasing = "corrected_concentration_decreases_with_dilution"
    if len(in_range) >= 3:
        fit = linregress(np.log10(in_range), np.log10(means))
        slope, p_value = float(fit.slope), float(fit.pvalue)
        if p_value < .05:
            diagnostics.append(increasing if slope > 0 else decreasing)
    # A failed series whose means move monotonically with dilution is flagged by pattern, not by a test.
    if not passed and not diagnostics and np.all(np.diff(means) > 0):
        diagnostics.append(increasing)
    elif not passed and not diagnostics and np.all(np.diff(means) < 0):
        diagnostics.append(decreasing)
    summary.update(linearity_status="passed" if passed else "failed", dilution_cv_percent=cv,
                   dilution_recovery_min_percent=float(recoveries.min()), dilution_recovery_max_percent=float(recoveries.max()),
                   dilution_trend_log_slope=slope, dilution_trend_p=p_value, linearity_diagnostics=";".join(diagnostics))
    return summary, rows
