"""Deterministic, source-linked interpretation facts. No fitting or inference here.

Only saved report-facing numbers are copied. Audit values stay in results.json.
Submodel gates take precedence over the top-level workflow's completion status.
"""
from copy import deepcopy
import json
from pathlib import Path

SUPPORTED = frozenset({"repeated_measures", "time_to_event", "tumor_growth", "nonparametric", "mmrm"})


def _ref(path):
    return "results.json#" + path


def _escape(key):
    return str(key).replace("~", "~0").replace("/", "~1")


def disclosures(value, path=""):
    """Inventory saved states, including new diagnostic codes, without interpreting them."""
    out = []
    if isinstance(value, dict):
        for key, child in value.items():
            here = path + "/" + _escape(key)
            if key in ("diagnostic", "diagnostics"):
                if isinstance(child, list):
                    out.extend({"kind": "diagnostic", "value": v, "source": _ref(here + f"/{i}")}
                               for i, v in enumerate(child) if v)
                elif child:
                    out.append({"kind": "diagnostic", "value": child, "source": _ref(here)})
            elif key == "status" or key.endswith(("_status", "reportable")) or key == "computable":
                out.append({"kind": "saved_state", "value": child, "source": _ref(here)})
            else:
                out.extend(disclosures(child, here))
    elif isinstance(value, list):
        for i, child in enumerate(value):
            out.extend(disclosures(child, path + f"/{i}"))
    return out


def _result(identifier, path, row, estimand, *, estimate=None, bounds=None, level=None,
            method=None, multiplicity="none", unit=None, tests=(), context=(),
            allowed=True, reasons=(), limited=False):
    reasons = list(dict.fromkeys(r for r in reasons if r))
    point = row.get(estimate) if estimate else None
    test = {key: row[key] for key in tests if key in row}
    usable = point is not None if estimate else any(row.get(k) is not None for k in tests if k.startswith("p_"))
    allowed = bool(allowed and usable)
    if not allowed and not reasons:
        reasons.append("estimate_or_test_not_available_in_saved_result")
    interval = {"status": "not_applicable", "level": None, "lower": None, "upper": None,
                "method": None, "multiplicity": "none", "sources": []}
    if bounds:
        lo, hi = (row.get(k) if allowed else None for k in bounds)
        interval.update(status="available" if lo is not None and hi is not None else
                        "partial" if lo is not None or hi is not None else "unavailable",
                        level=level, lower=lo, upper=hi, method=method, multiplicity=multiplicity,
                        sources=[_ref(path + "/" + k) for k in bounds if k in row])
        if allowed and interval["status"] != "available":
            reasons.append("interval_not_fully_available_in_saved_result")
            limited = True
    return {"id": identifier, "estimand": estimand, "context": {k: row[k] for k in context if k in row},
            "reportability": "withheld" if not allowed else "limited" if limited or reasons else "reportable",
            "reasons": list(dict.fromkeys(reasons)), "estimate": point if allowed else None, "unit": unit,
            "estimate_source": _ref(path + "/" + estimate) if estimate and estimate in row else None,
            "interval": interval, "test": test if allowed else {},
            "test_sources": {k: _ref(path + "/" + k) for k in test}, "source": _ref(path)}


def _repeated(result, cfg):
    q, fit = cfg["comparison"], result["fits"][0]
    two = q["design"].startswith("two_way")
    anova = q["design"].endswith("rm_anova")
    diagnostics = fit["diagnostics"]
    model = ("Split-plot ANOVA; GG-corrected within-unit tests" if two else "One-factor repeated-measures ANOVA; GG correction") if anova else (
        "Gaussian random-intercept REML; arm x categorical condition cell means" if two else
        "Gaussian random-intercept REML; categorical condition means")
    limitations = ["Independent units, not wells or measurements, define replication.",
                  "Conditions are categorical; these results do not estimate a continuous-time growth rate.",
                  "Non-significance does not establish equivalence or parallel profiles."]
    if not anova:
        limitations += ["The covariance model assumes a Gaussian unit intercept and equal, independent Gaussian residuals.",
                        "Available-case inference requires MAR conditional on the declared model; observed data cannot verify MAR."]
        if q["inference"] == "parametric_bootstrap":
            limitations.append("Bootstrap calibration is design-specific: consult the 0.8.0 B=199 and 0.8.2 B=999 records, including misses. Neither validates the current experiment.")
    else:
        limitations.append("The corrected omnibus test and the predeclared contrast family answer different hypotheses.")
    facts, primary = [], []
    rows = result.get("tests", []) if two else [fit]
    for i, row in enumerate(rows):
        effect = row.get("effect", "condition")
        identifier = "test:" + effect
        facts.append(_result(identifier, f"/tests/{i}" if two else "/fits/0", row,
            f"{effect} effect on {q['outcome']}", tests=("f_statistic", "wald_statistic", "df_numerator", "df_denominator", "p_value", "epsilon_gg"),
            allowed=fit["reportable"], reasons=diagnostics))
        if not two or effect == "arm_x_condition":
            primary.append(identifier)
    if not rows:
        facts.append(_result("test:withheld", "/fits/0", fit, "Declared repeated-measures tests",
                             allowed=False, reasons=diagnostics))
        primary.append("test:withheld")
    bootstrap = not anova and q["inference"] == "parametric_bootstrap"
    for i, row in enumerate(result["contrasts"]):
        reasons = diagnostics + ([row["diagnostic"]] if row.get("diagnostic") else [])
        if not row["reportable"] and not row.get("diagnostic"):
            reasons = reasons + ["contrast_inference_withheld_by_model"]
        identifier = f"contrast:{i}"
        facts.append(_result(identifier, f"/contrasts/{i}", row, f"Mean {q['outcome']} difference B minus A",
            estimate="difference_b_minus_a", bounds=("ci_low", "ci_high"), unit=q["unit"], level=q["confidence_level"],
            method=fit.get("contrast_method"), multiplicity="bootstrap max-abs-t family" if bootstrap else "Bonferroni family intervals; Holm p values",
            tests=("statistic", "df", "p_adjusted"), context=("arm", "condition", "group_a", "group_b"),
            allowed=fit["reportable"] and row["reportable"], reasons=reasons))
        primary.append(identifier)
    return model, "Condition mean differences; arm x condition interaction for two-way designs", facts, primary, limitations


def _survival(result, cfg):
    q, study, fit = cfg["comparison"], cfg["study"], result["fits"][0]
    diagnostics = fit["diagnostics"]
    facts, primary = [], []
    for i, row in enumerate(result["arm_summaries"]):
        path = f"/arm_summaries/{i}"
        identifier = f"median:{i}"
        reasons = ["median_not_reached"] if row["median_status"] == "not_reached" else []
        facts.append(_result(identifier, path, row, "Median time to declared event", estimate="median",
            bounds=("median_lower", "median_upper"), level=q["confidence_level"], unit=study["time_unit"],
            method=f"Inversion of KM {q['conf_type']} confidence limits", context=("arm", "n", "events", "censored", "max_follow_up"),
            reasons=reasons))
        # When the point median is not reached its lower confidence limit may still be estimable.
        # Keep it as a separate limited result; never use it as a median point estimate.
        if row["median"] is None and row["median_lower"] is not None:
            facts.append(_result(f"median_lower_limit:{i}", path, row, "Lower confidence limit for median time; point median not reached",
                estimate="median_lower", unit=study["time_unit"], context=("arm",), reasons=["median_not_reached"], limited=True))
        primary.append(identifier)
        for t in q["landmarks"]:
            key = f"survival_at_{t:g}"
            identifier = f"landmark:{i}:{t:g}"
            reasons = ["landmark_beyond_arm_follow_up"] if t > row["max_follow_up"] else []
            if row.get(key) == 1 and row.get(key + "_lower") == row.get(key + "_upper") == 1:
                reasons.append("no_events_before_landmark_greenwood_interval_degenerate_not_certainty")
            facts.append(_result(identifier, path, row, f"Survival probability at {t:g} {study['time_unit']}",
                estimate=key, bounds=(key + "_lower", key + "_upper"), level=q["confidence_level"], unit="probability",
                method=f"Kaplan-Meier Greenwood {q['conf_type']}", multiplicity="pointwise, not simultaneous across times or arms",
                context=("arm", "n", "events", "max_follow_up"),
                reasons=reasons))
            primary.append(identifier)
    lr = result["logrank"]
    reasons = diagnostics + (["logrank_has_zero_informative_degrees_of_freedom"] if lr["df"] == 0 else
                             ["logrank_reduced_rank_does_not_test_all_declared_arms"] if not lr["computable"] else [])
    facts.append(_result("logrank", "/logrank", lr, "Unadjusted equality of survival curves over observed follow-up",
        tests=("chisq", "df", "p_value", "p_resolution"), context=("inference",), allowed=lr["df"] > 0, reasons=reasons))
    primary.append("logrank")
    cox_ok = result["cox"]["reportable"]
    for i, row in enumerate(result["contrasts"]):
        facts.append(_result(f"pairwise_logrank:{i}", f"/contrasts/{i}", row, "Pairwise unadjusted survival-curve comparison; Holm family p",
            tests=("logrank_chisq", "p_adjusted"), context=("arm_a", "arm_b"), allowed=lr["computable"],
            reasons=diagnostics + ([] if lr["computable"] else ["pairwise_risk_set_rank_not_saved_withhold_until_review"])))
        facts.append(_result(f"pairwise_hr:{i}", f"/contrasts/{i}", row, "Cox hazard ratio B versus A, conditional on declared covariates",
            estimate="hazard_ratio_b_vs_a", bounds=("hr_lower", "hr_upper"), level=q["confidence_level"], unit="ratio",
            method="Cox Efron; Wald on log HR", multiplicity="Bonferroni simultaneous family intervals",
            context=("arm_a", "arm_b"), allowed=cox_ok, reasons=diagnostics))
    if not cox_ok:
        facts.append(_result("cox:withheld", "/cox", result["cox"], "Cox hazard ratios", allowed=False, reasons=diagnostics))
    for i, row in enumerate(result["cox_terms"]):
        facts.append(_result(f"cox_term:{i}", f"/cox_terms/{i}", row, "Conditional Cox hazard ratio per declared term",
            estimate="hazard_ratio", bounds=("hr_lower", "hr_upper"), level=q["confidence_level"], unit="ratio",
            method="Cox Efron; Wald on log HR", multiplicity="individual term interval, no family adjustment",
            tests=("z", "p_value"), context=("term", "adjusted_for"), allowed=cox_ok, reasons=diagnostics))
    for i, row in enumerate(result["ph_test"]):
        facts.append(_result(f"ph_test:{i}", f"/ph_test/{i}", row, "Proportional-hazards score diagnostic, KM time transform",
            tests=("chisq", "df", "p_value"), context=("term",), allowed=cox_ok))
    limitations = ["Endpoint, time zero and reasons for censoring are declarations requiring scientific review.",
                   "Independent, non-informative censoring is assumed and cannot be verified from these observations alone.",
                   "A hazard ratio requires the proportional-hazards assumption for a constant effect interpretation; inspect curves and PH diagnostics.",
                   "Log-rank is unadjusted even when the Cox model includes covariates; its p value is not a test of the adjusted HR.",
                   "No extrapolation, cure, equivalence, competing-risk or clustered-subject conclusion is supported."]
    return "Kaplan-Meier; log-rank; Cox proportional hazards with Efron ties", study["endpoint"], facts, primary, limitations


def _tumor(result, cfg):
    a, s, fit = cfg["analysis"], cfg["study"], result["fits"][0]
    diagnostics, level = fit["diagnostics"], a["confidence_level"]
    allowed = fit["model_reportable"]
    facts, primary = [], []
    if not allowed:
        facts.append(_result("growth_model:withheld", "/fits/0", fit, "Log-volume growth model inference", allowed=False, reasons=diagnostics))
        primary.append("growth_model:withheld")
    for i, row in enumerate(result["growth_tests"]):
        identifier = f"growth_test:{i}"
        facts.append(_result(identifier, f"/growth_tests/{i}", row, "Equality of arm-specific log growth rates",
            tests=("f_statistic", "df_numerator", "df_denominator", "p_value"), allowed=allowed, reasons=diagnostics))
        primary.append(identifier)
    for i, row in enumerate(result["growth_rates"]):
        for label, est, bounds, estimand, unit in (
            ("rate", "log_growth_rate_per_time", ("rate_lower", "rate_upper"), "Rate of change of log(V + offset)", "1/" + s["time_unit"]),
            ("doubling_time", "doubling_time", ("doubling_time_lower", "doubling_time_upper"), "Time for model V + offset to double", s["time_unit"])):
            reasons = diagnostics.copy()
            if label == "doubling_time" and row[est] is None:
                reasons.append("nonpositive_growth_rate_has_no_positive_doubling_time")
            facts.append(_result(f"{label}:{i}", f"/growth_rates/{i}", row, estimand, estimate=est, bounds=bounds,
                level=level, unit=unit, method="Satterthwaite t limits" if label == "rate" else "Transformed Satterthwaite rate limits when strictly positive",
                multiplicity="individual arm interval", context=("arm",), allowed=allowed, reasons=reasons))
    for i, row in enumerate(result["model_contrasts"]):
        for label, est, bounds, estimand, unit in (
            ("rate_difference", "log_rate_difference", ("rate_difference_lower", "rate_difference_upper"), "Log growth rate, arm minus control", "1/" + s["time_unit"]),
            ("model_tc", "model_tc_geometric_ratio_at_analysis_day", ("model_tc_lower", "model_tc_upper"),
             f"Geometric-mean (V + offset) arm/control ratio at {a['analysis_day']} {s['time_unit']}", "ratio")):
            identifier = f"{label}:{i}"
            facts.append(_result(identifier, f"/model_contrasts/{i}", row, estimand, estimate=est, bounds=bounds,
                level=level, unit=unit, method="Satterthwaite t; model T/C exponentiates log-ratio limits",
                multiplicity="Bonferroni family intervals; rate-difference p values use Holm",
                tests=("df", "p_adjusted") if label == "rate_difference" else (), context=("arm", "control_arm"), allowed=allowed, reasons=diagnostics))
            primary.append(identifier)
    for i, row in enumerate(result["observed_tgi"]):
        for label, estimand in (("tgi", "100 * (1 - mean change in arm / mean change in control)"),
                                ("tc", "100 * mean volume in arm / mean volume in control")):
            reasons = diagnostics.copy()
            if min(row["n_arm_on_day"], row["n_control_on_day"]) < 2:
                reasons.append("fewer_than_two_animals_measured_in_arm_or_control_on_analysis_day")
            elif row[label + "_percent"] is None:
                reasons.append("observed_ratio_not_defined_in_saved_result")
            facts.append(_result(f"observed_{label}:{i}", f"/observed_tgi/{i}", row,
                estimand + "; animals measured on the analysis day only", estimate=label + "_percent",
                bounds=(label + "_lower", label + "_upper"), level=level, unit="percent",
                method="Fieller with Welch-Satterthwaite df", multiplicity="Bonferroni family intervals",
                context=("arm", "control_arm", "n_arm_on_day", "n_control_on_day"), reasons=reasons))
    limitations = [f"Model is linear in log(V + {a['log_offset']}) over continuous {s['time_unit']}; results depend on the declared offset and fit range.",
                   "Random animal intercepts/slopes are Gaussian; residuals are independent with equal variance. Cage effects are not fitted.",
                   "Growth-model inference requires MAR conditional on the model; it does not correct MNAR removal.",
                   "Observed TGI and T/C describe measured animals only; growth-related removal can bias them toward survivors.",
                   "Model T/C concerns geometric means of V + offset; observed T/C concerns arithmetic means of V.",
                   "Missing Fieller limits do not supply a finite confidence interval; do not invent one or infer precision.",
                   "Tumor regression or TGI above 100% does not establish cure; time to humane endpoint is a separate analysis."]
    return "Log(V + offset) arm intercepts/slopes + correlated animal random intercept/slope; REML, Satterthwaite", "Growth-rate differences and model T/C; secondary observed TGI and T/C", facts, primary, limitations


def _nonparametric(result, cfg):
    q, fit = cfg["comparison"], result["fits"][0]
    shift = q["design"] in ("mann_whitney", "wilcoxon_signed_rank")
    estimand = ("Common location shift B minus A" if q["design"] == "mann_whitney" else
                "Pseudomedian of paired B minus A differences") if shift else "Differences in rank distributions"
    rows = [_result("primary", "/fits/0", fit, estimand,
        estimate="estimate" if shift else None, bounds=("ci_low", "ci_high") if shift else None,
        level=fit.get("achieved_confidence_level"), method=fit.get("interval_method"), unit=q["unit"] if shift else None,
        tests=("statistic", "df", "p_value"), allowed=fit["reportable"], reasons=fit["diagnostics"])]
    if shift:
        rows.append(_result("hodges_lehmann", "/fits/0", fit, estimand,
            estimate="hodges_lehmann", unit=q["unit"], allowed=fit["reportable"], reasons=fit["diagnostics"]))
    for i, row in enumerate(result["contrasts"]):
        rows.append(_result(f"dunn:{i}", f"/contrasts/{i}", row, "Pooled mean-rank difference B minus A",
            estimate="mean_rank_difference_b_minus_a", tests=("statistic", "p_adjusted"),
            context=("group_a", "group_b", "adjustment"), multiplicity=q["adjustment"], reasons=fit["diagnostics"]))
    return q["design"], estimand, rows, [r["id"] for r in rows], [
        "Rank significance does not establish biological importance, equivalence or a general difference in medians.",
        "HL requires a common location-shift model; signed ranks require symmetric paired differences.",
        "Exact inference is restricted to small samples without ties or zero differences; tied data use normal approximation.",
        "An unavailable or unbounded interval must not be presented as finite; use the saved confidence level. For normal approximation this is nominal, not empirical or exact attained coverage.",
        "Dunn adjusted p values do not supply location-shift confidence intervals; Friedman has no post-hoc family here."]


def _mmrm(result, cfg):
    q, fit = cfg["comparison"], result["fits"][0]
    rows = [_result("interaction", "/fits/0", fit, "Arm by categorical visit interaction",
        tests=("f_statistic", "df_numerator", "df_denominator", "p_value"), allowed=fit["reportable"], reasons=fit["diagnostics"])]
    for i, row in enumerate(result["contrasts"]):
        rows.append(_result(f"contrast:{i}", f"/contrasts/{i}", row, "Arm minus control mean difference at visit",
            estimate="estimate", bounds=("ci_low", "ci_high"), level=q["confidence_level"], unit=q["unit"],
            method=q["inference"], multiplicity="Bonferroni family intervals; Holm p values",
            tests=("statistic", "df", "p_adjusted"), context=("arm", "control_arm", "condition"),
            allowed=fit["reportable"], reasons=fit["diagnostics"]))
    return "Marginal REML " + q["covariance"], "Arm by visit interaction and visit-specific arm differences", rows, [r["id"] for r in rows], [
        "Covariance and contrast family must be predeclared; no data-driven model selection is performed.",
        "Available-case inference requires MAR conditional on this model; MNAR dropout is not corrected.",
        "AR(1) uses the declared visit order, not elapsed time; covariance is shared across arms.",
        "A non-significant interaction does not establish parallel profiles or equivalence.",
        "Kenward-Roger uses optional R mmrm; Satterthwaite uses Python. No Prism or SAS equivalence claim.",
        "0.8.2 calibration, 1000 Gaussian null datasets per row, 12 units per arm, 3 visits, 15% MCAR: interaction rejection US/Satterthwaite 7.0% and US/KR 7.1% exceeded the 6.3784% bound. US/KR family rejection 6.7% and coverage 93.3% also missed bounds. AR1 rows passed only this matching-covariance scenario. These results do not validate the current experiment."]


def build_facts(result, cfg):
    """Pure extraction from JSON-compatible saved values; does not mutate inputs."""
    kind = result["analysis_type"]
    if kind not in SUPPORTED or cfg["analysis_type"] != kind:
        raise ValueError("Interpretation facts are not implemented for this analysis type")
    model, estimand, rows, primary, limitations = {
        "repeated_measures": _repeated, "time_to_event": _survival, "tumor_growth": _tumor, "nonparametric": _nonparametric, "mmrm": _mmrm}[kind](result, cfg)
    return deepcopy({"schema_version": 1, "analysis_type": kind, "model": model, "estimand": estimand,
                     "declared_design": {k: v for k, v in cfg.items() if k in ("study", "comparison", "analysis")},
                     "sample_counts": {k: {"value": v, "source": _ref("/fits/0/" + k)}
                                       for k, v in result["fits"][0].items() if k.startswith("n_")},
                     "primary_result_ids": primary, "results": rows,
                     "required_disclosures": disclosures(result), "limitations": limitations})


def write_facts(run):
    """Write once, before manifest creation, from canonical saved JSON."""
    from .workflow import dump, sha
    run = Path(run)
    target = run / "interpretation_facts.json"
    if target.exists():
        raise FileExistsError(target)
    result = json.loads((run / "results.json").read_text())
    cfg = json.loads((run / "config.resolved.json").read_text())
    facts = build_facts(result, cfg)
    facts["source_sha256"] = {name: sha(run / name) for name in
                              ("results.json", "config.resolved.json", "diagnostics.json", "normalized_data.csv")}
    dump(target, facts)
