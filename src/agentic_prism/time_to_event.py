"""Time-to-event workflow: configuration, data contract and analysis (numerics in survival.py)."""
from copy import deepcopy
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import norm
from . import survival as sv

DEFAULTS = {
    "analysis_type": "time_to_event", "schema_version": 1, "input": None,
    "source": "User-supplied time-to-event records",
    "study": {"endpoint": None, "time_origin": None, "time_unit": None, "censoring_rationale": None, "rationale": None},
    "comparison": {"arms": None, "control_arm": None, "post_hoc": "arm_vs_control", "logrank_inference": "asymptotic",
                   "permutation_draws": 9999, "seed": 20260930, "covariates": [], "landmarks": [],
                   "conf_type": "log-log", "confidence_level": .95},
    "report": {"plot_style": "prism_like", "show_confidence_bands": True}}
TIME_UNITS = ("hour", "day", "week", "month", "year")


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported time-to-event config")
    cfg = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(cfg[key], dict):
            if not isinstance(value, dict) or set(value) - set(cfg[key]):
                raise ValueError(f"Unsupported {key} settings")
            cfg[key].update(value)
        else:
            cfg[key] = value
    s, q = cfg["study"], cfg["comparison"]
    if cfg["analysis_type"] != "time_to_event" or cfg["schema_version"] != 1 or not isinstance(cfg["input"], str) or not cfg["input"]:
        raise ValueError("Unsupported time-to-event schema or missing input")
    for key in ("endpoint", "time_origin", "censoring_rationale", "rationale"):
        if not isinstance(s[key], str) or not s[key].strip():
            raise ValueError(f"Declare study.{key}: the event definition, time zero, why censored subjects are censored, and the design rationale")
    if s["time_unit"] not in TIME_UNITS:
        raise ValueError(f"study.time_unit must be one of {TIME_UNITS}")
    arms = q["arms"]
    if not isinstance(arms, list) or len(arms) < 2 or any(not isinstance(a, str) or not a.strip() for a in arms) or len(set(arms)) != len(arms):
        raise ValueError("Declare at least two distinct arms")
    if q["post_hoc"] not in ("arm_vs_control", "all_pairs", "none"):
        raise ValueError("post_hoc must be arm_vs_control, all_pairs or none")
    if q["post_hoc"] == "arm_vs_control" and q["control_arm"] not in arms:
        raise ValueError("Declare control_arm as one of the arms")
    if q["post_hoc"] != "arm_vs_control" and q["control_arm"] is not None and q["control_arm"] not in arms:
        raise ValueError("control_arm must be one of the arms")
    if q["logrank_inference"] not in ("asymptotic", "permutation"):
        raise ValueError("logrank_inference must be asymptotic or permutation")
    if type(q["permutation_draws"]) is not int or not 999 <= q["permutation_draws"] <= 100000:
        raise ValueError("permutation_draws must be an integer from 999 to 100000")
    if type(q["seed"]) is not int or not 0 <= q["seed"] < 2**32:
        raise ValueError("seed must be an unsigned 32-bit integer")
    if not isinstance(q["covariates"], list) or any(not isinstance(c, str) or not c.strip() for c in q["covariates"]) \
            or len(set(q["covariates"])) != len(q["covariates"]) or set(q["covariates"]) & {"subject_id", "arm", "time", "event"}:
        raise ValueError("covariates must be distinct numeric column names")
    if not isinstance(q["landmarks"], list) or any(isinstance(t, bool) or not isinstance(t, (int, float)) or not t > 0 for t in q["landmarks"]):
        raise ValueError("landmarks must be positive times")
    if q["conf_type"] not in ("log", "log-log"):
        raise ValueError("conf_type must be log (R survfit default) or log-log (SAS/Prism asymmetrical)")
    if isinstance(q["confidence_level"], bool) or not isinstance(q["confidence_level"], (int, float)) or not .5 < q["confidence_level"] < 1:
        raise ValueError("Invalid confidence_level")
    if cfg["report"]["plot_style"] not in ("standard", "prism_like") or type(cfg["report"]["show_confidence_bands"]) is not bool:
        raise ValueError("Unsupported report settings")
    return cfg


def load_data(path, cfg):
    q = cfg["comparison"]
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"subject_id", "arm", "time", "event"} | set(q["covariates"])
    if d.empty or required - set(d) or d.columns.duplicated().any():
        raise ValueError(f"Invalid time-to-event table; missing {sorted(required - set(d))}")
    for key in ("subject_id", "arm"):
        if d[key].str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
    if d.subject_id.duplicated().any():
        raise ValueError("One row per subject: duplicate subject_id (recurrent events are not supported)")
    d["time"] = pd.to_numeric(d.time, errors="raise")
    if not np.isfinite(d.time).all() or (d.time < 0).any():
        raise ValueError("Times must be finite and nonnegative")
    if not d.event.isin(("0", "1")).all():
        raise ValueError("event must be 1 (event observed) or 0 (censored)")
    d["event"] = d.event.astype(int)
    for c in q["covariates"]:
        d[c] = pd.to_numeric(d[c], errors="raise")
        if not np.isfinite(d[c]).all():
            raise ValueError(f"Covariate {c} must be finite for every subject (no silent case deletion)")
    d["exclude"] = d.get("exclude", pd.Series("false", index=d.index)).str.lower()
    d["exclusion_reason"] = d.get("exclusion_reason", pd.Series("", index=d.index))
    if not d.exclude.isin(("true", "false")).all() or ((d.exclude == "true") & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Exclusions require true/false and a reason")
    d["exclude"] = d.exclude.eq("true")
    if set(d.arm) != set(q["arms"]):
        raise ValueError("Input arms do not match declared arms")
    used = d[~d.exclude]
    per_arm = used.groupby("arm").size().reindex(q["arms"]).fillna(0)
    if per_arm.min() < 2:
        raise ValueError("Every arm needs at least two included subjects")
    return d


def _holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    return adjusted


def _pairs(q):
    arms = q["arms"]
    if q["post_hoc"] == "none":
        return []
    if q["post_hoc"] == "arm_vs_control":
        c = arms.index(q["control_arm"])
        return [(c, i) for i in range(len(arms)) if i != c]
    return list(combinations(range(len(arms)), 2))


def compare(d, cfg):
    q, study = cfg["comparison"], cfg["study"]
    used = d[~d.exclude].reset_index(drop=True)
    arms = q["arms"]
    ai = used.arm.map({a: i for i, a in enumerate(arms)}).to_numpy()
    time, event = used.time.to_numpy(float), used.event.to_numpy(int)
    level, diagnostics = q["confidence_level"], []
    curves, summaries = [], []
    for i, arm in enumerate(arms):
        mask = ai == i
        km = sv.kaplan_meier(time[mask], event[mask], level, q["conf_type"])
        med = sv.median_survival(km)
        for j in range(len(km["time"])):
            curves.append({"arm": arm, "time": float(km["time"][j]), "n_risk": int(km["n_risk"][j]), "n_event": int(km["n_event"][j]),
                           "n_censor": int(km["n_censor"][j]), "survival": float(km["surv"][j]),
                           "std_err_cumhaz": float(km["std_err"][j]) if np.isfinite(km["std_err"][j]) else None,
                           "lower": None if np.isnan(km["lower"][j]) else float(km["lower"][j]),
                           "upper": None if np.isnan(km["upper"][j]) else float(km["upper"][j])})
        row = {"arm": arm, "n": int(mask.sum()), "events": int(event[mask].sum()), "censored": int((1 - event[mask]).sum()),
               "median": med[0], "median_lower": med[1], "median_upper": med[2],
               "median_status": "reached" if med[0] is not None else "not_reached",
               "max_follow_up": float(time[mask].max())}
        for t in q["landmarks"]:
            if t > time[mask].max() and km["surv"][-1] > 0:
                row[f"survival_at_{t:g}"] = None
                diagnostics.append(f"landmark_{t:g}_beyond_follow_up_in_{arm}")
            elif t > time[mask].max():
                row.update({f"survival_at_{t:g}": 0., f"survival_at_{t:g}_lower": None, f"survival_at_{t:g}_upper": None})
            else:
                s, lo, hi = sv.survival_at(km, t)
                row.update({f"survival_at_{t:g}": s, f"survival_at_{t:g}_lower": lo, f"survival_at_{t:g}_upper": hi})
        summaries.append(row)
        if row["events"] < 5:
            diagnostics.append(f"fewer_than_5_events_in_{arm}")
    overall = sv.logrank(time, event, ai, len(arms))
    if not overall["computable"]:
        diagnostics.append("logrank_reduced_or_undefined_arms_without_overlapping_risk_sets")
    overall.update(test="log-rank (all arms)", inference="asymptotic chi-square")
    if q["logrank_inference"] == "permutation":
        overall.update(p_value_asymptotic=overall["p_value"], inference=f"permutation ({q['permutation_draws']} draws, seed {q['seed']})",
                       p_value=sv.permutation_logrank(time, event, ai, len(arms), q["permutation_draws"], q["seed"]),
                       p_resolution=1 / (q["permutation_draws"] + 1))
    elif min(s["n"] for s in summaries) < 10:
        diagnostics.append("small_arms_asymptotic_logrank_consider_permutation")
    # Cox model: arm indicators against the reference (control arm, else the first declared arm) plus covariates.
    ref = arms.index(q["control_arm"]) if q["control_arm"] in arms else 0
    others = [i for i in range(len(arms)) if i != ref]
    x = np.column_stack([(ai == i).astype(float) for i in others] + [used[c].to_numpy(float) for c in q["covariates"]])
    names = [f"arm:{arms[i]}" for i in others] + q["covariates"]
    cox_rows, zph_rows, cox_summary = [], [], {"reportable": False}
    no_event_arms = [s["arm"] for s in summaries if s["events"] == 0]
    if no_event_arms:
        diagnostics.append("arm_without_events_cox_hazard_ratio_not_estimable")
    else:
        fit = sv.cox(time, event, x)
        estimable = fit["converged"] and np.all(np.isfinite(fit["covariance"])) and np.all(np.abs(fit["beta"]) < 20)
        if not estimable:
            diagnostics.append("cox_model_not_converged_or_monotone_likelihood")
        else:
            beta, cov = fit["beta"], fit["covariance"]
            cox_summary = {"reportable": True, "reference_arm": arms[ref], "likelihood_ratio": fit["lr"], "wald": fit["wald"],
                           "score": fit["score"], "df": fit["df"], "loglik_null": fit["loglik"][0], "loglik": fit["loglik"][1],
                           "events_per_parameter": float(event.sum() / x.shape[1])}
            if event.sum() / x.shape[1] < 10:
                diagnostics.append("fewer_than_10_events_per_cox_parameter")
            z1 = norm.ppf((1 + level) / 2)
            for name, b, se in zip(names, beta, np.sqrt(np.diag(cov))):
                cox_rows.append({"term": name, "log_hazard_ratio": float(b), "standard_error": float(se), "hazard_ratio": float(np.exp(b)),
                                 "hr_lower": float(np.exp(b - z1 * se)), "hr_upper": float(np.exp(b + z1 * se)),
                                 "z": float(b / se), "p_value": float(2 * norm.sf(abs(b / se))), "adjusted_for": ";".join(q["covariates"])})
            terms = [list(range(len(others)))] + [[len(others) + j] for j in range(len(q["covariates"]))]
            for name, t in zip(["arm"] + q["covariates"] + ["GLOBAL"], sv.cox_zph(time, event, fit, terms)):
                zph_rows.append({"term": name, **t})
            if zph_rows[-1]["p_value"] < .05:
                diagnostics.append("proportional_hazards_questionable_hazard_ratio_is_a_time_average")
            cox_summary["_beta"], cox_summary["_cov"] = beta, cov
    pairs = _pairs(q)
    contrast_rows = []
    zf = norm.ppf(1 - (1 - level) / (2 * max(len(pairs), 1)))
    for a, b in pairs:
        mask = (ai == a) | (ai == b)
        lr = sv.logrank(time[mask], event[mask], (ai[mask] == b).astype(int), 2)
        p = lr["p_value"]
        if q["logrank_inference"] == "permutation":
            p = sv.permutation_logrank(time[mask], event[mask], (ai[mask] == b).astype(int), 2, q["permutation_draws"], q["seed"])
        row = {"arm_a": arms[a], "arm_b": arms[b], "logrank_chisq": lr["chisq"], "p_unadjusted": p,
               "hazard_ratio_b_vs_a": None, "hr_lower": None, "hr_upper": None}
        if cox_summary["reportable"]:
            v = np.zeros(len(names))
            if b != ref:
                v[others.index(b)] += 1
            if a != ref:
                v[others.index(a)] -= 1
            est, se = float(v @ cox_summary["_beta"]), float(np.sqrt(v @ cox_summary["_cov"] @ v))
            row.update(hazard_ratio_b_vs_a=float(np.exp(est)), hr_lower=float(np.exp(est - zf * se)), hr_upper=float(np.exp(est + zf * se)))
        contrast_rows.append(row)
    for r, p in zip(contrast_rows, _holm([r["p_unadjusted"] for r in contrast_rows]) if contrast_rows else []):
        r["p_adjusted"] = float(p)
    cox_summary.pop("_beta", None), cox_summary.pop("_cov", None)
    fit_row = {"analysis_type": "time_to_event", "endpoint": study["endpoint"], "time_unit": study["time_unit"],
               "n_subjects": int(len(used)), "n_events": int(event.sum()), "n_arms": len(arms), "n_groups": len(arms),
               "logrank_chisq": overall["chisq"], "logrank_df": overall["df"], "p_value": overall["p_value"],
               "logrank_inference": overall["inference"], "cox_reportable": cox_summary["reportable"],
               "status": "estimated", "reportable": True, "diagnostics": list(dict.fromkeys(diagnostics)),
               "contrast_method": "pairwise log-rank (Holm p); Cox hazard ratios with Bonferroni simultaneous Wald limits"}
    return {"schema_version": 1, "analysis_type": "time_to_event", "fits": [fit_row], "arm_summaries": summaries,
            "km_curves": curves, "logrank": overall, "cox": cox_summary, "cox_terms": cox_rows, "ph_test": zph_rows,
            "contrasts": contrast_rows}
