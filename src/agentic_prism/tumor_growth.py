"""Tumor growth: log-linear random-slope mixed model, model-based T/C, and observed TGI% / T/C% with Fieller limits.

Model: log(V + offset) = a_arm + b_arm * day + u0_animal + u1_animal * day + e, with
(u0, u1) ~ N(0, G) unstructured and e ~ N(0, sigma^2) independent; REML. Arms have
their own intercepts and slopes (cell-means coding). Satterthwaite t/F as in lmerTest.
Animals removed at a humane endpoint leave the data; the mixed model remains valid
if removal depends only on observed volumes (MAR), which TGI% at a fixed day does not.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
from .mixed_inference import RandomSlopeModel, fit_reml, satterthwaite_unconstrained

DEFAULTS = {
    "analysis_type": "tumor_growth", "schema_version": 1, "input": None, "source": "User-supplied tumor measurements",
    "study": {"volume_unit": None, "time_unit": None, "time_origin": None, "removal_rule": None,
              "dropout_rationale": None, "rationale": None},
    "analysis": {"arms": None, "control_arm": None, "baseline_day": None, "analysis_day": None, "log_offset": None,
                 "confidence_level": .95},
    "report": {"plot_style": "prism_like"}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported tumor-growth config")
    cfg = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(cfg[key], dict):
            if not isinstance(value, dict) or set(value) - set(cfg[key]):
                raise ValueError(f"Unsupported {key} settings")
            cfg[key].update(value)
        else:
            cfg[key] = value
    s, a = cfg["study"], cfg["analysis"]
    if cfg["analysis_type"] != "tumor_growth" or cfg["schema_version"] != 1 or not isinstance(cfg["input"], str) or not cfg["input"]:
        raise ValueError("Unsupported tumor-growth schema or missing input")
    if s["volume_unit"] not in ("mm3", "cm3") or s["time_unit"] not in ("day", "week"):
        raise ValueError("study.volume_unit must be mm3 or cm3 and study.time_unit day or week")
    for key in ("time_origin", "removal_rule", "dropout_rationale", "rationale"):
        if not isinstance(s[key], str) or not s[key].strip():
            raise ValueError(f"Declare study.{key}: time zero, the humane-endpoint/removal rule, why dropout may be treated as MAR, and the design")
    arms = a["arms"]
    if not isinstance(arms, list) or len(arms) < 2 or any(not isinstance(x, str) or not x.strip() for x in arms) or len(set(arms)) != len(arms):
        raise ValueError("Declare at least two distinct arms")
    if a["control_arm"] not in arms:
        raise ValueError("control_arm must be one of the arms")
    for key in ("baseline_day", "analysis_day"):
        if isinstance(a[key], bool) or not isinstance(a[key], (int, float)):
            raise ValueError(f"Declare analysis.{key} (a scheduled measurement day)")
    if not a["analysis_day"] > a["baseline_day"]:
        raise ValueError("analysis_day must follow baseline_day")
    if isinstance(a["log_offset"], bool) or not isinstance(a["log_offset"], (int, float)) or a["log_offset"] < 0:
        raise ValueError("Declare analysis.log_offset (>= 0, in volume units) for log(V + offset); zero volumes need a positive offset")
    if isinstance(a["confidence_level"], bool) or not isinstance(a["confidence_level"], (int, float)) or not .5 < a["confidence_level"] < 1:
        raise ValueError("Invalid confidence_level")
    if cfg["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    return cfg


def load_data(path, cfg):
    a = cfg["analysis"]
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"animal_id", "arm", "day", "volume"}
    if d.empty or required - set(d) or d.columns.duplicated().any():
        raise ValueError(f"Invalid tumor table; missing {sorted(required - set(d))}")
    for key in ("animal_id", "arm"):
        if d[key].str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
    d["day"] = pd.to_numeric(d.day, errors="raise")
    d["volume"] = pd.to_numeric(d.volume, errors="raise")
    if not (np.isfinite(d.day).all() and np.isfinite(d.volume).all()) or (d.volume < 0).any():
        raise ValueError("Days and volumes must be finite; volumes nonnegative")
    d["exclude"] = d.get("exclude", pd.Series("false", index=d.index)).str.lower()
    d["exclusion_reason"] = d.get("exclusion_reason", pd.Series("", index=d.index))
    if not d.exclude.isin(("true", "false")).all() or ((d.exclude == "true") & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Exclusions require true/false and a reason")
    d["exclude"] = d.exclude.eq("true")
    if set(d.arm) != set(a["arms"]) or d.groupby("animal_id").arm.nunique().gt(1).any():
        raise ValueError("Arms must match the declaration and each animal belongs to one arm")
    used = d[~d.exclude]
    if used.duplicated(["animal_id", "day"]).any():
        raise ValueError("One volume per animal and day; average caliper replicates upstream with provenance")
    if a["log_offset"] == 0 and (used.volume <= 0).any():
        raise ValueError("Zero volumes present: declare a positive log_offset")
    base = used[used.day == a["baseline_day"]]
    if set(base.animal_id) != set(used.animal_id):
        raise ValueError("Every included animal needs a baseline_day measurement")
    per_arm = used.groupby("arm").animal_id.nunique().reindex(a["arms"]).fillna(0)
    if per_arm.min() < 3:
        raise ValueError("At least three animals per arm are required")
    if used.groupby("animal_id").size().lt(2).any():
        raise ValueError("Every included animal needs at least two measurement days")
    return d


def fieller(mean_t, var_t, n_t, mean_c, var_c, n_c, level):
    """Fieller limits for mean_t / mean_c from independent samples (Welch-Satterthwaite df at the estimate)."""
    ratio = mean_t / mean_c
    a, b = var_t / n_t, var_c / n_c
    df = (a + ratio ** 2 * b) ** 2 / (a ** 2 / (n_t - 1) + (ratio ** 2 * b) ** 2 / (n_c - 1))
    tq = t_dist.ppf((1 + level) / 2, df)
    g = tq ** 2 * b / mean_c ** 2
    if g >= 1:
        return ratio, None, None, df  # denominator not distinguishable from zero: unbounded set
    # Roots of (m_t - R m_c)^2 = t^2 (a + R^2 b): [r +/- (t/|m_c|) sqrt(a(1-g) + r^2 b)] / (1 - g).
    half = tq / abs(mean_c) * np.sqrt(a * (1 - g) + ratio ** 2 * b)
    return ratio, float((ratio - half) / (1 - g)), float((ratio + half) / (1 - g)), float(df)


def compare(d, cfg):
    a, s = cfg["analysis"], cfg["study"]
    arms, level = a["arms"], a["confidence_level"]
    used = d[~d.exclude].sort_values(["animal_id", "day"]).reset_index(drop=True)
    diagnostics = []
    days = sorted(used.day.unique())
    summaries = []
    for arm in arms:
        for day in days:
            v = used.volume[(used.arm == arm) & (used.day == day)]
            if len(v):
                summaries.append({"arm": arm, "day": float(day), "n": int(len(v)), "mean": float(v.mean()),
                                  "sem": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None, "median": float(v.median())})
    last_day = used.groupby("animal_id").day.max()
    removed = last_day[last_day < a["analysis_day"]]
    dropout = [{"animal_id": aid, "arm": used.arm[used.animal_id == aid].iloc[0], "last_day": float(day)} for aid, day in removed.items()]
    if dropout:
        diagnostics.append("animals_removed_before_analysis_day_mixed_model_assumes_MAR")
    # Growth model.
    ai = used.arm.map({x: i for i, x in enumerate(arms)}).to_numpy()
    day = used.day.to_numpy(float)
    y = np.log(used.volume.to_numpy(float) + a["log_offset"])
    k = len(arms)
    x = np.zeros((len(used), 2 * k))
    x[np.arange(len(used)), ai] = 1
    x[np.arange(len(used)), k + ai] = day
    z = np.column_stack([np.ones(len(used)), day])
    ids = used.animal_id.to_numpy(str)
    model = RandomSlopeModel(y, x, ids, z)
    vy, vt = float(np.var(y)), float(np.var(day)) or 1.
    starts = [model.pack(np.diag([vy * f, vy * f / vt * .1]), vy * (1 - f) + 1e-6) for f in (.3, .6, .9)]
    model_rows, tests, growth = [], [], []
    fit_row = {"analysis_type": "tumor_growth", "n_animals": int(used.animal_id.nunique()), "n_observations": int(len(used)),
               "n_groups": k, "n_arms": k, "log_offset": a["log_offset"], "status": "estimated", "reportable": True,
               "p_value": None, "model_reportable": False}
    try:
        res = fit_reml(model, starts)
        g, sigma = model.unpack(res.x)
        beta, cov = model.gls(res.x)[:2]
        corr = g[0, 1] / np.sqrt(g[0, 0] * g[1, 1])
        boundary = g[1, 1] < 1e-8 * max(vy / vt, 1e-12) or abs(corr) > .999
        fit_row.update(random_intercept_variance=float(g[0, 0]), random_slope_variance=float(g[1, 1]),
                       random_effect_correlation=float(corr), residual_variance=float(sigma), reml_deviance=float(res.fun))
        if boundary:
            diagnostics.append("random_slope_variance_or_correlation_at_boundary_satterthwaite_withheld")
            raise ArithmeticError("boundary")
        dev, covf = model.deviance, (lambda u: model.gls(u)[1])
        slope_l = [np.eye(2 * k)[k + i] for i in range(k)]
        srows, _, _ = satterthwaite_unconstrained(dev, covf, res.x, beta, slope_l, level=level, n_family=1)
        for arm, sr in zip(arms, srows):
            b, lo, hi = sr["estimate"], sr["ci_low"], sr["ci_high"]
            growth.append({"arm": arm, "log_growth_rate_per_time": b, "standard_error": sr["standard_error"], "df": sr["df"],
                           "rate_lower": lo, "rate_upper": hi,
                           "doubling_time": float(np.log(2) / b) if b > 0 else None,
                           "doubling_time_lower": float(np.log(2) / hi) if lo > 0 else None,
                           "doubling_time_upper": float(np.log(2) / lo) if lo > 0 else None})
        equal_slopes = np.array([np.eye(2 * k)[k + i] - np.eye(2 * k)[k] for i in range(1, k)])
        _, om, _ = satterthwaite_unconstrained(dev, covf, res.x, beta, [], joint=equal_slopes)
        tests.append({"test": "equal_growth_rates_across_arms", **{kk: om[kk] for kk in ("f_statistic", "df_numerator", "df_denominator", "p_value")}})
        fit_row["p_value"] = om["p_value"]
        c = arms.index(a["control_arm"])
        others = [i for i in range(k) if i != c]
        slope_diff = [np.eye(2 * k)[k + i] - np.eye(2 * k)[k + c] for i in others]
        tc_log = [np.eye(2 * k)[i] - np.eye(2 * k)[c] + a["analysis_day"] * (np.eye(2 * k)[k + i] - np.eye(2 * k)[k + c]) for i in others]
        drows, _, _ = satterthwaite_unconstrained(dev, covf, res.x, beta, slope_diff, level=level)
        trows, _, _ = satterthwaite_unconstrained(dev, covf, res.x, beta, tc_log, level=level)
        from .time_to_event import _holm
        pd_adj = _holm([r["p_unadjusted"] for r in drows])
        for i, dr, tr, p in zip(others, drows, trows, pd_adj):
            model_rows.append({"arm": arms[i], "control_arm": arms[c], "log_rate_difference": dr["estimate"], "rate_difference_lower": dr["ci_low"],
                               "rate_difference_upper": dr["ci_high"], "df": dr["df"], "p_adjusted": float(p),
                               "model_tc_geometric_ratio_at_analysis_day": float(np.exp(tr["estimate"])),
                               "model_tc_lower": float(np.exp(tr["ci_low"])), "model_tc_upper": float(np.exp(tr["ci_high"]))})
        fit_row["model_reportable"] = True
        fit_row["_fitted"] = {"beta": beta.tolist()}
    except ArithmeticError:
        if "random_slope_variance_or_correlation_at_boundary_satterthwaite_withheld" not in diagnostics:
            diagnostics.append("growth_model_failed_or_hessian_not_positive_definite")
    # Observed TGI% and T/C% at the analysis day (animals measured on that day only).
    at_day = used[used.day == a["analysis_day"]].set_index("animal_id")
    base = used[used.day == a["baseline_day"]].set_index("animal_id").volume
    obs_rows = []
    c_ids = at_day.index[at_day.arm == a["control_arm"]]
    delta_c = at_day.volume[c_ids] - base[c_ids]
    others = [x for x in arms if x != a["control_arm"]]
    fam_level = 1 - (1 - level) / max(len(others), 1)
    for arm in others:
        t_ids = at_day.index[at_day.arm == arm]
        row = {"arm": arm, "control_arm": a["control_arm"], "n_arm_on_day": int(len(t_ids)), "n_control_on_day": int(len(c_ids)),
               "tgi_percent": None, "tgi_lower": None, "tgi_upper": None, "tc_percent": None, "tc_lower": None, "tc_upper": None}
        if len(t_ids) >= 2 and len(c_ids) >= 2:
            delta_t = at_day.volume[t_ids] - base[t_ids]
            if delta_c.mean() > 0:
                r, lo, hi, _ = fieller(delta_t.mean(), delta_t.var(ddof=1), len(delta_t), delta_c.mean(), delta_c.var(ddof=1), len(delta_c), fam_level)
                row.update(tgi_percent=float(100 * (1 - r)), tgi_lower=None if hi is None else float(100 * (1 - hi)),
                           tgi_upper=None if lo is None else float(100 * (1 - lo)))
            else:
                diagnostics.append("control_did_not_grow_tgi_not_defined")
            vt_, vc_ = at_day.volume[t_ids], at_day.volume[c_ids]
            r, lo, hi, _ = fieller(vt_.mean(), vt_.var(ddof=1), len(vt_), vc_.mean(), vc_.var(ddof=1), len(vc_), fam_level)
            row.update(tc_percent=float(100 * r), tc_lower=None if lo is None else float(100 * lo), tc_upper=None if hi is None else float(100 * hi))
        obs_rows.append(row)
    if dropout:
        diagnostics.append("observed_tgi_uses_only_animals_measured_on_analysis_day_biased_if_removal_relates_to_growth")
    if len(days) < 4:
        diagnostics.append("fewer_than_4_measurement_days_growth_rate_imprecise")
    fit_row["diagnostics"] = list(dict.fromkeys(diagnostics))
    fitted = fit_row.pop("_fitted", None)
    return {"schema_version": 1, "analysis_type": "tumor_growth", "fits": [fit_row], "arm_day_summaries": summaries,
            "growth_rates": growth, "growth_tests": tests, "model_contrasts": model_rows, "observed_tgi": obs_rows,
            "dropout": dropout, "fitted_coefficients": fitted}
