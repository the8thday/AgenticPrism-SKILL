"""Declared ADA cut points on complete repeated drug-naive subject panels."""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats, special
from .variance_components import crossed_anova

DEFAULTS = {"analysis_type": "ada_cut_point", "schema_version": 1, "input": None, "source": None,
    "assay": {"tier": None, "signal_unit": None, "population": None, "drug_naive": None,
        "representative_negative_panel": None, "independent_subjects": None,
        "balanced_crossed_design": None, "single_reagent_lot": None, "rationale": None},
    "cut_point": {"false_positive_rate": None, "target_rationale": None,
        "method": "parametric", "normalization": "none", "transformation": "none",
        "boxcox_lambda": None, "selection": "auto", "comparison_alpha": .05,
        "nc_correlation_min": .7, "normality_alpha": .05, "skewness_limit": 1.,
        "variance_confidence": .95, "bound": "point", "bound_confidence": .90, "bound_applicable": False, "bound_rationale": None, "independent_pairs": [],
        "nonparametric_bound": "auto", "bootstrap_reps": 2000, "bootstrap_seed": None},
    "outliers": {"prespecified": None, "rationale": None, "analytical": "none", "biological": "none", "iqr_multiplier": 3.},
    "report": {"plot_style": "prism_like"}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw)-set(DEFAULTS):
        raise ValueError("Unsupported ADA configuration")
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v)-set(c[k]):
                raise ValueError(f"Unsupported {k} settings")
            c[k].update(v)
        else:
            c[k] = v
    if c["analysis_type"] != "ada_cut_point" or type(c["schema_version"]) is not int or c["schema_version"] != 1:
        raise ValueError("Unsupported schema")
    a, q, o = c["assay"], c["cut_point"], c["outliers"]
    for obj, keys in ((c, ("input", "source")), (a, ("signal_unit", "population", "rationale")),
                      (q, ("target_rationale",)), (o, ("rationale",))):
        if any(not isinstance(obj[k], str) or not obj[k].strip() for k in keys):
            raise ValueError("Inputs, units, population and all rationales must be nonempty strings")
    if any(a[k] is not True for k in ("drug_naive", "representative_negative_panel", "independent_subjects", "balanced_crossed_design", "single_reagent_lot")) or o["prespecified"] is not True:
        raise ValueError("ADA applicability declarations and prespecified outlier policy require literal true")
    if a["tier"] not in ("screening", "confirmatory", "titer"):
        raise ValueError("Declare screening, confirmatory or titer")
    if q['bound'] not in ('point','lower'):raise ValueError('bound must be point or lower')
    if type(q['bound_confidence']) not in (int,float) or not .5<q['bound_confidence']<1:raise ValueError('Invalid bound_confidence')
    if q['bound']=='lower':
        if q['bound_applicable'] is not True or not isinstance(q['bound_rationale'],str) or not q['bound_rationale'].strip():raise ValueError('Lower bound requires literal true applicability and rationale for marginal independent future subject/run inference')
        if o['biological']!='none' or o['analytical']!='none':raise ValueError('Lower bound does not support data-dependent IQR handling; declare none')
        if q['selection']=='dynamic':raise ValueError('Dynamic lower-bound deployment is unsupported')
        pairs=q['independent_pairs']
        if q['nonparametric_bound'] not in ('auto','independent_pairs','two_way_bootstrap'):raise ValueError('nonparametric_bound must be auto, independent_pairs or two_way_bootstrap')
        if q['method']=='nonparametric':
            if q['nonparametric_bound']=='auto':q['nonparametric_bound']='independent_pairs' if pairs else 'two_way_bootstrap'
            if q['nonparametric_bound']=='independent_pairs':
                if not isinstance(pairs,list) or len(pairs)<3 or any(not isinstance(p,list) or len(p)!=2 or any(not isinstance(v,str) or not v.strip() for v in p) for p in pairs):raise ValueError('Nonparametric lower bound needs >=3 prespecified [subject_id,run_id] independent pairs')
                if len({p[0] for p in pairs})!=len(pairs) or len({p[1] for p in pairs})!=len(pairs):raise ValueError('Independent pairs must use distinct subjects AND distinct runs')
                if q['bootstrap_seed'] is not None or q['bootstrap_reps']!=2000:raise ValueError('Bootstrap settings apply only to two_way_bootstrap')
            else:
                if pairs:raise ValueError('independent_pairs apply only to the independent_pairs bound')
                if type(q['bootstrap_reps']) is not int or not 999<=q['bootstrap_reps']<=100000:raise ValueError('bootstrap_reps must be an integer from 999 to 100000')
                if type(q['bootstrap_seed']) is not int or q['bootstrap_seed']<0:raise ValueError('two_way_bootstrap requires a declared nonnegative integer bootstrap_seed')
        elif pairs or q['nonparametric_bound']!='auto' or q['bootstrap_seed'] is not None or q['bootstrap_reps']!=2000:raise ValueError('independent_pairs and bootstrap settings apply only to nonparametric lower bounds')
    elif q['bound_applicable'] is not False or q['bound_rationale'] is not None or q['independent_pairs'] or q['nonparametric_bound']!='auto' or q['bootstrap_seed'] is not None or q['bootstrap_reps']!=2000:raise ValueError('Bound declarations apply only to lower bounds')
    if q["method"] not in ("parametric", "nonparametric") or q["normalization"] not in ("none", "ratio", "difference"):
        raise ValueError("Unsupported percentile method or normalization")
    if q["transformation"] not in ("none", "log", "boxcox") or q["selection"] not in ("auto", "fixed", "floating", "dynamic"):
        raise ValueError("Unsupported transformation or selection")
    if q["transformation"] == "boxcox":
        if type(q["boxcox_lambda"]) not in (float, int) or not np.isfinite(q["boxcox_lambda"]):
            raise ValueError("Box-Cox requires a prespecified finite lambda; no data-driven optimization")
    elif q["boxcox_lambda"] is not None:
        raise ValueError("boxcox_lambda applies only to Box-Cox")
    for k, low, high in (("false_positive_rate", 0, .5), ("comparison_alpha", 0, .5),
                         ("normality_alpha", 0, .5), ("nc_correlation_min", 0, 1), ("variance_confidence", .5, 1)):
        if type(q[k]) not in (float, int) or not low < q[k] < high:
            raise ValueError(f"Invalid {k}")
    if type(q["skewness_limit"]) not in (float, int) or not 0 < q["skewness_limit"] <= 2:
        raise ValueError("Invalid skewness_limit")
    if a["tier"] == "confirmatory" and (q["normalization"] != "none" or q["transformation"] != "none"):
        raise ValueError("Confirmatory currently supports raw percent inhibition only; no log of negative inhibition")
    if q["selection"] == "floating" and q["normalization"] == "none":
        raise ValueError("Floating cut points require declared NC normalization")
    if o["analytical"] not in ("none", "iqr_flag") or o["biological"] not in ("none", "iqr_flag", "iqr_exclude"):
        raise ValueError("Unsupported outlier policy; analytical flags require investigation, not automatic deletion")
    if type(o["iqr_multiplier"]) not in (float, int) or not .5 <= o["iqr_multiplier"] <= 10:
        raise ValueError("Invalid IQR multiplier")
    if c["report"]["plot_style"] not in ("prism_like", "standard"):
        raise ValueError("Invalid plot style")
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"observation_id", "subject_id", "run_id", "replicate_id", "signal", "nc_signal", "inhibited_signal", "lot", "population", "treatment_status", "exclude", "exclusion_reason"}
    if not required <= set(d) or not len(d):
        raise ValueError("Missing required ADA columns or empty data")
    for k in required-{"inhibited_signal", "exclusion_reason"}:
        if (d[k].str.strip() == "").any():
            raise ValueError(f"Empty {k}")
    if d.observation_id.duplicated().any() or d.duplicated(["subject_id", "run_id", "replicate_id"]).any():
        raise ValueError("Duplicate observation or subject/run/replicate")
    if set(d.exclude)-{"true", "false"}:
        raise ValueError("exclude must be literal true/false")
    d["exclude"] = d.exclude == "true"
    if (d.exclude & (d.exclusion_reason.str.strip() == "")).any():
        raise ValueError("Excluded rows require a reason")
    if set(d.treatment_status) != {"drug_naive"} or set(d.population) != {cfg["assay"]["population"]} or d.lot.nunique() != 1:
        raise ValueError("Drug-naive matching population and a single declared reagent lot required")
    for k in ("signal", "nc_signal") + (("inhibited_signal",) if cfg["assay"]["tier"] == "confirmatory" else ()):
        d[k] = pd.to_numeric(d[k], errors="raise")
        if not np.isfinite(d[k]).all() or (d[k] <= 0).any():
            raise ValueError(f"Finite positive {k} required")
    if cfg["assay"]["tier"] != "confirmatory" and (d.inhibited_signal != "").any():
        raise ValueError("inhibited_signal applies only to confirmatory")
    if (d.groupby("run_id").nc_signal.nunique() != 1).any():
        raise ValueError("nc_signal must be the declared aggregate NC for that run, constant across rows")
    _panel(d[~d.exclude])
    return d


def _panel(d):
    counts = d.groupby(["subject_id", "run_id"]).size().unstack()
    if counts.shape[0] < 50 or counts.shape[1] < 3 or counts.isna().any().any() or counts.to_numpy().min() < 2 or np.unique(counts.to_numpy()).size != 1:
        raise ValueError("Require >=50 subjects measured in every >=3 runs, equal >=2 technical replicates per cell after exclusions")


def transform(x, q):
    x = np.asarray(x, float)
    if q["transformation"] == "none":
        return x
    if (x <= 0).any():
        raise ValueError("Log/Box-Cox requires positive values; no automatic shifts")
    return np.log(x) if q["transformation"] == "log" else special.boxcox(x, q["boxcox_lambda"])


def inverse(x, q):
    v = x if q["transformation"] == "none" else np.exp(x) if q["transformation"] == "log" else special.inv_boxcox(x, q["boxcox_lambda"])
    if not np.isfinite(v):
        raise ValueError("Cut point outside inverse transformation domain")
    return float(v)


def percentile(values, probability, method):
    x = np.asarray(values, float).ravel()
    if len(x) < 3 or not np.isfinite(x).all() or not 0 < probability < 1:
        raise ValueError("Invalid percentile inputs")
    if method == "parametric":
        return float(x.mean()+stats.norm.ppf(probability)*x.std(ddof=1))
    if method == "nonparametric":
        return float(np.quantile(x, probability, method="linear"))
    raise ValueError("Unsupported method")


def iqr_flags(x, multiplier):
    x = np.asarray(x, float)
    lo, hi = np.quantile(x, [.25, .75], method="linear")
    return (x < lo-multiplier*(hi-lo)) | (x > hi+multiplier*(hi-lo))


def panel_diagnostics(y):
    y = np.asarray(y, float)
    vc = crossed_anova(y)
    lev = stats.levene(*y.T, center="median")
    flat = y.ravel()
    # Shapiro nominal p assumes independence: diagnostic only for repeated cells.
    sw = stats.shapiro(flat) if 3 <= len(flat) <= 5000 else None
    return {"run_mean_f": vc["factor_b_f"], "run_mean_p": vc["factor_b_p"],
            "run_variance_statistic": float(lev.statistic), "run_variance_p": float(lev.pvalue),
            "shapiro_w": float(sw.statistic) if sw else None, "shapiro_p": float(sw.pvalue) if sw else None,
            "skewness": float(stats.skew(flat, bias=False)),
            "test_caveat": "Run F uses subject-blocked ANOVA. Median Levene and pooled Shapiro p are descriptive with repeated subjects; non-significance is not equivalence."}


def compute(d, cfg):
    a, q, o = cfg["assay"], cfg["cut_point"], cfg["outliers"]
    used = d[~d.exclude].copy()
    analytical = []
    if o["analytical"] == "iqr_flag":
        residual = used.signal-used.groupby(["subject_id", "run_id"]).signal.transform("median")
        for _, indices in used.groupby("run_id").groups.items():
            flags = iqr_flags(residual.loc[indices], o["iqr_multiplier"])
            analytical.extend(used.loc[indices].loc[flags, "observation_id"].tolist())
    numeric = ["signal", "nc_signal"] + (["inhibited_signal"] if a["tier"] == "confirmatory" else [])
    cells = used.groupby(["subject_id", "run_id"], sort=True)[numeric].mean().reset_index()
    raw = cells.signal.to_numpy() if a["tier"] != "confirmatory" else 100*(1-cells.inhibited_signal.to_numpy()/cells.signal.to_numpy())
    norm = raw/cells.nc_signal.to_numpy() if q["normalization"] == "ratio" else raw-cells.nc_signal.to_numpy() if q["normalization"] == "difference" else raw
    cells["raw_value"], cells["normalized_value"] = raw, norm
    cells["transformed_value"] = transform(norm, q)
    biological = []
    means = cells.groupby("subject_id").transformed_value.mean()
    if o["biological"] != "none":
        biological = means.index[iqr_flags(means, o["iqr_multiplier"])].tolist()
    cells["biological_flag"] = cells.subject_id.isin(biological)
    cells["used"] = ~cells.biological_flag if o["biological"] == "iqr_exclude" else True
    clean = cells[cells.used]
    if clean.subject_id.nunique() < 50:
        raise ValueError("Fewer than 50 subjects remain after declared biological exclusion")
    y = clean.pivot(index="subject_id", columns="run_id", values="transformed_value")
    raw_y = clean.assign(raw_transformed=transform(clean.raw_value, q)).pivot(index="subject_id", columns="run_id", values="raw_transformed")
    diag, raw_diag = panel_diagnostics(y), panel_diagnostics(raw_y)
    alpha = q["comparison_alpha"]
    reasons, cautions = [], []
    if analytical:
        reasons.append("Analytical IQR flags require investigation and a documented rerun; no automatic well deletion.")
    if biological and o["biological"] == "iqr_flag":
        cautions.append("Biological IQR flags retained under the prespecified flag-only policy.")
    corr = None
    nc = clean.groupby("run_id").nc_signal.first().reindex(y.columns).to_numpy()
    if np.ptp(nc) > 0 and np.ptp(raw_y.mean(0)) > 0:
        corr = float(np.corrcoef(transform(nc, q), raw_y.mean(0))[0, 1])
    if q['bound']=='lower':
        selection='fixed' if q['normalization']=='none' else 'floating'
        cautions.append('Lower bound targets the marginal distribution over independent future subjects and random runs, not every realized run. Run diagnostics remain descriptive.')
        if q['normalization']!='none' and (corr is None or corr<q['nc_correlation_min']):reasons.append('Negative-control tracking fails the prespecified correlation criterion.')
    elif diag["run_variance_p"] < alpha:
        selection = "dynamic"
    elif diag["run_mean_p"] < alpha:
        selection = "dynamic"
        cautions.append("Residual run mean differences remain; a common normalized factor is unsupported.")
    elif q["normalization"] != "none" and raw_diag["run_mean_p"] < alpha:
        selection = "floating"
        if corr is None or corr < q["nc_correlation_min"]:
            reasons.append("Negative-control tracking fails the prespecified correlation criterion.")
    else:
        selection = "fixed" if q["normalization"] == "none" else "floating"
        if selection == "floating":
            reasons.append("Raw run mean shift not demonstrated; normalization benefit is unestablished in this scoped decision rule.")
    if q["selection"] != "auto" and q["selection"] != selection:
        reasons.append(f"Requested {q['selection']} conflicts with diagnostic selection {selection}.")
    if q["method"] == "parametric" and abs(diag["skewness"]) > q["skewness_limit"] and (diag["shapiro_p"] is None or diag["shapiro_p"] < q["normality_alpha"]):
        reasons.append("Prespecified normality/skewness gate fails; do not switch methods after seeing results without a revised protocol.")
    if selection == "dynamic":
        reasons.append("Heterogeneous runs: dynamic is a recommendation only; per-run audit cut points do not validate future-run deployment.")
    probability = 1-q["false_positive_rate"]
    audit_cp = inverse(percentile(y, probability, q["method"]), q)
    bound_result=None
    if q['bound']=='lower':
        if d.exclude.any():raise ValueError('Lower-bound scope requires the complete prespecified panel without outcome-dependent exclusions')
        from .ada_bounds import parametric_panel_lower,order_statistic_lower,two_way_bootstrap_lower
        if q['method']=='parametric':bound_result=parametric_panel_lower(y.to_numpy(),probability,q['bound_confidence'])
        elif q['nonparametric_bound']=='two_way_bootstrap':
            bound_result=two_way_bootstrap_lower(y.to_numpy(),probability,q['bound_confidence'],q['bootstrap_reps'],q['bootstrap_seed'])
            if not bound_result['tail_support_sufficient']:reasons.append('Too few subjects for a nonparametric bound on this percentile: n_subjects x FPR is below 3; collect more subjects or use a prespecified parametric bound.')
        else:
            pairs=q['independent_pairs']
            if any(s not in y.index or r not in y.columns for s,r in pairs):raise ValueError('A declared independent pair is absent from the panel')
            bound_result=order_statistic_lower([y.loc[s,r] for s,r in pairs],probability,q['bound_confidence'])
            bound_result['independent_pairs']=pairs
        audit_cp=inverse(bound_result['lower_bound'],q)
        bound_result['decision_scale_lower_bound']=audit_cp
    per_run = []
    for run in y:
        cp = inverse(percentile(y[run], probability, q["method"]), q)
        per_run.append({"run_id": run, "audit_cut_point": cp, "n_subjects": len(y), "reportable_for_future_run": False})
    vc = crossed_anova(y, q["variance_confidence"])
    vc["factor_labels"] = {"factor_a": "subject", "factor_b": "run", "residual": "subject_by_run_plus_cell_mean_error"}
    limitations = ["Cut point is an assay decision threshold, not ADA incidence, concentration or clinical risk.",
        "Single population and reagent lot; disease status alone does not exclude a drug-naive representative panel.",
        "Equal technical replicates are averaged before analysis; inference unit is the subject.",
        "Point mode has no cut-point confidence interval; lower mode is a separate opt-in estimand.",
        "Pooled percentiles estimate this validation panel's marginal response distribution; no validated future-run prediction guarantee.",
        "Variance decomposition separates subject/run/residual only; analyst/day/plate effects and interactions are not separately estimated.",
        "Normality, Levene and correlation diagnostics do not prove assay suitability; observed in-study false positives must be checked."]
    if len(y)*q["false_positive_rate"] < 1:
        cautions.append("Fewer than one expected subject beyond the requested percentile per run; tail resolution is inadequate for a reliable empirical extreme-tail claim.")
    if q["bound"] == "point" and q["method"] == "nonparametric" and len(y)*q["false_positive_rate"] < 1:
        reasons.append("Nonparametric extreme percentile lacks independent-subject tail support.")
    fit = {"status": "withheld" if reasons else "ok", "reportable": not reasons, "tier": a["tier"],
        "selected_type": selection, "requested_type": q["selection"], "false_positive_rate_target": q["false_positive_rate"],
        "method": q["method"], "cut_point": None if reasons else audit_cp, "audit_cut_point": audit_cp,
        "scale": "percent_inhibition" if a["tier"] == "confirmatory" else q["normalization"],
        "n_subjects": len(y), "n_runs": y.shape[1], "normalization": q["normalization"], "transformation": q["transformation"],
        "diagnostics": reasons+cautions, "withholding_reasons": reasons, "raw_run_diagnostics": raw_diag,
        "analysis_diagnostics": diag, "nc_run_correlation": corr, "variance_components": vc,
        "biological_flag_subjects": biological, "analytical_flag_observations": analytical,
        "deployment": "raw signal >= cut_point" if q["normalization"] == "none" else "signal >= run_NC * cut_point" if q["normalization"] == "ratio" else "signal >= run_NC + cut_point"}
    if bound_result is not None:
        fit['cut_point_bound']=bound_result
        limitations=[v for v in limitations if not v.startswith('Point mode has no cut-point confidence') and not v.startswith('Pooled percentiles')]+bound_result['limitations']
    if a["tier"] == "confirmatory":
        fit["deployment"] = "100*(1 - inhibited_signal/uninhibited_signal) >= cut_point"
    return {"schema_version": 1, "analysis_type": "ada_cut_point", "fits": [fit], "cells": cells.to_dict("records"),
            "per_run_audit": per_run, "limitations": limitations}
