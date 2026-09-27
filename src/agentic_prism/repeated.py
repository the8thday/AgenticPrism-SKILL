"""One-factor repeated measures and Gaussian random-intercept models.

MixedLM estimates variance components by REML. Fixed-effect covariance is the
plug-in GLS covariance (X' V^-1 X)^-1. Inference options: Satterthwaite t/F
(lmerTest-style, see mixed_inference.py); a parametric bootstrap that draws
whole subject vectors, refits REML and calibrates centered Wald and max-|t|
pivots; or asymptotic Wald. None is Kenward-Roger, a likelihood-ratio
bootstrap, or a distribution-free method.
"""
from copy import deepcopy
from itertools import combinations
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import f as f_dist, t as t_dist, norm, chi2
from statsmodels.regression.mixed_linear_model import MixedLM
from .groups import read_unit_table
from .mixed_inference import RandomInterceptModel, satterthwaite

DEFAULTS = {
    "analysis_type": "repeated_measures", "schema_version": 1, "input": None,
    "source": "User-supplied unit-level repeated measurements",
    "comparison": {"design": None, "groups": None, "outcome": None, "unit": None,
        "rationale": None, "post_hoc": "none", "control": None, "confidence_level": .95,
        "missing_policy": "require_complete", "missingness_rationale": None,
        "covariance_rationale": None, "inference": None, "arms": None, "control_arm": None,
        "bootstrap_reps": 999, "seed": 20260927},
    "report": {"plot_style": "prism_like"}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported repeated-measures config")
    cfg = deepcopy(DEFAULTS)
    for key, val in raw.items():
        if isinstance(cfg[key], dict):
            if not isinstance(val, dict) or set(val) - set(cfg[key]):
                raise ValueError(f"Unsupported {key} settings")
            cfg[key].update(val)
        else:
            cfg[key] = val
    q = cfg["comparison"]
    if cfg["analysis_type"] != "repeated_measures" or type(cfg["schema_version"]) is not int or cfg["schema_version"] != 1:
        raise ValueError("Unsupported repeated-measures schema")
    for val in [cfg["input"], cfg["source"], *[q[k] for k in ("outcome", "unit", "rationale")]]:
        if not isinstance(val, str) or not val.strip():
            raise ValueError("Declare input, source, outcome, unit and rationale")
    if q["design"] not in ("one_way_rm_anova", "random_intercept", "two_way_rm_anova", "two_way_mixed"):
        raise ValueError("Supported design: one_way_rm_anova, random_intercept, two_way_rm_anova or two_way_mixed")
    two_way = q["design"].startswith("two_way")
    groups = q["groups"]
    if (not isinstance(groups, list) or len(groups) < (2 if two_way else 3) or any(not isinstance(g, str) or not g.strip() for g in groups)
            or len(set(groups)) != len(groups)):
        raise ValueError("Declare distinct ordered within-unit conditions (at least three; two for between-by-within designs)")
    if two_way:
        arms = q["arms"]
        if not isinstance(arms, list) or len(arms) < 2 or any(not isinstance(a, str) or not a.strip() for a in arms) or len(set(arms)) != len(arms):
            raise ValueError("Between-by-within designs need at least two distinct declared arms")
        if q["post_hoc"] not in ("none", "arm_vs_control_each_condition", "condition_vs_first_each_arm"):
            raise ValueError("post_hoc for between-by-within designs: none, arm_vs_control_each_condition or condition_vs_first_each_arm")
        if q["post_hoc"] == "arm_vs_control_each_condition":
            if q["control_arm"] not in arms:
                raise ValueError("Declare control_arm as one of the arms")
        elif q["control_arm"] is not None:
            raise ValueError("control_arm applies only to arm_vs_control_each_condition")
        if q["control"] is not None:
            raise ValueError("Use control_arm, not control, in between-by-within designs")
    else:
        if q["arms"] is not None or q["control_arm"] is not None:
            raise ValueError("arms/control_arm apply only to between-by-within designs")
        if q["post_hoc"] not in ("none", "vs_control", "all_pairs"):
            raise ValueError("Unsupported post_hoc family")
        if q["post_hoc"] == "vs_control":
            if q["control"] not in groups:
                raise ValueError("Declare a control in groups")
        elif q["control"] is not None:
            raise ValueError("control applies only to vs_control")
    if type(q["confidence_level"]) not in (float, int) or not .5 < q["confidence_level"] < 1:
        raise ValueError("Invalid confidence_level")
    if q["missing_policy"] not in ("require_complete", "available_case"):
        raise ValueError("Unsupported missing_policy")
    if q["missing_policy"] == "available_case" and (not isinstance(q["missingness_rationale"], str) or not q["missingness_rationale"].strip()):
        raise ValueError("Available-case fitting requires a missingness_rationale supporting MAR conditional on the model")
    if q["design"] in ("one_way_rm_anova", "two_way_rm_anova"):
        if q["missing_policy"] != "require_complete" or q["inference"] not in (None, "greenhouse_geisser"):
            raise ValueError("RM ANOVA requires complete data and Greenhouse-Geisser inference")
        q["inference"] = "greenhouse_geisser"
    elif q["design"] == "two_way_mixed":
        if q["inference"] not in ("satterthwaite", "wald_asymptotic"):
            raise ValueError("two_way_mixed supports satterthwaite or wald_asymptotic inference")
        if not isinstance(q["covariance_rationale"], str) or not q["covariance_rationale"].strip():
            raise ValueError("Declare covariance_rationale for Gaussian random intercept + iid residuals")
    else:
        if q["inference"] not in ("satterthwaite", "parametric_bootstrap", "wald_asymptotic"):
            raise ValueError("Explicitly select satterthwaite, parametric_bootstrap or wald_asymptotic inference")
        if not isinstance(q["covariance_rationale"], str) or not q["covariance_rationale"].strip():
            raise ValueError("Declare covariance_rationale for Gaussian random intercept + iid residuals")
    if type(q["bootstrap_reps"]) is not int or not 199 <= q["bootstrap_reps"] <= 10000:
        raise ValueError("bootstrap_reps must be an integer from 199 to 10000")
    if type(q["seed"]) is not int or not 0 <= q["seed"] < 2**32:
        raise ValueError("seed must be an unsigned 32-bit integer")
    if cfg["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    return cfg


def load_data(path, cfg):
    d = read_unit_table(path)
    q = cfg["comparison"]
    if q["design"].startswith("two_way"):
        return _load_two_way(d, q)
    if set(d.group) != set(q["groups"]) or set(d.outcome) != {q["outcome"]} or set(d.unit) != {q["unit"]}:
        raise ValueError("Input group/outcome/unit does not match declaration")
    used = d[~d.exclude]
    if used.duplicated(["independent_unit_id", "group"]).any():
        raise ValueError("Duplicate unit/condition; aggregate technical replicates upstream with provenance")
    if set(used.group) != set(q["groups"]):
        raise ValueError("Every declared condition needs included observations")
    counts = used.groupby("independent_unit_id").size()
    if q["missing_policy"] == "require_complete" and not counts.eq(len(q["groups"])).all():
        raise ValueError("Incomplete repeated measurements; no silent complete-case deletion")
    minimum = 3 if q["design"] == "one_way_rm_anova" else 6
    if len(counts) < minimum or used.groupby("group").size().min() < minimum:
        raise ValueError(f"At least {minimum} units overall and per condition are required")
    if q["design"] == "random_intercept":
        if counts.ge(2).sum() < 6:
            raise ValueError("At least six units with repeated measurements are required")
        # Require an overlap-connected within-unit design. Disjoint condition
        # sets must not be reinterpreted as a repeated-measures experiment.
        reached = {q["groups"][0]}
        sets = list(used.groupby("independent_unit_id").group.apply(set))
        for _ in q["groups"]:
            for conditions in sets:
                if reached & conditions:
                    reached |= conditions
        if reached != set(q["groups"]):
            raise ValueError("Disconnected within-unit condition design")
    return d


def _load_two_way(d, q):
    if "arm" not in d or d.arm.str.strip().eq("").any():
        raise ValueError("Between-by-within designs need a nonempty arm column")
    if set(d.group) != set(q["groups"]) or set(d.arm) != set(q["arms"]) or set(d.outcome) != {q["outcome"]} or set(d.unit) != {q["unit"]}:
        raise ValueError("Input arm/group/outcome/unit does not match declaration")
    if d.groupby("independent_unit_id").arm.nunique().gt(1).any():
        raise ValueError("Each independent unit must belong to exactly one arm")
    used = d[~d.exclude]
    if used.duplicated(["independent_unit_id", "group"]).any():
        raise ValueError("Duplicate unit/condition; aggregate technical replicates upstream with provenance")
    counts = used.groupby("independent_unit_id").size()
    if q["missing_policy"] == "require_complete" and not counts.eq(len(q["groups"])).all():
        raise ValueError("Incomplete repeated measurements; no silent complete-case deletion")
    per_arm = used.groupby("arm").independent_unit_id.nunique().reindex(q["arms"]).fillna(0)
    if per_arm.min() < 3:
        raise ValueError("At least three units per arm are required")
    cell = used.groupby(["arm", "group"]).size().reindex(pd.MultiIndex.from_product([q["arms"], q["groups"]])).fillna(0)
    if cell.min() < 3:
        raise ValueError("Every arm x condition cell needs at least three observed units")
    return d


def _pairs(q):
    groups = q["groups"]
    if q["post_hoc"] == "none":
        return []
    if q["post_hoc"] == "vs_control":
        return [(groups.index(q["control"]), i) for i, g in enumerate(groups) if g != q["control"]]
    return list(combinations(range(len(groups)), 2))


def _holm(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    return adjusted


def rm_anova(matrix, q):
    """Balanced one-factor within-subject ANOVA; always apply GG correction."""
    y = np.asarray(matrix, float)
    n, k = y.shape
    centered = y - y.mean()
    residual = centered - centered.mean(axis=0) - centered.mean(axis=1)[:, None]
    ss_condition = n * np.sum(centered.mean(axis=0)**2)
    ss_error = np.sum(residual**2)
    if ss_error <= np.finfo(float).eps * max(np.sum(centered**2), np.finfo(float).tiny):
        raise ValueError("Zero residual variance; repeated-measures inference withheld")
    df1, df2 = k - 1, (n - 1) * (k - 1)
    stat = ss_condition / df1 / (ss_error / df2)
    h = np.eye(k) - np.ones((k, k)) / k
    s = h @ np.cov(y, rowvar=False, ddof=1) @ h
    epsilon = float(np.clip(np.trace(s)**2 / ((k-1)*np.sum(s*s)), 1/(k-1), 1))
    omnibus = {"f_statistic": float(stat), "df_numerator_uncorrected": df1,
        "df_denominator_uncorrected": df2, "epsilon_gg": epsilon,
        "df_numerator": epsilon*df1, "df_denominator": epsilon*df2,
        "p_uncorrected": float(f_dist.sf(stat, df1, df2)),
        "p_value": float(f_dist.sf(stat, epsilon*df1, epsilon*df2)),
        "ss_condition": float(ss_condition), "ss_error": float(ss_error),
        "partial_eta_squared": float(ss_condition/(ss_condition+ss_error)),
        "status": "estimated", "reportable": True, "diagnostics": [],
        "contrast_method": "paired_t_holm_p_bonferroni_ci"}
    rows = []
    pairs = _pairs(q)
    for a, b in pairs:
        delta = y[:, b] - y[:, a]
        diff, se = float(delta.mean()), float(delta.std(ddof=1)/np.sqrt(n))
        row = {"group_a": q["groups"][a], "group_b": q["groups"][b], "difference_b_minus_a": diff,
               "standard_error": se, "df": n-1, "reportable": se > 0}
        if se <= 0:
            row.update(statistic=None, p_unadjusted=None, p_adjusted=None, ci_low=None, ci_high=None,
                       diagnostic="zero_variance_paired_difference")
            omnibus["diagnostics"].append("some_contrasts_have_zero_variance")
        else:
            t = diff/se
            margin = t_dist.ppf(1-(1-q["confidence_level"])/(2*len(pairs)), n-1)*se
            row.update(statistic=float(t), p_unadjusted=float(2*t_dist.sf(abs(t), n-1)),
                       ci_low=float(diff-margin), ci_high=float(diff+margin), diagnostic="")
        rows.append(row)
    # Keep the original family size when a contrast is untestable.
    adjusted = _holm([r["p_unadjusted"] if r["p_unadjusted"] is not None else 1 for r in rows])
    for r, p in zip(rows, adjusted):
        r["p_adjusted"] = float(p) if r["reportable"] else None
    return omnibus, rows


def fit_random_intercept(y, x, groups, polish=True):
    """Fit a scaled Gaussian model, then compute explicit GLS covariance/BLUPs."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    offset, scale = float(y.mean()), float(y.std())
    if scale <= 0 or not np.isfinite(scale):
        raise ValueError("Zero outcome variance")
    z = (y-offset)/scale
    attempts, candidates = [], []
    for method in ("bfgs", "powell"):
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                fit = MixedLM(z, x, groups=groups).fit(reml=True, method=method, disp=False, maxiter=500)
            attempts.append({"method": method, "converged": bool(fit.converged),
                             "warnings": [str(w.message) for w in caught]})
            if fit.converged and np.isfinite(fit.llf) and np.isfinite(fit.scale) and fit.scale > 0 and np.isfinite(fit.cov_re).all():
                candidates.append(fit)
                if method == "bfgs":
                    break
        except (ValueError, np.linalg.LinAlgError) as exc:
            attempts.append({"method": method, "converged": False, "error": str(exc)})
    if not candidates:
        return {"ok": False, "attempts": attempts}
    fit = max(candidates, key=lambda f: f.llf)
    tau, sigma = float(fit.cov_re[0, 0])*scale**2, float(fit.scale)*scale**2
    if tau < 0 or sigma / scale**2 < 1e-12:
        return {"ok": False, "attempts": attempts}
    llf = float(fit.llf-(len(y)-x.shape[1])*np.log(scale))
    if polish and tau/sigma >= 1e-6:
        # Polish the interior REML optimum on the exact deviance (statsmodels stops early by ~1e-5 in theta);
        # this makes Satterthwaite df agree with lmerTest to ~1e-7. Accepted only if the deviance decreases.
        model = RandomInterceptModel(y, x, groups)
        deviance = lambda z: model.deviance(np.exp(z))
        start = np.log([tau, sigma])
        polished = minimize(deviance, start, method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 4000})
        if polished.success and polished.fun < deviance(start):
            llf -= (polished.fun - deviance(start)) / 2
            tau, sigma = map(float, np.exp(polished.x))
            attempts.append({"method": "nelder_mead_reml_polish", "converged": True, "warnings": []})
    info = np.zeros((x.shape[1], x.shape[1]))
    rhs = np.zeros(x.shape[1])
    blocks = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    for ix in blocks:
        v_inv = (np.eye(len(ix)) - tau/(sigma+len(ix)*tau)*np.ones((len(ix),len(ix))))/sigma
        info += x[ix].T @ v_inv @ x[ix]
        rhs += x[ix].T @ v_inv @ (y[ix]-offset)
    covariance = np.linalg.inv(info)
    beta = covariance @ rhs + offset  # cell-means coding: every column is one condition
    if not np.isfinite(covariance).all() or np.linalg.eigvalsh(covariance).min() <= 0:
        return {"ok": False, "attempts": attempts}
    marginal = x @ beta
    random = np.zeros(len(y))
    for ix in blocks:
        random[ix] = tau/(sigma+len(ix)*tau)*np.sum(y[ix]-marginal[ix])
    return {"ok": True, "beta": beta, "covariance": covariance, "tau": tau, "sigma": sigma,
            "boundary": tau/sigma < 1e-6, "marginal": marginal, "conditional": marginal+random,
            "random_intercept": random, "attempts": attempts,
            "reml_log_likelihood": llf}


def _wald(beta, covariance, contrast):
    delta = contrast @ beta
    cv = contrast @ covariance @ contrast.T
    return float(delta @ np.linalg.solve(cv, delta))


def mixed_model(y, x, groups, q):
    fit = fit_random_intercept(y, x, groups)
    omnibus = {"status": "failed", "reportable": False, "diagnostics": [],
                "p_value": None, "optimizer_attempts": fit["attempts"], "inference": q["inference"],
                "contrast_method": {"parametric_bootstrap": "bootstrap_max_abs_t",
                                    "satterthwaite": "satterthwaite_t_holm_p_bonferroni_ci"}.get(q["inference"], "wald_z_holm_p_bonferroni_ci")}
    if not fit["ok"]:
        omnibus["diagnostics"] = ["mixed_model_did_not_converge_or_invalid_covariance"]
        return omnibus, [], [], [], []
    beta, cov = fit["beta"], fit["covariance"]
    k = len(beta)
    contrast = np.eye(k)[1:] - np.eye(k)[0]
    w = _wald(beta, cov, contrast)
    omnibus.update(status="estimated", reportable=True, wald_statistic=w, df_numerator=k-1,
                    random_intercept_variance=fit["tau"], residual_variance=fit["sigma"],
                    icc=fit["tau"]/(fit["tau"]+fit["sigma"]), reml_log_likelihood=fit["reml_log_likelihood"])
    if fit["boundary"]:
        omnibus["diagnostics"].append("random_intercept_variance_near_zero_boundary")
    if len(np.unique(groups)) < 30:
        omnibus["diagnostics"].append("fewer_than_30_units_small_sample_inference_requires_caution")
    pairs = _pairs(q)
    c = np.array([np.eye(k)[b]-np.eye(k)[a] for a,b in pairs]).reshape(-1,k)
    diffs = c @ beta
    ses = np.sqrt(np.einsum('ij,jk,ik->i',c,cov,c))
    rows = [{"group_a":q["groups"][a], "group_b":q["groups"][b],
             "difference_b_minus_a":float(diff), "standard_error":float(se),
             "statistic":float(diff/se), "df":None, "reportable":True}
            for (a,b),diff,se in zip(pairs,diffs,ses)]
    samples = []
    if q["inference"] == "satterthwaite":
        try:
            if fit["boundary"]:
                raise ArithmeticError("variance boundary")
            model = RandomInterceptModel(y, x, groups)
            srows, som, _ = satterthwaite(model, [fit["tau"], fit["sigma"]], beta, c if len(rows) else [],
                                          joint=contrast, level=q["confidence_level"])
        except (ArithmeticError, np.linalg.LinAlgError):
            omnibus.update(status="limited", reportable=False)
            omnibus["diagnostics"].append("satterthwaite_unavailable_at_variance_boundary_or_singular_hessian_consider_bootstrap")
            for r in rows:
                r.update(reportable=False, p_unadjusted=None, p_adjusted=None, ci_low=None, ci_high=None)
        else:
            omnibus.update(f_statistic=som["f_statistic"], df_denominator=som["df_denominator"],
                           p_value=som["p_value"], satterthwaite_component_df=som["component_df"])
            adjusted = _holm([r["p_unadjusted"] for r in srows]) if srows else []
            for r, sr, pa in zip(rows, srows, adjusted):
                r.update(df=sr["df"], statistic=sr["statistic"], p_unadjusted=sr["p_unadjusted"], p_adjusted=float(pa),
                         ci_low=sr["ci_low"], ci_high=sr["ci_high"])
    elif q["inference"] == "wald_asymptotic":
        omnibus["p_value"] = float(chi2.sf(w,k-1))
        omnibus["diagnostics"].append("asymptotic_wald_no_small_sample_df_correction")
        if rows:
            adjusted = _holm([2*norm.sf(abs(r["statistic"])) for r in rows])
            for r,p in zip(rows, adjusted):
                margin = norm.ppf(1-(1-q["confidence_level"])/(2*len(rows)))*r["standard_error"]
                r.update(p_unadjusted=float(2*norm.sf(abs(r["statistic"]))), p_adjusted=float(p),
                         ci_low=r["difference_b_minus_a"]-margin, ci_high=r["difference_b_minus_a"]+margin)
    else:
        rng = np.random.default_rng(q["seed"])
        _, gi = np.unique(groups, return_inverse=True)
        draws, n_boundary = [], 0
        for b in range(q["bootstrap_reps"]):
            simulated = x@beta + rng.normal(0,np.sqrt(fit["tau"]),gi.max()+1)[gi] + rng.normal(0,np.sqrt(fit["sigma"]),len(y))
            bf = fit_random_intercept(simulated,x,groups,polish=False)
            record = {"draw":b, "success":bool(bf["ok"]), "wald_centered":None, "max_abs_t":None,
                      "boundary":None, "optimizer_attempts":bf["attempts"]}
            if bf["ok"]:
                delta = bf["beta"]-beta
                wb = _wald(delta,bf["covariance"],contrast)
                tb = (c@delta)/np.sqrt(np.einsum('ij,jk,ik->i',c,bf["covariance"],c))
                draws.append((wb,tb))
                n_boundary += int(bf["boundary"])
                record.update(wald_centered=wb, max_abs_t=float(np.max(abs(tb))) if len(tb) else None,
                              boundary=bool(bf["boundary"]))
                for j,t in enumerate(tb):
                    record[f"contrast_{j}_t"] = float(t)
            samples.append(record)
        n = len(draws)
        omnibus.update(bootstrap_requested=q["bootstrap_reps"], bootstrap_successful=n,
                       bootstrap_boundary_fits=n_boundary, seed=q["seed"], p_resolution=1/(n+1))
        omnibus["diagnostics"].append("parametric_bootstrap_assumes_gaussian_random_intercept_and_iid_residuals")
        if n < 199 or n < .95*q["bootstrap_reps"]:
            omnibus.update(status="limited", reportable=False)
            omnibus["diagnostics"].append("insufficient_successful_bootstrap_refits_inference_withheld")
            for r in rows:
                r.update(reportable=False,p_unadjusted=None,p_adjusted=None,ci_low=None,ci_high=None)
        else:
            bw = np.array([v[0] for v in draws])
            omnibus["p_value"] = float((1+np.sum(bw>=w))/(n+1))
            omnibus["p_mc_standard_error"] = float(np.sqrt(omnibus["p_value"]*(1-omnibus["p_value"])/(n+1)))
            if n < q["bootstrap_reps"]:
                omnibus["diagnostics"].append("some_bootstrap_refits_failed_see_saved_draws")
            if rows:
                bt = np.abs(np.array([v[1] for v in draws]))
                maxima = bt.max(axis=1)
                critical = float(np.quantile(maxima,q["confidence_level"],method="higher"))
                for j,r in enumerate(rows):
                    obs = abs(r["statistic"])
                    r.update(p_unadjusted=float((1+np.sum(bt[:,j]>=obs))/(n+1)),
                             p_adjusted=float((1+np.sum(maxima>=obs))/(n+1)),
                             ci_low=r["difference_b_minus_a"]-critical*r["standard_error"],
                             ci_high=r["difference_b_minus_a"]+critical*r["standard_error"])
    fixed = [{"group":g,"estimated_mean":float(v),"standard_error":float(np.sqrt(cov[i,i]))}
             for i,(g,v) in enumerate(zip(q["groups"],beta))]
    residuals = [{"marginal_fitted":float(m),"conditional_fitted":float(cy),
                  "random_intercept_blup":float(b),"marginal_residual":float(yy-m),
                  "conditional_residual":float(yy-cy)}
                 for yy,m,cy,b in zip(y,fit["marginal"],fit["conditional"],fit["random_intercept"])]
    return omnibus,rows,fixed,residuals,samples


def compare(d, cfg):
    q = cfg["comparison"]
    if q["design"].startswith("two_way"):
        return compare_two_way(d, cfg)
    used = d[~d.exclude].sort_values(["independent_unit_id", "group"], kind="stable")
    matrix = used.pivot(index="independent_unit_id", columns="group", values="value").reindex(columns=q["groups"])
    missing = [{"independent_unit_id":uid,"group":g,"observed":bool(pd.notna(matrix.loc[uid,g]))}
               for uid in matrix.index for g in matrix.columns]
    summaries = [{"group":g,"n":int(matrix[g].count()),"mean":float(matrix[g].mean()),
                  "sd":float(matrix[g].std(ddof=1))} for g in matrix.columns]
    residuals, samples, fixed = [], [], []
    if q["design"] == "one_way_rm_anova":
        omnibus,rows = rm_anova(matrix.to_numpy(float),q)
    else:
        groups = used.independent_unit_id.to_numpy(str)
        x = np.column_stack([(used.group==g).to_numpy(float) for g in q["groups"]])
        omnibus,rows,fixed,residuals,samples = mixed_model(used.value.to_numpy(float),x,groups,q)
        for r,(_,observation) in zip(residuals,used.iterrows()):
            r.update(observation_id=observation.observation_id,independent_unit_id=observation.independent_unit_id,
                     group=observation.group)
    omnibus.update(analysis_type="repeated_measures",design=q["design"],inference=q["inference"],
                    n_units=len(matrix),n_observations=len(used),n_groups=len(q["groups"]),
                    n_missing_cells=int(matrix.isna().sum().sum()),confidence_level=q["confidence_level"])
    if matrix.isna().any().any():
        omnibus["diagnostics"].append("missing_measurements_available_case_MAR_assumption_not_testable")
    return {"schema_version":1,"analysis_type":"repeated_measures","fits":[omnibus],
            "group_summaries":summaries,"contrasts":rows,"fixed_effects":fixed,
            "residuals":residuals,"missingness":missing,"bootstrap_samples":samples}


def compare_two_way(d, cfg):
    from . import repeated_two_way as tw
    q = cfg["comparison"]
    arms, conds = q["arms"], q["groups"]
    used = d[~d.exclude].copy()
    used["arm_index"] = used.arm.map({a: i for i, a in enumerate(arms)})
    used["condition_index"] = used.group.map({c: i for i, c in enumerate(conds)})
    used = used.sort_values(["independent_unit_id", "condition_index"], kind="stable")
    unit_arm = used.groupby("independent_unit_id").arm.first()
    matrix = used.pivot(index="independent_unit_id", columns="group", values="value").reindex(columns=conds)
    missing = [{"independent_unit_id": uid, "arm": unit_arm[uid], "group": g, "observed": bool(pd.notna(matrix.loc[uid, g]))}
               for uid in matrix.index for g in conds]
    summaries = [{"arm": a, "group": g, "n": int(len(v)), "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else None}
                 for a in arms for g in conds for v in [used.value[(used.arm == a) & (used.group == g)].astype(float)]]
    fit_row = {"analysis_type": "repeated_measures", "design": q["design"], "inference": q["inference"],
               "n_units": int(len(matrix)), "n_observations": int(len(used)), "n_arms": len(arms), "n_groups": len(conds),
               "n_missing_cells": int(matrix.isna().sum().sum()), "confidence_level": q["confidence_level"],
               "status": "estimated", "reportable": True, "diagnostics": [], "p_value": None}
    fixed, residuals, tests, rows = [], [], [], []
    if q["design"] == "two_way_rm_anova":
        arm_index = unit_arm.loc[matrix.index].map({a: i for i, a in enumerate(arms)}).to_numpy()
        tests = tw.split_plot_anova(matrix.to_numpy(float), arm_index, len(arms))
        rows = tw.anova_contrasts(matrix.to_numpy(float), arm_index, q, arms, conds)
        fit_row.update(epsilon_gg=tests[1]["epsilon_gg"],
                       contrast_method="welch_between_arms_paired_within_arm_holm_p_bonferroni_ci")
    else:
        x = tw.cell_design(used.arm_index.to_numpy(), used.condition_index.to_numpy(), len(arms), len(conds))
        groups = used.independent_unit_id.to_numpy(str)
        y = used.value.to_numpy(float)
        fit = fit_random_intercept(y, x, groups)
        fit_row["optimizer_attempts"] = fit["attempts"]
        if not fit["ok"]:
            fit_row.update(status="failed", reportable=False, diagnostics=["mixed_model_did_not_converge_or_invalid_covariance"])
        else:
            fit_row.update(random_intercept_variance=fit["tau"], residual_variance=fit["sigma"],
                           icc=fit["tau"] / (fit["tau"] + fit["sigma"]), reml_log_likelihood=fit["reml_log_likelihood"],
                           contrast_method="satterthwaite_t_holm_p_bonferroni_ci" if q["inference"] == "satterthwaite" else "wald_z_holm_p_bonferroni_ci")
            if fit["boundary"]:
                fit_row["diagnostics"].append("random_intercept_variance_near_zero_boundary")
            if len(matrix) < 30:
                fit_row["diagnostics"].append("fewer_than_30_units_small_sample_inference_requires_caution")
            tests, rows, diagnostics = tw.mixed_tests(y, x, groups, fit, q, arms, conds)
            fit_row["diagnostics"] += diagnostics
            if tests is None:
                fit_row.update(status="limited", reportable=False)
                tests, rows = [], []
            fixed = [{"arm": arms[i], "group": conds[j], "estimated_mean": float(fit["beta"][i * len(conds) + j]),
                      "standard_error": float(np.sqrt(fit["covariance"][i * len(conds) + j, i * len(conds) + j]))}
                     for i, j in tw.cells(len(arms), len(conds))]
            residuals = [{"observation_id": o.observation_id, "independent_unit_id": o.independent_unit_id, "arm": o.arm,
                          "group": o.group, "marginal_fitted": float(m), "conditional_fitted": float(c),
                          "random_intercept_blup": float(b), "marginal_residual": float(yy - m), "conditional_residual": float(yy - c)}
                         for o, yy, m, c, b in zip(used.itertuples(), y, fit["marginal"], fit["conditional"], fit["random_intercept"])]
    if matrix.isna().any().any():
        fit_row["diagnostics"].append("missing_measurements_available_case_MAR_assumption_not_testable")
    interaction = next((t for t in tests if t["effect"] == "arm_x_condition"), None)
    if interaction is not None:
        fit_row["p_value"] = interaction["p_value"]
        fit_row["primary_test"] = "arm_x_condition"
    return {"schema_version": 1, "analysis_type": "repeated_measures", "fits": [fit_row], "tests": tests,
            "group_summaries": summaries, "contrasts": rows, "fixed_effects": fixed, "residuals": residuals,
            "missingness": missing, "bootstrap_samples": []}
