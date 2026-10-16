"""Independent-unit multivariable Gaussian, binomial-logit and Poisson-log models.

Explicit design matrices (no formula evaluation), model-based or declared HC3/
Pearson-scale inference, and conservative reportability diagnostics. All fit
and prediction numbers are saved here; reporting and facts never refit.
"""
from copy import deepcopy
from itertools import product
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
from scipy import optimize, special, stats
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan

TYPES = {"generalized_regression"}
DEFAULTS = {
    "schema_version": 1, "analysis_type": "generalized_regression", "input": None, "source": None,
    "design": {"unit": None, "rationale": None, "outcome": None, "outcome_unit": None, "estimand": None},
    "model": {"family": None, "predictors": None, "interactions": [], "rationale": None,
              "covariance": "model_based", "trials": None, "exposure": None, "offset": None,
              "offset_rationale": None, "confidence_level": .95},
    "missing": {"policy": "error", "rationale": None},
    "prediction": {"rows": []}, "report": {"plot_style": "prism_like"}}
RESERVED = {"unit_id", "outcome", "exclude", "exclusion_reason", "prediction_id", "Intercept"}


def _text(v, name):
    if not isinstance(v, str) or not v.strip(): raise ValueError("Declare " + name)


def _column(v):
    return isinstance(v, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", v) is not None and v not in RESERVED


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw)-set(DEFAULTS): raise ValueError("Unsupported regression config keys")
    c = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(c[key], dict):
            if not isinstance(value, dict) or set(value)-set(c[key]): raise ValueError("Unsupported " + key + " keys")
            c[key].update(deepcopy(value))
        else: c[key] = deepcopy(value)
    if c["schema_version"] != 1 or c["analysis_type"] not in TYPES: raise ValueError("Unsupported regression schema")
    for key, value in [("input", c["input"]), ("source", c["source"]), *c["design"].items()]: _text(value, key)
    m = c["model"]; _text(m["rationale"], "model rationale")
    if m["family"] not in ("gaussian", "binomial", "poisson"): raise ValueError("Family must be gaussian, binomial or poisson")
    allowed = ("model_based", "HC3") if m["family"] == "gaussian" else ("model_based", "pearson") if m["family"] == "poisson" else ("model_based",)
    if m["covariance"] not in allowed: raise ValueError("Covariance unsupported for the declared family")
    pr = m["predictors"]
    if not isinstance(pr, dict) or not pr: raise ValueError("Declare named predictors")
    for name, spec in pr.items():
        if not _column(name) or not isinstance(spec, dict): raise ValueError("Invalid predictor declaration")
        if spec.get("kind") == "continuous":
            if set(spec) != {"kind", "unit"}: raise ValueError("Continuous predictors require kind and unit only")
            _text(spec["unit"], name + " unit")
        elif spec.get("kind") == "categorical":
            if set(spec) != {"kind", "levels", "reference"}: raise ValueError("Declare categorical levels and reference")
            levels = spec["levels"]
            if not isinstance(levels, list) or len(levels) < 2 or any(not isinstance(v, str) or not v.strip() for v in levels) or len(set(levels)) != len(levels) or spec["reference"] not in levels:
                raise ValueError("Categorical levels must be distinct strings with a declared reference")
        else: raise ValueError("Predictor kind must be continuous or categorical")
    interactions = m["interactions"]
    if not isinstance(interactions, list): raise ValueError("Interactions must be declared pairs")
    seen = set()
    for pair in interactions:
        if not isinstance(pair, list) or len(pair) != 2 or any(not isinstance(v, str) or v not in pr for v in pair) or pair[0] == pair[1]: raise ValueError("Interactions require two distinct declared main-effect predictors")
        key = tuple(sorted(pair))
        if key in seen: raise ValueError("Duplicate interaction")
        seen.add(key)
    if type(m["confidence_level"]) not in (int, float) or not .5 < m["confidence_level"] < 1: raise ValueError("Invalid confidence level")
    extra = [m[k] for k in ("trials", "exposure", "offset") if m[k] is not None]
    if any(not _column(v) or v in pr for v in extra) or len(extra) != len(set(extra)): raise ValueError("Trials/offset/exposure columns must be distinct from predictors")
    if m["trials"] is not None and m["family"] != "binomial": raise ValueError("Trials are binomial only")
    if m["offset"] is not None or m["exposure"] is not None:
        if m["family"] != "poisson" or (m["offset"] is not None and m["exposure"] is not None): raise ValueError("Declare either Poisson log offset or exposure")
        _text(m["offset_rationale"], "offset/exposure rationale and source")
    elif m["offset_rationale"] is not None: raise ValueError("Offset rationale requires an offset/exposure")
    if c["missing"]["policy"] not in ("error", "complete_case"): raise ValueError("Invalid missing policy")
    if c["missing"]["policy"] == "complete_case": _text(c["missing"]["rationale"], "complete-case rationale")
    rows = c["prediction"]["rows"]
    if not isinstance(rows, list): raise ValueError("Prediction rows must be a list")
    ids = []
    for row in rows:
        required = {"prediction_id", *pr, *[m[k] for k in ("exposure", "offset") if m[k] is not None]}
        if not isinstance(row, dict) or set(row) != required: raise ValueError("Prediction rows require ID and exactly the declared predictors/offset/exposure")
        _text(row["prediction_id"], "prediction_id"); ids.append(row["prediction_id"])
    if len(ids) != len(set(ids)): raise ValueError("Duplicate prediction IDs")
    if c["report"]["plot_style"] not in ("prism_like", "standard"): raise ValueError("Invalid plot style")
    return c


def load_data(path, cfg):
    # Keep raw qualifiers, category labels (including 'NA') and leading-zero IDs.
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _numeric(frame, name):
    text = frame[name].astype(str).str.strip()
    missing = frame[name].isna() | text.isin(("", "NA", "NaN", "nan"))
    values = pd.to_numeric(text.mask(missing), errors="coerce")
    if ((~missing) & (~np.isfinite(values))).any():
        raise ValueError(f"{name}: nonfinite value or nonnumeric qualifier; no zero/pseudocount substitution")
    return values


def _prepare(data, c):
    m = c["model"]; d = data.copy()
    numeric = ["outcome"] + [k for k, v in m["predictors"].items() if v["kind"] == "continuous"]
    numeric += [m[k] for k in ("trials", "exposure", "offset") if m[k] is not None]
    required = {"unit_id", *numeric, *m["predictors"]}
    if d.empty or required-set(d): raise ValueError("Missing regression columns: " + str(sorted(required-set(d))))
    d["unit_id"] = d.unit_id.astype(str)
    if d.unit_id.str.strip().eq("").any() or d.unit_id.duplicated().any(): raise ValueError("One nonempty unique ID per independent unit; repeated or technical rows need another model")
    excluded = d.get("exclude", pd.Series(False, index=d.index)).astype(str).str.lower()
    reason = d.get("exclusion_reason", pd.Series("", index=d.index)).fillna("").astype(str)
    if not excluded.isin(("true", "false")).all() or (excluded.eq("true") & reason.str.strip().eq("")).any(): raise ValueError("Explicit exclusions require true/false and reasons")
    log = {"input_rows": len(d), "excluded": [{"unit_id": str(d.loc[i, "unit_id"]), "reason": reason.loc[i]} for i in d.index[excluded.eq("true")]],
           "missing_policy": c["missing"], "missing_rows": [], "transformations": ["Declared reference-category indicator coding and interactions; intercept included"]}
    u = d.loc[excluded.eq("false")].copy()
    for name in numeric: u[name] = _numeric(u, name)
    missing = u[numeric].isna().any(axis=1)
    for name, spec in m["predictors"].items():
        if spec["kind"] != "categorical": continue
        col = u[name].fillna("").astype(str)
        blank = col.str.strip().eq("")
        if (~blank & ~col.isin(spec["levels"])).any(): raise ValueError("Undeclared category in " + name)
        u[name] = col; missing |= blank
    log["missing_rows"] = u.loc[missing, "unit_id"].tolist()
    if missing.any() and c["missing"]["policy"] == "error": raise ValueError("Missing model data; declare a justified complete_case policy instead of silent deletion")
    u = u.loc[~missing].reset_index(drop=True)
    if u.empty: raise ValueError("No complete included units")
    for name, spec in m["predictors"].items():
        if spec["kind"] == "categorical" and set(u[name]) != set(spec["levels"]): raise ValueError("A declared category is absent after exclusions/missingness: " + name)
    y = u.outcome.to_numpy(float)
    if m["family"] in ("binomial", "poisson") and ((y < 0).any() or (y != np.floor(y)).any()): raise ValueError("Counts/successes must be nonnegative integers, not normalized percentages")
    trials = u[m["trials"]].to_numpy(float) if m["trials"] else np.ones(len(u))
    if m["family"] == "binomial" and ((trials <= 0).any() or (trials != np.floor(trials)).any() or (y > trials).any()): raise ValueError("Binary outcome must be 0/1, or successes between zero and positive integer trials")
    offset = _offset(u, m)
    log["used_unit_ids"] = u.unit_id.tolist(); log["used_rows"] = len(u)
    if m["exposure"]: log["transformations"].append("Natural log of declared exposure used as fixed offset, coefficient exactly one")
    return u, y, trials, offset, log


def _offset(d, m):
    if m["exposure"]:
        v = np.asarray(d[m["exposure"]], float)
        if not np.isfinite(v).all() or (v <= 0).any(): raise ValueError("Exposure must be finite and strictly positive")
        return np.log(v)
    if m["offset"]:
        v = np.asarray(d[m["offset"]], float)
        if not np.isfinite(v).all(): raise ValueError("Log offset must be finite")
        return v
    return np.zeros(len(d))


def design_matrix(d, m):
    columns, names, groups, blocks = [np.ones(len(d))], ["Intercept"], {}, {}
    for name, spec in m["predictors"].items():
        if spec["kind"] == "continuous":
            a = np.asarray(d[name], float)
            if not np.isfinite(a).all(): raise ValueError("Predictor must be finite: " + name)
            block = [(name, a)]
        else:
            v = d[name].astype(str)
            if not v.isin(spec["levels"]).all(): raise ValueError("Unknown prediction category: " + name)
            block = [(name + "[" + level + "]", (v == level).to_numpy(float)) for level in spec["levels"] if level != spec["reference"]]
        blocks[name] = block; groups[name] = []
        for label, values in block:
            groups[name].append(len(names)); columns.append(values); names.append(label)
    for pair in m["interactions"]:
        term = ":".join(pair); groups[term] = []
        for a, b in product(blocks[pair[0]], blocks[pair[1]]):
            names.append(a[0] + ":" + b[0]); groups[term].append(len(columns)); columns.append(a[1]*b[1])
    if len(names) != len(set(names)): raise ValueError("Encoded coefficient names collide; rename categorical levels")
    X = np.column_stack(columns)
    if not np.isfinite(X).all(): raise ValueError("Design matrix overflow")
    return X, names, groups


def _separated(X, y, trials, family):
    """Find any likelihood recession direction, including quasi-separation.

    Bounded coefficients fix the arbitrary scale. Boundary binomial rows allow
    a signed nonnegative direction, interior rows require zero. Poisson positive
    rows require zero and zero rows allow negative directions. Maximizing the
    sum detects at least one strict inequality, not only complete separation.
    """
    if family == "gaussian": return False
    if family == "binomial":
        edge = (y == 0) | (y == trials)
        A = X[edge] * np.where(y[edge] == 0, -1, 1)[:, None]
        E = X[~edge]
    else:
        A = -X[y == 0]; E = X[y > 0]
    if len(A) == 0: return False
    fit = optimize.linprog(-A.sum(axis=0), A_ub=-A, b_ub=np.zeros(len(A)),
                           A_eq=E if len(E) else None, b_eq=np.zeros(len(E)) if len(E) else None,
                           bounds=[(-1, 1)]*X.shape[1], method="highs")
    if not fit.success: raise ValueError("Separation feasibility check failed; cannot establish finite MLE")
    return bool(-fit.fun > 1e-7)


def _finite(value):
    if isinstance(value, dict): return {k: _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_finite(v) for v in value]
    if isinstance(value, np.ndarray): return _finite(value.tolist())
    if isinstance(value, (float, np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer): return int(value)
    if isinstance(value, np.bool_): return bool(value)
    return value


def compute(data, c):
    m = c["model"]; family = m["family"]; covariance = m["covariance"]
    u, y, trials, offset, log = _prepare(data, c)
    X, names, groups = design_matrix(u, m); n, p = X.shape
    # Scaling changes neither coefficients nor predictions after back-transform.
    scales = np.linalg.norm(X, axis=0)/np.sqrt(n); scales[scales == 0] = 1
    Z = X/scales; rank = np.linalg.matrix_rank(Z); condition = float(np.linalg.cond(Z)); df = n-p
    flags, blockers = [], []
    notes = ["Independent unit: " + c["design"]["unit"] + ". Conditional associations do not establish causation or assay acceptance.",
             "Coefficients and mean-response intervals are pointwise, not a simultaneous family; no automatic variable/model selection or exclusions.",
             "Interactions make main effects conditional on other predictors at zero/reference; grouped Wald tests concern the declared coded coefficients, not Type II/III ANOVA."]
    if log["missing_rows"]: notes.append("Complete-case analysis removed the listed units; the table cannot verify its missingness assumptions and selection may bias results.")
    if family == "binomial": notes.append("The positive event is outcome=1 (or successes/trials). Exponentiated slopes are conditional odds ratios, not risk ratios; grouped binomial variance assumes the declared trials follow the binomial model.")
    if family == "poisson": notes.append("Exponentiated slopes are conditional expected-count ratios, or rate ratios when a declared exposure/log offset is used. They are not automatically risk ratios.")
    if family != "gaussian": notes.append("GLM Wald intervals are asymptotic. Sparse outcomes and small independent sample counts can yield poor coverage even without detected separation.")
    if covariance == "HC3": notes.append("HC3 with a residual-df t reference is approximate; it does not correct dependence, confounding, omitted predictors or a wrong mean model.")
    if covariance == "pearson": notes.append("Quasi-Poisson uses Pearson chi-square/residual df as the variance scale and a t reference; proportional variance is assumed. No likelihood/AIC or negative-binomial fit is claimed.")
    r = {"analysis_type": "generalized_regression", "primary": {"status": "withheld", "reportable": False,
         "family": family, "link": "identity" if family == "gaussian" else "logit" if family == "binomial" else "log",
         "covariance": covariance, "n_units": n, "coefficient_count": p, "df_residual": df,
         "confidence_level": m["confidence_level"], "estimand": c["design"]["estimand"],
         "coefficients": [], "term_tests": [], "predictions": []},
         "design": {"columns": names, "matrix": X.tolist(), "unit_ids": u.unit_id.tolist(),
                    "predictors": m["predictors"], "interactions": m["interactions"], "outcome": c["design"]["outcome"], "outcome_unit": c["design"]["outcome_unit"]},
         "preprocessing": log, "diagnostics": {"rank": int(rank), "scaled_condition_number": condition, "flags": flags},
         "must_mention": notes, "failing_items": [], "audit": {},
         "limitations": ["No causal, predictive-validation, simultaneous-coverage or Prism-equivalence claim.",
                         "No repeated/clustered outcomes, GLMM/GEE, penalization/Firth, imputation, weights, negative-binomial/Gamma models, or automatic model selection.",
                         "Mean prediction intervals are not individual prediction intervals. Support checks do not prove transportability.",
                         "Diagnostic gates and conditional coverage are design-dependent; passing checks does not prove assumptions."]}

    def finish():
        flags.extend(b for b in blockers if b not in flags)
        r["failing_items"] = [{"reason": b} for b in blockers]
        if blockers: notes.append("Reportable inference withheld: " + ", ".join(blockers) + ". Saved audit numbers must not be promoted.")
        return _finite(r)

    if rank < p: blockers.append("rank_deficient_design")
    if df < 2: blockers.append("insufficient_residual_degrees_of_freedom")
    if not np.isfinite(condition) or condition > 1e10: blockers.append("ill_conditioned_design")
    if blockers: return finish()
    if _separated(Z, y, trials, family):
        blockers.append("separation_or_nonfinite_mle"); return finish()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            if family == "gaussian":
                fit = sm.OLS(y, Z).fit()
                covz = fit.get_robustcov_results(cov_type="HC3").cov_params() if covariance == "HC3" else fit.cov_params()
                converged = True; phi = float(fit.scale)
            else:
                endog = np.column_stack([y, trials-y]) if family == "binomial" else y
                fam = sm.families.Binomial() if family == "binomial" else sm.families.Poisson()
                fit = sm.GLM(endog, Z, family=fam, offset=offset).fit(maxiter=200, tol=1e-10)
                converged = fit.converged; phi = float(fit.pearson_chi2/df) if covariance == "pearson" else 1.
                covz = fit.cov_params()*phi
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
            r["diagnostics"]["fit_error"] = type(e).__name__ + ": " + str(e)
            blockers.append("fit_failed"); return finish()
    r["diagnostics"]["fit_warnings"] = [str(w.message) for w in caught]
    beta = np.asarray(fit.params)/scales; cov = np.asarray(covz)/np.outer(scales, scales)
    if not converged: blockers.append("fit_not_converged")
    if not np.isfinite(beta).all() or not np.isfinite(cov).all() or (np.diag(cov) <= 0).any() or not np.isfinite(phi) or phi <= 0:
        blockers.append("invalid_parameter_uncertainty"); return finish()
    mu = np.asarray(fit.fittedvalues); residual = y/trials-mu if family == "binomial" else y-mu
    if family == "gaussian":
        pearson = residual/np.sqrt(phi)
        hat = fit.get_influence().hat_matrix_diag
        bp = float(het_breuschpagan(residual, Z, robust=True)[1])
        r["diagnostics"]["breusch_pagan_p"] = bp
        if bp < .01:
            flags.append("heteroscedasticity_signal")
            notes.append("Breusch-Pagan heteroscedasticity diagnostic p < .01; this approximate check does not identify the variance mechanism.")
            if covariance == "model_based": blockers.append("classical_variance_assumption_flagged")
    else:
        pearson = np.asarray(fit.resid_pearson)/np.sqrt(phi)
        hat = fit.get_hat_matrix_diag()
        dispersion = float(fit.pearson_chi2/df)
        r["diagnostics"].update(pearson_chi_square=float(fit.pearson_chi2), dispersion_ratio=dispersion)
        if family == "poisson" or m["trials"]:
            dp = float(stats.chi2.sf(fit.pearson_chi2, df))
            r["diagnostics"]["dispersion_upper_tail_p"] = dp
            notes.append("Pearson overdispersion diagnostic is an approximate upper-tail chi-square check; sparse fitted counts can make it unreliable.")
            if dp < .01:
                flags.append("overdispersion_signal")
                notes.append("Overdispersion signal: fixed Poisson/binomial variance may understate uncertainty; no automatic switch of distribution or deletion of observations.")
                if covariance == "model_based": blockers.append("fixed_dispersion_assumption_flagged")
        if family == "binomial":
            r["diagnostics"].update(successes=float(y.sum()), failures=float((trials-y).sum()))
            if min(y.sum(), (trials-y).sum()) < 10*p:
                flags.append("sparse_event_information"); notes.append("Fewer than ten successes or failures per coefficient: a heuristic sparse-information flag, not an adequacy threshold; Wald coverage may be poor.")
    hsafe = np.maximum(1-hat, np.finfo(float).eps)
    cooks = pearson**2*hat/(p*hsafe**2)
    influential = (hat > 2*p/n) | (cooks > 4/n)
    if influential.any():
        flags.append("influence_review")
        notes.append("Listed high leverage (>2p/n) or approximate Cook distance (>4/n) values are review heuristics, not exclusion rules.")
    if (hat >= 1-1e-10).any(): blockers.append("unit_leverage_one")
    r["diagnostics"].update(converged=bool(converged), scale=phi,
        influence_unit_ids=u.loc[influential, "unit_id"].tolist(),
        observations=[{"unit_id": uid, "observed": float(obs), "trials": float(nt) if family == "binomial" else None,
                       "fitted_mean": float(pred), "response_residual": float(res), "scaled_pearson_residual": float(pe),
                       "leverage": float(h), "approx_cooks_distance": float(cd)}
                      for uid, obs, nt, pred, res, pe, h, cd in zip(u.unit_id, y, trials, mu, residual, pearson, hat, cooks)])
    use_t = family == "gaussian" or covariance == "pearson"
    distribution = stats.t(df) if use_t else stats.norm()
    critical = distribution.ppf((1+m["confidence_level"])/2)
    rows = []
    for i, name in enumerate(names):
        se = float(np.sqrt(cov[i, i])); stat = float(beta[i]/se)
        ci = [float(beta[i]-critical*se), float(beta[i]+critical*se)]
        row = {"term": name, "estimate": float(beta[i]), "standard_error": se, "ci": ci,
               "statistic": stat, "reference": "t" if use_t else "normal", "df": df if use_t else None,
               "p_value": float(2*distribution.sf(abs(stat)))}
        if family != "gaussian":
            with np.errstate(over="ignore", under="ignore"): exponentiated = np.exp([beta[i], *ci])
            row.update(exponentiated_estimate=float(exponentiated[0]), exponentiated_ci=exponentiated[1:].tolist())
            if not np.isfinite(exponentiated).all() or (exponentiated == 0).any(): blockers.append("unbounded_exponentiated_uncertainty")
        rows.append(row)
    terms = []
    for term, indices in groups.items():
        b = beta[indices]; v = cov[np.ix_(indices, indices)]
        try: wald = float(b @ np.linalg.solve(v, b))
        except np.linalg.LinAlgError:
            blockers.append("singular_term_covariance"); continue
        k = len(indices)
        terms.append({"term": term, "coefficients": [names[i] for i in indices], "statistic": wald/k if use_t else wald,
                      "reference": "F" if use_t else "chi_square", "df_numerator": k,
                      "df_denominator": df if use_t else None,
                      "p_value": float(stats.f.sf(wald/k, k, df) if use_t else stats.chi2.sf(wald, k))})
    predictions = []
    if c["prediction"]["rows"]:
        new = pd.DataFrame(c["prediction"]["rows"])
        Xnew, _, _ = design_matrix(new, m); offsets = _offset(new, m)
        for row, x, off in zip(c["prediction"]["rows"], Xnew, offsets):
            support = np.vstack([Z.T, np.ones(n)])
            target = np.r_[x/scales, 1.]
            if m["exposure"] or m["offset"]:
                offset_scale = max(1., float(np.max(np.abs(offset))))
                support = np.vstack([support, offset/offset_scale]); target = np.r_[target, off/offset_scale]
            feasible = optimize.linprog(np.zeros(n), A_eq=support,
                b_eq=target, bounds=(0, None), method="highs")
            pr = {"prediction_id": row["prediction_id"], "input": row,
                  "reportable": bool(feasible.success) and not blockers, "interval_target": "conditional_mean_response"}
            if not feasible.success:
                pr.update(status="withheld", reason="outside_observed_design_convex_hull")
            elif blockers:
                pr.update(status="withheld", reason="parent_model_withheld")
            else:
                eta = float(x @ beta+off); se = float(np.sqrt(max(0., x @ cov @ x)))
                vals = np.array([eta, eta-critical*se, eta+critical*se])
                with np.errstate(over="ignore", under="ignore"):
                    vals = vals if family == "gaussian" else special.expit(vals) if family == "binomial" else np.exp(vals)
                if not np.isfinite(vals).all() or (family == "poisson" and (vals <= 0).any()):
                    pr.update(status="withheld", reportable=False, reason="prediction_numeric_range")
                else: pr.update(status="estimated", mean=float(vals[0]), ci=vals[1:].tolist(), link_standard_error=se)
            predictions.append(pr)
    if any(not q["reportable"] for q in predictions):
        notes.append("Some requested mean predictions are withheld; retain each row's stated support or parent-model reason.")
    fit_summary = {"scale": phi, "residual_df": df,
                   "r_squared": float(fit.rsquared) if family == "gaussian" else None,
                   "deviance": float(fit.deviance) if family != "gaussian" else float(fit.ssr),
                   "log_likelihood": float(fit.llf) if covariance != "pearson" else None,
                   "aic": float(fit.aic) if covariance != "pearson" else None}
    if blockers:
        r["audit"] = {"coefficients": rows, "term_tests": terms, "covariance_matrix": cov.tolist(), "fit_summary": fit_summary}
    else:
        r["primary"].update(status="estimated", reportable=True, coefficients=rows, term_tests=terms, fit_summary=fit_summary)
        r["covariance_matrix"] = cov.tolist()
    r["primary"]["predictions"] = predictions
    return finish()


def write_artifacts(out, result):
    from .workflow import dump
    out = Path(out)
    dump(out/"diagnostics.json", result["diagnostics"])
    dump(out/"preprocessing_log.json", result["preprocessing"])
    design = result["design"]
    matrix = pd.DataFrame(design["matrix"], columns=design["columns"])
    matrix.insert(0, "unit_id", design["unit_ids"]); matrix.to_csv(out/"design_matrix.csv", index=False)
    coefficients = result["primary"]["coefficients"]
    pd.DataFrame(coefficients, columns=["term", "estimate", "standard_error", "ci", "statistic", "reference", "df", "p_value", "exponentiated_estimate", "exponentiated_ci"]).to_csv(out/"coefficients.csv", index=False)
    pd.DataFrame(result["primary"]["predictions"], columns=["prediction_id", "status", "reportable", "mean", "ci", "reason", "interval_target"]).to_csv(out/"predictions.csv", index=False)


def render(run, result, style):
    """Plot saved coefficients and diagnostic observations; no fitting/inference."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/"figures"; out.mkdir(exist_ok=True)
    names = []; primary = result["primary"]
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        if primary["coefficients"]:
            rows = primary["coefficients"][1:]
            fig, ax = plt.subplots(figsize=(7, max(3, .3*len(rows)+1)), layout="constrained")
            estimates = np.array([x["estimate"] for x in rows]); ci = np.array([x["ci"] for x in rows])
            ax.errorbar(estimates, np.arange(len(rows)), xerr=np.array([estimates-ci[:, 0], ci[:, 1]-estimates]), fmt="o", color=pstyle.color(token, 0))
            ax.axvline(0, color=pstyle.MUTED, ls="--", lw=1)
            ax.set(yticks=np.arange(len(rows)), yticklabels=[x["term"] for x in rows], xlabel="Coefficient on linear-predictor scale (pointwise CI)", title=primary["family"] + " / " + primary["covariance"])
            pstyle.save(fig, out/"regression_coefficients", token, 300); names.append("regression_coefficients")
        obs = result["diagnostics"].get("observations", [])
        if obs:
            fig, ax = plt.subplots(figsize=(6, 3.5), layout="constrained")
            ax.scatter([x["fitted_mean"] for x in obs], [x["scaled_pearson_residual"] for x in obs], s=14, color=pstyle.color(token, 0), alpha=.7)
            ax.axhline(0, color=pstyle.MUTED, ls="--", lw=1)
            ax.set(xlabel="Saved fitted mean", ylabel="Scaled Pearson residual", title="Diagnostic only" + (" — parent inference withheld" if not primary["reportable"] else ""))
            pstyle.save(fig, out/"regression_residuals", token, 300); names.append("regression_residuals")
    return names
