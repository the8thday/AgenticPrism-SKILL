"""Between-arm by within-condition designs (for example treatment arm x time).

two_way_rm_anova: complete split-plot ANOVA. The arm effect is a one-way ANOVA on
unit means. Condition and arm x condition use the multivariate linear model on
orthonormal within-unit contrasts Z = Y C with type III (equal arm weight)
hypotheses; univariate F statistics use tr(H)/tr(E), with Greenhouse-Geisser
epsilon from the pooled within-arm covariance of Z. This matches car::Anova
type III univariate tests with GG correction (afex default).

two_way_mixed: cell means for every arm x condition plus a unit random
intercept, REML. Type III tests are equal-weight cell-mean hypotheses with
Satterthwaite F (lmerTest-style) or asymptotic Wald chi-square.
"""
from itertools import product
import numpy as np
from scipy.stats import chi2, f as f_dist, norm, t as t_dist, ttest_ind, ttest_rel
from .mixed_inference import RandomInterceptModel, satterthwaite

EFFECTS = ("arm", "condition", "arm_x_condition")


def _holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    return adjusted


def _orthonormal_contrasts(k):
    """k x (k-1) columns orthonormal and orthogonal to the unit vector."""
    basis = np.linalg.qr(np.column_stack([np.ones(k), np.eye(k)[:, :k - 1]]))[0]
    return basis[:, 1:]


def split_plot_anova(y, arm_index, n_arms):
    """y: units x conditions (complete). Returns the three type III tests."""
    y = np.asarray(y, float)
    n, k = y.shape
    d = np.eye(n_arms)[arm_index]
    counts = d.sum(0)
    if n - n_arms < 1 or counts.min() < 2:
        raise ValueError("Split-plot ANOVA needs at least two units per arm and more units than arms")
    tests = []
    # Arm: one-way ANOVA on unit means (the k factor cancels in F).
    m = y.mean(1)
    arm_means = (d.T @ m) / counts
    ss_a = k * float(np.sum(counts * (arm_means - m.mean()) ** 2))
    ss_s = k * float(np.sum((m - arm_means[arm_index]) ** 2))
    df_a, df_s = n_arms - 1, n - n_arms
    if ss_s <= 0:
        raise ValueError("Zero between-unit variance within arms")
    f_a = (ss_a / df_a) / (ss_s / df_s)
    tests.append({"effect": "arm", "f_statistic": f_a, "df_numerator": float(df_a), "df_denominator": float(df_s),
                  "p_value": float(f_dist.sf(f_a, df_a, df_s)), "epsilon_gg": None, "p_uncorrected": None,
                  "partial_eta_squared": ss_a / (ss_a + ss_s)})
    c = _orthonormal_contrasts(k)
    z = y @ c
    gamma = np.linalg.solve(d.T @ d, d.T @ z)
    resid = z - d @ gamma
    e = resid.T @ resid
    ss_e = float(np.trace(e))
    if ss_e <= 0:
        raise ValueError("Zero within-unit residual variance")
    df_e = (n - n_arms) * (k - 1)
    s = e / (n - n_arms)
    eps = float(np.clip(np.trace(s) ** 2 / ((k - 1) * np.sum(s * s)), 1 / (k - 1), 1))
    inv_dd = np.linalg.inv(d.T @ d)
    hypotheses = {"condition": np.full((1, n_arms), 1 / n_arms),
                  "arm_x_condition": np.column_stack([np.eye(n_arms - 1), -np.ones(n_arms - 1)])}
    for effect, lm in hypotheses.items():
        lg = lm @ gamma
        h = lg.T @ np.linalg.solve(lm @ inv_dd @ lm.T, lg)
        ss_h = float(np.trace(h))
        df_h = lm.shape[0] * (k - 1)
        fstat = (ss_h / df_h) / (ss_e / df_e)
        tests.append({"effect": effect, "f_statistic": fstat, "df_numerator": eps * df_h, "df_denominator": eps * df_e,
                      "p_value": float(f_dist.sf(fstat, eps * df_h, eps * df_e)), "epsilon_gg": eps,
                      "df_numerator_uncorrected": float(df_h), "df_denominator_uncorrected": float(df_e),
                      "p_uncorrected": float(f_dist.sf(fstat, df_h, df_e)), "partial_eta_squared": ss_h / (ss_h + ss_e)})
    return tests


def family(q, arms, conditions):
    """Declared contrast family as (label dict, arm index pair, condition index pair)."""
    kind = q["post_hoc"]
    if kind == "none":
        return []
    if kind == "arm_vs_control_each_condition":
        c = arms.index(q["control_arm"])
        return [({"condition": conditions[j], "group_a": arms[c], "group_b": arms[i]}, (c, j), (i, j))
                for j in range(len(conditions)) for i in range(len(arms)) if i != c]
    return [({"arm": arms[i], "group_a": conditions[0], "group_b": conditions[j]}, (i, 0), (i, j))
            for i in range(len(arms)) for j in range(1, len(conditions))]


def anova_contrasts(y, arm_index, q, arms, conditions):
    """Assumption-light pairwise tests: Welch between arms at a condition, paired t within an arm."""
    rows = []
    fam = family(q, arms, conditions)
    for label, (ia, ja), (ib, jb) in fam:
        if ja == jb:
            a, b = y[arm_index == ia, ja], y[arm_index == ib, jb]
            sa, sb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
            se = float(np.sqrt(sa + sb))
            df = float((sa + sb) ** 2 / (sa ** 2 / (len(a) - 1) + sb ** 2 / (len(b) - 1))) if se > 0 else None
            test = ttest_ind(b, a, equal_var=False) if se > 0 else None
            diff, method = float(b.mean() - a.mean()), "welch"
        else:
            delta = y[arm_index == ia, jb] - y[arm_index == ia, ja]
            se = float(delta.std(ddof=1) / np.sqrt(len(delta)))
            df = float(len(delta) - 1)
            test = ttest_rel(y[arm_index == ia, jb], y[arm_index == ia, ja]) if se > 0 else None
            diff, method = float(delta.mean()), "paired_t"
        row = {**label, "difference_b_minus_a": diff, "standard_error": se, "df": df, "method": method,
               "reportable": test is not None}
        if test is None:
            row.update(statistic=None, p_unadjusted=None, ci_low=None, ci_high=None)
        else:
            margin = float(t_dist.ppf(1 - (1 - q["confidence_level"]) / (2 * len(fam)), df) * se)
            row.update(statistic=float(test.statistic), p_unadjusted=float(test.pvalue), ci_low=diff - margin, ci_high=diff + margin)
        rows.append(row)
    adjusted = _holm([r["p_unadjusted"] if r["p_unadjusted"] is not None else 1 for r in rows]) if rows else []
    for r, p in zip(rows, adjusted):
        r["p_adjusted"] = float(p) if r["reportable"] else None
    return rows


def cell_design(arm_index, cond_index, n_arms, k):
    x = np.zeros((len(arm_index), n_arms * k))
    x[np.arange(len(arm_index)), np.asarray(arm_index) * k + np.asarray(cond_index)] = 1
    return x


def type3_matrices(n_arms, k):
    ca = np.column_stack([-np.ones(n_arms - 1), np.eye(n_arms - 1)])
    ck = np.column_stack([-np.ones(k - 1), np.eye(k - 1)])
    return {"arm": np.kron(ca, np.full((1, k), 1 / k)),
            "condition": np.kron(np.full((1, n_arms), 1 / n_arms), ck),
            "arm_x_condition": np.kron(ca, ck)}


def mixed_tests(y, x, groups, fit, q, arms, conditions):
    """Type III tests and the declared family for the cell-means random-intercept model."""
    n_arms, k = len(arms), len(conditions)
    beta, cov = fit["beta"], fit["covariance"]
    fam = family(q, arms, conditions)
    vectors = []
    for _, (ia, ja), (ib, jb) in fam:
        v = np.zeros(n_arms * k)
        v[ib * k + jb] += 1
        v[ia * k + ja] -= 1
        vectors.append(v)
    labels = [label for label, _, _ in fam]
    tests, rows, diagnostics = [], [], []
    mats = type3_matrices(n_arms, k)
    if q["inference"] == "satterthwaite":
        if fit["boundary"]:
            diagnostics.append("satterthwaite_unavailable_at_variance_boundary_or_singular_hessian")
            return None, None, diagnostics
        model = RandomInterceptModel(y, x, groups)
        try:
            srows, _, _ = satterthwaite(model, [fit["tau"], fit["sigma"]], beta, vectors, level=q["confidence_level"])
            for effect in EFFECTS:
                _, om, _ = satterthwaite(model, [fit["tau"], fit["sigma"]], beta, [], joint=mats[effect])
                tests.append({"effect": effect, "f_statistic": om["f_statistic"], "df_numerator": float(om["df_numerator"]),
                              "df_denominator": om["df_denominator"], "p_value": om["p_value"]})
        except (ArithmeticError, np.linalg.LinAlgError):
            diagnostics.append("satterthwaite_unavailable_at_variance_boundary_or_singular_hessian")
            return None, None, diagnostics
        for label, sr in zip(labels, srows):
            rows.append({**label, "difference_b_minus_a": sr["estimate"], "standard_error": sr["standard_error"],
                         "df": sr["df"], "statistic": sr["statistic"], "p_unadjusted": sr["p_unadjusted"],
                         "ci_low": sr["ci_low"], "ci_high": sr["ci_high"], "method": "satterthwaite_t", "reportable": True})
    else:
        for effect in EFFECTS:
            lm = mats[effect]
            d = lm @ beta
            w = float(d @ np.linalg.solve(lm @ cov @ lm.T, d))
            tests.append({"effect": effect, "wald_statistic": w, "df_numerator": float(lm.shape[0]), "df_denominator": None,
                          "p_value": float(chi2.sf(w, lm.shape[0]))})
        z = norm.ppf(1 - (1 - q["confidence_level"]) / (2 * max(len(vectors), 1)))
        for label, v in zip(labels, vectors):
            est, se = float(v @ beta), float(np.sqrt(v @ cov @ v))
            rows.append({**label, "difference_b_minus_a": est, "standard_error": se, "df": None, "statistic": est / se,
                         "p_unadjusted": float(2 * norm.sf(abs(est / se))), "ci_low": est - z * se, "ci_high": est + z * se,
                         "method": "wald_z", "reportable": True})
        diagnostics.append("asymptotic_wald_no_small_sample_df_correction")
    adjusted = _holm([r["p_unadjusted"] for r in rows]) if rows else []
    for r, p in zip(rows, adjusted):
        r["p_adjusted"] = float(p)
    return tests, rows, diagnostics


def cells(n_arms, k):
    return list(product(range(n_arms), range(k)))
