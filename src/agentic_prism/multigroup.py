"""One-way comparison of three or more independent groups with a predeclared family.

Classic one-way ANOVA pairs with Dunnett (vs one control) or Tukey-Kramer (all
pairs); Welch ANOVA pairs with Games-Howell (all pairs) or Holm-adjusted Welch t
tests against one control. Contrasts are always reported for the declared family
and are not gated on the omnibus test. No repeated-measures, covariate or
nonparametric method is included.
"""
from copy import deepcopy
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import dunnett, f as f_dist, t as t_dist, ttest_ind, tukey_hsd

DEFAULTS = {"analysis_type": "multi_group_comparison", "schema_version": 1, "input": None,
            "source": "User-supplied independent experimental units",
            "comparison": {"design": None, "groups": None, "control": None, "post_hoc": None,
                           "outcome": None, "unit": None, "rationale": None,
                           "confidence_level": .95, "sd_ratio_warning": 3.},
            "report": {"plot_style": "prism_like"}}
POST_HOC = {"one_way_anova": ("dunnett_vs_control", "tukey_all_pairs", "none"),
            "welch_anova": ("games_howell_all_pairs", "holm_welch_vs_control", "none")}
NEEDS_CONTROL = ("dunnett_vs_control", "holm_welch_vs_control")
# Dunnett's multivariate-t probabilities use randomized quasi-Monte Carlo in SciPy.
DUNNETT_SEED = 20260926


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported multi-group comparison config")
    c = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(c[key], dict):
            if not isinstance(value, dict) or set(value) - set(c[key]):
                raise ValueError(f"Unsupported {key} settings")
            c[key].update(value)
        else:
            c[key] = value
    q = c["comparison"]
    if c["analysis_type"] != "multi_group_comparison" or c["schema_version"] != 1 or not isinstance(c["input"], str) or not c["input"]:
        raise ValueError("Unsupported multi-group schema or missing input")
    if q["design"] not in POST_HOC:
        raise ValueError("comparison.design must be one_way_anova or welch_anova")
    if q["post_hoc"] not in POST_HOC[q["design"]]:
        raise ValueError(f"post_hoc for {q['design']} must be one of {POST_HOC[q['design']]}")
    groups = q["groups"]
    if (not isinstance(groups, list) or len(groups) < 3 or any(not isinstance(g, str) or not g.strip() for g in groups)
            or len(set(groups)) != len(groups)):
        raise ValueError("Declare at least three distinct group names; use group_comparison for two groups")
    if q["post_hoc"] in NEEDS_CONTROL:
        if q["control"] not in groups:
            raise ValueError("Declare the control as one of the groups for a vs-control family")
    elif q["control"] is not None:
        raise ValueError("control applies only to a vs-control family")
    if any(not isinstance(q[k], str) or not q[k].strip() for k in ("outcome", "unit", "rationale")):
        raise ValueError("Declare outcome, unit and design rationale")
    for key, low, high in (("confidence_level", .5, 1), ("sd_ratio_warning", 1, np.inf)):
        if isinstance(q[key], bool) or not isinstance(q[key], (int, float)) or not low < q[key] < high:
            raise ValueError(f"Invalid {key}")
    if c["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    return c


def load_data(path, cfg):
    from .groups import read_unit_table
    q = cfg["comparison"]
    d = read_unit_table(path)
    if set(d.group) != set(q["groups"]) or set(d.outcome) != {q["outcome"]} or set(d.unit) != {q["unit"]}:
        raise ValueError("Input group/outcome/unit does not match declared comparison")
    used = d[~d.exclude]
    if used.independent_unit_id.duplicated().any():
        raise ValueError("Multiple measurements of one independent unit; this design needs disjoint units, "
                         "one prespecified value each (aggregate technical replicates upstream with provenance)")
    return d


def _summaries(values, groups):
    rows = []
    for g in groups:
        v = values[g]
        rows.append({"group": g, "n": len(v), "mean": float(np.mean(v)), "sd": float(np.std(v, ddof=1)),
                     "sem": float(np.std(v, ddof=1) / np.sqrt(len(v))), "min": float(np.min(v)), "max": float(np.max(v))})
    return rows


def oneway(values, groups):
    """Classic equal-variance one-way ANOVA by sums of squares."""
    arrays = [np.asarray(values[g], float) for g in groups]
    n = np.array([len(a) for a in arrays])
    grand = np.concatenate(arrays).mean()
    ssb = float(sum(len(a) * (a.mean() - grand) ** 2 for a in arrays))
    ssw = float(sum(((a - a.mean()) ** 2).sum() for a in arrays))
    df1, df2 = len(arrays) - 1, int(n.sum() - len(arrays))
    if ssw <= 1e-24:
        raise ValueError("Within-group variance is zero; ANOVA inference withheld")
    stat = (ssb / df1) / (ssw / df2)
    return {"f_statistic": float(stat), "df_numerator": float(df1), "df_denominator": float(df2),
            "p_value": float(f_dist.sf(stat, df1, df2)), "ss_between": ssb, "ss_within": ssw,
            "mse": ssw / df2, "eta_squared": ssb / (ssb + ssw)}


def welch(values, groups):
    """Welch (1951) heteroscedastic one-way ANOVA."""
    arrays = [np.asarray(values[g], float) for g in groups]
    k = len(arrays)
    n = np.array([len(a) for a in arrays], float)
    var = np.array([a.var(ddof=1) for a in arrays])
    if (var <= 1e-24).any():
        raise ValueError("A group has zero variance; Welch ANOVA inference withheld")
    w = n / var
    mean = np.array([a.mean() for a in arrays])
    centre = (w * mean).sum() / w.sum()
    tmp = ((1 - w / w.sum()) ** 2 / (n - 1)).sum()
    a = (w * (mean - centre) ** 2).sum() / (k - 1)
    b = 1 + 2 * (k - 2) / (k ** 2 - 1) * tmp
    stat, df2 = a / b, (k ** 2 - 1) / (3 * tmp)
    return {"f_statistic": float(stat), "df_numerator": float(k - 1), "df_denominator": float(df2),
            "p_value": float(f_dist.sf(stat, k - 1, df2))}


def holm(p):
    """Holm step-down adjusted p values, in input order."""
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., (len(p) - rank) * p[i]))
        adjusted[i] = running
    return adjusted


def contrasts(values, groups, q):
    """Declared contrast family; differences are later group minus earlier group/control."""
    kind, level = q["post_hoc"], q["confidence_level"]
    rows = []
    if kind == "none":
        return rows, None
    if kind == "dunnett_vs_control":
        others = [g for g in groups if g != q["control"]]
        res = dunnett(*[values[g] for g in others], control=values[q["control"]], rng=DUNNETT_SEED)
        ci = res.confidence_interval(level)
        for i, g in enumerate(others):
            rows.append({"group_a": q["control"], "group_b": g, "difference_b_minus_a": float(np.mean(values[g]) - np.mean(values[q["control"]])),
                         "statistic": float(res.statistic[i]), "p_adjusted": float(res.pvalue[i]),
                         "ci_low": float(ci.low[i]), "ci_high": float(ci.high[i])})
        return rows, f"Dunnett many-to-one, pooled variance; multivariate-t probabilities by SciPy randomized QMC (seed {DUNNETT_SEED})"
    if kind in ("tukey_all_pairs", "games_howell_all_pairs"):
        res = tukey_hsd(*[values[g] for g in groups], equal_var=kind == "tukey_all_pairs")
        ci = res.confidence_interval(level)
        for i, j in combinations(range(len(groups)), 2):
            # SciPy reports mean_i - mean_j (its "statistic" is that difference, not q); flip to later-minus-earlier order.
            rows.append({"group_a": groups[i], "group_b": groups[j], "difference_b_minus_a": float(np.mean(values[groups[j]]) - np.mean(values[groups[i]])),
                         "statistic": None, "p_adjusted": float(res.pvalue[i, j]),
                         "ci_low": float(-ci.high[i, j]), "ci_high": float(-ci.low[i, j])})
        return rows, ("Tukey-Kramer studentized range, pooled variance" if kind == "tukey_all_pairs"
                      else "Games-Howell studentized range with per-pair Welch df")
    control = values[q["control"]]
    others = [g for g in groups if g != q["control"]]
    m = len(others)
    raw = []
    for g in others:
        b, a = np.asarray(values[g], float), np.asarray(control, float)
        sa, sb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = np.sqrt(sa + sb)
        df = (sa + sb) ** 2 / (sa ** 2 / (len(a) - 1) + sb ** 2 / (len(b) - 1))
        test = ttest_ind(b, a, equal_var=False)
        margin = t_dist.ppf(1 - (1 - level) / (2 * m), df) * se
        diff = float(b.mean() - a.mean())
        raw.append(float(test.pvalue))
        rows.append({"group_a": q["control"], "group_b": g, "difference_b_minus_a": diff, "statistic": float(test.statistic),
                     "df": float(df), "p_unadjusted": float(test.pvalue), "ci_low": diff - margin, "ci_high": diff + margin})
    for row, p in zip(rows, holm(raw)):
        row["p_adjusted"] = float(p)
    return rows, "Welch t vs control, Holm step-down p; Bonferroni simultaneous intervals"


def compare(d, cfg):
    q = cfg["comparison"]
    used = d[~d.exclude]
    values = {g: used.loc[used.group == g, "value"].to_numpy(float) for g in q["groups"]}
    small = [g for g, v in values.items() if len(v) < 3]
    if small:
        raise ValueError(f"At least three independent units per group are required: {small}")
    summaries = _summaries(values, q["groups"])
    omnibus = oneway(values, q["groups"]) if q["design"] == "one_way_anova" else welch(values, q["groups"])
    family, method = contrasts(values, q["groups"], q)
    sds = [s["sd"] for s in summaries]
    ratio = max(sds) / min(sds) if min(sds) > 0 else float("inf")
    diagnostics = []
    if q["design"] == "one_way_anova" and ratio > q["sd_ratio_warning"]:
        diagnostics.append("unequal_group_sd_consider_welch_design")
    if len(set(len(v) for v in values.values())) > 1:
        diagnostics.append("unequal_group_sizes")
    # validation/multigroup_error_rates.json: with an n=4 high-variance group, Welch ANOVA and
    # Games-Howell exceeded nominal error (0.069 and 0.077 at alpha 0.05).
    if q["design"] == "welch_anova" and min(len(v) for v in values.values()) < 6:
        diagnostics.append("small_group_welch_procedures_may_be_liberal")
    return {"analysis_type": "multi_group_comparison", "design": q["design"], "post_hoc": q["post_hoc"],
            "groups": q["groups"], "control": q["control"], "outcome": q["outcome"], "unit": q["unit"],
            "n_groups": len(q["groups"]), "n_total": int(sum(len(v) for v in values.values())),
            **omnibus, "max_min_sd_ratio": float(ratio), "confidence_level": q["confidence_level"],
            "contrast_method": method, "n_contrasts": len(family),
            "status": "estimated", "reportable": True, "diagnostics": diagnostics}, summaries, family
