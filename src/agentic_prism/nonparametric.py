"""Unit-level rank tests and location-shift intervals; Python runtime only.

Exact untied distributions are counted by integer dynamic programming. Normal
intervals invert tie-corrected rank statistics (R wilcox.test, correct=TRUE).
R 4.6 conditional exact inference with ties is deliberately not selected.
"""
from copy import deepcopy
from functools import lru_cache
from itertools import combinations
import math
import numpy as np
from scipy import stats
from scipy.optimize import brentq
from .groups import read_unit_table

DEFAULTS = {"analysis_type": "nonparametric", "schema_version": 1, "input": None,
    "source": "User-supplied unit-level observations",
    "comparison": {"design": None, "groups": [], "outcome": None, "unit": None,
        "rationale": None, "independent_units": None, "rank_model_appropriate": None,
        "location_shift_or_symmetry": None, "inference": "auto", "confidence_level": .95,
        "post_hoc": "none", "adjustment": "holm"},
    "report": {"plot_style": "prism_like"}}
DESIGNS = ("mann_whitney", "wilcoxon_signed_rank", "kruskal_wallis", "friedman")


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported nonparametric config")
    cfg = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(cfg[key], dict):
            if not isinstance(value, dict) or set(value) - set(cfg[key]):
                raise ValueError(f"Unsupported {key} settings")
            cfg[key].update(value)
        else:
            cfg[key] = value
    q = cfg["comparison"]
    if cfg["analysis_type"] != "nonparametric" or type(cfg["schema_version"]) is not int or cfg["schema_version"] != 1:
        raise ValueError("Unsupported schema")
    if any(not isinstance(cfg[k], str) or not cfg[k].strip() for k in ("input", "source")):
        raise ValueError("input and source required")
    if q["design"] not in DESIGNS:
        raise ValueError("Declare a supported rank-test design")
    if not isinstance(q["groups"], list) or any(not isinstance(g, str) or not g.strip() for g in q["groups"]):
        raise ValueError("groups must be a list of names in contrast order")
    n = len(q["groups"])
    if len(set(q["groups"])) != n or (n != 2 if q["design"] in DESIGNS[:2] else n < 3):
        raise ValueError("Two groups for MW/Wilcoxon; at least three for KW/Friedman")
    if any(not isinstance(q[k], str) or not q[k].strip() for k in ("outcome", "unit", "rationale")):
        raise ValueError("Declare outcome, unit and design rationale")
    if any(q[k] is not True for k in ("independent_units", "rank_model_appropriate")):
        raise ValueError("Applicability requires literal true: independent_units and rank_model_appropriate")
    if q["design"] in DESIGNS[:2] and q["location_shift_or_symmetry"] is not True:
        raise ValueError("HL inference requires literal true location_shift_or_symmetry and rationale")
    if q["design"] in DESIGNS[2:] and q["location_shift_or_symmetry"] is not None:
        raise ValueError("location_shift_or_symmetry applies only to MW/Wilcoxon")
    if q["inference"] not in ("auto", "exact", "asymptotic") or (q["design"] in DESIGNS[2:] and q["inference"] != "auto"):
        raise ValueError("inference applies only to MW/Wilcoxon; KW/Friedman use chi-square")
    if type(q["confidence_level"]) not in (int, float) or not .5 < q["confidence_level"] < 1:
        raise ValueError("Invalid confidence_level")
    if q["design"] in DESIGNS[2:] and q["confidence_level"] != .95:
        raise ValueError("confidence_level applies only to MW/Wilcoxon intervals")
    if q["post_hoc"] not in ("none", "dunn_all_pairs") or (q["post_hoc"] != "none" and q["design"] != "kruskal_wallis"):
        raise ValueError("Dunn is available only after declared Kruskal-Wallis")
    if q["adjustment"] not in ("holm", "bonferroni") or (q["post_hoc"] == "none" and q["adjustment"] != "holm"):
        raise ValueError("Adjustment applies to Dunn only")
    if cfg["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    return cfg


def load_data(path, cfg):
    d = read_unit_table(path)
    q = cfg["comparison"]
    if set(d.group) != set(q["groups"]) or set(d.outcome) != {q["outcome"]} or set(d.unit) != {q["unit"]}:
        raise ValueError("Data group/outcome/unit differ from configuration")
    u = d[~d.exclude]
    if u.duplicated(["independent_unit_id", "group"]).any():
        raise ValueError("One observation per independent unit per group required; aggregate upstream with provenance")
    if set(u.group) != set(q["groups"]) or u.groupby("group").size().min() < 3:
        raise ValueError("At least three units in each group required")
    if q["design"] in ("mann_whitney", "kruskal_wallis"):
        if u.independent_unit_id.duplicated().any():
            raise ValueError("Independent groups must not share unit IDs")
    elif not (u.groupby("independent_unit_id").size() == len(q["groups"])).all():
        raise ValueError("Complete matched units/blocks required; no silent incomplete-block deletion")
    return d


@lru_cache(maxsize=128)
def exact_cdf(n, m=None):
    """Distribution of positive-rank sum or Mann-Whitney U, without ties."""
    if m is None:
        counts = [1] + [0] * (n * (n + 1) // 2)
        for r in range(1, n + 1):
            for s in range(len(counts)-1, r-1, -1):
                counts[s] += counts[s-r]
        total = 2**n
    else:
        # Count selections of n ranks out of n+m, then remove the minimum sum.
        size = n*(2*(n+m)-n+1)//2
        dp = [[0]*(size+1) for _ in range(n+1)]
        dp[0][0] = 1
        for r in range(1, n+m+1):
            for j in range(min(n, r), 0, -1):
                for s in range(size, r-1, -1):
                    dp[j][s] += dp[j-1][s-r]
        shift = n*(n+1)//2
        counts = dp[n][shift:shift+n*m+1]
        total = math.comb(n+m, n)
    return np.cumsum(np.array([v/total for v in counts], dtype=float))


def rank_stat(x, y=None, correct=True):
    if y is None:
        z = x[x != 0]
        n = len(z)
        r = stats.rankdata(abs(z))
        _, counts = np.unique(r, return_counts=True)
        val, mean = float(sum(r[z > 0])), n*(n+1)/4
        var = n*(n+1)*(2*n+1)/24 - sum(counts**3-counts)/48
    else:
        n, m = len(x), len(y)
        r = stats.rankdata(np.r_[x, y])
        _, counts = np.unique(r, return_counts=True)
        val, mean = float(sum(r[:n])-n*(n+1)/2), n*m/2
        var = n*m/12 * (n+m+1 - sum(counts**3-counts)/((n+m)*(n+m-1)))
    delta = val-mean
    return val, (delta-(.5*np.sign(delta) if correct else 0))/np.sqrt(var) if var > 0 else np.nan


def rank_shift(x, y=None, inference="auto", level=.95):
    """x minus y (or paired differences x). Two-sided zero-shift test."""
    x = np.asarray(x, float)
    y = None if y is None else np.asarray(y, float)
    n, m = len(x), None if y is None else len(y)
    tied = len(np.unique(abs(x))) != n or np.any(x == 0) if y is None else len(np.unique(np.r_[x,y])) != n+m
    small = n < 50 and (m is None or m < 50)
    if inference == "exact" and (tied or not small):
        raise ValueError("Exact mode requires fewer than 50 observations per sample, without ties or zero paired differences")
    exact = inference == "exact" or (inference == "auto" and small and not tied)
    diffs = np.sort((x[:,None]+x[None,:])[np.triu_indices(n)]/2 if y is None else (x[:,None]-y).ravel())
    val, z = rank_stat(x, y)
    diagnostics = []
    if tied:
        diagnostics.append("ties_or_zero_differences_normal_approximation_selected")
    estimate = float(np.median(diffs))
    lo, hi, achieved = None, None, level
    if exact:
        cdf = exact_cdf(n, m)
        k = int(val)
        # Both exact null distributions are symmetric. Evaluate the opposite
        # tail directly: 1-CDF loses tiny probabilities for separated samples.
        p = min(1., 2*min(cdf[k], cdf[len(cdf)-1-k]))
        qu = int(np.searchsorted(cdf, (1-level)/2))
        if cdf[qu] <= (1-level)/2 + 10*np.finfo(float).eps:
            qu += 1
        if qu:
            lo, hi = float(diffs[qu-1]), float(diffs[len(diffs)-qu])
            achieved = float(1-2*cdf[qu-1])
        else:
            achieved = 1.
            diagnostics.append("exact_interval_unbounded_at_requested_confidence")
    else:
        p = float(2*stats.norm.sf(abs(z))) if np.isfinite(z) else None
        diagnostics.append("normal_approximation_small_samples_require_caution")
        left, right = (float(min(x)), float(max(x))) if y is None else (float(diffs[0]), float(diffs[-1]))
        def w(d, correct=True):
            return rank_stat(x-d, y, correct)[1]
        wl, wr = w(left), w(right)
        alpha = 1-level
        if np.isfinite(wl) and np.isfinite(wr) and left < right:
            if y is None:
                while alpha < 1 and (wl < stats.norm.isf(alpha/2) or wr > stats.norm.ppf(alpha/2)):
                    alpha *= 2
                if alpha > 1-level+1e-12:
                    achieved = 1-min(1,alpha)
                    diagnostics.append("requested_confidence_not_achievable")
            def root(target, correction=True):
                a, b = w(left, correction)-target, w(right, correction)-target
                if a <= 0: return left
                if b >= 0: return right
                return float(brentq(lambda v: w(v, correction)-target, left, right, xtol=1e-4))
            if alpha < 1:
                lo, hi = root(stats.norm.isf(alpha/2)), root(stats.norm.ppf(alpha/2))
            else:
                lo = hi = float(np.median(x))
            # Match R's normal-approximation location estimate; preserve exact HL too.
            estimate = root(0, False)
        else:
            diagnostics.append("rank_interval_not_computable_all_tied_or_zero")
    return {"statistic": val, "p_value": p, "estimate": estimate,
        "hodges_lehmann": float(np.median(diffs)), "ci_low": lo, "ci_high": hi,
        "confidence_level": level, "achieved_confidence_level": achieved,
        "confidence_level_basis": "attained_exact" if exact else "nominal_normal_approximation_not_empirical_coverage",
        "interval_method": "exact_rank_inversion" if exact else "normal_rank_inversion_continuity_corrected",
        "inference": "exact" if exact else "asymptotic", "reportable": p is not None,
        "status": "estimated" if p is not None else "withheld", "diagnostics": diagnostics}


def adjust(p, method):
    p = np.asarray(p, float)
    if method == "bonferroni": return np.minimum(1, len(p)*p)
    order = np.argsort(p)
    out = np.empty(len(p))
    out[order] = np.minimum(1, np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
    return out


def omnibus(samples, design, adjustment="holm", post_hoc="none"):
    """Kruskal-Wallis + Dunn or complete-block Friedman; chi-square asymptotics."""
    diagnostics = ["chi_square_approximation_small_samples_require_caution"]
    if design == "friedman":
        if all(np.ptp(row) == 0 for row in np.array(samples).T):
            stat, p = None, None
        else:
            r = stats.friedmanchisquare(*samples)
            stat, p = float(r.statistic), float(r.pvalue)
        rows = []
    else:
        values = np.concatenate(samples)
        ranks = stats.rankdata(values)
        _, counts = np.unique(values, return_counts=True)
        total = len(values)
        variance = total*(total+1)/12-sum(counts**3-counts)/(12*(total-1))
        rows = []
        if variance <= 0:
            stat, p = None, None
        else:
            r = stats.kruskal(*samples)
            stat, p = float(r.statistic), float(r.pvalue)
            means = [a.mean() for a in np.split(ranks, np.cumsum([len(a) for a in samples])[:-1])]
            if post_hoc != "none":
                for a, b in combinations(range(len(samples)), 2):
                    z = (means[b]-means[a])/np.sqrt(variance*(1/len(samples[a])+1/len(samples[b])))
                    rows.append({"a": a, "b": b, "mean_rank_difference_b_minus_a": float(means[b]-means[a]),
                        "statistic": float(z), "p_unadjusted": float(2*stats.norm.sf(abs(z)))})
                for row, adj in zip(rows, adjust([r["p_unadjusted"] for r in rows], adjustment)):
                    row.update(p_adjusted=float(adj), adjustment=adjustment)
    if p is None: diagnostics.append("zero_rank_variance_inference_withheld")
    return {"statistic": stat, "df": len(samples)-1, "p_value": p, "reportable": p is not None,
        "status": "estimated" if p is not None else "withheld", "diagnostics": diagnostics}, rows


def compare(d, cfg):
    q = cfg["comparison"]
    u = d[~d.exclude]
    samples = [u[u.group == g].set_index("independent_unit_id").value.sort_index() for g in q["groups"]]
    arrays = [v.to_numpy(float) for v in samples]
    rows = []
    if q["design"] in DESIGNS[:2]:
        if q["design"] == "wilcoxon_signed_rank":
            delta = samples[1].reindex(samples[0].index).to_numpy(float)-arrays[0]
            fit = rank_shift(delta, inference=q["inference"], level=q["confidence_level"])
        else:
            fit = rank_shift(arrays[1], arrays[0], q["inference"], q["confidence_level"])
    else:
        fit, rows = omnibus(arrays, q["design"], q["adjustment"], q["post_hoc"])
        for row in rows:
            row["group_a"], row["group_b"] = q["groups"][row.pop("a")], q["groups"][row.pop("b")]
    fit.update(design=q["design"], n_total=len(u), n_units=u.independent_unit_id.nunique(), n_groups=len(samples))
    summaries = [{"group": g, "n": len(v), "median": float(np.median(v))} for g,v in zip(q["groups"],arrays)]
    return {"schema_version": 1, "analysis_type": "nonparametric", "fits": [fit], "contrasts": rows, "group_summaries": summaries}
