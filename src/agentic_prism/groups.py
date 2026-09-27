"""Two-group independent Welch and explicitly paired t comparisons."""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist, ttest_ind, ttest_rel

DEFAULTS = {"analysis_type": "group_comparison", "schema_version": 1, "input": None,
            "source": "User-supplied independent experimental units",
            "comparison": {"design": None, "group_a": None, "group_b": None, "outcome": None,
                           "unit": None, "rationale": None, "confidence_level": .95},
            "report": {"plot_style": "prism_like"}}


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError("Unsupported group comparison config")
    c = deepcopy(DEFAULTS)
    for key, value in raw.items():
        if isinstance(c[key], dict):
            if not isinstance(value, dict) or set(value) - set(c[key]):
                raise ValueError(f"Unsupported {key} settings")
            c[key].update(value)
        else:
            c[key] = value
    q = c["comparison"]
    if c["analysis_type"] != "group_comparison" or c["schema_version"] != 1 or not isinstance(c["input"], str) or not c["input"]:
        raise ValueError("Unsupported group schema or missing input")
    if q["design"] not in ("independent_welch", "paired_t") or any(not isinstance(q[k], str) or not q[k].strip() for k in
        ("group_a", "group_b", "outcome", "unit", "rationale")) or q["group_a"] == q["group_b"]:
        raise ValueError("Declare design, two groups, outcome, unit and design rationale")
    if isinstance(q["confidence_level"], bool) or not isinstance(q["confidence_level"], (int, float)) or not .5 < q["confidence_level"] < 1:
        raise ValueError("Invalid confidence level")
    if c["report"]["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Unsupported plot style")
    return c


def read_unit_table(path):
    """Columns, IDs, finite values and explicit exclusions shared by group designs."""
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"observation_id", "independent_unit_id", "group", "outcome", "value", "unit"}
    if d.empty or required - set(d) or d.columns.duplicated().any():
        raise ValueError(f"Invalid group table; missing {sorted(required-set(d))}")
    for key in required - {"value"}:
        if d[key].str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
    if d.observation_id.duplicated().any():
        raise ValueError("Duplicate observation ID")
    d["value"] = pd.to_numeric(d.value, errors="raise")
    if not np.isfinite(d.value).all():
        raise ValueError("Nonfinite value")
    d["exclude"] = d.get("exclude", pd.Series("false", index=d.index)).str.lower()
    d["exclusion_reason"] = d.get("exclusion_reason", pd.Series("", index=d.index))
    if not d.exclude.isin(("true", "false")).all() or ((d.exclude == "true") & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Exclusions require true/false and a reason")
    d["exclude"] = d.exclude.eq("true")
    return d


def load_data(path, cfg):
    d = read_unit_table(path)
    q = cfg["comparison"]
    if set(d.group) != {q["group_a"], q["group_b"]} or set(d.outcome) != {q["outcome"]} or set(d.unit) != {q["unit"]}:
        raise ValueError("Input group/outcome/unit does not match declared comparison")
    used = d[~d.exclude]
    if used.duplicated(["independent_unit_id", "group"]).any():
        raise ValueError("Multiple measurements of one independent unit/group; aggregate upstream with provenance")
    if q["design"] == "independent_welch":
        if used.independent_unit_id.duplicated().any():
            raise ValueError("Independent groups share unit IDs")
    else:
        if set(used.loc[used.group == q["group_a"], "independent_unit_id"]) != set(used.loc[used.group == q["group_b"], "independent_unit_id"]):
            raise ValueError("Paired groups need complete matching independent_unit_id values")
    return d


def compare(d, cfg):
    q = cfg["comparison"]
    used = d[~d.exclude]
    a = used[used.group == q["group_a"]].set_index("independent_unit_id").value
    b = used[used.group == q["group_b"]].set_index("independent_unit_id").value
    if min(len(a), len(b)) < 3:
        raise ValueError("At least three independent units per group or three complete pairs are required")
    if q["design"] == "paired_t":
        b = b.reindex(a.index)
        delta = b.to_numpy(float) - a.to_numpy(float)
        if np.var(delta, ddof=1) <= 1e-24:
            raise ValueError("Paired differences have zero variance; t inference withheld")
        stat = ttest_rel(b, a)
        diff, se, df = float(np.mean(delta)), float(np.std(delta, ddof=1) / np.sqrt(len(delta))), len(delta)-1
    else:
        av, bv = a.to_numpy(float), b.to_numpy(float)
        va, vb = float(np.var(av, ddof=1)), float(np.var(bv, ddof=1))
        if va + vb <= 1e-24:
            raise ValueError("Both groups have zero variance; t inference withheld")
        sa, sb = va/len(av), vb/len(bv)
        df = float((sa+sb)**2/(sa**2/(len(av)-1) + sb**2/(len(bv)-1)))
        diff, se = float(np.mean(bv)-np.mean(av)), float(np.sqrt(sa+sb))
        stat = ttest_ind(bv, av, equal_var=False)
    margin = float(t_dist.ppf((1+q["confidence_level"])/2, df)*se)
    return {"analysis_type": "group_comparison", "design": q["design"], "group_a": q["group_a"],
            "group_b": q["group_b"], "outcome": q["outcome"], "unit": q["unit"],
            "n_a": len(a), "n_b": len(b), "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "mean_difference_b_minus_a": diff, "standard_error": se, "df": float(df),
            "t_statistic": float(stat.statistic), "p_two_sided": float(stat.pvalue),
            "confidence_level": q["confidence_level"], "ci_difference": [diff-margin, diff+margin],
            "status": "estimated", "diagnostics": []}
