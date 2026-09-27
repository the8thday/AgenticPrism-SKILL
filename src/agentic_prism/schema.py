"""Strict configuration and observation contracts for the first supported module."""
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import pandas as pd

UNITS = {"M": 1., "mM": 1e-3, "uM": 1e-6, "µM": 1e-6, "μM": 1e-6, "nM": 1e-9, "pM": 1e-12}
DEFAULTS = {
    "schema_version": 1,
    "analysis_type": "equilibrium_binding",
    "model": "one_site_with_baseline",
    "input": None,
    "column_map": {},
    "source": "User-supplied data; see input snapshot and hash.",
    "assay": {"equilibrium_supported": None, "signal_proportional": None,
              "single_site_supported": None, "concentration_basis": None,
              "free_approximation_supported": None, "interpretation": "KD", "rationale": ""},
    "fit": {"baseline_mode": "fitted", "baseline_value": None,
            "residual_scale": "linear", "weighting": "unweighted", "multistart": 5,
            "kd_bounds_M": [1e-15, 1.], "max_nfev": 3000},
    "uncertainty": {"parameter_ci": "profile_f", "level": .95},
    "replicates": {"independent_unit": "none", "conditions_comparable": False},
    "sensitivity": [],
    "report": {"language": "zh-CN", "plot_style": "prism_like", "allow_style_switch": True,
               "export_formats": ["svg", "pdf", "png"], "figure_width_mm": 85,
               "png_dpi": 300},
}


def resolve_config(raw):
    def merge(base, given, path=""):
        if not isinstance(given, dict):
            raise ValueError(f"{path or 'config'} must be an object")
        unknown = set(given) - set(base)
        if unknown:
            raise ValueError(f"Unsupported setting at {path}: {sorted(unknown)}")
        out = deepcopy(base)
        for k, v in given.items():
            out[k] = merge(base[k], v, path + k + ".") if isinstance(base[k], dict) and k != "column_map" else v
        return out
    c = merge(DEFAULTS, raw)
    if c["schema_version"] != 1 or c["analysis_type"] != "equilibrium_binding" or c["model"] != "one_site_with_baseline":
        raise ValueError("Only schema 1 / equilibrium_binding / one_site_with_baseline is supported")
    a = c["assay"]
    if any(a[k] is not True for k in ("equilibrium_supported", "signal_proportional", "single_site_supported")):
        raise ValueError("Assay applicability unresolved: equilibrium, proportional signal and single-site support must be documented")
    if a["concentration_basis"] not in ("free", "total") or (not isinstance(a["rationale"], str) or not a["rationale"].strip()):
        raise ValueError("Specify free/total concentration and an assay rationale")
    if a["concentration_basis"] == "total" and a["free_approximation_supported"] is not True:
        raise ValueError("Total concentration requires a justified free-concentration approximation; depletion model unavailable")
    if a["interpretation"] not in ("KD", "apparent_KD"):
        raise ValueError("interpretation must be KD or a justified apparent_KD")
    f = c["fit"]
    if f["baseline_mode"] not in ("fitted", "fixed", "zero_control") or f["residual_scale"] not in ("linear", "log"):
        raise ValueError("Unsupported baseline or residual mode")
    if f["weighting"] not in ("unweighted", "inverse_sd"):
        raise ValueError("weighting must be unweighted or inverse_sd")
    if f["weighting"] == "inverse_sd" and f["residual_scale"] != "linear":
        raise ValueError("inverse_sd is supported only for linear residuals")
    if f["baseline_mode"] == "fixed" and (not isinstance(f["baseline_value"], (float, int)) or not np.isfinite(f["baseline_value"])):
        raise ValueError("A finite baseline_value is required for fixed mode")
    if type(f["multistart"]) is not int or not 1 <= f["multistart"] <= 20:
        raise ValueError("multistart must be an integer between 1 and 20")
    if type(f["max_nfev"]) is not int or f["max_nfev"] < 100:
        raise ValueError("max_nfev must be an integer >=100")
    b = f["kd_bounds_M"]
    if not isinstance(b, list) or len(b) != 2 or not all(np.isfinite(b)) or not 0 < b[0] < b[1]:
        raise ValueError("kd_bounds_M must contain two increasing positive finite values")
    if c["uncertainty"]["parameter_ci"] not in ("profile_f", "none") or not .5 < c["uncertainty"]["level"] < 1:
        raise ValueError("Unsupported interval method/level")
    if c["replicates"]["independent_unit"] not in ("none", "experiment_id"):
        raise ValueError("Independent unit must be none or experiment_id")
    if type(c["replicates"]["conditions_comparable"]) is not bool:
        raise ValueError("conditions_comparable must be boolean")
    if not isinstance(c["sensitivity"], list) or set(c["sensitivity"]) - {"alternate_loss", "exclude_highest"}:
        raise ValueError("Supported sensitivity scenarios: alternate_loss, exclude_highest")
    r = c["report"]
    if r["language"] != "zh-CN" or r["plot_style"] not in ("standard", "prism_like"):
        raise ValueError("Supported report language zh-CN; styles standard/prism_like")
    if set(r["export_formats"]) != {"svg", "pdf", "png"} or len(r["export_formats"]) != 3:
        raise ValueError("This release exports all of SVG, PDF and PNG")
    if r["figure_width_mm"] not in (85, 180) or r["png_dpi"] not in (300, 600) or type(r["allow_style_switch"]) is not bool:
        raise ValueError("Invalid report size, DPI or switching option")
    if not isinstance(c["column_map"], dict) or not all(isinstance(k, str) and isinstance(v, str) for k,v in c["column_map"].items()):
        raise ValueError("column_map must map original column names to canonical names")
    json.dumps(c, allow_nan=False)
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False).rename(columns=cfg["column_map"])
    if d.columns.duplicated().any():
        raise ValueError("Duplicate columns after mapping")
    required = {"sample_id", "experiment_id", "curve_id", "concentration", "concentration_unit", "response", "response_unit"}
    if required - set(d):
        raise ValueError(f"Missing columns: {sorted(required - set(d))}")
    if d.empty:
        raise ValueError("Empty input")
    d["source_row"] = np.arange(2, len(d) + 2)
    for k in ("sample_id", "curve_id", "response_unit", "concentration_unit"):
        if d[k].str.strip().eq("").any():
            raise ValueError(f"Missing {k}")
    if cfg["replicates"]["independent_unit"] == "experiment_id" and d.experiment_id.str.strip().eq("").any():
        raise ValueError("Independent experiment_id is missing; do not infer it from date")
    for k in ("concentration", "response"):
        d[k] = pd.to_numeric(d[k], errors="raise")
        if not np.isfinite(d[k]).all():
            raise ValueError(f"Non-finite {k}")
    if (d.concentration < 0).any() or not d.concentration_unit.isin(UNITS).all():
        raise ValueError("Concentrations must be nonnegative with supported molar units; no implicit mass conversion")
    d["concentration_M"] = d.concentration * d.concentration_unit.map(UNITS)
    if "observation_id" not in d:
        d["observation_id"] = [f"row-{i}" for i in d.source_row]
    if d.observation_id.str.strip().eq("").any() or d.observation_id.duplicated().any():
        raise ValueError("observation_id must be nonempty and unique")
    if "exclude" not in d:
        d["exclude"] = "false"
    if not d.exclude.str.lower().isin(["true", "false"]).all():
        raise ValueError("exclude must be true or false")
    d["exclude"] = d.exclude.str.lower().eq("true")
    if "exclusion_reason" not in d:
        d["exclusion_reason"] = ""
    if (d.exclude & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Every excluded observation needs an exclusion_reason")
    if "response_sd" in d:
        d["response_sd"] = pd.to_numeric(d.response_sd.replace("", np.nan), errors="raise")
    if cfg["fit"]["weighting"] == "inverse_sd":
        if "response_sd" not in d or not np.isfinite(d.loc[~d.exclude, "response_sd"]).all() or (d.loc[~d.exclude, "response_sd"] <= 0).any():
            raise ValueError("inverse_sd requires positive finite response_sd per included observation")
    for cid, g in d.groupby("curve_id", sort=False):
        if any(g[k].nunique(dropna=False) != 1 for k in ("sample_id", "experiment_id", "response_unit")):
            raise ValueError(f"Curve {cid} maps to inconsistent sample, experiment or response unit")
    return d
