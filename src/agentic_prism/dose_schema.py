"""Strict input contract for relative EC50/IC50 curves: symmetric 4PL, or opt-in 5PL / bell-shaped (0.13.1)."""
from copy import deepcopy
import json
import numpy as np
import pandas as pd
from .schema import UNITS

DOSE_UNITS = {**UNITS, "mg/mL": 1., "ug/mL": 1e-3, "µg/mL": 1e-3,
              "μg/mL": 1e-3, "ng/mL": 1e-6, "source_unit": 1.}
MASS_UNITS = {"mg/mL", "ug/mL", "µg/mL", "μg/mL", "ng/mL"}

MODELS = ("relative_four_parameter_logistic", "relative_five_parameter_logistic", "bell_shaped")
DEFAULTS = {
    "schema_version": 1, "analysis_type": "dose_response_4pl", "model": "relative_four_parameter_logistic",
    # model_rationale is required for the opt-in models and removed from 4PL configs (legacy bytes).
    "model_rationale": "",
    "input": None, "source": "User-supplied dose-response observations", "column_map": {},
    "provenance": None,
    "effect_levels": {"percentages": [], "rationale": ""},
    "dose_scale": "linear",  # "log10" means log10 of values in the supplied concentration_unit.
    "assay": {"endpoint": None, "direction": None, "tested_agent": "", "response_definition": "",
              "relative_half_response_supported": None, "rationale": ""},
    # fixed_bottom/fixed_top constrain the lower/higher plateau (Prism naming), in response units.
    "fit": {"weighting": "unweighted", "fixed_bottom": None, "fixed_top": None,
            "hill_bounds": [0.05, 8.], "log50_bounds": None, "multistart": 9, "max_nfev": 2000},
    "uncertainty": {"method": "profile_f", "level": .95},
    "replicates": {"independent_unit": "none", "conditions_comparable": False},
    # Each item: {"id", "reference_curve", "test_curve"}; curves share one canonical unit.
    "comparisons": [],
    # parallelism_method "equivalence" replaces the F-test gate with predeclared margins (USP <1032> style).
    "comparison_settings": {"parallelism_alpha": .05, "parallelism_method": "f_test",
                            "equivalence": {"confidence_level": .90, "hill_ratio_limits": None,
                                            "bottom_difference_limits": None, "top_difference_limits": None,
                                            "rationale": ""},
                            "rp_acceptance_limits": None, "rp_acceptance_rationale": ""},
    "report": {"language": "zh-CN", "plot_style": "prism_like", "figure_width_mm": 180,
               "png_dpi": 300},
}


def resolve_dose_config(raw):
    def merge(base, given, path="config"):
        if not isinstance(given, dict):
            raise ValueError(f"{path} must be an object")
        extra = set(given) - set(base)
        if extra:
            raise ValueError(f"Unsupported setting at {path}: {sorted(extra)}")
        result = deepcopy(base)
        for key, value in given.items():
            result[key] = merge(base[key], value, path + "." + key) if isinstance(base[key], dict) and key != "column_map" else value
        return result

    c = merge(DEFAULTS, raw)
    levels = c['effect_levels']
    values = levels['percentages']
    if not isinstance(values, list) or any(type(v) not in (int, float) or not np.isfinite(v) or not 0 < v < 100 for v in values) or len(set(values)) != len(values):
        raise ValueError('effect_levels.percentages must be distinct finite percentages strictly between 0 and 100')
    if values:
        if c['model'] != MODELS[0]:
            raise ValueError('ECx/ICx effect_levels currently require the symmetric 4PL; no 5PL or bell-phase endpoints')
        if not isinstance(levels['rationale'], str) or not levels['rationale'].strip():
            raise ValueError('Declare effect_levels.rationale and the requested relative-effect targets before fitting')
    else:
        if levels['rationale']:
            raise ValueError('effect_levels.rationale requires percentages')
        c.pop('effect_levels')
    if c["schema_version"] != 1 or c["analysis_type"] != "dose_response_4pl" or c["model"] not in MODELS:
        raise ValueError("Unsupported dose-response schema or model; models are " + ", ".join(MODELS))
    if c["model"] == MODELS[0]:
        if c["model_rationale"]:
            raise ValueError("model_rationale applies only to the opt-in 5PL and bell-shaped models")
        c.pop("model_rationale")
    else:
        if not isinstance(c["model_rationale"], str) or not c["model_rationale"].strip():
            raise ValueError("Declare model_rationale: why this shape was chosen before seeing the fit (assay biology or prior runs)")
        if c["comparisons"]:
            raise ValueError("Relative-potency comparisons use the parallel symmetric 4PL only")
        if c["model"] == "bell_shaped" and (c["fit"]["weighting"] != "unweighted" or c["fit"]["fixed_bottom"] is not None
                                            or c["fit"]["fixed_top"] is not None):
            raise ValueError("bell_shaped supports unweighted fits with three free plateaus only")
    a = c["assay"]
    if a["endpoint"] not in ("EC50", "IC50") or a["direction"] not in ("increasing", "decreasing"):
        raise ValueError("Specify the biological endpoint EC50/IC50 and observed response direction independently")
    if a["relative_half_response_supported"] is not True or any(not isinstance(a[k], str) or not a[k].strip()
                                                         for k in ("tested_agent", "response_definition", "rationale")):
        raise ValueError("Document the tested agent, response meaning and relative-half-response applicability")
    if c["dose_scale"] not in ("linear", "log10") or not isinstance(c["input"], str) or not c["input"]:
        raise ValueError("Specify an input file and explicit linear/log10 dose scale")
    if c["provenance"] is not None and (not isinstance(c["provenance"], str) or not c["provenance"]):
        raise ValueError("provenance must be null or a JSON file path")
    if not isinstance(c["column_map"], dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in c["column_map"].items()):
        raise ValueError("column_map must map source names to canonical names")
    f = c["fit"]
    if f["weighting"] not in ("unweighted", "relative"):
        raise ValueError("weighting must be unweighted or relative (1/Ycurve^2)")
    for key in ("fixed_bottom", "fixed_top"):
        if f[key] is not None and (isinstance(f[key], bool) or not isinstance(f[key], (int, float)) or not np.isfinite(f[key])):
            raise ValueError(f"{key} must be null or a finite response value")
    if f["fixed_bottom"] is not None and f["fixed_top"] is not None and not f["fixed_top"] > f["fixed_bottom"]:
        raise ValueError("fixed_top must exceed fixed_bottom; Top is the higher plateau for either direction")
    if f["weighting"] == "relative" and any(f[k] is not None and f[k] <= 0 for k in ("fixed_bottom", "fixed_top")):
        raise ValueError("Relative weighting divides by the fitted curve; fixed plateaus must be positive")
    h = f["hill_bounds"]
    if not isinstance(h, list) or len(h) != 2 or not all(np.isfinite(h)) or not 0 < h[0] < h[1]:
        raise ValueError("hill_bounds must be increasing positive finite values")
    b = f["log50_bounds"]
    if b is not None and (not isinstance(b, list) or len(b) != 2 or not all(np.isfinite(b)) or b[0] >= b[1]):
        raise ValueError("log50_bounds must be null or two increasing finite values in canonical units")
    if type(f["multistart"]) is not int or not 1 <= f["multistart"] <= 25 or type(f["max_nfev"]) is not int or f["max_nfev"] < 100:
        raise ValueError("Invalid multistart or maximum function evaluations")
    u = c["uncertainty"]
    if u["method"] not in ("profile_f", "none") or not isinstance(u["level"], (int, float)) or not .5 < u["level"] < 1:
        raise ValueError("Unsupported uncertainty method or level")
    rep = c["replicates"]
    if rep["independent_unit"] not in ("none", "experiment_id") or type(rep["conditions_comparable"]) is not bool:
        raise ValueError("replicates.independent_unit must be none/experiment_id; conditions_comparable boolean")
    if not isinstance(c["comparisons"], list):
        raise ValueError("comparisons must be a list")
    ids = set()
    for item in c["comparisons"]:
        if not isinstance(item, dict) or set(item) != {"id", "reference_curve", "test_curve"} or \
                not all(isinstance(v, str) and v.strip() for v in item.values()):
            raise ValueError("Each comparison needs exactly id, reference_curve and test_curve strings")
        if item["id"] in ids or item["reference_curve"] == item["test_curve"]:
            raise ValueError("Comparison ids must be unique and compare two different curves")
        ids.add(item["id"])
    alpha = c["comparison_settings"]["parallelism_alpha"]
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not 0 < alpha < .5:
        raise ValueError("parallelism_alpha must be in (0, 0.5)")
    cs = c["comparison_settings"]
    eq = cs["equivalence"]

    def interval(value, around):
        return (isinstance(value, list) and len(value) == 2 and all(type(v) in (int, float) and np.isfinite(v) for v in value)
                and value[0] < around < value[1])
    if cs["parallelism_method"] not in ("f_test", "equivalence"):
        raise ValueError("parallelism_method must be f_test or equivalence")
    margins = ("hill_ratio_limits", "bottom_difference_limits", "top_difference_limits")
    if cs["parallelism_method"] == "equivalence":
        if not interval(eq["hill_ratio_limits"], 1.) or eq["hill_ratio_limits"][0] <= 0:
            raise ValueError("equivalence.hill_ratio_limits must be [low, high] with 0 < low < 1 < high")
        for key in margins[1:]:
            if eq[key] is not None and not interval(eq[key], 0.):
                raise ValueError(f"equivalence.{key} must be null or [low, high] around 0 in response units")
        if not isinstance(eq["rationale"], str) or not eq["rationale"].strip():
            raise ValueError("Equivalence margins need a rationale naming their historical or predeclared source")
        if isinstance(eq["confidence_level"], bool) or not isinstance(eq["confidence_level"], (int, float)) \
                or not .5 < eq["confidence_level"] < 1:
            raise ValueError("equivalence.confidence_level must be in (0.5, 1)")
    elif any(eq[k] is not None for k in margins) or eq["rationale"]:
        raise ValueError("Equivalence margins are set but parallelism_method is f_test")
    if cs["rp_acceptance_limits"] is not None:
        if not interval(cs["rp_acceptance_limits"], 1.) or cs["rp_acceptance_limits"][0] <= 0:
            raise ValueError("rp_acceptance_limits must be null or [low, high] with 0 < low < 1 < high")
        if not isinstance(cs["rp_acceptance_rationale"], str) or not cs["rp_acceptance_rationale"].strip():
            raise ValueError("rp_acceptance_limits need a rationale naming their source")
    elif cs["rp_acceptance_rationale"]:
        raise ValueError("rp_acceptance_rationale given without rp_acceptance_limits")
    r = c["report"]
    if r["language"] != "zh-CN" or r["plot_style"] not in ("prism_like", "standard") or r["figure_width_mm"] not in (85, 180) or r["png_dpi"] not in (300, 600):
        raise ValueError("Unsupported report settings")
    json.dumps(c, allow_nan=False)
    return c


def load_dose_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False).rename(columns=cfg["column_map"])
    required = {"sample_id", "experiment_id", "curve_id", "concentration", "concentration_unit",
                "response", "response_unit"}
    if d.empty or required - set(d) or d.columns.duplicated().any():
        raise ValueError(f"Invalid 4PL table; missing columns {sorted(required - set(d))}")
    for key in ("sample_id", "experiment_id", "curve_id", "concentration_unit", "response_unit"):
        if d[key].str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
    for key in ("concentration", "response"):
        d[key] = pd.to_numeric(d[key], errors="raise")
        if not np.isfinite(d[key]).all():
            raise ValueError(f"Nonfinite {key}")
    if not d.concentration_unit.isin(DOSE_UNITS).all():
        raise ValueError("Unsupported concentration unit")
    if cfg["dose_scale"] == "linear":
        if (d.concentration < 0).any():
            raise ValueError("Linear concentration cannot be negative")
        d["concentration_linear"] = d.concentration.astype(float)
    else:
        d["concentration_linear"] = 10. ** d.concentration.to_numpy(float)
        if not np.isfinite(d.concentration_linear).all() or (d.concentration_linear <= 0).any():
            raise ValueError("Invalid log10 concentration conversion")
    d["concentration_canonical"] = d.concentration_linear * d.concentration_unit.map(DOSE_UNITS)
    if not np.isfinite(d.concentration_canonical).all():
        raise ValueError("Nonfinite canonical concentration")
    d["canonical_unit"] = d.concentration_unit.map(lambda x: "mg/mL" if x in MASS_UNITS else "source_unit" if x == "source_unit" else "M")
    d["source_row"] = np.arange(2, len(d) + 2)
    if "observation_id" not in d:
        d["observation_id"] = [f"row-{x}" for x in d.source_row]
    if d.observation_id.str.strip().eq("").any() or d.observation_id.duplicated().any():
        raise ValueError("Missing or duplicate observation_id")
    if "exclude" not in d:
        d["exclude"] = "false"
    if not d.exclude.str.lower().isin(("true", "false")).all():
        raise ValueError("exclude must be true or false")
    d["exclude"] = d.exclude.str.lower().eq("true")
    if "exclusion_reason" not in d:
        d["exclusion_reason"] = ""
    if (d.exclude & d.exclusion_reason.str.strip().eq("")).any():
        raise ValueError("Every exclusion needs a reason")
    for curve, group in d.groupby("curve_id", sort=False):
        for key in ("sample_id", "experiment_id", "concentration_unit", "canonical_unit", "response_unit"):
            if group[key].nunique() != 1:
                raise ValueError(f"Curve {curve} has inconsistent {key}")
    curves = set(d.curve_id)
    for item in cfg["comparisons"]:
        missing = {item["reference_curve"], item["test_curve"]} - curves
        if missing:
            raise ValueError(f"Comparison {item['id']} names unknown curves: {sorted(missing)}")
        units = d.loc[d.curve_id.isin([item["reference_curve"], item["test_curve"]])]
        if units.canonical_unit.nunique() != 1 or units.response_unit.nunique() != 1:
            raise ValueError(f"Comparison {item['id']} mixes dose-unit families or response units")
    return d
