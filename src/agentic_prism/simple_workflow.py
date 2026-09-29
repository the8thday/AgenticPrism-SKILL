"""Provenance-preserving runs for ELISA calibration and group comparisons."""
from datetime import datetime, timezone
from pathlib import Path
import importlib.metadata
import json
import platform
import shutil
import pandas as pd
from . import __version__
from .workflow import dump, sha, implementation_hash


def analyze_simple(config_path, output, render=True):
    src, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(f"Output already exists: {out}")
    out.mkdir(parents=True)
    try:
        raw = json.loads(src.read_text())
        kind = raw.get("analysis_type")
        if kind == "elisa_quantification":
            from . import elisa as module
        elif kind == "group_comparison":
            from . import groups as module
        elif kind == "multi_group_comparison":
            from . import multigroup as module
        elif kind == "mmrm":
            from . import mmrm as module
        elif kind == "nonparametric":
            from . import nonparametric as module
        elif kind == "repeated_measures":
            from . import repeated as module
        elif kind == "time_to_event":
            from . import time_to_event as module
        elif kind == "tumor_growth":
            from . import tumor_growth as module
        else:
            raise ValueError("Unknown analysis type")
        cfg = module.resolve_config(raw)
        inp = (src.parent / cfg["input"]).resolve()
        d = module.load_data(inp) if kind == "elisa_quantification" else module.load_data(inp, cfg)
        shutil.copyfile(inp, out / "input.csv")
        cfg["input"] = "input.csv"
        dump(out / "config.resolved.json", cfg)
        d.to_csv(out / "normalized_data.csv", index=False)
        dump(out / "preprocessing_log.json", {
            "automatic_outlier_removal": False, "aggregation_before_fit": "none",
            "excluded_observations": d.loc[d.exclude, ["exclusion_reason"] +
                 (["well_id"] if kind == "elisa_quantification" else ["subject_id"] if kind == "time_to_event" else ["animal_id", "day"] if kind == "tumor_growth" else ["observation_id"])].to_dict("records"),
            "elisa_blank_subtraction": "none" if kind == "elisa_quantification" else None})
        if kind == "elisa_quantification":
            fits, wells, summaries = module.quantify(d, cfg)
            result = {"schema_version": 1, "analysis_type": kind, "fits": fits, "wells": wells, "summaries": summaries}
            pd.DataFrame([{"plate_id": f["plate_id"], **row} for f in fits for row in f["dilution_linearity"]], columns=["plate_id", "sample_id", "dilution_factor", "n_wells", "in_range",
                                             "mean_corrected_concentration", "recovery_vs_dilution_mean_percent"]).to_csv(
                out / "dilution_linearity.csv", index=False)
            from .elisa_validation import across_plates
            result["cross_plate_qc"] = across_plates(fits)
            pd.DataFrame([{ "plate_id":f["plate_id"], "status":f["independent_qc"]["status"], **row}
                          for f in fits for row in (f["independent_qc"]["levels"] or [{}])]).to_csv(out / "independent_qc.csv", index=False)
            pd.DataFrame(result["cross_plate_qc"]).to_csv(out / "cross_plate_qc.csv", index=False)
            pd.DataFrame(wells).assign(diagnostics=lambda x: x.diagnostics.map(";".join)).to_csv(out / "unknown_wells.csv", index=False)
            pd.DataFrame(summaries).to_csv(out / "sample_summary.csv", index=False)
            pd.DataFrame([{"plate_id": f["plate_id"], **{k: v for k, v in level.items() if k != "back_calculated_input_unit"},
                           "back_calculated_input_unit": ";".join("" if b is None else repr(b) for b in level["back_calculated_input_unit"]),
                           "limits_percent": "-".join(map(str, level["limits_percent"])),
                           "lloq_input_unit": f["qc"]["lloq_input_unit"], "uloq_input_unit": f["qc"]["uloq_input_unit"],
                           "calibration_accepted": f["qc"]["accepted"]}
                          for f in fits for level in f["qc"]["levels"]]).to_csv(out / "standards_qc.csv", index=False)
            dump(out / "diagnostics.json", {f["plate_id"]: {"status": f["status"], "reportable": f["reportable"],
                                                "calibration_qc_accepted": f["qc"]["accepted"],
                                                "independent_qc_status": f["independent_qc"]["status"],
                                                "diagnostics": f["diagnostics"]} for f in fits})
        elif kind in ("nonparametric", "mmrm"):
            result = module.compare(d, cfg)
            for key, name in (("fits", "omnibus"), ("contrasts", "contrasts"), ("group_summaries", "group_summary")):
                pd.DataFrame(result.get(key, [])).to_csv(out / f"{name}.csv", index=False)
            dump(out / "diagnostics.json", {k: result["fits"][0][k] for k in ("status", "reportable", "diagnostics")})
        elif kind == "repeated_measures":
            result = module.compare(d, cfg)
            columns = {
                "contrasts": ["arm", "condition", "group_a", "group_b", "difference_b_minus_a", "standard_error", "statistic", "df",
                              "p_unadjusted", "p_adjusted", "ci_low", "ci_high", "method", "reportable", "diagnostic"],
                "fixed_effects": ["arm", "group", "estimated_mean", "standard_error"],
                "residuals": ["observation_id", "independent_unit_id", "arm", "group", "marginal_fitted",
                              "conditional_fitted", "random_intercept_blup", "marginal_residual", "conditional_residual"],
                "missingness": ["independent_unit_id", "arm", "group", "observed"],
                "tests": ["effect", "f_statistic", "wald_statistic", "df_numerator", "df_denominator", "p_value",
                          "epsilon_gg", "p_uncorrected", "partial_eta_squared"]}
            two_way = "tests" in result
            if not two_way:
                for key in ("contrasts", "fixed_effects", "residuals", "missingness"):
                    columns[key] = [c for c in columns[key] if c not in ("arm", "condition", "method")]
            tables = (("fits", "omnibus"), ("group_summaries", "group_summary"), ("contrasts", "contrasts"),
                      ("fixed_effects", "fixed_effects"), ("residuals", "residuals"), ("missingness", "missingness"))
            for key, name in tables + ((("tests", "tests"),) if two_way else ()):
                table = pd.DataFrame(result[key], columns=columns.get(key))
                if "diagnostics" in table:
                    table["diagnostics"] = table.diagnostics.map(";".join)
                if "optimizer_attempts" in table:
                    table["optimizer_attempts"] = table.optimizer_attempts.map(json.dumps)
                table.to_csv(out / f"{name}.csv", index=False)
            dump(out / "bootstrap_samples.json", result.pop("bootstrap_samples"))
            dump(out / "diagnostics.json", {k: result["fits"][0][k] for k in ("status", "reportable", "diagnostics")})
        elif kind == "time_to_event":
            result = module.compare(d, cfg)
            for key, name in (("fits", "omnibus"), ("arm_summaries", "arm_summary"), ("km_curves", "km_curves"),
                              ("cox_terms", "cox_terms"), ("ph_test", "ph_test"), ("contrasts", "contrasts")):
                table = pd.DataFrame(result[key])
                if "diagnostics" in table:
                    table["diagnostics"] = table.diagnostics.map(";".join)
                table.to_csv(out / f"{name}.csv", index=False)
            pd.DataFrame([{k: (";".join(map(str, v)) if isinstance(v, list) else v) for k, v in result["logrank"].items()}]).to_csv(
                out / "logrank.csv", index=False)
            dump(out / "diagnostics.json", {k: result["fits"][0][k] for k in ("status", "reportable", "diagnostics")})
        elif kind == "tumor_growth":
            result = module.compare(d, cfg)
            for key, name in (("fits", "omnibus"), ("arm_day_summaries", "arm_day_summary"), ("growth_rates", "growth_rates"),
                              ("growth_tests", "growth_tests"), ("model_contrasts", "model_contrasts"), ("observed_tgi", "observed_tgi"),
                              ("dropout", "dropout")):
                table = pd.DataFrame(result[key])
                if "diagnostics" in table:
                    table["diagnostics"] = table.diagnostics.map(";".join)
                table.to_csv(out / f"{name}.csv", index=False)
            result["contrasts"] = result["model_contrasts"]
            dump(out / "diagnostics.json", {k: result["fits"][0][k] for k in ("status", "reportable", "diagnostics")})
        elif kind == "multi_group_comparison":
            omnibus, summaries, family = module.compare(d, cfg)
            result = {"schema_version": 1, "analysis_type": kind, "fits": [omnibus],
                      "group_summaries": summaries, "contrasts": family}
            pd.DataFrame([omnibus]).assign(groups=lambda x: x.groups.map(";".join),
                                           diagnostics=lambda x: x.diagnostics.map(";".join)).to_csv(out / "omnibus.csv", index=False)
            pd.DataFrame(summaries).to_csv(out / "group_summary.csv", index=False)
            pd.DataFrame(family, columns=["group_a", "group_b", "difference_b_minus_a", "statistic", "df", "p_unadjusted",
                                          "p_adjusted", "ci_low", "ci_high"]).to_csv(out / "contrasts.csv", index=False)
            dump(out / "diagnostics.json", {"status": omnibus["status"], "diagnostics": omnibus["diagnostics"]})
        else:
            comparison = module.compare(d, cfg)
            result = {"schema_version": 1, "analysis_type": kind, "fits": [comparison]}
            pd.DataFrame([comparison]).assign(ci_low=lambda x: x.ci_difference.map(lambda z: z[0]),
                                              ci_high=lambda x: x.ci_difference.map(lambda z: z[1]),
                                              diagnostics=lambda x: x.diagnostics.map(";".join)).drop(columns="ci_difference").to_csv(
                                                  out / "comparison.csv", index=False)
            dump(out / "diagnostics.json", {"status": comparison["status"], "diagnostics": comparison["diagnostics"]})
        dump(out / "results.json", result)
        from .interpretation import SUPPORTED, write_facts
        if kind in SUPPORTED:
            write_facts(out)
        elif kind == "elisa_quantification":
            from .legacy_facts import write_facts as legacy_facts
            legacy_facts(out, all_states=True)
        elif kind in ("group_comparison", "multi_group_comparison"):
            from .legacy_facts import write_facts as legacy_facts
            legacy_facts(out)
        (out / "rerun.txt").write_text("agentic-prism analyze --config config.resolved.json --output ../rerun-new\n"
                                       "Run from this directory; choose a new output directory.\n")
        dump(out / "manifest.json", {"schema_version": 1, "analysis_type": kind,
             "package_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
             "input_original_path": str(inp), "input_sha256": sha(inp), "source": cfg["source"],
             "python": platform.python_version(), "platform": platform.platform(),
             "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas", "matplotlib") + (("statsmodels",) if kind == "repeated_measures" else ())},
             **({"r_environment": result["r_environment"]} if kind == "mmrm" and result.get("r_environment") else {}),
             "implementation_sha256": implementation_hash(),
             **({"random_seed": module.DUNNETT_SEED} if kind == "multi_group_comparison" else
                {"random_seed": cfg["comparison"]["seed"]} if kind == "time_to_event" and cfg["comparison"]["logrank_inference"] == "permutation" else {}),
             "scientific_artifacts_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
             "prism_numerical_equivalence": "not_claimed"})
    except Exception as exc:
        dump(out / "failure.json", {"status": "failed", "exception": type(exc).__name__, "message": str(exc)})
        raise
    if render:
        from .simple_report import render_simple
        render_simple(out)
    return out
