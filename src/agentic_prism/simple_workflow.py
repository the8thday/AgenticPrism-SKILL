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
                 (["well_id"] if kind == "elisa_quantification" else ["observation_id"])].to_dict("records"),
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
        (out / "rerun.txt").write_text("agentic-prism analyze --config config.resolved.json --output ../rerun-new\n"
                                       "Run from this directory; choose a new output directory.\n")
        dump(out / "manifest.json", {"schema_version": 1, "analysis_type": kind,
             "package_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
             "input_original_path": str(inp), "input_sha256": sha(inp), "source": cfg["source"],
             "python": platform.python_version(), "platform": platform.platform(),
             "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas", "matplotlib")},
             "implementation_sha256": implementation_hash(),
             **({"random_seed": module.DUNNETT_SEED} if kind == "multi_group_comparison" else {}),
             "scientific_artifacts_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
             "prism_numerical_equivalence": "not_claimed"})
    except Exception as exc:
        dump(out / "failure.json", {"status": "failed", "exception": type(exc).__name__, "message": str(exc)})
        raise
    if render:
        from .simple_report import render_simple
        render_simple(out)
    return out
