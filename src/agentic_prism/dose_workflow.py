"""Versioned, auditable 4PL analysis independent of rendering."""
from datetime import datetime, timezone
from pathlib import Path
import importlib.metadata
import json
import platform
import shutil
import pandas as pd
from . import __version__
from .dose_schema import resolve_dose_config, load_dose_data
from .dose_fit import fit_dose_curve, summarize_dose, compare_curves, summarize_potency
from .workflow import dump, sha, implementation_hash


def analyze_dose(config_path, output, render=True):
    src, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError("Choose a new dose-response run directory")
    out.mkdir(parents=True)
    try:
        cfg = resolve_dose_config(json.loads(src.read_text()))
        inp = (src.parent / cfg["input"]).resolve()
        d = load_dose_data(inp, cfg)
        shutil.copyfile(inp, out / "input.csv")
        if cfg["provenance"]:
            provenance = (src.parent / cfg["provenance"]).resolve()
            provenance_data = json.loads(provenance.read_text())
            if provenance_data.get("data_sha256") != sha(inp):
                raise ValueError("Input no longer matches source provenance")
            shutil.copyfile(provenance, out / "source_provenance.json")
            cfg["provenance"] = "source_provenance.json"
        cfg["input"] = "input.csv"
        dump(out / "config.resolved.json", cfg)
        d.to_csv(out / "normalized_data.csv", index=False)
        dump(out / "preprocessing_log.json", {"dose_scale": cfg["dose_scale"],
             "unit_conversion": "M, mg/mL, or source_unit only within its unit family; original values retained",
             "normalization": "none", "baseline_subtraction": "none", "response_aggregation": "none",
             "automatic_outlier_removal": False,
             "excluded_observations": d.loc[d.exclude, ["observation_id", "exclusion_reason"]].to_dict("records")})
        fits, grids, observations, groups = [], [], [], {}
        for curve_id, group in d.groupby("curve_id", sort=False):
            fit, grid, rows = fit_dose_curve(group, cfg)
            fits.append(fit)
            grids.extend(grid)
            observations.extend(rows)
            groups[curve_id] = group
        by_curve = {fit["curve_id"]: fit for fit in fits}
        summaries = summarize_dose(fits, cfg)
        comparisons, comparison_grid = [], []
        for item in cfg["comparisons"]:
            comparison, grid = compare_curves(item, groups, by_curve, cfg)
            comparisons.append(comparison)
            comparison_grid.extend(grid)
        potency = summarize_potency(comparisons, cfg)
        dump(out / "results.json", {"schema_version": 1, "analysis_type": "dose_response_4pl", "fits": fits,
                                    "summaries": summaries, "comparisons": comparisons,
                                    "potency_summaries": potency})
        pd.DataFrame(summaries).to_csv(out / "sample_summary.csv", index=False)
        if comparisons:
            pd.DataFrame([{k: v for k, v in c.items() if k not in ("rp_ci", "diagnostics", "parallel_model",
                                                                     "parallelism_f_test", "shared_c50_f_test",
                                                                     "parallelism_equivalence", "rp_acceptance")} |
                          {"rp_ci_low": c["rp_ci"][0], "rp_ci_high": c["rp_ci"][1],
                           "parallelism_p": (c["parallelism_f_test"] or {}).get("p_value"),
                           "shared_c50_p": (c["shared_c50_f_test"] or {}).get("p_value"),
                           "diagnostics": ";".join(c["diagnostics"])} |
                          ({"parallelism_equivalent": c["parallelism_equivalence"]["equivalent"]} if c["parallelism_equivalence"] else {}) |
                          ({"rp_acceptance_within": c["rp_acceptance"]["within"]} if c["rp_acceptance"] else {})
                          for c in comparisons]).to_csv(out / "comparisons.csv", index=False)
            equivalence_rows = [{"comparison_id": c["comparison_id"], **{k: v for k, v in t.items() if k not in ("ci", "limits")},
                                 "ci_low": t["ci"][0], "ci_high": t["ci"][1], "limit_low": t["limits"][0], "limit_high": t["limits"][1],
                                 "confidence_level": c["parallelism_equivalence"]["confidence_level"]}
                                for c in comparisons if c["parallelism_equivalence"] for t in c["parallelism_equivalence"]["tests"]]
            if cfg["comparison_settings"]["parallelism_method"] == "equivalence":
                pd.DataFrame(equivalence_rows).to_csv(out / "parallelism_equivalence.csv", index=False)
            pd.DataFrame(potency).to_csv(out / "potency_summary.csv", index=False)
            pd.DataFrame(comparison_grid, columns=["comparison_id", "curve_id", "concentration_canonical",
                                                   "predicted_response"]).to_csv(out / "comparison_grid.csv", index=False)
        pd.DataFrame([{k: v for k, v in fit.items() if k not in ("diagnostics", "ci_canonical")} |
                      {"diagnostics": ";".join(fit["diagnostics"]), "ci_low_canonical": fit["ci_canonical"][0],
                       "ci_high_canonical": fit["ci_canonical"][1]} for fit in fits]).to_csv(out / "fit_results.csv", index=False)
        pd.DataFrame(grids, columns=["curve_id", "concentration_canonical", "predicted_response"]).to_csv(out / "curve_grid.csv", index=False)
        pd.DataFrame(observations, columns=["observation_id", "curve_id", "concentration_canonical",
                   "observed_response", "predicted_response", "residual", "used_in_fit"]).to_csv(out / "predictions.csv", index=False)
        dump(out / "diagnostics.json", {fit["curve_id"]: {key: fit[key] for key in
             ("status", "reportable", "range_status", "ci_status", "diagnostics")} for fit in fits})
        (out / "rerun.txt").write_text("agentic-prism analyze --config config.resolved.json --output ../dose-rerun-new\n"
                                        "Run from this directory in the pinned environment; choose a fresh output directory.\n")
        dump(out / "manifest.json", {"schema_version": 1, "analysis_type": "dose_response_4pl",
             "package_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
             "input_original_path": str(inp), "input_sha256": sha(inp), "source": cfg["source"],
             "python": platform.python_version(), "platform": platform.platform(),
             "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas", "matplotlib")},
             "implementation_sha256": implementation_hash(),
             "scientific_artifacts_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
             "n_curves": len(fits), "n_reportable": sum(f["reportable"] for f in fits),
             "n_comparisons": len(comparisons),
             "prism_numerical_equivalence": "not_claimed"})
    except Exception as exc:
        dump(out / "failure.json", {"status": "failed", "exception": type(exc).__name__, "message": str(exc)})
        raise
    if render:
        from .dose_report import render_dose
        render_dose(out)
    return out
