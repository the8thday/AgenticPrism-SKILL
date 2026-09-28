"""Analysis artifacts are independent of render artifacts and Agent prose."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import shutil
import numpy as np
import pandas as pd
from . import __version__
from .schema import resolve_config, load_data
from .fit import fit_curve, response, summarize, sensitivity


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def implementation_hash():
    root = Path(__file__).parent
    return {p.name: sha(p) for p in sorted(root.glob("*.py"))}


def analyze(config_path, output, render=True):
    try:
        requested = json.loads(Path(config_path).read_text())
    except (OSError, ValueError):
        requested = None
    if isinstance(requested, dict) and requested.get("analysis_type") in ("stability", "potency_assay", "comparability", "specification"):
        from .cmc_workflow import analyze_cmc
        return analyze_cmc(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") == "variance_components":
        from .precision_workflow import analyze_precision
        return analyze_precision(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") == "ada_cut_point":
        from .ada_workflow import analyze_ada
        return analyze_ada(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") == "method_validation":
        from .method_validation_workflow import analyze_method_validation
        return analyze_method_validation(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") in ("ada_sensitivity", "ada_drug_tolerance"):
        from .ada_performance import analyze_ada_performance
        return analyze_ada_performance(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") == "binding_kinetics":
        from .kinetics_workflow import analyze_kinetics
        return analyze_kinetics(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") == "dose_response_4pl":
        from .dose_workflow import analyze_dose
        return analyze_dose(config_path, output, render)
    if isinstance(requested, dict) and requested.get("analysis_type") in ("elisa_quantification", "group_comparison", "multi_group_comparison", "nonparametric", "mmrm", "repeated_measures", "time_to_event", "tumor_growth"):
        from .simple_workflow import analyze_simple
        return analyze_simple(config_path, output, render)
    config_path, out = Path(config_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError(f"Output already exists; choose a new run directory: {out}")
    out.mkdir(parents=True)
    try:
        raw = json.loads(config_path.read_text())
        cfg = resolve_config(raw)
        inp = Path(cfg["input"])
        if not inp.is_absolute():
            inp = config_path.parent / inp
        inp = inp.resolve()
        d = load_data(inp, cfg)
        shutil.copyfile(inp, out / "input.csv")
        cfg["input"] = "input.csv"
        dump(out / "config.resolved.json", cfg)
        d.to_csv(out / "normalized_data.csv", index=False)
        log = {"concentration_conversion": "Original concentration and unit retained; internal M.",
               "excluded_observations": d.loc[d.exclude, ["observation_id", "exclusion_reason"]].to_dict("records"),
               "baseline": cfg["fit"]["baseline_mode"],
               "automatic_outlier_removal": False,
               "response_aggregation_before_fit": "none; supplied rows fit as supplied",
               "replicate_policy": cfg["replicates"], "source": cfg["source"]}
        dump(out / "preprocessing_log.json", log)
        fits, preds, residuals, sens = [], [], [], []
        for cid, g in d.groupby("curve_id", sort=False):
            ft = fit_curve(g, cfg)
            fits.append(ft)
            if not g.loc[~g.exclude].empty:
                sens.extend(sensitivity(g, cfg, ft))
            if ft["kd_M"] is None:
                continue
            positive = g.loc[g.concentration_M > 0, "concentration_M"]
            grid = np.r_[0., np.geomspace(positive.min(), positive.max(), 160)]
            for c, y in zip(grid, response(grid, ft["kd_M"], ft["amplitude"], ft["baseline"])):
                preds.append({"curve_id": cid, "concentration_M": float(c), "predicted_response": float(y)})
            pred = response(g.concentration_M.to_numpy(), ft["kd_M"], ft["amplitude"], ft["baseline"])
            for (_, row), yy in zip(g.iterrows(), pred):
                valid = row.response > 0 and yy > 0
                included = not row.exclude and not (cfg["fit"]["baseline_mode"] == "zero_control" and row.concentration_M == 0)
                residuals.append({"observation_id": row.observation_id, "curve_id": cid,
                                  "concentration_M": float(row.concentration_M), "predicted_response": float(yy),
                                  "residual_linear": float(row.response - yy),
                                  "residual_log10": float(np.log10(row.response) - np.log10(yy)) if valid else None,
                                  "used_in_objective": included})
        dump(out / "results.json", {"schema_version": 1, "fits": fits, "summaries": summarize(fits, cfg)})
        pd.DataFrame([{**f, "diagnostics": ";".join(f["diagnostics"])} for f in fits]).to_csv(out / "fit_results.csv", index=False)
        pd.DataFrame(summarize(fits, cfg)).to_csv(out / "sample_summary.csv", index=False)
        pd.DataFrame(preds, columns=["curve_id", "concentration_M", "predicted_response"]).to_csv(out / "predictions.csv", index=False)
        pd.DataFrame(residuals, columns=["observation_id", "curve_id", "concentration_M", "predicted_response", "residual_linear", "residual_log10", "used_in_objective"]).to_csv(out / "residuals.csv", index=False)
        pd.DataFrame(sens, columns=["curve_id", "scenario", "status", "kd_M", "ratio_to_primary"]).to_csv(out / "sensitivity.csv", index=False)
        dump(out / "diagnostics.json", {f["curve_id"]: {k:f[k] for k in ("status", "range_status", "identifiability", "ci_status", "diagnostics")} for f in fits})
        (out / "rerun.txt").write_text('agentic-prism analyze --config config.resolved.json --output ../rerun-new\n\nRun from this run directory in the pinned environment. Choose a fresh output directory.\n')
        dump(out / "manifest.json", {"schema_version": 1, "package_version": __version__, "created_utc": datetime.now(timezone.utc).isoformat(),
             "input_original_path": str(inp), "input_sha256": sha(out / "input.csv"), "source": cfg["source"],
             "python": platform.python_version(), "platform": platform.platform(),
             "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pandas", "matplotlib")},
             "implementation_sha256": implementation_hash(), "random_seed": None,
             "randomness": "Deterministic multistart grid; no stochastic inference in production.",
             "scientific_artifacts_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
             "n_curves": len(fits), "n_reportable": sum(f["reportable"] for f in fits),
             "prism_numerical_benchmark": "not_completed"})
    except Exception as exc:
        dump(out / "failure.json", {"status": "failed", "exception": type(exc).__name__, "message": str(exc)})
        raise
    if render:
        from .report import render_report
        render_report(out)
    return out


def verify_run(run):
    run = Path(run)
    manifest = json.loads((run / "manifest.json").read_text())
    for name, expected in manifest["scientific_artifacts_sha256"].items():
        if sha(run / name) != expected:
            raise ValueError(f"Scientific artifact changed: {name}")
    return manifest
