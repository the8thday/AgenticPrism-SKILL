from pathlib import Path
from datetime import datetime,timezone
import json
import shutil
import platform
import importlib.metadata
import pandas as pd
from . import __version__
from .workflow import dump,sha,implementation_hash
from .kinetics_schema import resolve_kinetic_config,load_kinetic_data
from .kinetics_fit import fit_kinetic_group


def analyze_kinetics(config_path,output,render=True):
    raw=json.loads(Path(config_path).read_text())
    if raw.get("model", "one_to_one")!="one_to_one":
        from .advanced_kinetics_workflow import analyze_advanced
        return analyze_advanced(config_path,output,render)
    src=Path(config_path).resolve();out=Path(output).resolve()
    if out.exists(): raise FileExistsError("Choose a new kinetic run directory")
    out.mkdir(parents=True)
    try:
        cfg=resolve_kinetic_config(json.loads(src.read_text()))
        inp=(src.parent/cfg["input"]).resolve()
        d,pre=load_kinetic_data(inp,cfg)
        shutil.copyfile(inp,out/"input.csv")
        if cfg["provenance"]:
            provenance=(src.parent/cfg["provenance"]).resolve()
            provenance_data=json.loads(provenance.read_text())
            if provenance_data.get("data_sha256")!=sha(inp): raise ValueError("Imported data no longer matches importer provenance")
            shutil.copyfile(provenance,out/"source_provenance.json")
            cfg["provenance"]="source_provenance.json"
        cfg["input"]="input.csv"
        dump(out/"config.resolved.json",cfg)
        d.to_csv(out/"normalized_data.csv",index=False)
        dump(out/"preprocessing_log.json",{"configuration":cfg["preprocessing"],"curves":pre,"fit_window":cfg["fit"],
             "automatic_smoothing":False,"automatic_outlier_removal":False,"phase_boundaries":"Input metadata, not inferred from trace shape"})
        fits=[];preds=[];boot=[]
        for _,group in d.groupby("fit_group_id",sort=False):
            fit,pp,bb=fit_kinetic_group(group,cfg)
            fits.append(fit);preds.extend(pp);boot.extend(bb)
        dump(out/"results.json",{"schema_version":1,"analysis_type":"binding_kinetics","fits":fits})
        pd.DataFrame([{k:v for k,v in f.items() if k not in ("curve_parameters","intervals","diagnostics","audit_intervals","window_sensitivity","reliability_blockers")}|
                      {"diagnostics":";".join(f["diagnostics"]),**{k+"_ci_"+s:bounds[i] for k,bounds in f["intervals"].items() for i,s in enumerate(("low","high"))}} for f in fits]).to_csv(out/"fit_results.csv",index=False)
        pd.DataFrame([dict(fit_group_id=f["fit_group_id"],**c) for f in fits for c in f["curve_parameters"]]).to_csv(out/"curve_parameters.csv",index=False)
        pd.DataFrame([{"fit_group_id":f["fit_group_id"], "dissociation_fraction":w["dissociation_fraction"],
                       "status":w["status"], **w.get("ratios_to_primary",{})}
                      for f in fits for w in f["window_sensitivity"]]).to_csv(out/"window_sensitivity.csv",index=False)
        pd.DataFrame(preds,columns=["observation_id","fit_group_id","curve_id","time_s","time_since_association_s","phase","concentration_M","observed_response","predicted_response","residual","used_in_fit"]).to_csv(out/"predictions.csv",index=False)
        pd.DataFrame(boot,columns=["fit_group_id","iteration","kon_M_inv_s_inv","koff_s_inv","kd_M","accepted"]).to_csv(out/"bootstrap.csv",index=False)
        dump(out/"diagnostics.json",{f["fit_group_id"]:{k:f[k] for k in ("status","reportable","ci_status","diagnostics")} for f in fits})
        (out/"rerun.txt").write_text('agentic-prism analyze --config config.resolved.json --output ../kinetic-rerun-new\nRun from this run directory in the pinned environment.\n')
        if cfg.get("steady_state", {}).get("enabled"):
            from .steady_state import compute
            dump(out/"steady_state.json", compute(d, fits, cfg))
        from .legacy_facts import write_facts
        write_facts(out)
        dump(out/"manifest.json",{"schema_version":1,"analysis_type":"binding_kinetics","package_version":__version__,
              "created_utc":datetime.now(timezone.utc).isoformat(),"source":cfg["source"],"input_original_path":str(inp),"input_sha256":sha(inp),
              "python":platform.python_version(),"platform":platform.platform(),"random_seed":cfg["uncertainty"]["seed"],
              "implementation_sha256":implementation_hash(),"dependencies":{p:importlib.metadata.version(p) for p in ("numpy","scipy","pandas","matplotlib")},
              "scientific_artifacts_sha256":{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},
              "vendor_software_equivalence":"not_claimed; see scoped public reference benchmark"})
    except Exception as exc:
        dump(out/"failure.json",{"status":"failed","exception":type(exc).__name__,"message":str(exc)})
        raise
    if render:
        from .kinetics_report import render_kinetics
        render_kinetics(out)
    return out
