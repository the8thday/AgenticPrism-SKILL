"""Adapter for an explicitly identified Octet Time1/Data1 Results.txt export.

This does not parse .frd files, native projects, or infer phase boundaries.
"""
from pathlib import Path
import csv
import json
import shutil
import numpy as np
import pandas as pd
from .workflow import dump,sha
from .schema import UNITS


def read_octet_results(path):
    lines=Path(path).read_text(encoding="utf-8-sig").splitlines()
    header={}
    for i,line in enumerate(lines):
        fields=line.split("\t")
        if fields[:2]==["Time1","Data1"]:
            break
        if len(fields)>=2: header[fields[0]]=fields[1]
    else: raise ValueError(f"Unsupported Octet text layout: {path}")
    if not all(k in header for k in ("Conc1","Start1","Stop1")):
        raise ValueError("Octet export must contain Conc1, Start1, Stop1")
    values=[]
    for line in lines[i+1:]:
        if not line.strip(): continue
        cols=line.split("\t")
        if len(cols)<2: raise ValueError("Malformed Octet observation")
        values.append((float(cols[0]),float(cols[1])))
    d=pd.DataFrame(values,columns=["time","response"])
    if d.empty or not np.isfinite(d.to_numpy()).all() or not np.all(np.diff(d.time)>0):
        raise ValueError("Octet observations must be finite and time ordered")
    if not np.isfinite(float(header["Conc1"])) or float(header["Conc1"])<0:
        raise ValueError("Invalid Octet concentration header")
    if not np.isclose(d.time.iloc[0],float(header["Start1"])) or not np.isclose(d.time.iloc[-1],float(header["Stop1"])):
        raise ValueError("Octet start/stop headers disagree with recorded time coverage")
    return d,header


def import_octet(manifest_path,output):
    src=Path(manifest_path).resolve();out=Path(output).resolve()
    if out.exists(): raise FileExistsError("Choose a new import directory")
    raw=json.loads(src.read_text())
    allowed={"format","source","time_unit","concentration_unit","response_unit","traces"}
    if set(raw)-allowed or raw.get("format")!="octet_results_txt": raise ValueError("Unsupported Octet manifest")
    if raw.get("time_unit") not in ("s","min") or raw.get("concentration_unit") not in UNITS or not raw.get("response_unit") or not raw.get("source"):
        raise ValueError("Manifest must explicitly supply time, concentration, response units and source")
    if not isinstance(raw.get("traces"),list) or not raw["traces"]: raise ValueError("No traces")
    required={"file","fit_group_id","sample_id","experiment_id","curve_id","sensor_id","cycle_id","association_start","dissociation_start"}
    rows=[];files=[];ids=set()
    for trace in raw["traces"]:
        optional={"reference_file","blank_file","blank_reference_file"}
        if required-set(trace) or set(trace)-(required|optional): raise ValueError("Invalid trace metadata")
        if ("blank_file" in trace)!=("blank_reference_file" in trace) or ("blank_file" in trace and "reference_file" not in trace):
            raise ValueError("Double referencing needs reference_file, blank_file and blank_reference_file together")
        if trace["curve_id"] in ids: raise ValueError("Duplicate curve_id in import manifest")
        ids.add(trace["curve_id"])
        a,b=float(trace["association_start"]),float(trace["dissociation_start"])
        if not np.isfinite([a,b]).all() or b<=a: raise ValueError("Invalid stage boundaries")
        path=(src.parent/trace["file"]).resolve()
        d,header=read_octet_results(path)
        if not (d.time.min()<=a<b<d.time.max()): raise ValueError("Export does not cover the declared association and dissociation")
        matched={}
        for key,column,role in (("reference_file","reference_response","reference"),("blank_file","blank_response","blank"),
                                ("blank_reference_file","blank_reference_response","blank_reference")):
            if key not in trace: continue
            rp=(src.parent/trace[key]).resolve()
            ref,_=read_octet_results(rp)
            if len(ref)!=len(d) or not np.allclose(ref.time,d.time,rtol=0,atol=1e-9):
                raise ValueError(f"{role} times must match exactly; no implicit interpolation")
            matched[column]=ref.response.to_numpy()
            files.append({"original_path":str(rp),"sha256":sha(rp),"role":role})
        for i,row in d.iterrows():
            rec={k:trace[k] for k in required-{"file"}}
            rec.update(observation_id=f"{trace['curve_id']}-{i+1}",time=float(row.time),response=float(row.response),
                       concentration=float(header["Conc1"]),phase="baseline" if row.time<a else "association" if row.time<b else "dissociation",
                       time_unit=raw["time_unit"],concentration_unit=raw["concentration_unit"],response_unit=raw["response_unit"],
                       source_file=path.name,source_data_row=i+1)
            for column,values in matched.items(): rec[column]=float(values[i])
            rows.append(rec)
        files.append({"original_path":str(path),"sha256":sha(path),"role":"sample","export_header":header})
    out.mkdir(parents=True)
    (out/"source_files").mkdir()
    for i,file in enumerate(files):
        dest=out/"source_files"/f"{i:03d}_{Path(file['original_path']).name}"
        shutil.copyfile(file["original_path"],dest)
        file["snapshot_path"]=str(dest.relative_to(out))
    pd.DataFrame(rows).to_csv(out/"data.csv",index=False)
    dump(out/"import_manifest.json",{"adapter":"octet_results_txt_v1","source":raw["source"],"configuration":raw,
          "input_files":files,"data_sha256":sha(out/"data.csv"),"fit_columns_ignored":True,
          "preprocessing":"No reference/blank/baseline subtraction, smoothing or fitting applied by importer; original Data1 retained. "
                         "Matched reference/blank traces are carried as separate columns for explicit reference_mode selection."})
    dump(out/"config.template.json",{"analysis_type":"binding_kinetics","input":"data.csv","provenance":"import_manifest.json",
          "source":raw["source"],"assay":{"one_to_one_supported":None,"independent_cycles":None,"concentration_known":True,"rationale":"Supply actual experimental justification before fitting."}})
    return out
