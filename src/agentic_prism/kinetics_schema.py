"""Explicit contracts for single-site kinetic sensorgrams (independent cycles or single-cycle series)."""
from copy import deepcopy
import json
import numpy as np
import pandas as pd
from .schema import UNITS

DEFAULTS = {
    "schema_version": 1, "analysis_type": "binding_kinetics", "model": "one_to_one",
    "input": None, "source": "User-supplied sensorgrams", "column_map": {}, "provenance": None,
    # injection_design single_cycle: sequential increasing injections without regeneration on one surface.
    "assay": {"one_to_one_supported": None, "independent_cycles": None,
              "concentration_known": None, "rationale": "", "injection_design": "multi_cycle"},
    "preprocessing": {"reference_mode": "none", "baseline_subtraction": "none", "stride": 1},
    "fit": {"rmax": "per_curve", "offset": "fixed_zero", "weighting": "unweighted",
            "kon_bounds": [1., 1e9], "koff_bounds": [1e-8, 10.], "multistart": 9,
            "association_skip_s": 0., "dissociation_skip_s": 0., "dissociation_duration_s": None},
    "uncertainty": {"method": "block_bootstrap", "level": .95, "replicates": 200,
                    "block_length": 20, "seed": 20260922},
    "report": {"language": "zh-CN", "plot_style": "prism_like", "figure_width_mm": 180, "png_dpi": 300},
}


def resolve_kinetic_config(raw):
    def merge(base, given, path=""):
        if not isinstance(given, dict):
            raise ValueError(f"{path} must be an object")
        if set(given) - set(base):
            raise ValueError(f"Unsupported kinetic settings at {path}: {sorted(set(given)-set(base))}")
        out = deepcopy(base)
        for k,v in given.items():
            out[k] = merge(base[k],v,path+k+".") if isinstance(base[k],dict) and k != "column_map" else v
        return out
    c = merge(DEFAULTS,raw)
    if c["schema_version"] != 1 or c["analysis_type"] != "binding_kinetics" or c["model"] != "one_to_one":
        raise ValueError("Unsupported kinetics schema or model")
    a = c["assay"]
    if a["injection_design"] not in ("multi_cycle", "single_cycle"):
        raise ValueError("assay.injection_design must be multi_cycle or single_cycle")
    # Multi-cycle requires regenerated independent cycles; single-cycle must declare that it has none.
    cycles_ok = a["independent_cycles"] is True if a["injection_design"] == "multi_cycle" else a["independent_cycles"] is False
    if any(a[k] is not True for k in ("one_to_one_supported", "concentration_known")) or not cycles_ok or (not isinstance(a["rationale"], str) or not a["rationale"].strip()):
        raise ValueError("Document 1:1 applicability, cycle design (independent_cycles true for multi_cycle, "
                         "false for single_cycle), known concentration and rationale")
    pre = c["preprocessing"]
    if pre["reference_mode"] not in ("none", "reference_column", "double_reference") or pre["baseline_subtraction"] not in ("none", "baseline_mean"):
        raise ValueError("Unsupported reference/baseline preprocessing")
    if type(pre["stride"]) is not int or pre["stride"] < 1:
        raise ValueError("stride must be a positive integer")
    f = c["fit"]
    if f["rmax"] not in ("shared", "per_curve") or f["offset"] not in ("fixed_zero", "fitted_per_curve") or f["weighting"] != "unweighted":
        raise ValueError("Unsupported kinetic parameter sharing or weighting")
    for key in ("kon_bounds", "koff_bounds"):
        b = f[key]
        if not isinstance(b,list) or len(b)!=2 or not all(np.isfinite(b)) or not 0 < b[0] < b[1]:
            raise ValueError(f"Invalid {key}")
    if type(f["multistart"]) is not int or not 1 <= f["multistart"] <= 25:
        raise ValueError("multistart must be 1..25")
    for k in ("association_skip_s", "dissociation_skip_s"):
        if not isinstance(f[k],(int,float)) or not np.isfinite(f[k]) or f[k]<0:
            raise ValueError(f"Invalid {k}")
    dur=f["dissociation_duration_s"]
    if dur is not None and (not isinstance(dur,(int,float)) or not np.isfinite(dur) or dur<=f["dissociation_skip_s"]):
        raise ValueError("Dissociation duration must exceed skip")
    u=c["uncertainty"]
    if u["method"] not in ("none", "block_bootstrap") or not .5<u["level"]<1:
        raise ValueError("Unsupported uncertainty method/level")
    if type(u["replicates"]) is not int or not 50<=u["replicates"]<=2000 or type(u["block_length"]) is not int or u["block_length"]<2 or type(u["seed"]) is not int or u["seed"]<0:
        raise ValueError("Bootstrap requires 50..2000 replicates, block_length >=2 and nonnegative integer seed")
    r=c["report"]
    if r["language"]!="zh-CN" or r["plot_style"] not in ("standard","prism_like") or r["figure_width_mm"] not in (85,180) or r["png_dpi"] not in (300,600):
        raise ValueError("Unsupported kinetic report option")
    if not isinstance(c["input"],str) or not c["input"] or not isinstance(c["column_map"],dict):
        raise ValueError("input path and column_map object required")
    json.dumps(c, allow_nan=False)
    return c


def load_kinetic_data(path, cfg):
    d=pd.read_csv(path,dtype=str,keep_default_na=False).rename(columns=cfg["column_map"])
    req={"fit_group_id","sample_id","experiment_id","curve_id","time","time_unit","response","response_unit",
         "concentration","concentration_unit","phase","association_start","dissociation_start"}
    if req-set(d) or d.empty or d.columns.duplicated().any():
        raise ValueError(f"Invalid kinetic table; missing columns {sorted(req-set(d))}")
    for k in ("fit_group_id","sample_id","experiment_id","curve_id","time_unit","response_unit","concentration_unit","phase"):
        if d[k].str.strip().eq("").any(): raise ValueError(f"Missing {k}")
    for k in ("time","response","concentration","association_start","dissociation_start"):
        d[k]=pd.to_numeric(d[k],errors="raise")
        if not np.isfinite(d[k]).all(): raise ValueError(f"Nonfinite {k}")
    if not d.time_unit.isin(["s","min"]).all() or not d.concentration_unit.isin(UNITS).all():
        raise ValueError("Supported time units s/min; concentrations must use explicit molar units")
    if (d.concentration<=0).any():
        raise ValueError("Only positive analyte traces are fitted; use a matched reference column for blank subtraction")
    if not d.phase.isin(["baseline","association","dissociation"]).all():
        raise ValueError("Unknown phase; use baseline, association or dissociation")
    single=cfg["assay"]["injection_design"]=="single_cycle"
    fac=d.time_unit.map({"s":1.,"min":60.})
    d["time_s"]=d.time*fac
    d["association_start_s"]=d.association_start*fac
    d["dissociation_start_s"]=d.dissociation_start*fac
    d["concentration_M"]=d.concentration*d.concentration_unit.map(UNITS)
    d["source_row"]=np.arange(2,len(d)+2)
    if "observation_id" not in d: d["observation_id"]=[f"row-{i}" for i in d.source_row]
    if d.observation_id.duplicated().any() or d.observation_id.str.strip().eq("").any(): raise ValueError("Nonunique/missing observation_id")
    if "exclude" not in d: d["exclude"]="false"
    if not d.exclude.str.lower().isin(["true","false"]).all(): raise ValueError("exclude must be true/false")
    d["exclude"]=d.exclude.str.lower().eq("true")
    if "exclusion_reason" not in d: d["exclusion_reason"]=""
    if (d.exclude & d.exclusion_reason.str.strip().eq("")).any(): raise ValueError("Missing exclusion reason")
    d["response_processed"]=d.response.astype(float)
    mode=cfg["preprocessing"]["reference_mode"]
    columns={"none":[],"reference_column":["reference_response"],
             "double_reference":["reference_response","blank_response","blank_reference_response"]}[mode]
    for k in columns:
        if k not in d: raise ValueError(f"{k} column is missing for reference_mode={mode}")
        d[k]=pd.to_numeric(d[k],errors="raise")
        if not np.isfinite(d[k]).all(): raise ValueError(f"Nonfinite {k}")
    if mode!="none":
        d["response_processed"]-=d.reference_response
    if mode=="double_reference":
        # Buffer-cycle (blank) trace on the same active/reference pair, matched row by row.
        d["response_processed"]-=d.blank_response-d.blank_reference_response
    d["baseline_subtracted"]=0.
    d["used_in_fit"]=False
    log=[]
    for cid,g in d.groupby("curve_id",sort=False):
        per_curve=("fit_group_id","sample_id","experiment_id","response_unit")+(() if single else ("concentration_M","association_start_s","dissociation_start_s"))
        for k in per_curve:
            if g[k].nunique()!=1: raise ValueError(f"Curve {cid} has inconsistent {k}")
        if g.time_s.duplicated().any() or (np.diff(g.time_s.to_numpy())<=0).any():
            raise ValueError(f"Curve {cid}: times must strictly increase with no duplicates")
        if single:
            check_single_cycle_steps(cid,g)
        else:
            start,end=g.association_start_s.iloc[0],g.dissociation_start_s.iloc[0]
            if end<=start: raise ValueError("Dissociation start must follow association start")
            expected=np.where(g.time_s<start,"baseline",np.where(g.time_s<end,"association","dissociation"))
            if not np.array_equal(expected,g.phase.to_numpy()): raise ValueError(f"Curve {cid}: phases disagree with explicit boundaries")
        start,end=g.association_start_s.to_numpy(),g.dissociation_start_s.to_numpy()
        if cfg["preprocessing"]["baseline_subtraction"]=="baseline_mean":
            base=g[(g.phase=="baseline")&~g.exclude].response_processed
            if len(base)<2: raise ValueError("baseline_mean requires >=2 included pre-association observations")
            d.loc[g.index,"response_processed"]-=float(base.mean())
            d.loc[g.index,"baseline_subtracted"]=float(base.mean())
        # Stride is explicit, anchored to the first non-baseline observation, before fit windows.
        idx=g.loc[g.phase!="baseline"].index[::cfg["preprocessing"]["stride"]]
        d.loc[idx,"used_in_fit"]=True
        f=cfg["fit"]
        eligible=((g.phase=="association") & (g.time_s>=start+f["association_skip_s"])) | ((g.phase=="dissociation") & (g.time_s>=end+f["dissociation_skip_s"]))
        if f["dissociation_duration_s"] is not None: eligible &= g.time_s<=end+f["dissociation_duration_s"]
        d.loc[g.index,"used_in_fit"] &= eligible & ~g.exclude
        log.append({"curve_id":cid,"baseline_subtracted":float(d.loc[g.index,"baseline_subtracted"].iloc[0]),"input_rows":len(g),"fitted_rows":int(d.loc[g.index,"used_in_fit"].sum())})
    for gid,g in d.groupby("fit_group_id",sort=False):
        if any(g[k].nunique()!=1 for k in ("sample_id","experiment_id","response_unit")):
            raise ValueError(f"Fit group {gid} cannot pool samples, independent experiments or response units")
    return d,log


def single_cycle_steps(g):
    """Ordered (association_start_s, dissociation_start_s, concentration_M) steps of one curve."""
    steps=g[["association_start_s","dissociation_start_s","concentration_M"]].drop_duplicates()
    return [tuple(map(float,row)) for row in steps.sort_values("association_start_s").itertuples(index=False)]


def check_single_cycle_steps(cid,g):
    """Each row carries its own step's boundaries and concentration; phases must follow from time."""
    steps=single_cycle_steps(g)
    starts=[s[0] for s in steps]
    if len(set(starts))!=len(steps):
        raise ValueError(f"Curve {cid}: one association start maps to several dissociation starts or concentrations")
    if len(steps)<2 or len({s[2] for s in steps})<2:
        raise ValueError(f"Curve {cid}: single-cycle series needs at least two injections with distinct concentrations")
    for (a,b,_),nxt in zip(steps,starts[1:]+[np.inf]):
        if not a<b<=nxt: raise ValueError(f"Curve {cid}: injection steps must not overlap (association < dissociation <= next association)")
    t=g.time_s.to_numpy();a=g.association_start_s.to_numpy();b=g.dissociation_start_s.to_numpy()
    following=np.array([dict(zip(starts,starts[1:]+[np.inf]))[x] for x in a])
    first=starts[0]
    expected=np.where(t<first,"baseline",np.where(t<a,"mismatch",np.where(t<b,"association",np.where(t<following,"dissociation","mismatch"))))
    if (expected=="baseline").any() and (a[expected=="baseline"]!=first).any():
        raise ValueError(f"Curve {cid}: baseline rows must carry the first injection's metadata")
    if not np.array_equal(expected,g.phase.to_numpy()):
        raise ValueError(f"Curve {cid}: phases or step metadata disagree with explicit injection boundaries")
