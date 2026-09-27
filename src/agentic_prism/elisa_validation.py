"""Independent control recovery and descriptive cross-plate validation evidence."""
import numpy as np


def independent_controls(group, fit, cfg, factor):
    from .elisa import inverse
    controls=group[(group.role=="qc") & ~group.exclude]
    rule=cfg["independent_qc"]
    result={"status":"not_assessed", "levels":[], "criteria":rule,
            "range_evidence":"No assay-wide LLOQ/ULOQ established by this check"}
    if controls.empty:
        if (group.role=="qc").any():
            result["status"]="failed"
        return result
    if rule["preparation_independent"] is not True or not rule["rationale"].strip():
        result["status"]="independence_unconfirmed"
        return result
    for nominal,rows in controls.groupby("concentration",sort=True):
        vals=[]
        for _,row in rows.iterrows():
            value=inverse(float(row.response),fit)
            q=fit["qc"]
            if value is None or not q["accepted"] or not q["lloq_input_unit"]<=value/factor<=q["uloq_input_unit"]:
                vals.append(None)
            else:
                vals.append(value/factor*float(row.dilution_factor))
        complete=all(v is not None for v in vals)
        mean=float(np.mean(vals)) if complete else None
        recovery=100*mean/nominal if mean is not None else None
        cv=float(100*np.std(vals,ddof=1)/mean) if complete and len(vals)>1 and mean>0 else None
        passed=bool(fit["reportable"] and len(vals)>=rule["min_replicates"] and recovery is not None
                    and rule["recovery_limits_percent"][0]<=recovery<=rule["recovery_limits_percent"][1]
                    and cv is not None and cv<=rule["max_cv_percent"])
        result["levels"].append({"nominal":float(nominal),"unit":str(rows.concentration_unit.iloc[0]),
            "n":len(vals),"mean":mean,"recovery_percent":recovery,"cv_percent":cv,"passed":passed})
    result["status"]="passed" if len(result["levels"])>=rule["min_levels"] and all(x["passed"] for x in result["levels"]) else "failed"
    return result


def across_plates(fits):
    """Equal weight to each plate, no automatic assay validation claim."""
    groups={}
    for fit in fits:
        for row in fit.get("independent_qc",{}).get("levels",[]):
            groups.setdefault((row["nominal"],row["unit"]),[]).append(row)
    summaries=[]
    for (nominal,unit),rows in groups.items():
        vals=[r["mean"] for r in rows if r["passed"]]
        complete=len(vals)==len(rows) and len(vals)>=3
        summaries.append({"nominal":nominal,"unit":unit,"n_plates":len(rows),
            "status":"descriptive_only" if complete else "insufficient_or_failed_controls",
            "mean_recovery_percent":float(100*np.mean(vals)/nominal) if complete else None,
            "between_plate_mean_cv_percent":float(100*np.std(vals,ddof=1)/np.mean(vals)) if complete else None})
    return summaries
