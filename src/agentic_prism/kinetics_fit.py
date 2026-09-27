"""1:1 association/dissociation global fit with variable-projected Rmax/offset.

Rates are shared within one declared experiment. Linear nuisance parameters are
solved exactly at each rate proposal. Bootstrap resamples contiguous residual
blocks separately within each trace segment (one association or dissociation
period), and derives KD from each joint refit. A single-cycle curve carries
several injection steps; the analytic 1:1 state is propagated across them.
"""
import hashlib
import numpy as np
from scipy.optimize import least_squares


def unit_response(t, concentration, association_duration, kon, koff):
    t=np.asarray(t,dtype=float)
    rate=kon*concentration+koff
    occupancy=kon*concentration/rate
    bound=occupancy*(-np.expm1(-rate*np.minimum(np.maximum(t,0),association_duration)))
    return bound*np.exp(-koff*np.maximum(t-association_duration,0))


def schedule_response(t, steps, kon, koff):
    """Unit-Rmax 1:1 response for injection steps [(start, end, C)] relative to the first start."""
    if len(steps)==1:
        start,end,c=steps[0]
        return unit_response(np.asarray(t,dtype=float)-start,c,end-start,kon,koff)
    t=np.asarray(t,dtype=float)
    out=np.zeros_like(t)
    state=0.
    for k,(start,end,c) in enumerate(steps):
        nxt=steps[k+1][0] if k+1<len(steps) else np.inf
        rate=kon*c+koff
        eq=kon*c/rate
        mask=(t>=start)&(t<end)
        out[mask]=eq+(state-eq)*np.exp(-rate*(t[mask]-start))
        state=eq+(state-eq)*np.exp(-rate*(end-start))
        mask=(t>=end)&(t<nxt)
        out[mask]=state*np.exp(-koff*(t[mask]-end))
        if np.isfinite(nxt): state*=np.exp(-koff*(nxt-end))
    return out


def _steps(g):
    """Injection steps relative to the curve's first association start."""
    rows=g[["association_start_s","dissociation_start_s","concentration_M"]].drop_duplicates().sort_values("association_start_s")
    first=float(rows.association_start_s.iloc[0])
    return first,[(float(a)-first,float(b)-first,float(c)) for a,b,c in rows.itertuples(index=False)]


def _arrays(group):
    curves=[]
    for cid,g in group.groupby("curve_id",sort=False):
        h=g[g.used_in_fit]
        first,steps=_steps(g)
        t=h.time_s.to_numpy()-first
        # Segment 2k is step k's association and 2k+1 its dissociation; the last is the final dissociation.
        step=np.searchsorted([s[0] for s in steps],t,side="right")-1
        segment=2*step+(h.phase.to_numpy()=="dissociation")
        curves.append({"curve_id":cid,"t":t,"steps":steps,"c":float(g.concentration_M.iloc[0]),
                       "duration":steps[0][1]-steps[0][0],
                       "y":h.response_processed.to_numpy(float),"phase":h.phase.to_numpy(),
                       "segment":segment,"final_dissociation":t>=steps[-1][1],
                       "indices":h.index.to_numpy()})
    return curves


def _project(curves,z,config,ys=None):
    hs=[schedule_response(c["t"],c["steps"],10.**z[0],10.**z[1]) for c in curves]
    ys=[c["y"] for c in curves] if ys is None else ys
    offset=config["fit"]["offset"]=="fitted_per_curve"
    centered_h=[h-h.mean() if offset else h for h in hs]
    centered_y=[y-y.mean() if offset else y for y in ys]
    if config["fit"]["rmax"]=="shared":
        denom=sum(float(h@h) for h in centered_h)
        r=max(0.,sum(float(h@y) for h,y in zip(centered_h,centered_y))/max(denom,1e-300))
        rmax=[r]*len(curves)
    else:
        rmax=[max(0.,float(h@y)/max(float(h@h),1e-300)) for h,y in zip(centered_h,centered_y)]
    offsets=[float(y.mean()-r*h.mean()) if offset else 0. for h,y,r in zip(hs,ys,rmax)]
    predictions=[r*h+b for h,r,b in zip(hs,rmax,offsets)]
    return predictions,rmax,offsets


def _solve(curves,cfg,ys=None,initial=None):
    ys=[c["y"] for c in curves] if ys is None else ys
    scale=max(max(float(np.max(np.abs(y))) for y in ys),1e-12)
    lower=np.log10([cfg["fit"]["kon_bounds"][0],cfg["fit"]["koff_bounds"][0]])
    upper=np.log10([cfg["fit"]["kon_bounds"][1],cfg["fit"]["koff_bounds"][1]])
    def fun(z):
        pp,_,_=_project(curves,z,cfg,ys)
        return np.concatenate([(p-y)/scale for p,y in zip(pp,ys)])
    if initial is None:
        m=int(np.ceil(np.sqrt(cfg["fit"]["multistart"])))
        starts=[np.array([a,b]) for a in np.linspace(3,7,m) for b in np.linspace(-5,-1,m)][:cfg["fit"]["multistart"]]
    else: starts=[initial]
    best=None
    for start in starts:
        sol=least_squares(fun,np.clip(start,lower+1e-7,upper-1e-7),bounds=(lower,upper),
                          ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=1200)
        if sol.success and np.isfinite(sol.fun).all() and (best is None or sol.fun@sol.fun<best.fun@best.fun): best=sol
    if best is None: raise ArithmeticError("Kinetic optimization failed")
    pp,rmax,offsets=_project(curves,best.x,cfg,ys)
    boundary=bool(np.any(np.minimum(best.x-lower,upper-best.x)<1e-4) or any(r<=scale*1e-10 for r in rmax))
    return best,pp,rmax,offsets,boundary


def _resample_blocks(residual,length,rng):
    centered=residual-residual.mean()
    starts=rng.integers(0,len(centered),size=int(np.ceil(len(centered)/length)))
    indices=(starts[:,None]+np.arange(length)[None,:])%len(centered)
    return centered[indices.ravel()[:len(centered)]]


def fit_kinetic_group(group,cfg):
    gid=str(group.fit_group_id.iloc[0])
    result={"fit_group_id":gid,"sample_id":str(group.sample_id.iloc[0]),"experiment_id":str(group.experiment_id.iloc[0]),
            "status":"failed","reportable":False,"model":"one_to_one","kon_M_inv_s_inv":None,"koff_s_inv":None,"kd_M":None,
            "n_curves":int(group.curve_id.nunique()),"n_input":len(group),"n_fit":int(group.used_in_fit.sum()),
            "diagnostics":[],"ci_method":cfg["uncertainty"]["method"],"ci_status":"not_computed", "ci_level":cfg["uncertainty"]["level"],
            "intervals":{},"audit_intervals":{},"ci_reportable":False,"window_sensitivity":[],
            "reliability_blockers":[],"curve_parameters":[],"bootstrap_successes":0}
    curves=_arrays(group)
    def fail(reason):
        result["diagnostics"].append(reason)
        return result,[],[]
    single=cfg["assay"].get("injection_design","multi_cycle")=="single_cycle"
    concentrations={s[2] for c in curves for s in c["steps"]}
    if (not single and len(curves)<2) or len(concentrations)<2:
        return fail("requires_two_distinct_positive_concentrations")
    for c in curves:
        if any(np.count_nonzero(c["phase"]==p)<5 for p in ("association","dissociation")):
            return fail("insufficient_association_or_dissociation_points")
        if np.ptp(c["y"])<=max(float(np.max(np.abs(c["y"]))),1e-12)*1e-9:
            return fail("flat_trace")
    try: sol,pp,rmax,offsets,boundary=_solve(curves,cfg)
    except (ValueError,ArithmeticError): return fail("optimizer_failed")
    kon,koff=map(float,10.**sol.x)
    p=2+(1 if cfg["fit"]["rmax"]=="shared" else len(curves))+(len(curves) if cfg["fit"]["offset"]=="fitted_per_curve" else 0)
    rank=int(np.linalg.matrix_rank(sol.jac))
    sse=sum(float((y-c["y"])@(y-c["y"])) for c,y in zip(curves,pp))
    cond=float(np.linalg.cond(sol.jac))
    result.update(kon_M_inv_s_inv=kon,koff_s_inv=koff,kd_M=koff/kon,optimizer_success=True,
                  objective_sse=sse,residual_df=result["n_fit"]-p,numerical_boundary_hit=boundary,
                  rate_jacobian_rank=rank,rate_jacobian_condition=cond if np.isfinite(cond) else None)
    if boundary: result["diagnostics"].append("numerical_boundary_hit")
    if rank<2: result["diagnostics"].append("rates_not_locally_identifiable")
    if cfg["preprocessing"]["baseline_subtraction"]=="baseline_mean":
        result["diagnostics"].append("intervals_conditional_on_subtracted_baseline")
    if cfg["preprocessing"]["reference_mode"] in ("reference_column","double_reference"):
        result["diagnostics"].append("reference_uncertainty_not_separately_propagated")
    predictions=[]
    for c,pred,rm,b in zip(curves,pp,rmax,offsets):
        acfs=[]
        for phase in ("association","dissociation"):
            # Lag-1 pairs are taken within segments only (one segment per phase for multi-cycle curves).
            pairs=[(rr[:-1],rr[1:]) for seg in np.unique(c["segment"][c["phase"]==phase])
                   for rr in [c["y"][c["segment"]==seg]-pred[c["segment"]==seg]] if len(rr)>1]
            lead=np.concatenate([a for a,_ in pairs]) if pairs else np.zeros(0)
            lag=np.concatenate([b for _,b in pairs]) if pairs else np.zeros(0)
            acf=float(np.corrcoef(lead,lag)[0,1]) if len(lead)>1 and np.std(lead)>1e-14 and np.std(lag)>1e-14 else 0.
            acfs.append(acf if np.isfinite(acf) else 0.)
        diss_t=c["t"][c["final_dissociation"]]
        fractional_decay=float(-np.expm1(-koff*(diss_t.max()-diss_t.min())))
        if fractional_decay<.02: result["diagnostics"].append("little_dissociation_in_observation_window")
        if max(map(abs,acfs))>.5: result["diagnostics"].append("serially_correlated_residuals_review_model_and_block_length")
        for seg in np.unique(c["segment"]):
            dt=np.diff(c["t"][c["segment"]==seg])
            if len(dt)>1 and np.std(dt)>np.mean(dt)*.05:
                result["diagnostics"].append("irregular_sampling_block_duration_varies")
        result["curve_parameters"].append({"curve_id":c["curve_id"],"concentration_M":c["c"],"rmax":float(rm),"offset":float(b),
                                          "association_residual_lag1":acfs[0],"dissociation_residual_lag1":acfs[1],
                                          "fractional_decay_in_fit_window":fractional_decay,
                                          **({"concentration_M":None,"injection_concentrations_M":[s[2] for s in c["steps"]]}
                                             if len(c["steps"])>1 else {})})
        raw=group[group.curve_id==c["curve_id"]]
        first,_=_steps(raw)
        tt=raw.time_s.to_numpy()-first
        full=rm*schedule_response(tt,c["steps"],kon,koff)+b
        for (_,row),yh in zip(raw.iterrows(),full):
            predictions.append({"observation_id":row.observation_id,"fit_group_id":gid,"curve_id":c["curve_id"],
                                "time_s":float(row.time_s),"time_since_association_s":float(row.time_s-first),
                                "phase":row.phase,"concentration_M":float(row.concentration_M),
                                "observed_response":float(row.response_processed),"predicted_response":float(yh),
                                "residual":float(row.response_processed-yh),"used_in_fit":bool(row.used_in_fit)})
    boot=[]
    u=cfg["uncertainty"]
    if u["method"]=="block_bootstrap":
        if sse<=1e-22*max(1.,sum(float(c["y"]@c["y"]) for c in curves)):
            result["ci_status"]="noise_scale_not_estimable"
        elif any(np.count_nonzero(c["segment"]==seg)<2*u["block_length"] for c in curves for seg in np.unique(c["segment"])):
            result["ci_status"]="insufficient_points_for_block_length"
            result["diagnostics"].append("bootstrap_block_length_exceeds_half_phase_length")
        else:
            salt=int(hashlib.sha256(gid.encode()).hexdigest()[:8],16)
            rng=np.random.default_rng(np.random.SeedSequence([u["seed"],salt]))
            for iteration in range(u["replicates"]):
                ys=[]
                for c,pred in zip(curves,pp):
                    y=pred.copy()
                    for seg in np.unique(c["segment"]):
                        mask=c["segment"]==seg
                        y[mask]+=_resample_blocks(c["y"][mask]-pred[mask],u["block_length"],rng)
                    ys.append(y)
                try:
                    bs,_,_,_,bd=_solve(curves,cfg,ys,sol.x)
                    a,b=map(float,10.**bs.x)
                    boot.append({"fit_group_id":gid,"iteration":iteration,"kon_M_inv_s_inv":a,"koff_s_inv":b,"kd_M":b/a,"accepted":not bd})
                except (ValueError,ArithmeticError):
                    boot.append({"fit_group_id":gid,"iteration":iteration,"kon_M_inv_s_inv":None,"koff_s_inv":None,"kd_M":None,"accepted":False})
            accepted=[b for b in boot if b["accepted"]]
            result["bootstrap_successes"]=len(accepted)
            if len(accepted)>=.9*u["replicates"]:
                q=[(1-u["level"])/2,(1+u["level"])/2]
                result["intervals"]={k:list(map(float,np.quantile([b[k] for b in accepted],q))) for k in ("kon_M_inv_s_inv","koff_s_inv","kd_M")}
                result["ci_status"]="conditional_block_bootstrap"
            else:
                result["ci_status"]="insufficient_valid_bootstrap_refits"
    else: result["ci_status"]="not_requested"
    # Deterministic sensitivity checks use earlier portions of the declared
    # dissociation window. These are diagnostics, not additional replicates.
    windows=[]
    for fraction in (.75, .5):
        shortened=[]
        for c in curves:
            times=c["t"][c["final_dissociation"]]
            cutoff=times.min()+fraction*(times.max()-times.min())
            keep=~c["final_dissociation"] | (c["t"]<=cutoff)
            shortened.append({**c, **{key:c[key][keep] for key in ("t","y","phase","indices","segment","final_dissociation")}})
        item={"dissociation_fraction":fraction,"status":"insufficient_points"}
        if all(np.count_nonzero(c["final_dissociation"])>=5 for c in shortened):
            try:
                ws,_,_,_,wb=_solve(shortened,cfg)
                wa,wk=map(float,10.**ws.x)
                ratios={"kon":wa/kon,"koff":wk/koff,"kd":(wk/wa)/(koff/kon)}
                item.update(status="limited" if wb else "estimated", ratios_to_primary=ratios)
                if wb or max(max(v,1/v) for v in ratios.values())>2:
                    result["diagnostics"].append("fit_window_sensitive")
            except (ValueError,ArithmeticError):
                item["status"]="optimizer_failed"
                result["diagnostics"].append("fit_window_check_failed")
        else:
            result["diagnostics"].append("fit_window_check_insufficient_points")
        windows.append(item)
    result["window_sensitivity"]=windows
    blockers={"serially_correlated_residuals_review_model_and_block_length",
              "irregular_sampling_block_duration_varies","little_dissociation_in_observation_window",
              "fit_window_sensitive","fit_window_check_failed","fit_window_check_insufficient_points"}
    result["reliability_blockers"]=sorted(blockers.intersection(result["diagnostics"]))
    result["audit_intervals"]=result["intervals"].copy()
    if result["reliability_blockers"]:
        result["intervals"]={}
        result["ci_status"]="withheld_reliability_check"
    result["ci_reportable"]=bool(result["intervals"]) and not boundary and rank==2
    result["diagnostics"]=list(dict.fromkeys(result["diagnostics"]))
    result["reportable"]=not boundary and rank==2 and result["ci_status"] in ("conditional_block_bootstrap","not_requested")
    result["status"]="estimated_with_diagnostics" if result["reportable"] and result["diagnostics"] else "estimated" if result["reportable"] else "limited"
    return result,predictions,boot
