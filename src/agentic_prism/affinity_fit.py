"""Opt-in affinity fitting with projected linear nuisance terms and profile-F."""
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, brentq
from scipy.stats import f as f_dist
from .binding_models import complex_concentration, free_fraction, competition_fraction, regime


def fit_affinity(g,cfg):
    h=g.loc[~g.exclude].copy(); model=cfg['model']; cell=cfg['analysis_type']=='cell_binding'
    result={'curve_id':str(g.fit_group_id.iloc[0]),'sample_id':str(g.sample_id.iloc[0]),
        'experiment_id':str(g.experiment_id.iloc[0]),'response_unit':str(g.response_unit.iloc[0]),
        'model':model,'interpretation':'apparent KD' if cell else 'Ki' if model=='competition_exact' else 'KD',
        'status':'failed','reportable':False,'kd_M':None,'ki_M':None,'audit_estimate_M':None,
        'ci_low_M':None,'ci_high_M':None,'supported_upper_bound_M':None,'ci_status':'not_computed',
        'ci_method':'profile_f','ci_level':cfg['uncertainty']['level'],'diagnostics':[],'design_warnings':[],
        'n_total':len(g),'n_fit':len(h),'profile':[],'parameters':[], 'predictions':[]}
    def fail(reason): result['diagnostics'].append(reason); return result
    x=h.concentration_M.to_numpy(float); y=h.response.to_numpy(float)
    if not len(y) or not np.any(x>0): return fail('no_positive_concentration_or_observations')
    if np.ptp(y)<=max(np.max(np.abs(y)),1e-300)*1e-10: return fail('flat_response')
    curves=list(dict.fromkeys(h.curve_id)); groups=[h.curve_id.to_numpy()==c for c in curves]
    scale=max(float(np.max(np.abs(y))),float(np.ptp(y)),1e-300); ys=y/scale
    fitted_pt=not cell and cfg['active_sites']['mode']=='fitted'
    separate=cell and cfg['nonspecific']['baseline']=='separate'
    n_nonlin=1+int(fitted_pt); n_linear=len(curves)*((4 if separate else 3) if cell else 2); df=len(x)-n_linear-n_nonlin
    if df<=0 or any(len(np.unique(x[m]))<4 for m in groups): return fail('insufficient_distinct_concentrations_or_df')
    if cell:
        s=cfg['receptors']; pt=s['cells_per_mL']*1000*s['receptors_per_cell']/6.02214076e23 if s['depletion']=='quadratic' else 0.
        total=(h.series.to_numpy()=='total').astype(float)
    else: pt=cfg['active_sites']['concentration_M']
    weights=h.response_sd.to_numpy(float)/scale if cfg['fit']['weighting']=='inverse_sd' else np.ones(len(y))
    klo,khi=np.log10(cfg['fit']['kd_bounds_M']); bounds=[(klo,khi)]
    if fitted_pt: bounds.append(tuple(np.log10(cfg['active_sites']['bounds_M'])))
    xscale=float(x.max())
    def evaluate(z, full=False):
        kd=10.**float(z[0]); p=10.**float(z[1]) if fitted_pt else pt
        pred=np.zeros(len(y)); params=[]
        for cid,mask in zip(curves,groups):
            xx=x[mask]
            if cell:
                occupancy=complex_concentration(p,xx,kd)/p if p else xx/(kd+xx)
                occ=occupancy*total[mask]
            elif model=='solution_equilibrium_titration': occ=free_fraction(h.loc[mask,'constant_species_M'].to_numpy(float),xx,kd)
            elif model=='competition_exact':
                q=cfg['competition']; occ=competition_fraction(xx,p,q['tracer_total_M'],q['tracer_kd_M'],kd)
            else: occ=complex_concentration(p,xx,kd)/p
            design=np.column_stack([np.ones(len(xx)),occ]+([xx/xscale] if cell else [])+([1.-total[mask]] if separate else []))
            yy=ys[mask]; w=weights[mask]
            beta=np.linalg.lstsq(design/w[:,None],yy/w,rcond=None)[0]
            if cfg['fit']['weighting']=='relative':
                if np.any(yy<=0): raise ValueError('Relative weighting requires positive background-handled responses')
                # Same predicted-response denominator convention as dose-response.
                def rr(b):
                    pp=design@b
                    return (pp-yy)/np.maximum(pp,1e-12)
                opt=least_squares(rr,beta,max_nfev=cfg['fit']['max_nfev'],ftol=1e-11,xtol=1e-11,gtol=1e-11)
                beta=opt.x
            pred[mask]=design@beta
            params.append({'curve_id':str(cid),'baseline':float(beta[0]*scale),'amplitude':float(beta[1]*scale),
                'nonspecific_slope_per_M':float(beta[2]*scale/xscale) if cell else None,
                **({'control_baseline_offset':float(beta[3]*scale) if separate else None} if cell else {})})
        residual=(pred-ys)/np.maximum(pred,1e-12) if cfg['fit']['weighting']=='relative' else (pred-ys)/weights
        return (residual,pred*scale,params,p) if full else residual
    lo=np.array([b[0] for b in bounds]); hi=np.array([b[1] for b in bounds])
    starts=np.linspace(max(klo+.01,np.log10(x[x>0].min())-2),min(khi-.01,np.log10(x.max())+2),cfg['fit']['multistart'])
    sols=[]
    for k in starts:
        pts=[np.log10(max(x.max()/10,1e-30))] if fitted_pt else [None]
        if fitted_pt: pts+=[np.log10(max(np.median(x[x>0]),1e-30))]
        for p0 in pts:
            z=np.array([k]+([p0] if fitted_pt else [])); z=np.clip(z,lo+1e-6,hi-1e-6)
            opt=least_squares(evaluate,z,bounds=(lo,hi),max_nfev=cfg['fit']['max_nfev'],ftol=1e-11,xtol=1e-11,gtol=1e-11)
            if opt.success and np.isfinite(opt.fun).all(): sols.append(opt)
    if not sols: return fail('optimizer_failed')
    best=min(sols,key=lambda s:float(s.fun@s.fun)); ss=float(best.fun@best.fun)
    residual,pred,params,p=evaluate(best.x,True); estimate=float(10.**best.x[0])
    result.update(audit_estimate_M=estimate,objective_sse=ss,residual_df=df,parameters=params,
        objective_scale=scale,weighting=cfg['fit']['weighting'],pt_M=float(p) if p is not None else None,
        min_positive_M=float(x[x>0].min()),max_concentration_M=float(x.max()),optimizer_success=True)
    if model!='competition_exact':
        ratios=[regime(v,estimate) for v in sorted(set(h.constant_species_M))] if model=='solution_equilibrium_titration' else [regime(p,estimate)]
        result['regime_diagnostics']=ratios
    # Non-blocking design check: a titration that never exceeds the constant
    # species cannot show the equivalence break, so KD rests on sub-saturation
    # curvature alone. Operational factor 2; reportability is unchanged.
    for cid,mask in zip(curves,groups):
        const=float(h.loc[mask,'constant_species_M'].iloc[0]) if model=='solution_equilibrium_titration' else \
            None if model=='competition_exact' or not p else float(p)
        if const and float(x[mask].max())<2*const:
            result['design_warnings'].append({'code':'titrant_max_below_twice_constant_species','curve_id':str(cid),
                'max_titrant_M':float(x[mask].max()),'constant_species_M':const,'operational_factor':2.})
    for row,yy in zip(h.to_dict('records'),pred):
        result['predictions'].append({'observation_id':row['observation_id'],'curve_id':row['curve_id'],
            'concentration_M':row['concentration_M'],'series':row.get('series','binding'),'response':row['response'],'predicted_response':float(yy),'residual':float(row['response']-yy)})
    threshold=ss*(1+float(f_dist.ppf(cfg['uncertainty']['level'],1,df))/df)
    cache={}
    def profile(k):
        key=float(k)
        if key not in cache:
            if fitted_pt:
                # Global scalar scan brackets the nuisance minimum; do not rely on one local Pt start.
                grid=np.linspace(lo[1],hi[1],31)
                vals=[float(evaluate([k,v])@evaluate([k,v])) for v in grid]
                j=int(np.argmin(vals)); left=grid[max(0,j-1)]; right=grid[min(len(grid)-1,j+1)]
                opt=minimize_scalar(lambda v:float(evaluate([k,v])@evaluate([k,v])),bounds=(left,right),method='bounded',options={'xatol':1e-9})
                val=float(opt.fun); nuisance=float(10.**opt.x)
            else:
                rr=evaluate([k]); val=float(rr@rr); nuisance=float(p) if p is not None else None
            cache[key]={'log10_parameter_M':key,'sse':val,'profiled_pt_M':nuisance}
        return cache[key]['sse']-threshold
    def endpoint(direction):
        prev=float(best.x[0]); boundary=klo if direction<0 else khi
        for distance in np.arange(.25,khi-klo+.25,.25):
            now=max(boundary,prev-distance) if direction<0 else min(boundary,prev+distance)
            if profile(now)>=0:
                return float(10.**brentq(profile,min(now,prev),max(now,prev),xtol=1e-8))
            if now==boundary: return None
        return None
    if ss<1e-24:
        result['ci_status']='noise_scale_not_estimable'; result['diagnostics'].append('noise_scale_not_estimable')
    else:
        try:
            low,high=endpoint(-1),endpoint(1)
            result.update(ci_low_M=low,ci_high_M=high,ci_status='closed' if low and high else 'lower_open' if high else 'upper_open' if low else 'both_open')
        except (ValueError,RuntimeError,FloatingPointError):
            result['ci_status']='profile_failed'; result['diagnostics'].append('profile_failed')
    result['profile']=sorted(cache.values(),key=lambda a:a['log10_parameter_M']); result['profile_threshold_sse']=threshold
    if any(v['amplitude']<=0 for v in params): result['diagnostics'].append('nonpositive_specific_amplitude')
    if cell and any(v['nonspecific_slope_per_M']<0 for v in params): result['diagnostics'].append('negative_nonspecific_slope')
    if np.any(np.minimum(best.x-lo,hi-best.x)<1e-4): result['diagnostics'].append('numerical_boundary_hit')
    if np.linalg.matrix_rank(best.jac)<n_nonlin: result['diagnostics'].append('rank_deficient')
    if not x[x>0].min()<=estimate<=x.max(): result['diagnostics'].append('parameter_outside_titrated_range')
    if result['ci_status'] not in ('closed','noise_scale_not_estimable'): result['diagnostics'].append('interval_not_closed')
    if fitted_pt and (p/estimate>=10 or result['ci_low_M'] is None):
        result['diagnostics'].append('fitted_Pt_titration_or_flat_profile_point_withheld')
        if result['ci_status'] in ('closed','lower_open') and not any(s in result['diagnostics'] for s in ('nonpositive_specific_amplitude','rank_deficient')):
            result['supported_upper_bound_M']=result['ci_high_M']
    blocking=[v for v in result['diagnostics'] if v!='noise_scale_not_estimable']
    result['reportable']=not blocking; result['status']='estimated' if result['reportable'] else 'limited'
    if result['reportable']:
        result['ki_M' if model=='competition_exact' else 'kd_M']=estimate
    return result
