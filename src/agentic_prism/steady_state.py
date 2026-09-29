"""Opt-in SPR/BLI plateau affinity, preserving default kinetic artifacts."""
import numpy as np
import pandas as pd
from .affinity_schema import positive,text
from .schema import resolve_config
from .fit import fit_curve

DEFAULTS={'enabled':False,'windows':[],'response_scales_comparable':None,'rationale':'',
          'flatness_max_fraction':None,'ratio_band':[.5,2.]}


def validate(s,cfg):
    if type(s['enabled']) is not bool: raise ValueError('steady_state.enabled must be boolean')
    if not s['enabled']:
        if s != DEFAULTS: raise ValueError('Disabled steady_state must not contain active options')
        return
    if cfg['assay']['injection_design']!='multi_cycle': raise ValueError('Steady-state currently requires independent cycles')
    if cfg['fit']['rmax']!='shared':
        raise ValueError('Steady-state affinity requires fit.rmax=shared: per-curve Rmax means plateau responses are not on one scale')
    if s['response_scales_comparable'] is not True: raise ValueError('Declare comparable steady-state response scales across sensors')
    text(s['rationale'],'steady-state rationale')
    if not isinstance(s['windows'],list) or not s['windows']: raise ValueError('Declare plateau windows per curve')
    ids=[]
    for w in s['windows']:
        if not isinstance(w,dict) or set(w)!={'curve_id','start_s','end_s'}: raise ValueError('Each plateau window requires curve_id/start_s/end_s')
        text(w['curve_id'],'window curve_id');positive(w['start_s'],'start_s');positive(w['end_s'],'end_s')
        if w['start_s']>=w['end_s']: raise ValueError('Increasing plateau window required')
        ids.append(w['curve_id'])
    if len(set(ids))!=len(ids): raise ValueError('Duplicate plateau windows')
    if s['flatness_max_fraction'] is not None:
        positive(s['flatness_max_fraction'],'flatness_max_fraction')
        if s['flatness_max_fraction']>.05: raise ValueError('Flatness proxy must be at most 5% change across the declared window')
    band=s['ratio_band']
    if not isinstance(band,list) or len(band)!=2: raise ValueError('Declare ratio interpretation band')
    for v in band: positive(v,'ratio band')
    if band[0]>=band[1]: raise ValueError('Increasing ratio band required')


def compute(d,fits,cfg):
    s=cfg['steady_state'];windows={w['curve_id']:w for w in s['windows']}
    if set(windows)!=set(d.curve_id): raise ValueError('Plateau windows must cover every supplied curve exactly')
    groups=[]
    for ft in fits:
        g=d[d.fit_group_id==ft['fit_group_id']];rows=[];obs=[]
        for cid,h in g.groupby('curve_id',sort=False):
            w=windows[cid];a=float(h.association_start_s.iloc[0]);end=float(h.dissociation_start_s.iloc[0])-a
            if w['end_s']>end: raise ValueError('Plateau window exceeds association phase')
            t=h.time_s-a;use=h[(h.phase=='association') & h.used_in_fit & (t>=w['start_s']) & (t<=w['end_s'])]
            row={'curve_id':str(cid),'concentration_M':float(h.concentration_M.iloc[0]),'window':w,'included':False,
                 'fraction_equilibrium_end':None,'fraction_equilibrium_window_start':None,'flatness_fraction':None,'reason':None,'req':None}
            if len(use)<3: row['reason']='insufficient_plateau_points'
            else:
                yy=use.response_processed.to_numpy(float);tt=use.time_s.to_numpy(float)
                offset=0.
                if ft['reportable']:
                    rate=ft['kon_M_inv_s_inv']*row['concentration_M']+ft['koff_s_inv']
                    row['fraction_equilibrium_end']=float(-np.expm1(-rate*end))
                    row['fraction_equilibrium_window_start']=float(-np.expm1(-rate*w['start_s']))
                    offset=next(v['offset'] for v in ft['curve_parameters'] if v['curve_id']==cid)
                    row['reason']='below_95_percent_equilibrium' if row['fraction_equilibrium_end']<.95 or row['fraction_equilibrium_window_start']<.95 else None
                    row['equilibrium_assessment']='reportable_1_to_1_fit'
                elif s['flatness_max_fraction'] is None: row['reason']='kinetics_unreportable_no_declared_flatness_criterion'
                elif cfg['fit']['offset']!='fixed_zero': row['reason']='unreportable_kinetic_offset_cannot_correct_req'
                else:
                    change=abs(float(np.polyfit(tt-tt.min(),yy,1)[0])*(tt.max()-tt.min()))
                    row['flatness_fraction']=change/max(abs(float(yy.mean())),1e-300)
                    row['reason']='plateau_not_flat' if row['flatness_fraction']>s['flatness_max_fraction'] else None
                    row['equilibrium_assessment']='declared_flatness_proxy_not_proven_95_percent'
                row['req']=float(yy.mean()-offset)
                if row['req']<=0: row['reason']='nonpositive_req'
                if row['reason'] is None:
                    row['included']=True
                    obs.append(dict(sample_id=str(h.sample_id.iloc[0]),experiment_id=str(h.experiment_id.iloc[0]),curve_id=str(ft['fit_group_id']),
                        concentration_M=row['concentration_M'],response=row['req'],response_unit=str(h.response_unit.iloc[0]),exclude=False))
            rows.append(row)
        kd_fit=None;ratio=None
        if len(obs)>=4 and len(set(r['concentration_M'] for r in obs))>=4:
            cc=resolve_config({'input':'saved plateau means','assay':{'equilibrium_supported':True,'signal_proportional':True,
                'single_site_supported':True,'concentration_basis':'free','rationale':s['rationale']}})
            kd_fit=fit_curve(pd.DataFrame(obs),cc)
            if kd_fit['reportable'] and ft['reportable']: ratio=kd_fit['kd_M']/ft['kd_M']
        groups.append({'fit_group_id':ft['fit_group_id'],'curves':rows,'fit':kd_fit,'kd_ss_over_kd_kin':ratio,'ratio_band':s['ratio_band'],
            'ratio_status':'unavailable' if ratio is None else 'within_declared_band' if s['ratio_band'][0]<=ratio<=s['ratio_band'][1] else 'discordant_diagnostic',
            'must_mention':['KD_ss and KD_kin are never averaged.','Plateau flatness alone does not prove equilibrium.',
                            'Mean-window responses and their profile-F intervals condition on preprocessing and kinetic offset.',
                            'At least four distinct eligible concentrations required; no imputation of missing plateaus.']})
    from .evidence_0110 import EVIDENCE
    return {'validation_evidence':EVIDENCE['steady_state'],'groups':groups,'method':'existing hyperbola on declared plateau means','new_published_example_gate':'UNMET'}
