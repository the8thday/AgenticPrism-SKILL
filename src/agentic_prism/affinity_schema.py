"""Opt-in equilibrium affinity contracts, separate from legacy resolved defaults."""
from copy import deepcopy
import json
import numpy as np
from .schema import DEFAULTS as LEGACY_DEFAULTS, load_data as load_equilibrium_data

MODELS = ('one_site_depletion', 'solution_equilibrium_titration', 'competition_exact')
DEFAULTS = deepcopy(LEGACY_DEFAULTS)
DEFAULTS.update(model=None,
    active_sites={'mode': None, 'concentration_M': None, 'basis': None, 'source': '', 'activity_basis': '', 'bounds_M': [1e-15, 1.]},
    equilibration={'incubation_s': None, 'method': None, 'koff_per_s': None, 'time_to_95_s': None, 'source': ''},
    competition={'tracer_total_M': None, 'tracer_kd_M': None, 'source': '', 'cheng_prusoff': False, 'ic50_M': None},
    valency={'constant_species': None, 'constant_concentration_basis': None, 'detected_signal': None,
             'titrant_concentration_basis': None, 'source': ''},
    model_rationale='')


def merge(base, raw, path=''):
    if not isinstance(raw, dict): raise ValueError(f'{path} must be an object')
    if set(raw)-set(base): raise ValueError(f'Unsupported settings {path}: {sorted(set(raw)-set(base))}')
    out = deepcopy(base)
    for k,v in raw.items():
        out[k] = merge(base[k],v,path+k+'.') if isinstance(base[k],dict) and k!='column_map' else v
    return out


def positive(v, name):
    if isinstance(v,bool) or not isinstance(v,(float,int)) or not np.isfinite(v) or v<=0:
        raise ValueError(f'{name} must be positive and finite')


def text(v, name):
    if not isinstance(v,str) or not v.strip(): raise ValueError(f'Declare {name} with rationale/source')


def validate_common(c):
    if c['schema_version'] != 1: raise ValueError('Only schema 1 supported')
    text(c['input'],'input'); text(c['source'],'source'); text(c['model_rationale'],'model_rationale')
    a=c['assay']
    for key in ('equilibrium_supported','signal_proportional','single_site_supported'):
        if a[key] is not True: raise ValueError(f'Assay gate: {key} must be literally true with evidence')
    text(a['rationale'],'assay rationale')
    f=c['fit']
    if f['baseline_mode']!='fitted' or f['residual_scale']!='linear':
        raise ValueError('New affinity models require jointly fitted baseline and linear residual scale')
    if f['weighting'] not in ('unweighted','inverse_sd','relative'): raise ValueError('Unsupported weighting')
    for key in ('kd_bounds_M',):
        b=f[key]
        if not isinstance(b,list) or len(b)!=2: raise ValueError('Invalid KD bounds')
        for v in b: positive(v,'KD bound')
        if b[0]>=b[1]: raise ValueError('Increasing KD bounds required')
    if type(f['multistart']) is not int or not 1<=f['multistart']<=20: raise ValueError('Invalid multistart')
    if type(f['max_nfev']) is not int or f['max_nfev']<100: raise ValueError('Invalid max_nfev')
    u=c['uncertainty']
    if u['parameter_ci']!='profile_f' or not .5<u['level']<1: raise ValueError('New models require profile_f')
    if c['sensitivity']: raise ValueError('Legacy sensitivity options do not apply to new models')
    if c['replicates']['independent_unit'] not in ('none','experiment_id') or type(c['replicates']['conditions_comparable']) is not bool:
        raise ValueError('Declare independence and comparability')
    r=c['report']
    if r['language']!='zh-CN' or r['plot_style'] not in ('standard','prism_like') or r['figure_width_mm'] not in (85,180) or r['png_dpi'] not in (300,600):
        raise ValueError('Unsupported report settings')
    if sorted(r['export_formats'])!=['pdf','png','svg'] or type(r['allow_style_switch']) is not bool: raise ValueError('Invalid exports')
    if not isinstance(c['column_map'],dict) or not all(isinstance(k,str) and isinstance(v,str) for k,v in c['column_map'].items()): raise ValueError('Invalid column_map')
    json.dumps(c,allow_nan=False)


def validate_active(s, allow_fit=True):
    if s['mode'] not in (('declared','fitted') if allow_fit else ('declared',)): raise ValueError('Declare active-site concentration mode')
    text(s['source'],'active-site source'); text(s['activity_basis'],'active-site activity basis')
    if s['basis']!='active_sites': raise ValueError('Nominal protein is not active Pt; specify active_sites and activity basis')
    if s['mode']=='declared': positive(s['concentration_M'],'active-site concentration_M')
    b=s['bounds_M']
    if not isinstance(b,list) or len(b)!=2: raise ValueError('Invalid Pt bounds')
    for v in b: positive(v,'Pt bound')
    if b[0]>=b[1]: raise ValueError('Increasing Pt bounds required')


def validate_set_valency(v):
    """The SET signal model is the free-site fraction of the constant species.

    A capture step that detects every bivalent molecule with at least one free
    site is not proportional to that fraction, so it is refused, not corrected.
    """
    text(v['source'],'SET valency and detection source')
    if v['constant_species'] not in ('monovalent','bivalent'): raise ValueError('Declare SET constant-species valency: monovalent or bivalent')
    if v['detected_signal'] not in ('free_sites','molecules_with_any_free_site'):
        raise ValueError('Declare what the SET signal detects: free_sites or molecules_with_any_free_site')
    if v['constant_concentration_basis'] not in ('binding_sites','molecules'): raise ValueError('Declare constant-species concentration basis')
    if v['titrant_concentration_basis'] not in ('binding_sites','monovalent_molecules'):
        raise ValueError('Declare titrant concentration as binding_sites (e.g. 2x IgG) or monovalent_molecules')
    if v['constant_species']=='bivalent':
        if v['detected_signal']=='molecules_with_any_free_site':
            raise ValueError('Bivalent constant species detected as molecules with any free site is not proportional to the '
                             'free-site fraction; no bivalent SET model is implemented. Use a monovalent format (e.g. Fab) or a free-site readout')
        if v['constant_concentration_basis']!='binding_sites':
            raise ValueError('Bivalent constant species must be supplied as binding-site concentrations')


def resolve_config(raw):
    c=merge(DEFAULTS,raw); validate_common(c)
    if c['analysis_type']!='equilibrium_binding' or c['model'] not in MODELS: raise ValueError('Unsupported affinity model')
    if c['assay']['concentration_basis']!='total' or c['assay']['interpretation']!='KD': raise ValueError('New solution models require total concentrations and KD interpretation (Ki for competition)')
    if c['fit']['weighting']=='relative': raise ValueError('Relative weighting is scoped to cell_binding')
    validate_active(c['active_sites'],allow_fit=c['model']=='one_site_depletion')
    if c['model']=='solution_equilibrium_titration':
        e=c['equilibration']; positive(e['incubation_s'],'incubation_s'); text(e['source'],'equilibration source')
        if e['method']=='koff':
            positive(e['koff_per_s'],'koff_per_s'); t95=-np.log(.05)/e['koff_per_s']
        elif e['method']=='measured_time_course': positive(e['time_to_95_s'],'time_to_95_s'); t95=e['time_to_95_s']
        else: raise ValueError('Declare koff or measured_time_course equilibration method')
        if e['incubation_s']<t95: raise ValueError('Incubation is shorter than the declared conservative time to 95% equilibrium')
        validate_set_valency(c['valency'])
    elif c['valency']!=DEFAULTS['valency']: raise ValueError('valency declarations apply only to solution_equilibrium_titration')
    q=c['competition']
    if type(q['cheng_prusoff']) is not bool: raise ValueError('cheng_prusoff must be boolean')
    if c['model']=='competition_exact':
        for key in ('tracer_total_M','tracer_kd_M'): positive(q[key],key)
        text(q['source'],'tracer concentrations and KD source')
        if q['cheng_prusoff']:
            from .binding_models import cheng_prusoff
            cheng_prusoff(q['ic50_M'],q['tracer_total_M'],q['tracer_kd_M'],c['active_sites']['concentration_M'])
    elif q['cheng_prusoff']: raise ValueError('Cheng-Prusoff is only a diagnostic for competition')
    return c


def load_data(path,cfg):
    d=load_equilibrium_data(path,cfg)
    if cfg['model']=='solution_equilibrium_titration':
        for k in ('fit_group_id','constant_species_M','constant_species_source'):
            if k not in d or d[k].astype(str).str.strip().eq('').any(): raise ValueError(f'SET requires {k}')
        import pandas as pd
        d['constant_species_M']=pd.to_numeric(d.constant_species_M,errors='raise')
        if not np.isfinite(d.constant_species_M).all() or (d.constant_species_M<=0).any(): raise ValueError('Positive constant species required')
        for _,g in d.groupby('curve_id',sort=False):
            if g.constant_species_M.nunique()!=1 or g.fit_group_id.nunique()!=1: raise ValueError('Constant species and fit group must be fixed per curve')
        for _,g in d.groupby('fit_group_id',sort=False):
            if any(g[k].nunique()!=1 for k in ('sample_id','experiment_id','response_unit')): raise ValueError('SET cannot pool independent experiments or samples')
            if g.loc[~g.exclude,'constant_species_M'].nunique()<2: raise ValueError('SET n-curve requires at least two included distinct constant-species concentrations')
    else: d['fit_group_id']=d.curve_id
    return d
