"""Cell binding is its own assay: apparent KD, never intrinsic KD or EC50."""
from copy import deepcopy
import numpy as np
from .affinity_schema import merge,validate_common,positive,text
from .schema import DEFAULTS as BASE, load_data as load_equilibrium_data

DEFAULTS=deepcopy(BASE)
DEFAULTS.update(analysis_type='cell_binding',model='cell_binding_apparent',model_rationale='',
    receptors={'depletion':None,'cells_per_mL':None,'receptors_per_cell':None,'source':'','density_id':''},
    nonspecific={'control':None,'baseline':None,'source':''},
    background={'mode':None,'value':None,'source':''},
    cell_assay={'cell_line':'','equilibrium_incubation_supported':None,'incubation_s':None,
        'temperature_C':None,'internalization_controlled':None,'wash_protocol':'','wash_dissociation_supported':None,
        'detection':None,'valency':None,'rationale':''})
DEFAULTS['assay']['interpretation']='apparent_KD'


def resolve_config(raw):
    c=merge(DEFAULTS,raw);validate_common(c)
    if c['analysis_type']!='cell_binding' or c['model']!='cell_binding_apparent': raise ValueError('Only cell_binding_apparent supported')
    if c['assay']['interpretation']!='apparent_KD': raise ValueError('Cell-binding EC50 cannot be labelled intrinsic KD; always apparent_KD')
    if c['assay']['concentration_basis']!='total': raise ValueError('Declare total applied ligand concentration')
    if c['fit']['weighting'] not in ('unweighted','relative'): raise ValueError('Cell weighting must be unweighted or relative')
    a=c['cell_assay']
    for k in ('equilibrium_incubation_supported','internalization_controlled','wash_dissociation_supported'):
        if a[k] is not True: raise ValueError(f'Cell applicability unresolved: {k}')
    for k in ('cell_line','wash_protocol','rationale'): text(a[k],k)
    positive(a['incubation_s'],'cell incubation_s')
    if isinstance(a['temperature_C'],bool) or not isinstance(a['temperature_C'],(float,int)) or not np.isfinite(a['temperature_C']): raise ValueError('Declare temperature_C and internalization rationale')
    if a['detection'] not in ('direct_label','secondary_antibody') or a['valency'] not in ('monovalent','bivalent_IgG','multivalent'):
        raise ValueError('Declare detection and valency')
    s=c['receptors'];text(s['source'],'receptor depletion source');text(s['density_id'],'density_id')
    if s['depletion']=='quadratic':
        positive(s['cells_per_mL'],'cells_per_mL');positive(s['receptors_per_cell'],'receptors_per_cell')
    elif s['depletion']=='negligible':
        if c['assay']['free_approximation_supported'] is not True: raise ValueError('Justify negligible receptor depletion')
    else: raise ValueError('Declare quadratic or negligible receptor depletion')
    n=c['nonspecific']
    if n['control'] not in ('isotype','antigen_negative','excess_competitor'): raise ValueError('Declare measured nonspecific control series')
    text(n['source'],'nonspecific control source')
    if n['baseline'] not in ('shared','separate'): raise ValueError('Declare nonspecific.baseline: shared or separate')
    if n['control']=='antigen_negative' and n['baseline']!='separate':
        raise ValueError('Antigen-negative cells can differ in autofluorescence; their control series needs a separate baseline')
    b=c['background'];text(b['source'],'unstained/secondary-only background handling source')
    if b['mode'] not in ('joint_baseline','subtract_declared'): raise ValueError('Explicit background mode required')
    if b['mode']=='subtract_declared' and (isinstance(b['value'],bool) or not isinstance(b['value'],(int,float)) or not np.isfinite(b['value'])): raise ValueError('Finite declared background required')
    return c


def load_data(path,cfg):
    d=load_equilibrium_data(path,cfg)
    for key in ('cell_line','density_id','series'):
        if key not in d or d[key].str.strip().eq('').any(): raise ValueError(f'Cell input requires {key}')
    if set(d.cell_line)!={cfg['cell_assay']['cell_line']} or set(d.density_id)!={cfg['receptors']['density_id']}:
        raise ValueError('Cannot pool different cell lines or receptor densities in one run/fit')
    if not d.series.isin(['total','nonspecific']).all(): raise ValueError('series must be total or nonspecific')
    for _,g in d.loc[~d.exclude].groupby('curve_id',sort=False):
        if set(g.series)!={'total','nonspecific'}: raise ValueError('Each curve needs total and matched nonspecific series')
        if g.loc[g.series=='total','concentration_M'].nunique()<4: raise ValueError('At least four total-series concentrations required')
        ns=g[g.series=='nonspecific']
        if ns.concentration_M.nunique()<3: raise ValueError('At least three control concentrations needed to fit nonspecific slope')
        if ns.concentration_M.min()>g.concentration_M.min() or ns.concentration_M.max()<g.concentration_M.max(): raise ValueError('Control series must cover total concentration range')
    if cfg['background']['mode']=='subtract_declared':
        d['response_before_background']=d.response
        d['response']=d.response-cfg['background']['value']
    d['fit_group_id']=d.curve_id
    return d
