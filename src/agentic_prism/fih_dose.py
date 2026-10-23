"""First-in-human dose calculations from declared inputs (0.15.2); design-only, no data file.

- NOAEL -> HED -> MRSD (FDA 2005 guidance): body-surface-area divisors (Table 1, divisor = 37/Km),
  the body-weight exponent form HED = NOAEL (W_animal/W_human)^(1 - 0.67), or mg/kg scaling
  (required for proteins > 100 kDa given parenterally); MRSD = HED / safety factor.
- MABEL / pharmacologically active dose: target concentration from receptor occupancy
  (free drug KD RO/(1 - RO), plus RO x total target under 1:1 depletion) or a declared effect
  concentration, times the human central volume (Cmax after an IV bolus).
- Allometric human PK: animal 2-compartment parameters per kg scaled by declared exponents on the
  absolute scale; candidate human doses simulated (IV bolus or infusion).
Every output is a calculation from declared inputs, never a recommended or approved FIH dose.
"""
from copy import deepcopy
import math
import numpy as np
from scipy.integrate import solve_ivp

TYPES = {'fih_dose'}
# FDA 2005 guidance Table 1: reference body weight (kg), Km, divisor (animal mg/kg -> HED mg/kg).
TABLE1 = {'mouse': (.02, 3, 12.3), 'hamster': (.08, 5, 7.4), 'rat': (.15, 6, 6.2), 'ferret': (.30, 7, 5.3),
          'guinea_pig': (.40, 8, 4.6), 'rabbit': (1.8, 12, 3.1), 'dog': (10., 20, 1.8), 'monkey': (3., 12, 3.1),
          'marmoset': (.35, 6, 6.2), 'squirrel_monkey': (.60, 7, 5.3), 'baboon': (12., 20, 1.8),
          'micro_pig': (20., 27, 1.4), 'mini_pig': (40., 35, 1.1)}
HUMAN_KM = 37
MODALITIES = ('small_molecule', 'protein_le_100kda', 'protein_gt_100kda')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'fih_dose', 'input': None, 'source': None,
            'molecule': None, 'human': None, 'noael': None, 'mabel': None, 'allometric_pk': None,
            'report': {'plot_style': 'prism_like'}}


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def _pos(v):
    return type(v) in (int, float) and math.isfinite(v) and v > 0


def _keys(obj, need, opt=(), name=''):
    if not isinstance(obj, dict) or not set(need) <= set(obj) <= set(need) | set(opt):
        raise ValueError(f'{name} needs ' + ', '.join(need) + (' (optional ' + ', '.join(opt) + ')' if opt else ''))


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported fih_dose configuration keys')
    c = deepcopy(DEFAULTS); c.update({k: deepcopy(v) for k, v in raw.items()})
    if c['schema_version'] != 1 or c['analysis_type'] != 'fih_dose': raise ValueError('Unsupported fih_dose schema')
    if c['input'] is not None: raise ValueError('fih_dose takes declared inputs, not a data file')
    if not _text(c['source']): raise ValueError('Declare source')
    m = c['molecule']; _keys(m, ('name', 'modality', 'molecular_weight_kda'), name='molecule')
    if m['modality'] not in MODALITIES: raise ValueError('molecule.modality must be ' + ', '.join(MODALITIES))
    if not _pos(m['molecular_weight_kda']): raise ValueError('molecule.molecular_weight_kda must be positive (needed for molar conversion)')
    h = c['human']; _keys(h, ('body_weight_kg', 'source'), name='human')
    if not _pos(h['body_weight_kg']) or not _text(h['source']): raise ValueError('human.body_weight_kg and source required')
    if not any(c[k] for k in ('noael', 'mabel', 'allometric_pk')): raise ValueError('Declare at least one of noael, mabel, allometric_pk')
    n = c['noael']
    if n is None: c.pop('noael')
    else:
        _keys(n, ('studies', 'scaling', 'species_selection', 'safety_factor', 'safety_factor_rationale'), ('selection_rationale',), 'noael')
        if n['scaling'] not in ('bsa_table', 'body_weight_exponent', 'mg_per_kg'): raise ValueError('noael.scaling must be bsa_table, body_weight_exponent or mg_per_kg')
        if m['modality'] == 'protein_gt_100kda' and n['scaling'] != 'mg_per_kg':
            raise ValueError('Proteins > 100 kDa given parenterally are scaled by mg/kg, not body surface area (FDA 2005 guidance)')
        if not isinstance(n['studies'], list) or not n['studies']: raise ValueError('noael.studies must be a non-empty list')
        for s in n['studies']:
            _keys(s, ('species', 'noael_mg_per_kg', 'source'), ('animal_body_weight_kg',), 'each NOAEL study')
            if s['species'] not in TABLE1: raise ValueError('species must be one of ' + ', '.join(TABLE1))
            if not _pos(s['noael_mg_per_kg']) or not _text(s['source']): raise ValueError('NOAEL must be positive with a source')
            if n['scaling'] == 'body_weight_exponent' and not _pos(s.get('animal_body_weight_kg')): raise ValueError('body_weight_exponent needs animal_body_weight_kg')
        sel = n['species_selection']
        if sel != 'most_sensitive':
            if sel not in {s['species'] for s in n['studies']}: raise ValueError('species_selection must be most_sensitive or a species in studies')
            if not _text(n.get('selection_rationale')): raise ValueError('A named most-appropriate species needs selection_rationale')
        if not _pos(n['safety_factor']) or n['safety_factor'] < 1 or not _text(n['safety_factor_rationale']): raise ValueError('safety_factor must be >= 1 with a rationale')
    mb = c['mabel']
    if mb is None: c.pop('mabel')
    else:
        _keys(mb, ('targets', 'central_volume_l_per_kg', 'central_volume_source'), name='mabel')
        if not _pos(mb['central_volume_l_per_kg']) or not _text(mb['central_volume_source']): raise ValueError('mabel needs central_volume_l_per_kg with a source')
        if not isinstance(mb['targets'], list) or not mb['targets']: raise ValueError('mabel.targets must be a non-empty list')
        for t in mb['targets']:
            if not isinstance(t, dict) or t.get('basis') not in ('receptor_occupancy', 'effect_concentration'): raise ValueError('Each MABEL target needs basis receptor_occupancy or effect_concentration')
            if t['basis'] == 'receptor_occupancy':
                _keys(t, ('basis', 'name', 'kd_nM', 'target_occupancy_percent', 'source'), ('total_target_nM',), 'receptor_occupancy target')
                if not _pos(t['kd_nM']) or not (type(t['target_occupancy_percent']) in (int, float) and 0 < t['target_occupancy_percent'] < 100):
                    raise ValueError('kd_nM must be positive and target_occupancy_percent in (0, 100)')
                if 'total_target_nM' in t and not _pos(t['total_target_nM']): raise ValueError('total_target_nM must be positive')
            else:
                _keys(t, ('basis', 'name', 'value', 'unit', 'source'), name='effect_concentration target')
                if not _pos(t['value']) or t['unit'] not in ('nM', 'ug/mL'): raise ValueError('effect_concentration needs a positive value with unit nM or ug/mL')
            if not _text(t['name']) or not _text(t['source']): raise ValueError('Each MABEL target needs a name and a source')
    a = c['allometric_pk']
    if a is None: c.pop('allometric_pk')
    else:
        _keys(a, ('species', 'animal_body_weight_kg', 'parameters_per_kg', 'units', 'exponents', 'exponents_source', 'candidate_doses_mg_per_kg', 'route', 'horizon_h', 'source'),
              ('infusion_duration_h', 'kd_nM'), 'allometric_pk')
        p = a['parameters_per_kg']
        if not isinstance(p, dict) or not {'cl', 'v1', 'q', 'v2'} <= set(p) <= {'cl', 'v1', 'q', 'v2', 'vmax', 'km'} or not all(_pos(v) for v in p.values()) or (('vmax' in p) != ('km' in p)):
            raise ValueError('parameters_per_kg needs positive cl, v1, q, v2 (optionally vmax with km)')
        if a['units'] != {'cl': 'L/h/kg', 'v': 'L/kg', 'vmax': 'mg/h/kg', 'km': 'ug/mL'}: raise ValueError('units must be {"cl": "L/h/kg", "v": "L/kg", "vmax": "mg/h/kg", "km": "ug/mL"} (convert first; no silent unit conversion)')
        e = a['exponents']; need = {'cl', 'q', 'v1', 'v2'} | ({'vmax'} if 'vmax' in p else set())
        if not isinstance(e, dict) or set(e) != need or not all(_pos(v) and v <= 1.2 for v in e.values()) or not _text(a['exponents_source']):
            raise ValueError('exponents must give ' + ', '.join(sorted(need)) + ' (0 < b <= 1.2) with exponents_source')
        if not _pos(a['animal_body_weight_kg']) or not _text(a['source']) or not _text(a['species']): raise ValueError('allometric_pk needs species, animal_body_weight_kg and source')
        if not isinstance(a['candidate_doses_mg_per_kg'], list) or not a['candidate_doses_mg_per_kg'] or not all(_pos(x) for x in a['candidate_doses_mg_per_kg']): raise ValueError('candidate_doses_mg_per_kg must be positive numbers')
        if a['route'] not in ('iv_bolus', 'iv_infusion'): raise ValueError('allometric_pk.route must be iv_bolus or iv_infusion')
        if a['route'] == 'iv_infusion' and not _pos(a.get('infusion_duration_h')): raise ValueError('iv_infusion needs infusion_duration_h')
        if not _pos(a['horizon_h']): raise ValueError('horizon_h must be positive')
        if 'kd_nM' in a and not _pos(a['kd_nM']): raise ValueError('kd_nM must be positive')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def hed(noael, species, scaling, animal_bw=None, human_bw=60.):
    if scaling == 'mg_per_kg': return float(noael)
    if scaling == 'bsa_table': return float(noael/TABLE1[species][2])
    return float(noael*(animal_bw/human_bw)**(1 - .67))


def nm_to_ug_ml(nm, mw_kda):
    return nm*mw_kda*1e-3


def target_concentration(t, mw):
    """Total drug concentration (ug/mL) for a MABEL target."""
    if t['basis'] == 'effect_concentration':
        return t['value'] if t['unit'] == 'ug/mL' else nm_to_ug_ml(t['value'], mw)
    ro = t['target_occupancy_percent']/100; free = t['kd_nM']*ro/(1 - ro)
    return nm_to_ug_ml(free + ro*t.get('total_target_nM', 0.), mw)


def human_parameters(a, human_bw, exponents=None):
    e = exponents or a['exponents']; wa = a['animal_body_weight_kg']; out = {}
    for k, v in a['parameters_per_kg'].items():
        if k == 'km': out[k] = v; continue
        b = e['cl' if k == 'cl' else k] if k in e else e['cl']
        out[k] = v*wa*(human_bw/wa)**b/human_bw  # absolute scaling, back to per kg
    return out


def simulate(p, dose_mg_kg, route, horizon, duration=None, n=2001):
    """Central concentration (mg/L = ug/mL) for a per-kg dose; returns times and concentrations."""
    def rhs(t, y):
        c = y[0]/p['v1']; inf = dose_mg_kg/duration if route == 'iv_infusion' and t < duration else 0.
        el = p['cl']*c + (p['vmax']*c/(p['km'] + c) if 'vmax' in p else 0.)
        return [inf - el - p['q']*c + p['q']*y[1]/p['v2'], p['q']*c - p['q']*y[1]/p['v2']]
    t = np.linspace(0, horizon, n)
    if route == 'iv_infusion' and duration not in t: t = np.unique(np.r_[t, duration])
    s = solve_ivp(rhs, (0, horizon), [dose_mg_kg if route == 'iv_bolus' else 0., 0.], t_eval=t, method='LSODA', rtol=1e-10, atol=1e-12, max_step=horizon/2000)
    return s.t, s.y[0]/p['v1']


def _beta(p):
    k10, k12, k21 = p['cl']/p['v1'], p['q']/p['v1'], p['q']/p['v2']; s = k10 + k12 + k21
    return (s - math.sqrt(s*s - 4*k10*k21))/2


def _time_above(t, c, level):
    above = c >= level
    return float(np.trapezoid(above.astype(float), t)) if above.any() else 0.


def project(a, human_bw, mw, targets):
    p = human_parameters(a, human_bw); rows = []
    for dose in a['candidate_doses_mg_per_kg']:
        t, c = simulate(p, dose, a['route'], a['horizon_h'], a.get('infusion_duration_h'))
        r = {'dose_mg_per_kg': dose, 'dose_mg': dose*human_bw, 'cmax_ug_per_ml': float(c.max()), 'c_horizon_ug_per_ml': float(c[-1]),
             'auc_to_horizon_ug_h_per_ml': float(np.trapezoid(c, t))}
        if 'vmax' not in p:
            r['aucinf_ug_h_per_ml'] = float(dose/p['cl']); r['terminal_half_life_h'] = float(math.log(2)/_beta(p))
        r['hours_above_target'] = {x['name']: _time_above(t, c, x['concentration_ug_per_ml']) for x in targets}
        if a.get('kd_nM'):
            kd = nm_to_ug_ml(a['kd_nM'], mw); r['ro_percent_at_cmax'] = float(100*c.max()/(c.max() + kd)); r['ro_percent_at_horizon'] = float(100*c[-1]/(c[-1] + kd))
        rows.append(r)
    return p, rows


def compute(_, c):
    mw = c['molecule']['molecular_weight_kda']; bw = c['human']['body_weight_kg']; primary = {}; failing = []
    if c.get('noael'):
        n = c['noael']; studies = []
        for s in n['studies']:
            row = {**s, 'hed_mg_per_kg': hed(s['noael_mg_per_kg'], s['species'], n['scaling'], s.get('animal_body_weight_kg'), bw),
                   'divisor': TABLE1[s['species']][2] if n['scaling'] == 'bsa_table' else None,
                   'sensitivity': {k: hed(s['noael_mg_per_kg'], s['species'], k, s.get('animal_body_weight_kg') or TABLE1[s['species']][0], bw)
                                   for k in ('bsa_table', 'body_weight_exponent', 'mg_per_kg')}}
            studies.append(row)
        chosen = (min(studies, key=lambda r: r['hed_mg_per_kg']) if n['species_selection'] == 'most_sensitive'
                  else min((r for r in studies if r['species'] == n['species_selection']), key=lambda r: r['hed_mg_per_kg']))
        lowest = min(studies, key=lambda r: r['hed_mg_per_kg'])
        mrsd = chosen['hed_mg_per_kg']/n['safety_factor']
        primary['noael'] = {'scaling': n['scaling'], 'studies': studies, 'selected_species': chosen['species'], 'selection': n['species_selection'],
                            'hed_mg_per_kg': chosen['hed_mg_per_kg'], 'safety_factor': n['safety_factor'],
                            'mrsd_mg_per_kg': mrsd, 'mrsd_mg': mrsd*bw, 'human_body_weight_kg': bw,
                            'selected_is_not_most_sensitive': chosen is not lowest and chosen['hed_mg_per_kg'] > lowest['hed_mg_per_kg']}
        if primary['noael']['selected_is_not_most_sensitive']: failing.append({'item': 'noael', 'reason': f"selected species {chosen['species']} is not the most sensitive ({lowest['species']} HED {lowest['hed_mg_per_kg']:.3g} mg/kg)"})
        if n['safety_factor'] < 10: failing.append({'item': 'noael', 'reason': f"safety factor {n['safety_factor']:g} is below the guidance default of 10 (declared reason: {n['safety_factor_rationale']})"})
    targets = []
    if c.get('mabel'):
        mb = c['mabel']
        for t in mb['targets']:
            conc = target_concentration(t, mw); dose = conc*mb['central_volume_l_per_kg']
            row = {'name': t['name'], 'basis': t['basis'], 'concentration_ug_per_ml': conc, 'concentration_nM': conc/(mw*1e-3),
                   'dose_mg_per_kg': dose, 'dose_mg': dose*bw, 'source': t['source']}
            if t['basis'] == 'receptor_occupancy':
                row['sensitivity_kd_div3_x3_mg_per_kg'] = [target_concentration({**t, 'kd_nM': t['kd_nM']*f}, mw)*mb['central_volume_l_per_kg'] for f in (1/3, 3)]
            targets.append(row)
        best = min(targets, key=lambda r: r['dose_mg_per_kg'])
        primary['mabel'] = {'targets': targets, 'central_volume_l_per_kg': mb['central_volume_l_per_kg'], 'lowest_target': best['name'],
                            'mabel_mg_per_kg': best['dose_mg_per_kg'], 'mabel_mg': best['dose_mg']}
        if 'noael' in primary: primary['mabel']['mrsd_to_mabel_ratio'] = primary['noael']['mrsd_mg_per_kg']/best['dose_mg_per_kg']
    if c.get('allometric_pk'):
        a = c['allometric_pk']; p, rows = project(a, bw, mw, targets)
        sens = {f'cl_exponent_{b:g}': human_parameters({**a, 'parameters_per_kg': {'cl': a['parameters_per_kg']['cl'], 'v1': 1, 'q': 1, 'v2': 1}}, bw, {**a['exponents'], 'cl': b})['cl'] for b in (.75, .85, 1.)}
        primary['allometric_pk'] = {'species': a['species'], 'animal_body_weight_kg': a['animal_body_weight_kg'], 'human_parameters_per_kg': p,
                                    'exponents': a['exponents'], 'route': a['route'], 'horizon_h': a['horizon_h'], 'doses': rows,
                                    'sensitivity_human_cl_l_h_kg': sens, 'nonlinear': 'vmax' in p}
    return {'analysis_type': 'fih_dose', 'primary': primary, 'must_mention': _notes(c, primary), 'failing_items': failing,
            'limitations': ['Calculations from declared inputs only; none of these numbers is a recommended or approved first-in-human dose. Dose selection integrates toxicology, pharmacology, PK/PD and clinical considerations under FDA (2005) and EMA (2017) guidance.',
                            'MABEL uses Cmax after an IV bolus in the declared central volume; it ignores distribution, target turnover and receptor reserve unless they are built into the declared target concentration.',
                            'Allometric projection assumes the declared exponents and that target-mediated parameters (Km) do not scale; species differences in target expression, FcRn binding or immunogenicity are not modelled.',
                            'No published worked example was reproduced for MABEL or antibody allometry (see the 0.15.2 record).']}


def _notes(c, p):
    m = c['molecule']; notes = [f"First-in-human dose calculations for {m['name']} ({m['modality']}, {m['molecular_weight_kda']:g} kDa), human body weight {c['human']['body_weight_kg']:g} kg ({c['human']['source']}). These are calculations, not a recommended or approved dose."]
    if 'noael' in p:
        n = p['noael']; s = {'bsa_table': 'body-surface-area divisors (FDA 2005 Table 1)', 'body_weight_exponent': 'body-weight exponent (W_animal/W_human)^0.33', 'mg_per_kg': 'mg/kg (no BSA scaling)'}[n['scaling']]
        notes.append(f"NOAEL to HED by {s}; selected species {n['selected_species']} ({n['selection']}): HED {n['hed_mg_per_kg']:.4g} mg/kg; safety factor {n['safety_factor']:g} ({c['noael']['safety_factor_rationale']}); MRSD {n['mrsd_mg_per_kg']:.4g} mg/kg ({n['mrsd_mg']:.4g} mg).")
        if c['molecule']['modality'] == 'protein_gt_100kda': notes.append('For proteins > 100 kDa given parenterally the guidance scales by mg/kg; BSA scaling would not be appropriate.')
    if 'mabel' in p:
        mb = p['mabel']
        notes.append('MABEL targets (dose = target concentration x human central volume ' + f"{mb['central_volume_l_per_kg']:g} L/kg, Cmax after an IV bolus): " +
                     '; '.join(f"{t['name']} {t['concentration_ug_per_ml']:.4g} ug/mL ({t['concentration_nM']:.4g} nM) -> {t['dose_mg_per_kg']:.4g} mg/kg" for t in mb['targets']) +
                     f". Lowest: {mb['lowest_target']} {mb['mabel_mg_per_kg']:.4g} mg/kg ({mb['mabel_mg']:.4g} mg)." + (f" MRSD/MABEL ratio {mb['mrsd_to_mabel_ratio']:.3g}." if 'mrsd_to_mabel_ratio' in mb else ''))
        notes.append('Receptor-occupancy targets assume 1:1 binding at equilibrium with the declared KD (nM converted with the declared molecular weight); KD from a different species, format or assay changes the result (sensitivity: KD /3 and x3 reported).')
    if 'allometric_pk' in p:
        a = p['allometric_pk']; e = a['exponents']
        notes.append(f"Human PK projected from {a['species']} ({a['animal_body_weight_kg']:g} kg) by absolute allometry with declared exponents " + ', '.join(f'{k} {v:g}' for k, v in e.items()) +
                     f" ({c['allometric_pk']['exponents_source']}); human CL for exponents 0.75/0.85/1.0: " + ', '.join(f'{v:.3g}' for v in a['sensitivity_human_cl_l_h_kg'].values()) + ' L/h/kg.' +
                     (' Michaelis-Menten Vmax is scaled with its exponent and Km is carried over unchanged (an assumption).' if a['nonlinear'] else ''))
    return notes


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    a = result['primary'].get('allometric_pk')
    if not a: return []
    cfg = __import__('json').loads((Path(run)/'config.resolved.json').read_text())['allometric_pk']
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.2, 3.4), layout='constrained')
        for k, d in enumerate(a['doses']):
            t, c = simulate(a['human_parameters_per_kg'], d['dose_mg_per_kg'], a['route'], a['horizon_h'], cfg.get('infusion_duration_h'), 400)
            ax.plot(t, c, color=pstyle.color(token, k), lw=1, label=f"{d['dose_mg_per_kg']:g} mg/kg")
        for x in (result['primary'].get('mabel') or {}).get('targets', []):
            ax.axhline(x['concentration_ug_per_ml'], ls=':', color=pstyle.MUTED, lw=.8)
        ax.set_yscale('log'); ax.set(xlabel='Time (h)', ylabel='Projected concentration (ug/mL)', title='Projected human profiles (declared inputs)')
        ax.legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'fih_projection', token, 300)
    return ['fih_projection']
