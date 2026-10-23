"""Single-dose noncompartmental PK analysis computed by R PKNCA (0.15.0).

PKNCA performs every NCA calculation through the fixed bridge script
`r_scripts_pk/pk_runner.R`; this module only checks the declared contract, maps
the data, applies the declared reportability criteria to PKNCA's values,
summarises reportable subject values and fits the Smith et al. (2000) power
model for dose proportionality. No NCA formula is evaluated here.

IV bolus uses two PKNCA passes: pass 1 on the observed samples (time-0 predose
placeholder) gives C0 by PKNCA's back-extrapolation and the observed Cmax/Tmax;
pass 2 inserts that C0 at time 0 so PKNCA's AUC, AUMC, CL, Vz, MRT and Vss
include the back-extrapolated area (as WinNonlin does). When both passes select
the same lambda-z, pass-2 AUCinf must equal pass-1 aucivinf.obs.
"""
from copy import deepcopy
import math
import numpy as np
import pandas as pd
from scipy import stats

TYPES = {'pk_nca'}
ROUTES = ('intravascular_bolus', 'intravascular_infusion', 'extravascular')
BLQ_RULES = ('keep_as_zero', 'drop')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'pk_nca', 'input': None, 'source': None,
            'design': {'analyte': None, 'matrix': None, 'species': None, 'subject_unit': None, 'route': None,
                       'dose_unit': None, 'concentration_unit': None, 'time_unit': None, 'time_basis': None},
            'time_zero': None,
            'blq': None,
            'nca': {'auc_method': None, 'auc_method_rationale': None,
                    'lambda_z': {'selection': None, 'min_points': 3, 'adj_r2_factor': 1e-4}},
            'reportability': None,
            'dose_proportionality': None,
            'report': {'plot_style': 'prism_like'}}
AUC_METHODS = {'linear': 'linear', 'lin_up_log_down': 'lin up/log down'}
TERMINAL = ('lambda_z', 'half_life', 'r_squared', 'adj_r_squared', 'span_ratio', 'lambda_z_n_points',
            'lambda_z_time_first', 'lambda_z_time_last', 'clast_pred')
INF = ('aucinf', 'aucinf_extrap_percent', 'aumcinf', 'cl', 'vz', 'mrt', 'vss')
OBSERVED = ('c0', 'cmax', 'tmax', 'tlast', 'clast', 'auclast', 'aumclast')
PKNCA_NAMES = {'c0': 'c0', 'cmax': 'cmax', 'tmax': 'tmax', 'tlast': 'tlast', 'clast': 'clast.obs',
               'auclast': 'auclast', 'aumclast': 'aumclast', 'lambda_z': 'lambda.z', 'half_life': 'half.life',
               'r_squared': 'r.squared', 'adj_r_squared': 'adj.r.squared', 'span_ratio': 'span.ratio',
               'lambda_z_n_points': 'lambda.z.n.points', 'lambda_z_time_first': 'lambda.z.time.first',
               'lambda_z_time_last': 'lambda.z.time.last', 'clast_pred': 'clast.pred',
               'aucinf': 'aucinf.obs', 'aucinf_extrap_percent': 'aucpext.obs', 'aumcinf': 'aumcinf.obs',
               'cl': 'cl.obs', 'vz': 'vz.obs'}
DP_PARAMETERS = ('cmax', 'auclast', 'aucinf')


def _text(v):
    return isinstance(v, str) and bool(v.strip())


def _num(v):
    return type(v) in (int, float) and math.isfinite(v)


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS): raise ValueError('Unsupported pk_nca configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if k == 'nca':
            if not isinstance(v, dict) or set(v) - set(c['nca']): raise ValueError('Unsupported nca keys')
            lz = v.get('lambda_z', {})
            if not isinstance(lz, dict) or set(lz) - set(c['nca']['lambda_z']): raise ValueError('Unsupported nca.lambda_z keys')
            c['nca'].update({x: y for x, y in v.items() if x != 'lambda_z'}); c['nca']['lambda_z'].update(lz)
        elif isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]): raise ValueError(f'Unsupported {k} keys')
            c[k].update(v)
        else: c[k] = v
    if c['schema_version'] != 1 or c['analysis_type'] != 'pk_nca': raise ValueError('Unsupported pk_nca schema')
    for name, v in [('input', c['input']), ('source', c['source']), *c['design'].items()]:
        if not _text(v): raise ValueError('Declare ' + (name if name in ('input', 'source') else 'design.' + name))
    d = c['design']
    if d['route'] not in ROUTES: raise ValueError('design.route must be intravascular_bolus, intravascular_infusion or extravascular')
    if d['time_basis'] not in ('actual', 'nominal'): raise ValueError('design.time_basis must be actual or nominal')
    tz = c['time_zero']
    if d['route'] == 'intravascular_bolus':
        if tz is not None: raise ValueError('time_zero does not apply to intravascular_bolus: C0 is always back-extrapolated by PKNCA')
        c.pop('time_zero')
    else:
        if not isinstance(tz, dict) or set(tz) != {'rule', 'rationale'} or tz['rule'] not in ('observed', 'insert_zero') or not _text(tz['rationale']):
            raise ValueError('Declare time_zero {rule: observed (a predose time-0 row exists) or insert_zero, rationale}')
    b = c['blq']
    if b is not None:
        if not isinstance(b, dict) or set(b) != {'first', 'middle', 'last', 'rationale'} or not _text(b['rationale']) \
                or any(b[k] not in BLQ_RULES for k in ('first', 'middle', 'last')):
            raise ValueError('blq needs first, middle and last (keep_as_zero or drop) and a rationale')
    n = c['nca']
    if n['auc_method'] not in AUC_METHODS or not _text(n['auc_method_rationale']):
        raise ValueError('Declare nca.auc_method (linear or lin_up_log_down) and auc_method_rationale')
    lz = n['lambda_z']
    if lz['selection'] not in ('best_adjusted_r2', 'declared_points'): raise ValueError('nca.lambda_z.selection must be best_adjusted_r2 or declared_points')
    if type(lz['min_points']) is not int or lz['min_points'] < 3: raise ValueError('nca.lambda_z.min_points must be an integer of at least 3')
    if not _num(lz['adj_r2_factor']) or not 0 <= lz['adj_r2_factor'] < .01: raise ValueError('nca.lambda_z.adj_r2_factor must be in [0, 0.01)')
    r = c['reportability']
    if not isinstance(r, dict): raise ValueError('Declare reportability: criteria with a source, or descriptive_only with a rationale')
    if r.get('descriptive_only') is True:
        if set(r) != {'descriptive_only', 'rationale'} or not _text(r['rationale']): raise ValueError('reportability.descriptive_only needs only a rationale')
    else:
        keys = {'min_adj_r2', 'min_span_ratio', 'max_aucinf_extrap_percent', 'source'}
        if set(r) != keys or not _text(r['source']): raise ValueError('reportability needs min_adj_r2, min_span_ratio, max_aucinf_extrap_percent (each a number or null) and a source')
        for k in keys - {'source'}:
            if r[k] is not None and (not _num(r[k]) or r[k] < 0): raise ValueError(f'reportability.{k} must be a non-negative number or null')
        if r['min_adj_r2'] is not None and r['min_adj_r2'] > 1: raise ValueError('reportability.min_adj_r2 cannot exceed 1')
    dp = c['dose_proportionality']
    if dp is None:
        c.pop('dose_proportionality')
    else:
        keys = {'parameters', 'acceptance', 'confidence_level', 'source', 'rationale'}
        if not isinstance(dp, dict) or not {'parameters', 'acceptance', 'source', 'rationale'} <= set(dp) <= keys:
            raise ValueError('dose_proportionality needs parameters, acceptance [theta_low, theta_high], source and rationale (optional confidence_level)')
        dp.setdefault('confidence_level', .9)
        if not isinstance(dp['parameters'], list) or not dp['parameters'] or set(dp['parameters']) - set(DP_PARAMETERS) or len(set(dp['parameters'])) != len(dp['parameters']):
            raise ValueError('dose_proportionality.parameters must be distinct values from cmax, auclast, aucinf')
        a = dp['acceptance']
        if not isinstance(a, list) or len(a) != 2 or not all(_num(x) for x in a) or not 0 < a[0] < 1 < a[1]:
            raise ValueError('dose_proportionality.acceptance must be [theta_low < 1 < theta_high]')
        if not _num(dp['confidence_level']) or not .5 < dp['confidence_level'] < 1: raise ValueError('Invalid dose_proportionality.confidence_level')
        if not _text(dp['source']) or not _text(dp['rationale']): raise ValueError('Declare dose_proportionality source and rationale')
        if 'aucinf' in dp['parameters'] and r.get('descriptive_only'): raise ValueError('aucinf dose proportionality needs declared reportability criteria')
    if c['report']['plot_style'] not in ('standard', 'prism_like'): raise ValueError('Invalid plot style')
    return c


def _bool_column(d, col):
    v = d[col].astype(str).str.strip().str.lower()
    if not v.isin(['true', 'false']).all(): raise ValueError(f'{col} must be true or false')
    return v.eq('true')


def load_data(path, c):
    d = pd.read_csv(path, dtype={'subject': str, 'group': str})
    need = {'subject', 'time', 'concentration', 'dose'}
    route = c['design']['route']
    if route == 'intravascular_infusion': need.add('infusion_duration')
    if c['nca']['lambda_z']['selection'] == 'declared_points': need.add('lambda_z_point')
    if need - set(d.columns): raise ValueError('Missing columns: ' + ', '.join(sorted(need - set(d.columns))))
    extra = set(d.columns) - need - {'group', 'blq', 'exclude', 'exclusion_reason'}
    if extra: raise ValueError('Undeclared columns: ' + ', '.join(sorted(extra)))
    if d.subject.isna().any(): raise ValueError('subject is required on every row')
    for col in ('time', 'dose') + (('infusion_duration',) if 'infusion_duration' in need else ()):
        v = pd.to_numeric(d[col], errors='coerce')
        if v.isna().any() or not np.isfinite(v).all(): raise ValueError(f'{col} must be numeric on every row')
        d[col] = v
    if (d.time < 0).any(): raise ValueError('Single-dose analysis: time is measured from the dose and cannot be negative')
    d['blq'] = _bool_column(d, 'blq') if 'blq' in d else False
    conc = pd.to_numeric(d.concentration, errors='coerce')
    if (d.blq & d.concentration.notna()).any(): raise ValueError('A BLQ row must have an empty concentration (the LLOQ belongs in the source, not the value)')
    if (~d.blq & (conc.isna() | ~np.isfinite(conc))).any(): raise ValueError('concentration must be numeric unless blq is true; use blq for values below the limit')
    if (conc < 0).any(): raise ValueError('concentration cannot be negative')
    d['concentration'] = conc
    if 'exclude' in d:
        d['exclude'] = _bool_column(d, 'exclude')
        if d.exclude.any() and ('exclusion_reason' not in d or d.loc[d.exclude, 'exclusion_reason'].isna().any()):
            raise ValueError('Every excluded sample needs an exclusion_reason')
    else: d['exclude'] = False
    if 'lambda_z_point' in d: d['lambda_z_point'] = _bool_column(d, 'lambda_z_point')
    if d.duplicated(['subject', 'time']).any(): raise ValueError('One sample per subject and time; resolve replicate or duplicate records first')
    for col in ('dose', 'group', 'infusion_duration'):
        if col in d and d.groupby('subject')[col].nunique(dropna=False).max() > 1:
            raise ValueError(f'{col} must be constant within subject (single dose per subject)')
    if (d.dose <= 0).any(): raise ValueError('dose must be positive')
    if 'infusion_duration' in d and (d.infusion_duration <= 0).any(): raise ValueError('infusion_duration must be positive')
    if (d.blq & ~d.exclude & (d.time > 0)).any() and c['blq'] is None:
        raise ValueError('Post-dose BLQ samples present: declare blq first/middle/last rules with a rationale')
    return d


def _payload(d, c):
    """Map the analysis rows to PKNCA input; returns the payload and a preprocessing log."""
    route = c['design']['route']; bolus = route == 'intravascular_bolus'; log = []
    use = d[~d.exclude].copy()
    rows = []
    for s, part in use.groupby('subject', sort=True):
        part = part.sort_values('time')
        z = part[part.time == 0]
        if bolus:
            if len(z) and not (bool(z.blq.iloc[0]) or z.concentration.iloc[0] == 0):
                raise ValueError(f'{s}: an IV bolus time-0 row must be the predose sample (BLQ or 0); C0 is back-extrapolated')
            post = part[part.time > 0]
            if post.empty or post.blq.iloc[0]: raise ValueError(f'{s}: the first post-dose sample after an IV bolus must be measurable')
            if not len(z): log.append(f'{s}: predose time-0 placeholder added for PKNCA pass 1 (replaced by C0 in pass 2)')
            rows.append({'subject': s, 'time': 0., 'concentration': 0., 'include_hl': False})
            part = post
        elif not len(z):
            if c['time_zero']['rule'] == 'observed': raise ValueError(f'{s}: no time-0 predose row; time_zero.rule is observed')
            log.append(f'{s}: time-0 concentration 0 inserted (time_zero insert_zero: {c["time_zero"]["rationale"]})')
            rows.append({'subject': s, 'time': 0., 'concentration': 0., 'include_hl': False})
        for _, x in part.iterrows():
            rows.append({'subject': s, 'time': float(x.time), 'concentration': 0. if x.blq else float(x.concentration),
                         'include_hl': bool(x.lambda_z_point) if 'lambda_z_point' in x else False})
    p = pd.DataFrame(rows)
    blq = c['blq'] or {'first': 'keep_as_zero', 'middle': 'drop', 'last': 'drop'}
    pk = {k: 'keep' if blq[k] == 'keep_as_zero' else 'drop' for k in ('first', 'middle', 'last')}
    if bolus: pk['first'] = 'keep'  # Only the predose placeholder precedes the first measurable sample.
    elif pk['first'] == 'drop' and (p[p.time == 0].concentration == 0).any():
        raise ValueError('blq.first drop would remove the time-0 predose zero, so AUC from 0 cannot be computed; declare keep_as_zero')
    subj = use.groupby('subject', sort=True).first()
    common = ['cmax', 'tmax', 'tlast', 'clast.obs', 'auclast', 'aumclast', 'half.life', 'aucinf.obs', 'aucpext.obs', 'aumcinf.obs', 'cl.obs', 'vz.obs']
    mrt = ['mrt.obs'] if route == 'extravascular' else ['mrt.iv.obs', 'vss.iv.obs']
    payload = {'subject': p.subject.tolist(), 'time': p.time.tolist(), 'concentration': p.concentration.tolist(),
               'declared_points': c['nca']['lambda_z']['selection'] == 'declared_points', 'include_half_life': p.include_hl.tolist(),
               'dose_subject': subj.index.tolist(), 'dose': subj.dose.astype(float).tolist(),
               'duration': (subj.infusion_duration.astype(float).tolist() if route == 'intravascular_infusion' else [0.]*len(subj)),
               'route': 'extravascular' if route == 'extravascular' else 'intravascular', 'bolus': bolus,
               'auc_method': AUC_METHODS[c['nca']['auc_method']], 'blq': pk,
               'min_points': c['nca']['lambda_z']['min_points'], 'adj_r2_factor': c['nca']['lambda_z']['adj_r2_factor'],
               'parameters_pass1': (['c0', 'cmax', 'tmax', 'tlast', 'clast.obs', 'aucivlast', 'half.life', 'aucivinf.obs'] if bolus else common + mrt),
               'parameters_pass2': [x for x in common if x not in ('cmax', 'tmax')] + mrt if bolus else []}
    return payload, log, p


def _table(rows):
    t = pd.DataFrame({'subject': rows['subject'], 'parameter': rows['parameter'], 'value': rows['value'], 'exclude': rows['exclude']})
    out = {}
    for _, x in t.iterrows():
        v = x.value
        out.setdefault(x.subject, {})[x.parameter] = (None if v is None or (isinstance(v, float) and not math.isfinite(v)) else float(v),
                                                      None if x.exclude is None or (isinstance(x.exclude, float) and math.isnan(x.exclude)) else str(x.exclude))
    return out


def _units(d):
    uc, ut, ud = d['concentration_unit'], d['time_unit'], d['dose_unit']
    return {'c0': uc, 'cmax': uc, 'clast': uc, 'clast_pred': uc, 'tmax': ut, 'tlast': ut, 'half_life': ut, 'mrt': ut,
            'lambda_z_time_first': ut, 'lambda_z_time_last': ut, 'lambda_z': '1/' + ut,
            'auclast': f'{uc}*{ut}', 'aucinf': f'{uc}*{ut}', 'aumclast': f'{uc}*{ut}^2', 'aumcinf': f'{uc}*{ut}^2',
            'cl': f'({ud})/({uc}*{ut})', 'vz': f'({ud})/({uc})', 'vss': f'({ud})/({uc})', 'aucinf_extrap_percent': '%',
            'r_squared': '1', 'adj_r_squared': '1', 'span_ratio': '1', 'lambda_z_n_points': 'count'}


def _labels(route):
    ev = route == 'extravascular'
    return {'cl': 'CL/F' if ev else 'CL', 'vz': 'Vz/F' if ev else 'Vz', 'mrt': 'MRT (including absorption)' if ev else 'MRT_iv', 'vss': 'Vss'}


def _fmt(v):
    return 'not available' if v is None else f'{v:.4g}'


def _subject_rows(d, c, res):
    route = c['design']['route']; bolus = route == 'intravascular_bolus'
    p1 = _table(res['pass1']['rows']); p2 = _table(res['pass2']['rows']) if bolus else {}
    rep = c['reportability']; desc = rep.get('descriptive_only') is True
    subjects = []; failing = []; checks = []
    info = d[~d.exclude].groupby('subject', sort=True).first()
    for s in info.index:
        a = p1.get(s, {}); b = p2.get(s, {}) if bolus else a
        val = {}; why = {}
        for ours, theirs in PKNCA_NAMES.items():
            src = a if (not bolus or ours in ('c0', 'cmax', 'tmax')) else b
            v, e = src.get(theirs, (None, None)); val[ours] = v
            if e: why[ours] = e
        mk = ('mrt.obs',) if route == 'extravascular' else ('mrt.iv.obs',)
        val['mrt'], e = (b.get(mk[0], (None, None))); why.update({'mrt': e} if e else {})
        if route != 'extravascular':
            val['vss'], e = b.get('vss.iv.obs', (None, None)); why.update({'vss': e} if e else {})
        else: val['vss'] = None
        diag = []
        if bolus:
            first_obs = d[(d.subject == s) & ~d.exclude & (d.time > 0)].sort_values('time').concentration.iloc[0]
            if val['c0'] is not None and abs(val['c0'] - first_obs) <= 1e-12*max(1, abs(first_obs)): diag.append('c0_not_log_back_extrapolated')
            l1, _ = a.get('lambda.z', (None, None)); auciv, _ = a.get('aucivinf.obs', (None, None))
            if l1 is not None and val['lambda_z'] is not None and abs(l1 - val['lambda_z']) <= 1e-12*abs(l1):
                if auciv is None or val['aucinf'] is None or abs(auciv - val['aucinf']) > 1e-9*abs(auciv):
                    raise ArithmeticError(f'{s}: pass-2 AUCinf does not reproduce PKNCA aucivinf.obs')
                checks.append({'subject': s, 'aucivinf_obs_pass1': auciv, 'aucinf_obs_pass2': val['aucinf'], 'equal': True})
            elif l1 is not None and val['lambda_z'] is not None:
                diag.append('lambda_z_selection_differs_between_passes'); checks.append({'subject': s, 'equal': None, 'reason': 'different lambda-z selections'})
        # Reportability of terminal and AUCinf-dependent parameters.
        status = {}
        for k in OBSERVED: status[k] = 'reportable' if val.get(k) is not None else ('not_applicable' if k == 'c0' and not bolus else 'not_estimated')
        declared_short = (c['nca']['lambda_z']['selection'] == 'declared_points' and val['lambda_z_n_points'] is not None
                          and val['lambda_z_n_points'] < c['nca']['lambda_z']['min_points'])
        if val['lambda_z'] is None:
            term = 'not_estimated'; reasons = [why.get('lambda_z') or why.get('half_life') or 'lambda_z not estimated by PKNCA']
        elif declared_short:
            # PKNCA fits any number of manually selected points; the declared minimum still applies.
            term = 'withheld'; reasons = [f"{int(val['lambda_z_n_points'])} declared lambda-z points < min_points {c['nca']['lambda_z']['min_points']}"]
        elif desc:
            term = 'not_assessed'; reasons = []
        else:
            reasons = []
            if rep['min_adj_r2'] is not None and (val['adj_r_squared'] is None or val['adj_r_squared'] < rep['min_adj_r2']): reasons.append(f"adjusted R2 {_fmt(val['adj_r_squared'])} < {rep['min_adj_r2']}")
            if rep['min_span_ratio'] is not None and (val['span_ratio'] is None or val['span_ratio'] < rep['min_span_ratio']): reasons.append(f"span ratio {_fmt(val['span_ratio'])} < {rep['min_span_ratio']}")
            term = 'withheld' if reasons else 'reportable'
        for k in TERMINAL: status[k] = term if val.get(k) is not None or term != 'reportable' else 'not_estimated'
        inf = term; inf_reasons = list(reasons)
        if term == 'reportable' and rep.get('max_aucinf_extrap_percent') is not None and val['aucinf_extrap_percent'] is not None \
                and val['aucinf_extrap_percent'] > rep['max_aucinf_extrap_percent']:
            inf = 'withheld'; inf_reasons.append(f"AUCinf extrapolated {val['aucinf_extrap_percent']:.3g}% > {rep['max_aucinf_extrap_percent']}%")
        for k in INF:
            if k == 'vss' and route == 'extravascular': status[k] = 'not_applicable'; continue
            status[k] = inf if (val.get(k) is not None or inf != 'reportable') else 'not_estimated'
        for k, st in status.items():
            if st in ('withheld', 'not_estimated') and k in ('half_life', 'aucinf', 'cmax', 'auclast'):
                failing.append({'subject': s, 'parameter': k, 'status': st, 'reason': '; '.join(inf_reasons if k == 'aucinf' else reasons) or why.get(k) or st})
        row = {'subject': s, 'group': (str(info.loc[s, 'group']) if 'group' in info.columns else None), 'dose': float(info.loc[s, 'dose']),
               'values': val, 'status': status, 'terminal_reasons': reasons, 'aucinf_reasons': inf_reasons, 'pknca_exclusions': why, 'diagnostics': diag}
        if route == 'intravascular_infusion': row['infusion_duration'] = float(info.loc[s, 'infusion_duration'])
        subjects.append(row)
    return subjects, failing, checks


def _gstats(x):
    x = np.asarray(x, float); n = len(x); out = {'n': n, 'mean': float(x.mean()) if n else None, 'sd': float(x.std(ddof=1)) if n > 1 else None}
    if n and (x > 0).all():
        lx = np.log(x); out['geometric_mean'] = float(np.exp(lx.mean()))
        out['geometric_cv_percent'] = float(100*math.sqrt(math.expm1(lx.var(ddof=1)))) if n > 1 else None
    else: out['geometric_mean'] = None; out['geometric_cv_percent'] = None
    return out


def _summaries(subjects, route):
    keys = [k for k in OBSERVED + ('lambda_z', 'half_life') + INF if not (k == 'vss' and route == 'extravascular') and not (k == 'c0' and route != 'intravascular_bolus')]
    by = {}
    for r in subjects: by.setdefault(r['group'] if r['group'] is not None else f"dose {r['dose']:g}", []).append(r)
    out = []
    for g, rows in by.items():
        item = {'group': g, 'doses': sorted({r['dose'] for r in rows}), 'n_subjects': len(rows), 'parameters': {}}
        for k in keys:
            vals = [r['values'][k] for r in rows if r['status'][k] == 'reportable']
            st = {'n_not_reportable': sum(r['status'][k] != 'reportable' for r in rows)}
            if k in ('tmax', 'tlast'):
                st.update(n=len(vals), median=float(np.median(vals)) if vals else None, min=min(vals) if vals else None, max=max(vals) if vals else None)
            else: st.update(_gstats(vals))
            item['parameters'][k] = st
        out.append(item)
    return out


def dose_proportionality(dose, value, acceptance, level):
    """Smith et al. (2000) power model on one value per subject."""
    dose = np.asarray(dose, float); y = np.asarray(value, float); n = len(y)
    X = np.column_stack([np.ones(n), np.log(dose)]); beta, *_ = np.linalg.lstsq(X, np.log(y), rcond=None)
    res = np.log(y) - X @ beta; df = n - 2; s2 = float(res @ res)/df; cov = s2*np.linalg.inv(X.T @ X)
    se = math.sqrt(cov[1, 1]); q = stats.t.ppf(.5 + level/2, df); lo, hi = beta[1] - q*se, beta[1] + q*se
    r = float(dose.max()/dose.min()); acc = [1 + math.log(acceptance[0])/math.log(r), 1 + math.log(acceptance[1])/math.log(r)]
    concl = ('proportional_within_acceptance' if acc[0] <= lo and hi <= acc[1] else
             'disproportionate' if hi < acc[0] or lo > acc[1] else 'inconclusive')
    return {'n_subjects': n, 'n_dose_levels': int(len(np.unique(dose))), 'intercept': float(beta[0]), 'slope_b': float(beta[1]),
            'slope_se': se, 'df': df, 'slope_ci': [float(lo), float(hi)], 'confidence_level': level, 'dose_ratio_r': r,
            'acceptance_ratio': list(acceptance), 'acceptance_interval_b': acc, 'conclusion': concl,
            'dose_normalised_ratio_over_range': float(r**(beta[1] - 1)),
            'dose_normalised_ratio_ci': [float(r**(lo - 1)), float(r**(hi - 1))]}


def compute(d, c):
    from .r_bridge import call_pk
    route = c['design']['route']
    payload, log, mapped = _payload(d, c)
    res = call_pk('pknca', payload)
    subjects, failing, checks = _subject_rows(d, c, res)
    labels = _labels(route)
    out_dp = None
    if c.get('dose_proportionality'):
        dp = c['dose_proportionality']; out_dp = {'parameters': {}, 'source': dp['source'], 'rationale': dp['rationale']}
        for k in dp['parameters']:
            rows = [r for r in subjects if r['status'][k] == 'reportable' and r['values'][k] and r['values'][k] > 0]
            ndose = len({r['dose'] for r in rows})
            if len(rows) < 4 or ndose < 2:
                out_dp['parameters'][k] = {'status': 'withheld', 'reason': 'fewer_than_4_reportable_subjects_or_2_dose_levels', 'n_subjects': len(rows), 'n_dose_levels': ndose}
                failing.append({'parameter': 'dose_proportionality_' + k, 'status': 'withheld', 'reason': out_dp['parameters'][k]['reason']}); continue
            fit = dose_proportionality([r['dose'] for r in rows], [r['values'][k] for r in rows], dp['acceptance'], dp['confidence_level'])
            fit['status'] = 'estimated'; fit['n_excluded_not_reportable'] = len(subjects) - len(rows)
            if ndose == 2: fit['diagnostics'] = ['two_dose_levels_only']
            out_dp['parameters'][k] = fit
    excluded = d[d.exclude][['subject', 'time', 'exclusion_reason']].to_dict('records') if d.exclude.any() else []
    primary = {'route': route, 'parameter_labels': labels, 'units': _units(c['design']), 'subjects': subjects,
               'summaries': _summaries(subjects, route), 'excluded_samples': excluded, 'preprocessing_log': log,
               'pknca': {'version': res['r_environment']['packages'].get('PKNCA'), 'auc_method': payload['auc_method'],
                         'conc_blq': payload['blq'], 'min_hl_points': payload['min_points'], 'adj_r2_factor': payload['adj_r2_factor'],
                         'allow_tmax_in_half_life': False, 'c0_method_order': res.get('pknca_c0_method'),
                         'warnings': sorted(set(res['pass1']['warnings']) | set(res.get('pass2', {}).get('warnings', [])))},
               'bolus_two_pass_checks': checks if route == 'intravascular_bolus' else None,
               'counts': {'subjects': len(subjects), 'terminal_reportable': sum(r['status']['lambda_z'] == 'reportable' for r in subjects),
                          'terminal_withheld': sum(r['status']['lambda_z'] == 'withheld' for r in subjects),
                          'terminal_not_estimated': sum(r['status']['lambda_z'] == 'not_estimated' for r in subjects),
                          'aucinf_reportable': sum(r['status']['aucinf'] == 'reportable' for r in subjects)}}
    if out_dp is not None: primary['dose_proportionality'] = out_dp
    return {'analysis_type': 'pk_nca', 'primary': primary, 'mapped_input': mapped.to_dict('list'),
            'r_environment': res['r_environment'], 'must_mention': _notes(c, primary, labels, failing),
            'failing_items': failing,
            'limitations': ['Single dose at time 0 per subject only: no multiple-dose, steady-state, partial-interval, urine or sparse/destructive-sampling NCA.',
                            'Units are composed labels; no unit conversion is performed. mg and mg/kg doses give different CL and V units.',
                            'NCA describes the observed profiles; it does not identify linear PK, TMDD, absorption mechanisms or bioavailability.',
                            'Calculations by R PKNCA; no WinNonlin, Phoenix or Prism equivalence beyond the recorded published example.']}


def _notes(c, p, labels, failing):
    d = c['design']; rep = c['reportability']; lz = c['nca']['lambda_z']
    notes = [f"NCA by R PKNCA {p['pknca']['version']} for {d['analyte']} in {d['matrix']} ({d['species']}; unit: {d['subject_unit']}), route {d['route']}, {d['time_basis']} times; the subject is the unit and summaries use only reportable subject values.",
             f"AUC method {c['nca']['auc_method']} ({c['nca']['auc_method_rationale']}). Terminal phase: " +
             ('PKNCA best adjusted R2 rule over points after Tmax (at least ' + str(lz['min_points']) + ' points, factor ' + str(lz['adj_r2_factor']) + ').' if lz['selection'] == 'best_adjusted_r2'
              else 'declared lambda-z points (lambda_z_point column).')]
    if d['route'] == 'extravascular':
        notes.append('Extravascular dosing with unknown bioavailability: clearance and volume are CL/F and Vz/F; MRT includes absorption time. No Vss is reported.')
    elif d['route'] == 'intravascular_bolus':
        notes.append('IV bolus: C0 back-extrapolated by PKNCA (method order ' + ', '.join(p['pknca']['c0_method_order'] or []) + ') and inserted at time 0 for AUC, AUMC, CL, Vz, MRT and Vss; Cmax and Tmax are from observed samples.')
    else:
        notes.append('IV infusion: MRT_iv subtracts half the infusion duration; Vss = CL x MRT_iv.')
    if c['blq']: notes.append(f"BLQ rule (declared): first {c['blq']['first']}, middle {c['blq']['middle']}, last {c['blq']['last']} ({c['blq']['rationale']}).")
    if rep.get('descriptive_only'):
        notes.append('Terminal-phase parameters are descriptive only (no declared reportability criteria: ' + rep['rationale'] + '); they are excluded from summaries.')
    else:
        notes.append(f"Reportability criteria declared with source ({rep['source']}): adjusted R2 >= {rep['min_adj_r2']}, span ratio >= {rep['min_span_ratio']}, AUCinf extrapolated <= {rep['max_aucinf_extrap_percent']}%. Withheld lambda-z withholds t1/2, AUCinf, {labels['cl']}, {labels['vz']}, MRT and Vss; Cmax, Tmax and AUClast remain.")
    k = p['counts']
    notes.append(f"Terminal phase: {k['terminal_reportable']} reportable, {k['terminal_withheld']} withheld, {k['terminal_not_estimated']} not estimated of {k['subjects']} subjects.")
    if p['preprocessing_log']: notes.append('Preprocessing: ' + '; '.join(p['preprocessing_log'][:6]) + ('; ...' if len(p['preprocessing_log']) > 6 else '') + '.')
    if p['excluded_samples']: notes.append(f"{len(p['excluded_samples'])} sample(s) excluded by declared exclusion with reasons (listed in results).")
    if p.get('dose_proportionality'):
        for k2, f in p['dose_proportionality']['parameters'].items():
            if f['status'] == 'estimated':
                notes.append(f"Dose proportionality ({k2}, power model ln P = a + b ln dose, {f['n_subjects']} subjects, {f['n_dose_levels']} dose levels): b = {f['slope_b']:.3f}, {int(100*f['confidence_level'])}% CI {f['slope_ci'][0]:.3f} to {f['slope_ci'][1]:.3f}; acceptance for b {f['acceptance_interval_b'][0]:.3f} to {f['acceptance_interval_b'][1]:.3f} (ratio {f['acceptance_ratio'][0]}-{f['acceptance_ratio'][1]} over a {f['dose_ratio_r']:g}-fold range): {f['conclusion']}. This describes exposure scaling over the tested range; it does not identify a mechanism, and inconclusive is not evidence of proportionality.")
            else:
                notes.append(f"Dose proportionality ({k2}) withheld: {f['reason']}.")
    if p['pknca']['warnings']: notes.append('PKNCA warnings retained: ' + '; '.join(p['pknca']['warnings'][:4]) + '.')
    return notes


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    p = result['primary']; token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True)
    m = pd.DataFrame(result['mapped_input']); names = []
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), layout='constrained')
        groups = list(dict.fromkeys(r['group'] or f"dose {r['dose']:g}" for r in p['subjects']))
        for r in p['subjects']:
            g = groups.index(r['group'] or f"dose {r['dose']:g}"); part = m[(m.subject == r['subject']) & (m.time > 0)]
            pos = part[part.concentration > 0]
            for ax in axes: ax.plot(pos.time, pos.concentration, marker='o', ms=2.5, lw=.8, color=pstyle.color(token, g), label=groups[g] if r is next(x for x in p['subjects'] if (x['group'] or f"dose {x['dose']:g}") == groups[g]) else None)
            v = r['values']
            if v['lambda_z'] is not None and r['status']['lambda_z'] in ('reportable', 'not_assessed', 'withheld'):
                t = np.linspace(v['lambda_z_time_first'], v['lambda_z_time_last'], 20)
                ax = axes[1]; ax.plot(t, v['clast_pred']*np.exp(v['lambda_z']*(v['lambda_z_time_last'] - t)), ls='--', lw=.8,
                                      color=pstyle.MUTED if r['status']['lambda_z'] == 'withheld' else pstyle.color(token, g))
        axes[0].set(xlabel=f"Time ({result['primary']['units']['tmax']})", ylabel=f"Concentration ({p['units']['cmax']})", title='Observed profiles (linear)')
        axes[1].set_yscale('log'); axes[1].set(xlabel=f"Time ({p['units']['tmax']})", title='Log scale; dashed = lambda-z fit (grey = withheld)')
        axes[0].legend(frameon=False, fontsize=7)
        pstyle.save(fig, out/'pk_profiles', token, 300); names.append('pk_profiles')
        dp = p.get('dose_proportionality')
        if dp and any(f['status'] == 'estimated' for f in dp['parameters'].values()):
            est = {k: f for k, f in dp['parameters'].items() if f['status'] == 'estimated'}
            fig, axes = plt.subplots(1, len(est), figsize=(3.6*len(est), 3.2), layout='constrained', squeeze=False)
            for ax, (k, f) in zip(axes[0], est.items()):
                rows = [r for r in p['subjects'] if r['status'][k] == 'reportable' and r['values'][k]]
                x = np.array([r['dose'] for r in rows]); y = np.array([r['values'][k] for r in rows])
                ax.scatter(x, y, s=14, color=pstyle.color(token, 0)); xx = np.geomspace(x.min(), x.max(), 50)
                ax.plot(xx, np.exp(f['intercept'])*xx**f['slope_b'], color=pstyle.color(token, 1), lw=1)
                ax.set_xscale('log'); ax.set_yscale('log'); ax.set(xlabel='Dose', ylabel=k, title=f"b = {f['slope_b']:.2f} ({f['conclusion']})")
            pstyle.save(fig, out/'dose_proportionality', token, 300); names.append('dose_proportionality')
    return names
