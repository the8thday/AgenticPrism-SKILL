"""Linear and quadratic standard curves with inverse prediction (0.13.2).

Per plate: least squares of response on concentration (declared linear or quadratic; weights none,
1/x or 1/x^2, declared), standard back-calculation recovery against declared limits, pure-error
lack-of-fit F when standards are replicated, and inverse prediction of unknowns as the mean of their
replicate responses. Inverse intervals use the inversion method (as investr::invest with
interval = "inversion"): all x with |ybar0 - f(x)| <= t * sqrt(s^2 / m + Var f(x)), m the number of
unknown replicates, with s^2 pooled with the unknown's replicate variance (df + m - 1; unweighted fits). Concentrations outside the standard range are not extrapolated.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import optimize, stats

TYPES = {'standard_curve'}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'standard_curve', 'input': None, 'source': None,
            'assay': {'analyte': None, 'readout': None, 'concentration_unit': None, 'rationale': None},
            'model': {'form': None, 'form_rationale': None, 'weighting': 'none'},
            'acceptance': {'recovery_percent': None, 'recovery_percent_lowest': None, 'min_passing_fraction': None, 'source': None,
                           'lack_of_fit_alpha': .05},
            'uncertainty': {'level': .95}, 'report': {'plot_style': 'prism_like'}}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def _num(v, name, lo, hi):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not lo < v <= hi:
        raise ValueError(f'{name} out of range')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported standard-curve configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'standard_curve' or c['schema_version'] != 1:
        raise ValueError('Unsupported standard-curve schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    for k in ('analyte', 'readout', 'concentration_unit', 'rationale'):
        _text(c['assay'][k], 'assay.' + k)
    m = c['model']
    if m['form'] not in ('linear', 'quadratic'):
        raise ValueError('model.form must be linear or quadratic (sigmoid curves belong to elisa-quantification)')
    _text(m['form_rationale'], 'model.form_rationale (validated range or SOP; never chosen from this plate)')
    if m['weighting'] not in ('none', '1/x', '1/x^2'):
        raise ValueError('weighting must be none, 1/x or 1/x^2')
    a = c['acceptance']
    _num(a['recovery_percent'], 'recovery_percent', 0, 100); _num(a['recovery_percent_lowest'], 'recovery_percent_lowest', 0, 100)
    _num(a['min_passing_fraction'], 'min_passing_fraction', 0, 1); _text(a['source'], 'acceptance.source')
    _num(a['lack_of_fit_alpha'], 'lack_of_fit_alpha', 0, .5); _num(c['uncertainty']['level'], 'uncertainty.level', .5, .999)
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'plate_id', 'role', 'sample_id', 'concentration', 'response'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    if set(d.role) - {'standard', 'unknown'}:
        raise ValueError('role must be standard or unknown')
    d['response'] = pd.to_numeric(d.response, errors='raise')
    if not np.isfinite(d.response).all():
        raise ValueError('Nonfinite response')
    std = d.role == 'standard'
    if (d.concentration[std].str.strip() == '').any():
        raise ValueError('Standards need a nominal concentration')
    d['concentration'] = pd.to_numeric(d.concentration.where(std, None), errors='raise')
    if (d.concentration[std] < 0).any():
        raise ValueError('Negative standard concentration')
    d['dilution_factor'] = pd.to_numeric(d.get('dilution_factor', pd.Series('1', index=d.index)).replace('', '1'), errors='raise')
    if (d.dilution_factor < 1).any():
        raise ValueError('dilution_factor must be >= 1')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    for pid, g in d.groupby('plate_id'):
        s = g[(g.role == 'standard') & ~g.exclude]
        if s.concentration.nunique() < 4:
            raise ValueError(f'Plate {pid}: at least four distinct standard levels are required')
        if (g.role == 'unknown').any() and (g[g.role == 'unknown'].groupby('sample_id').dilution_factor.nunique() > 1).any():
            raise ValueError(f'Plate {pid}: one dilution per unknown sample_id; give each dilution its own sample_id')
    return d


def _design(x, form):
    x = np.asarray(x, float)
    return np.column_stack([np.ones_like(x), x] + ([x ** 2] if form == 'quadratic' else []))


def _weights(x, mode):
    x = np.asarray(x, float)
    if mode == 'none':
        return np.ones_like(x)
    if np.any(x <= 0):
        raise ValueError('1/x weighting needs positive standard concentrations (no zero standard)')
    return 1 / x if mode == '1/x' else 1 / x ** 2


def fit_plate(s, cfg):
    form, mode = cfg['model']['form'], cfg['model']['weighting']
    x, y = s.concentration.to_numpy(float), s.response.to_numpy(float); w = _weights(x, mode)
    X = _design(x, form); sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    resid = y - X @ beta; df = len(y) - X.shape[1]; s2 = float(np.sum(w * resid ** 2) / df)
    cov = s2 * np.linalg.inv((X * w[:, None]).T @ X)
    out = {'coefficients': beta.tolist(), 'covariance': cov.tolist(), 'residual_variance': s2, 'residual_df': int(df),
           'r_squared': float(1 - np.sum(w * resid ** 2) / np.sum(w * (y - np.average(y, weights=w)) ** 2)), 'lack_of_fit': None}
    levels = pd.Series(y).groupby(x).agg(['count', 'var']); pe_df = int((levels['count'] - 1).sum())
    if pe_df > 0 and len(levels) > X.shape[1]:
        sse_pe = float(sum(np.sum(w[x == lv] * (y[x == lv] - np.average(y[x == lv], weights=w[x == lv])) ** 2) for lv in levels.index))
        sse = float(np.sum(w * resid ** 2)); lof_df = len(levels) - X.shape[1]
        f = ((sse - sse_pe) / lof_df) / (sse_pe / pe_df) if sse_pe > 0 else np.inf
        out['lack_of_fit'] = {'f_statistic': float(f), 'df': [lof_df, pe_df], 'p_value': float(stats.f.sf(f, lof_df, pe_df))}
    return out, beta, cov, s2, df


def predict(xs, beta):
    return _design(xs, 'quadratic' if len(beta) == 3 else 'linear') @ beta


def inverse(y0, beta, lo, hi):
    """Concentration in [lo, hi] whose fitted response equals y0, or None (monotone fits only)."""
    f = lambda v: float(predict([v], beta)[0]) - y0
    a, b = f(lo), f(hi)
    if a == 0:
        return lo
    if b == 0:
        return hi
    if a * b > 0:
        return None
    return float(optimize.brentq(f, lo, hi, xtol=1e-12))


def inversion_interval(ybar, m, beta, cov, s2, df, lo, hi, level, y0=None):
    """With replicate unknown responses y0, sigma^2 is pooled with their within-sample variance and the df grow
    by m - 1 (Graybill 1976; investr::invest)."""
    if y0 is not None and m > 1:
        y0 = np.asarray(y0, float); pooled = (s2 * df + float(np.sum((y0 - y0.mean()) ** 2))) / (df + m - 1)
        cov = np.asarray(cov) * pooled / s2; s2, df = pooled, df + m - 1
    t = stats.t.ppf((1 + level) / 2, df)
    g = lambda v: (float(predict([v], beta)[0]) - ybar) ** 2 - t ** 2 * (s2 / m + float(_design([v], 'quadratic' if len(beta) == 3 else 'linear')[0] @ cov @ _design([v], 'quadratic' if len(beta) == 3 else 'linear')[0]))
    x0 = inverse(ybar, beta, lo, hi)
    if x0 is None:
        return [None, None], 'point_outside_range'
    span = hi - lo; ends = []
    for direction in (-1, 1):
        prev, root = x0, None
        for cur in np.linspace(x0, x0 + direction * 2 * span, 401)[1:]:
            if g(cur) > 0:
                root = float(optimize.brentq(g, min(prev, cur), max(prev, cur), xtol=1e-12)); break
            prev = cur
        ends.append(root)
    if None in ends:
        return ends, 'open'
    return ends, 'closed' if lo <= ends[0] and ends[1] <= hi else 'extends_beyond_standards'


def compute(d, cfg):
    acc, lv = cfg['acceptance'], cfg['uncertainty']['level']; plates, unknowns, must, failing = [], [], [], []
    for pid, g in d.groupby('plate_id', sort=False):
        s = g[(g.role == 'standard') & ~g.exclude]; fit, beta, cov, s2, df = fit_plate(s, cfg)
        lo, hi = float(s.concentration.min()), float(s.concentration.max()); diag = []
        grid = np.linspace(lo, hi, 400); slope = np.diff(predict(grid, beta))
        if not (np.all(slope > 0) or np.all(slope < 0)):
            diag.append('curve_not_monotone_over_standard_range')
        if fit['lack_of_fit'] and fit['lack_of_fit']['p_value'] < acc['lack_of_fit_alpha']:
            diag.append('lack_of_fit')
        back = []
        for conc, ss in s.groupby('concentration'):
            ybar = float(ss.response.mean()); xb = inverse(ybar, beta, lo, hi) if conc > 0 else None
            limit = acc['recovery_percent_lowest'] if conc == s.concentration[s.concentration > 0].min() else acc['recovery_percent']
            rec = None if xb is None or conc == 0 else 100 * xb / conc
            back.append({'nominal': float(conc), 'n': int(len(ss)), 'mean_response': ybar, 'back_calculated': xb, 'recovery_percent': rec,
                         'limit_percent': limit, 'passes': None if conc == 0 else bool(rec is not None and abs(rec - 100) <= limit)})
        scored = [b for b in back if b['passes'] is not None]; frac = sum(b['passes'] for b in scored) / len(scored)
        if frac < acc['min_passing_fraction']:
            diag.append('standards_fail_recovery_rule')
        accepted = not diag
        plates.append({'plate_id': str(pid), 'accepted': accepted, 'diagnostics': diag, 'fit': fit, 'standards': back, 'passing_fraction': frac,
                       'standard_range': [lo, hi], 'form': cfg['model']['form'], 'weighting': cfg['model']['weighting']})
        if not accepted:
            failing.append({'plate_id': str(pid), 'diagnostics': diag})
        for sid, u in g[(g.role == 'unknown') & ~g.exclude].groupby('sample_id', sort=False):
            ybar = float(u.response.mean()); m = len(u); dil = float(u.dilution_factor.iloc[0]); xw = inverse(ybar, beta, lo, hi)
            row = {'plate_id': str(pid), 'sample_id': sid, 'replicates': m, 'mean_response': ybar, 'dilution_factor': dil, 'reportable': False,
                   'in_well_concentration': None, 'concentration': None, 'interval': [None, None]}
            if not accepted:
                row['status'] = 'withheld_calibration_not_accepted'
            elif xw is None:
                low_side = float(predict([lo], beta)[0]); high_side = float(predict([hi], beta)[0]); below = (ybar < min(low_side, high_side)) == (low_side < high_side)
                row['status'] = 'below_standard_range' if below else 'above_standard_range'
            else:
                row.update(in_well_concentration=xw, concentration=xw * dil, status='quantified', reportable=True)
                if cfg['model']['weighting'] == 'none':
                    iv, st = inversion_interval(ybar, m, beta, cov, s2, df, lo, hi, lv, u.response.to_numpy(float))
                    row['interval'] = [None if v is None else v * dil for v in iv]; row['interval_status'] = st
                else:
                    row['interval_status'] = 'not_computed_for_weighted_fit'
            unknowns.append(row)
    must.append('Concentrations are interpolated only inside the standard range and multiplied by the declared dilution; intervals describe calibration and replicate noise on this plate only.')
    if cfg['model']['weighting'] != 'none':
        must.append('Weighted fits report point concentrations only; inverse intervals are computed for unweighted fits.')
    if any(u['status'].endswith('standard_range') for u in unknowns):
        must.append('Some unknowns fall outside the standard range; re-assay at another dilution instead of extrapolating.')
    return {'analysis_type': 'standard_curve', 'primary': {'plates': plates, 'unknowns': unknowns, 'assay': cfg['assay'], 'model': cfg['model'], 'acceptance': acc, 'level': lv},
            'must_mention': must, 'failing_items': failing,
            'limitations': ['Linear or quadratic only; sigmoid immunoassay curves belong to elisa-quantification.',
                            'Intervals assume independent Gaussian errors with constant variance (unweighted); no between-plate or matrix error.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); names = []
    d = pd.read_csv(Path(run)/'input.csv', dtype={'plate_id': str}, keep_default_na=False); p = result['primary']
    for i, pl in enumerate(p['plates']):
        s = d[(d.plate_id == pl['plate_id']) & (d.role == 'standard')]; lo, hi = pl['standard_range']; xs = np.linspace(lo, hi, 200)
        with plt.rc_context(pstyle.rc(token, *font_setup())):
            fig, ax = plt.subplots(figsize=(5.4, 3.8), layout='constrained')
            ax.plot(pd.to_numeric(s.concentration), s.response, 'o', ms=5, color=pstyle.color(token, 0), label='Standards')
            ax.plot(xs, predict(xs, np.array(pl['fit']['coefficients'])), '-', color=pstyle.color(token, 0), lw=1.4, label=f"{pl['form']} fit")
            u = [r for r in p['unknowns'] if r['plate_id'] == pl['plate_id'] and r['in_well_concentration'] is not None]
            ax.plot([r['in_well_concentration'] for r in u], [r['mean_response'] for r in u], 's', ms=5, mfc='white', color=pstyle.color(token, 1), label='Unknowns (in well)')
            ax.set(xlabel=f"Concentration ({p['assay']['concentration_unit']})", ylabel=p['assay']['readout'], title=f"Plate {pl['plate_id']}: {'accepted' if pl['accepted'] else 'not accepted'}")
            ax.legend(frameon=False, fontsize=7); name = f'standard_curve_{i+1}'; pstyle.save(fig, out/name, token, 300); names.append(name)
    return names
