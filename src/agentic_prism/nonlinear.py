"""Library of common nonlinear models fitted per curve by least squares (0.13.3).

Models (natural parameters; rates are optimized on the log scale and must be positive):
  one_phase_decay        y = plateau + (y0 - plateau) exp(-k x)            half_life = ln 2 / k
  one_phase_association  y = y0 + (plateau - y0)(1 - exp(-k x))            half_time = ln 2 / k
  two_phase_decay        y = plateau + span_fast exp(-k_fast x) + span_slow exp(-k_slow x), k_fast > k_slow
  exponential_growth     y = y0 exp(k x)                                    doubling_time = ln 2 / k
  michaelis_menten       y = vmax x / (km + x)
  logistic_growth        y = asym / (1 + exp((xmid - x) / scal))            rate = 1 / scal
Weights none or 1/y^2 (fixed, from observed y), declared. Each parameter gets a profile-F interval
(all other parameters re-optimized; S <= S_min (1 + F(1, df) / df), the same set as R confint.nls).
Model-specific reportability gates: plateau support for one-phase fits, a two- vs one-phase F test,
saturation for Michaelis-Menten and asymptote support for logistic growth.
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq, least_squares

TYPES = {'nonlinear_fit'}
LN2 = np.log(2.)
MODELS = {
    'one_phase_decay': (['y0', 'plateau', 'k'], ['k']),
    'one_phase_association': (['y0', 'plateau', 'k'], ['k']),
    'two_phase_decay': (['plateau', 'span_fast', 'span_slow', 'k_fast', 'k_slow'], ['k_fast', 'k_slow']),
    'exponential_growth': (['y0', 'k'], ['k']),
    'michaelis_menten': (['vmax', 'km'], ['vmax', 'km']),
    'logistic_growth': (['asym', 'xmid', 'scal'], ['scal']),
}
DEFAULTS = {'schema_version': 1, 'analysis_type': 'nonlinear_fit', 'input': None, 'source': None,
            'model': {'name': None, 'rationale': None, 'weighting': 'none'},
            'design': {'x_label': None, 'x_unit': None, 'y_label': None, 'y_unit': None, 'rationale': None},
            'replicates': {'independent_unit': 'none', 'independence_source': None},
            'uncertainty': {'level': .95}, 'report': {'plot_style': 'prism_like'}}
STRUCTURE_ALPHA = .01


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'Declare {name}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported nonlinear-fit configuration keys')
    c = deepcopy(DEFAULTS)
    for k, v in raw.items():
        if isinstance(c[k], dict):
            if not isinstance(v, dict) or set(v) - set(c[k]):
                raise ValueError(f'Unsupported {k} settings')
            c[k].update(v)
        else:
            c[k] = v
    if c['analysis_type'] != 'nonlinear_fit' or c['schema_version'] != 1:
        raise ValueError('Unsupported nonlinear-fit schema')
    _text(c['input'], 'input'); _text(c['source'], 'source')
    m = c['model']
    if m['name'] not in MODELS:
        raise ValueError('model.name must be one of ' + ', '.join(MODELS) + ' (sigmoid dose-response belongs to dose-response)')
    _text(m['rationale'], 'model.rationale (mechanism or prior data; never chosen by comparing fits on these data)')
    if m['weighting'] not in ('none', '1/y^2'):
        raise ValueError('model.weighting must be none or 1/y^2')
    for k in ('x_label', 'x_unit', 'y_label', 'y_unit', 'rationale'):
        _text(c['design'][k], 'design.' + k)
    r = c['replicates']
    if r['independent_unit'] not in ('none', 'experiment_id'):
        raise ValueError('replicates.independent_unit must be none or experiment_id')
    if r['independent_unit'] == 'experiment_id':
        _text(r['independence_source'], 'replicates.independence_source')
    lv = c['uncertainty']['level']
    if isinstance(lv, bool) or not isinstance(lv, (int, float)) or not .5 < lv < 1:
        raise ValueError('uncertainty.level must lie in (0.5, 1)')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    d = pd.read_csv(path, dtype=str, keep_default_na=False)
    req = {'sample_id', 'experiment_id', 'curve_id', 'x', 'y'}
    if d.empty or req - set(d):
        raise ValueError(f'Missing columns: {sorted(req - set(d))}')
    for k in ('x', 'y'):
        d[k] = pd.to_numeric(d[k], errors='raise')
        if not np.isfinite(d[k]).all():
            raise ValueError(f'Nonfinite {k}')
    d['exclude'] = d.get('exclude', pd.Series('false', index=d.index)).str.lower()
    d['exclusion_reason'] = d.get('exclusion_reason', pd.Series('', index=d.index))
    if not d.exclude.isin(('true', 'false')).all() or ((d.exclude == 'true') & d.exclusion_reason.str.strip().eq('')).any():
        raise ValueError('Exclusions need true/false and a reason')
    d['exclude'] = d.exclude.eq('true')
    for cid, g in d.groupby('curve_id'):
        if g.sample_id.nunique() != 1 or g.experiment_id.nunique() != 1:
            raise ValueError(f'Curve {cid} mixes samples or experiments')
    if cfg['model']['weighting'] == '1/y^2' and (d.y[~d.exclude] <= 0).any():
        raise ValueError('1/y^2 weighting needs positive y')
    if cfg['model']['name'] == 'michaelis_menten' and (d.x < 0).any():
        raise ValueError('Michaelis-Menten needs nonnegative substrate concentrations')
    return d


# ---------- model functions on natural parameters ----------

def predict(name, p, x):
    x = np.asarray(x, float)
    if name == 'one_phase_decay':
        return p['plateau'] + (p['y0'] - p['plateau']) * np.exp(-p['k'] * x)
    if name == 'one_phase_association':
        return p['y0'] + (p['plateau'] - p['y0']) * (1 - np.exp(-p['k'] * x))
    if name == 'two_phase_decay':
        return p['plateau'] + p['span_fast'] * np.exp(-p['k_fast'] * x) + p['span_slow'] * np.exp(-p['k_slow'] * x)
    if name == 'exponential_growth':
        return p['y0'] * np.exp(p['k'] * x)
    if name == 'michaelis_menten':
        return p['vmax'] * x / (p['km'] + x)
    return p['asym'] / (1 + np.exp((p['xmid'] - x) / p['scal']))


def _to_natural(name, z):
    names, logs = MODELS[name]; p = {}
    with np.errstate(over='ignore'):  # an overflowing rate becomes inf and is rejected by the residual guard
        for n, v in zip(names, z):
            p[n] = float(np.exp(v)) if n in logs else float(v)
    return p


def _to_internal(name, p):
    names, logs = MODELS[name]; z = []
    for n in names:
        v = p[n]
        z.append(np.log(v) if n in logs else v)
    return np.array(z, float)


def _starts(name, x, y):
    o = np.argsort(x); x, y = x[o], y[o]; span = max(float(np.ptp(x)), 1e-12); yr = max(float(np.ptp(y)), 1e-12)
    rates = [1 / span, 3 / span, 10 / span, .3 / span]
    if name == 'one_phase_decay':
        return [{'y0': y[0], 'plateau': y[-1], 'k': k} for k in rates]
    if name == 'one_phase_association':
        return [{'y0': y[0], 'plateau': y[-1], 'k': k} for k in rates]
    if name == 'two_phase_decay':
        return [{'plateau': y[-1], 'span_fast': (y[0] - y[-1]) * f, 'span_slow': (y[0] - y[-1]) * (1 - f), 'k_fast': kf, 'k_slow': kf / r}
                for f in (.3, .6) for kf in (5 / span, 15 / span) for r in (5, 20)]
    if name == 'exponential_growth':
        pos = y > 0
        k0 = np.polyfit(x[pos], np.log(y[pos]), 1)[0] if pos.sum() > 1 else 1 / span
        return [{'y0': max(y[0], 1e-12) if y[0] > 0 else yr / 10, 'k': max(k0, 1e-6 / span) * f} for f in (.5, 1, 2)]
    if name == 'michaelis_menten':
        return [{'vmax': y.max() * f, 'km': km} for f in (1, 1.5, 3) for km in (np.median(x[x > 0]) if (x > 0).any() else span, span / 10)]
    return [{'asym': y.max() * f, 'xmid': xm, 'scal': span / s} for f in (1.05, 1.5) for xm in (np.median(x), x[np.argmin(abs(y - y.max() / 2))]) for s in (4, 10)]


def fit_curve(g, cfg, _profile=True, model=None):
    name = model or cfg['model']['name']; names = MODELS[name][0]; lv = cfg['uncertainty']['level']
    u = g[~g.exclude].sort_values('x'); x, y = u.x.to_numpy(float), u.y.to_numpy(float)
    w = np.ones_like(y) if cfg['model']['weighting'] == 'none' else 1 / y ** 2; sw = np.sqrt(w)
    res = {'curve_id': str(g.curve_id.iloc[0]), 'sample_id': str(g.sample_id.iloc[0]), 'experiment_id': str(g.experiment_id.iloc[0]),
           'model': name, 'n_points': int(len(x)), 'status': 'failed', 'reportable': False, 'parameters': {}, 'derived': {}, 'diagnostics': []}
    df = len(x) - len(names)
    if df < 3 or len(np.unique(x)) < len(names) + 2:
        res['diagnostics'].append('too_few_points_for_model'); return res, [], []
    scale = max(float(np.max(np.abs(y))), 1e-12)
    def resid(z):
        q = _to_natural(name, z)
        if name == 'two_phase_decay' and q['k_fast'] <= q['k_slow'] * 1.0001:  # keep the phases ordered
            return np.full(len(y), 1e3)
        with np.errstate(over='ignore', invalid='ignore'):
            r = sw * (predict(name, q, x) - y) / scale
        return np.where(np.isfinite(r), r, 1e3)
    best = None
    for s in _starts(name, x, y):
        try:
            z0 = _to_internal(name, s)
        except (ValueError, FloatingPointError):
            continue
        if not np.all(np.isfinite(z0)):
            continue
        try:
            sol = least_squares(resid, z0, method='lm', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        except (ValueError, FloatingPointError):
            continue
        if np.all(np.isfinite(sol.fun)) and (best is None or sol.cost < best.cost):
            best = sol
    if best is None:
        res['diagnostics'].append('optimizer_failed'); return res, [], []
    p = _to_natural(name, best.x); fitted = predict(name, p, x); sse = float(np.sum(w * (fitted - y) ** 2))
    res.update(sse=sse, residual_df=int(df), rmse=float(np.sqrt(sse / df)), parameters={n: {'estimate': p[n]} for n in names}, status='estimated')
    if not _profile:
        return res, [], []
    d = res['diagnostics']
    rank = int(np.linalg.matrix_rank(best.jac)); jac_ok = rank == len(names)
    if not jac_ok:
        d.append('parameters_not_locally_identifiable')
    # runs test on residual signs (Wald-Wolfowitz), descriptive
    sgn = np.sign(y - fitted); sgn = sgn[sgn != 0]; n1, n2 = (sgn > 0).sum(), (sgn < 0).sum(); runs = 1 + int(np.sum(sgn[1:] != sgn[:-1]))
    if n1 > 0 and n2 > 0:
        mu = 2 * n1 * n2 / (n1 + n2) + 1; var = (mu - 1) * (mu - 2) / (n1 + n2 - 1)
        p_runs = float(stats.norm.cdf((runs - mu) / np.sqrt(var))) if var > 0 else 1.
        res['runs_test'] = {'runs': runs, 'expected': float(mu), 'p_too_few_runs': p_runs}
        if p_runs < .05:
            d.append('systematic_residual_pattern')
    # profile-F interval for each natural parameter
    threshold = sse * (1 + float(stats.f.ppf(lv, 1, df)) / df); z_hat = best.x.copy()
    for i, n in enumerate(names):
        def prof(value, i=i):
            def r(zr):
                z = np.insert(zr, i, value); return resid(z)
            zr0 = np.delete(z_hat, i); best_c = np.inf
            for start in (zr0, zr0 * 1.05):
                try:
                    s = least_squares(r, start, method='lm', xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=5000)
                    best_c = min(best_c, float(np.sum((s.fun * scale) ** 2)))
                except (ValueError, FloatingPointError):
                    pass
            return best_c
        se = np.sqrt(max(np.diag(np.linalg.pinv(best.jac.T @ best.jac))[i] * (sse / scale ** 2) / df, 1e-30)) if jac_ok else 1.
        ends = []
        for direction in (-1, 1):
            prev, root = z_hat[i], None
            for step in range(1, 81):
                cur = z_hat[i] + direction * step * .25 * se
                if prof(cur) > threshold:
                    try:
                        root = float(brentq(lambda v: prof(v) - threshold, min(prev, cur), max(prev, cur), xtol=1e-12))
                    except ValueError:
                        root = None
                    break
                prev = cur
            ends.append(root)
        log_scale = n in MODELS[name][1]
        tf = (lambda v: float(np.exp(v))) if log_scale else float
        ci = [None if e is None else tf(e) for e in ends]
        res['parameters'][n].update(ci=ci, ci_status='two_sided_profile_f' if None not in ci else 'open')
        if None in ci:
            d.append(f'profile_interval_open:{n}')
    # model-specific gates
    xmax = float(x.max()); withheld = set()
    if name in ('one_phase_decay', 'one_phase_association') and xmax * p['k'] < np.log(8):  # fewer than 3 half-lives observed
        d.append('plateau_not_approached'); withheld.add('plateau')
    if name == 'two_phase_decay':
        one, _, _ = fit_curve(g, cfg, _profile=False, model='one_phase_decay')
        if one.get('sse') is not None and sse > 0:
            fst = max(0., (one['sse'] - sse) / 2) / (sse / df); pv = float(stats.f.sf(fst, 2, df))
            res['structure_test'] = {'test': 'extra_sum_of_squares_F_two_vs_one_phase', 'f_statistic': float(fst), 'df': [2, df], 'p_value': pv, 'alpha': STRUCTURE_ALPHA}
            if pv >= STRUCTURE_ALPHA:
                d.append('two_phases_not_supported_by_data'); withheld.update(names)
        if xmax * p['k_slow'] < np.log(8):
            d.append('plateau_not_approached'); withheld.add('plateau')
    if name == 'michaelis_menten' and xmax < p['km']:
        d.append('saturation_not_approached'); withheld.update(['vmax', 'km'])
    if name == 'logistic_growth' and xmax < p['xmid'] + 2 * p['scal']:
        d.append('asymptote_not_approached'); withheld.add('asym')
    if name == 'exponential_growth' and p['k'] <= 0:
        withheld.add('k')
    for n in names:
        pr = res['parameters'][n]
        pr['reportable'] = bool(jac_ok and pr['ci_status'] == 'two_sided_profile_f' and n not in withheld)
    rate_names = {'one_phase_decay': ('half_life', 'k'), 'one_phase_association': ('half_time', 'k'), 'exponential_growth': ('doubling_time', 'k')}
    if name in rate_names:
        dn, kn = rate_names[name]; kp = res['parameters'][kn]
        res['derived'][dn] = {'estimate': LN2 / p[kn], 'ci': [None if kp['ci'][1] is None else LN2 / kp['ci'][1], None if kp['ci'][0] is None else LN2 / kp['ci'][0]],
                              'reportable': kp['reportable'], 'from': f'ln 2 / {kn}'}
    if name == 'two_phase_decay':
        for kn, dn in (('k_fast', 'half_life_fast'), ('k_slow', 'half_life_slow')):
            kp = res['parameters'][kn]
            res['derived'][dn] = {'estimate': LN2 / p[kn], 'ci': [None if kp['ci'][1] is None else LN2 / kp['ci'][1], None if kp['ci'][0] is None else LN2 / kp['ci'][0]], 'reportable': kp['reportable'], 'from': f'ln 2 / {kn}'}
        res['derived']['fraction_fast'] = {'estimate': p['span_fast'] / (p['span_fast'] + p['span_slow']), 'reportable': False, 'note': 'descriptive; no interval'}
    if name == 'logistic_growth':
        sp = res['parameters']['scal']
        res['derived']['rate'] = {'estimate': 1 / p['scal'], 'ci': [None if sp['ci'][1] is None else 1 / sp['ci'][1], None if sp['ci'][0] is None else 1 / sp['ci'][0]], 'reportable': sp['reportable'], 'from': '1 / scal'}
    res['reportable'] = bool(all(v['reportable'] for v in res['parameters'].values()))
    res['status'] = 'estimated' if res['reportable'] else 'limited'
    res['diagnostics'] = list(dict.fromkeys(d))
    grid_x = np.linspace(float(x.min()), xmax, 200)
    grid = [{'curve_id': res['curve_id'], 'x': float(a), 'predicted': float(b)} for a, b in zip(grid_x, predict(name, p, grid_x))]
    obs = [{'curve_id': res['curve_id'], 'x': float(a), 'y': float(b), 'predicted': float(c), 'residual': float(b - c)} for a, b, c in zip(x, y, fitted)]
    return res, grid, obs


def compute(d, cfg):
    fits, grid, obs = [], [], []
    for _, g in d.groupby('curve_id', sort=False):
        f, gr, ob = fit_curve(g, cfg); fits.append(f); grid += gr; obs += ob
    name = cfg['model']['name']; lv = cfg['uncertainty']['level']; names = MODELS[name][0]
    summaries = []
    key = {'one_phase_decay': 'half_life', 'one_phase_association': 'half_time', 'exponential_growth': 'doubling_time'}.get(name)
    for sample in dict.fromkeys(f['sample_id'] for f in fits):
        fs = [f for f in fits if f['sample_id'] == sample]
        for pn in ([key] if key else []) + [n for n in names if n not in MODELS[name][1]] + [n for n in names if n in MODELS[name][1]]:
            vals = [(f['derived'] if pn == key else f['parameters']).get(pn) for f in fs]
            row = {'sample_id': sample, 'quantity': pn, 'n_curves': len(fs), 'n_experiments': len({f['experiment_id'] for f in fs}), 'estimate': None, 'ci': [None, None]}
            if cfg['replicates']['independent_unit'] != 'experiment_id':
                row['status'] = 'independent_experiments_not_declared'
            elif any(v is None or not v['reportable'] for v in vals):
                row['status'] = 'withheld_incomplete_estimates'
            else:
                log_scale = pn == key or pn in MODELS[name][1]
                per = pd.Series([v['estimate'] for v in vals]).groupby([f['experiment_id'] for f in fs]).apply(lambda s: np.mean(np.log(s)) if log_scale else s.mean())
                m = float(per.mean()); row['estimate'] = float(np.exp(m)) if log_scale else m; row['scale'] = 'geometric mean' if log_scale else 'arithmetic mean'
                if len(per) > 1:
                    h = stats.t.ppf((1 + lv) / 2, len(per) - 1) * per.std(ddof=1) / np.sqrt(len(per))
                    row['ci'] = [float(np.exp(m - h)), float(np.exp(m + h))] if log_scale else [m - h, m + h]; row['status'] = 'summarized_t'
                else:
                    row['status'] = 'one_experiment_no_interval'
            summaries.append(row)
    must = [f"The {name.replace('_', ' ')} model was declared from prior knowledge; intervals are profile-F conditional on that model and on independent Gaussian errors"
            + (' with 1/y^2 weights' if cfg['model']['weighting'] == '1/y^2' else '') + '.']
    if any('plateau_not_approached' in f['diagnostics'] for f in fits):
        must.append('Some curves were observed for fewer than three half-lives; their plateaus are extrapolated and withheld.')
    if any('saturation_not_approached' in f['diagnostics'] for f in fits):
        must.append('Substrate concentrations stayed below Km in some curves; Vmax and Km are extrapolated and withheld (extend the concentration range).')
    if any('two_phases_not_supported_by_data' in f['diagnostics'] for f in fits):
        must.append('Two phases did not fit better than one (F test, p >= 0.01) in some curves; their two-phase parameters are withheld.')
    if any('systematic_residual_pattern' in f['diagnostics'] for f in fits):
        must.append('Residuals show a systematic pattern (runs test p < 0.05) in some curves; the model may not describe these data.')
    return {'analysis_type': 'nonlinear_fit', 'primary': {'model': cfg['model'], 'fits': fits, 'summaries': summaries, 'design': cfg['design'], 'level': lv},
            'curve_grid': grid, 'observations': obs, 'must_mention': must,
            'failing_items': [{'curve_id': f['curve_id'], 'diagnostics': f['diagnostics']} for f in fits if not f['reportable']],
            'limitations': ['One declared model per run; model choice is not data-driven. Compare models only for a predeclared question.',
                            'Profile intervals assume the model and error structure are right; replicate-experiment summaries carry between-experiment variability.']}


def render(run, result, style):
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); names = []
    obs, grid, p = pd.DataFrame(result['observations']), pd.DataFrame(result['curve_grid']), result['primary']
    if obs.empty:
        return []
    for n, sample in enumerate(dict.fromkeys(f['sample_id'] for f in p['fits'])):
        fs = [f for f in p['fits'] if f['sample_id'] == sample]
        with plt.rc_context(pstyle.rc(token, *font_setup())):
            fig, (ax, res) = plt.subplots(2, 1, figsize=(5.8, 4.6), sharex=True, height_ratios=[2.4, 1], layout='constrained')
            for j, f in enumerate(fs):
                o, gg = obs[obs.curve_id == f['curve_id']], grid[grid.curve_id == f['curve_id']]; col = pstyle.color(token, j)
                ax.plot(o.x, o.y, 'o', ms=4, color=col); ax.plot(gg.x, gg.predicted, '-', color=col, lw=1.4, label=f"{f['curve_id']} ({'reportable' if f['reportable'] else 'limited'})")
                res.plot(o.x, o.residual, 'o', ms=3, color=col, mfc='white')
            res.axhline(0, color=pstyle.MUTED, lw=.8)
            dz = p['design']; ax.set(ylabel=f"{dz['y_label']} ({dz['y_unit']})", title=f"{sample}: {p['model']['name'].replace('_', ' ')}")
            res.set(xlabel=f"{dz['x_label']} ({dz['x_unit']})", ylabel='Residual'); ax.legend(frameon=False, fontsize=7)
            name = f'nonlinear_{n+1}'; pstyle.save(fig, out/name, token, 300); names.append(name)
    return names
