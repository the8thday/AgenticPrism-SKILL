"""Prospective power and sample size for declared designs (0.13.1).

Every assumption (difference, SD, proportions, hazard ratio, event probability, margins)
must be declared with its source; nothing is estimated from the study being planned.

Exact methods: two-sample and paired t (noncentral t, both rejection tails, as
pwr.t.test and power.t.test(strict = TRUE)); one-way ANOVA (noncentral F with
lambda = k n f^2, f from the declared group means, as pwr.anova.test); TOST for two
independent means (exact integral over the sample-SD distribution, the quantity
PowerTOST computes with Owen's Q). Approximations, labelled as such: two proportions
by Cohen's arcsine h (pwr.2p2n.test) or the pooled normal approximation
(power.prop.test), and the Schoenfeld log-rank events formula.
"""
from copy import deepcopy
import numpy as np
from scipy import integrate, stats

TYPES = {'sample_size'}
DESIGNS = ('two_sample_t', 'paired_t', 'one_way_anova', 'two_proportions', 'logrank', 'tost_two_means')
DEFAULTS = {'schema_version': 1, 'analysis_type': 'sample_size', 'input': None, 'source': None, 'purpose': None,
            'design': None, 'solve_for': None, 'alpha': .05, 'sidedness': 'two_sided', 'target_power': None, 'n': None,
            'assumptions': {}, 'sensitivity': None, 'report': {'plot_style': 'prism_like'}}
# Required assumption keys per design; each needs '<key>_source' too unless listed in NO_SOURCE.
NEEDED = {'two_sample_t': ('difference', 'sd'), 'paired_t': ('difference', 'sd_difference'),
          'one_way_anova': ('group_means', 'sd'), 'two_proportions': ('p1', 'p2'),
          'logrank': ('hazard_ratio', 'event_probability'), 'tost_two_means': ('difference', 'sd', 'margin')}
OPTIONAL = {'allocation_ratio': 1., 'method': None}
N_MAX = 100000


def _number(v, name, low=None, high=None, open_low=True):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v):
        raise ValueError(f'{name} must be a finite number')
    if low is not None and (v <= low if open_low else v < low):
        raise ValueError(f'{name} must exceed {low}')
    if high is not None and v >= high:
        raise ValueError(f'{name} must be below {high}')


def resolve_config(raw):
    if not isinstance(raw, dict) or set(raw) - set(DEFAULTS):
        raise ValueError('Unsupported sample-size configuration keys')
    c = deepcopy(DEFAULTS)
    c.update({k: deepcopy(v) for k, v in raw.items()})
    if c['analysis_type'] != 'sample_size' or c['schema_version'] != 1:
        raise ValueError('Unsupported sample-size schema')
    if c['input'] is not None:
        raise ValueError('Sample-size planning takes declared assumptions, not an input file')
    for k in ('source', 'purpose'):
        if not isinstance(c[k], str) or not c[k].strip():
            raise ValueError(f'Declare {k}: the study being planned and where the assumptions come from')
    d = c['design']
    if d not in DESIGNS:
        raise ValueError('design must be one of ' + ', '.join(DESIGNS))
    if c['solve_for'] not in ('n', 'power'):
        raise ValueError('solve_for must be n or power')
    _number(c['alpha'], 'alpha', 0, .5)
    if c['sidedness'] not in ('two_sided', 'one_sided'):
        raise ValueError('sidedness must be two_sided or one_sided')
    if d in ('one_way_anova', 'tost_two_means') and c['sidedness'] != 'two_sided':
        raise ValueError(f'{d} has no sidedness choice; leave two_sided (TOST uses alpha for each one-sided test)')
    if c['solve_for'] == 'n':
        _number(c['target_power'], 'target_power', 0, 1)
        if c['target_power'] <= c['alpha'] or c['n'] is not None:
            raise ValueError('Solving for n needs target_power above alpha and no n')
    else:
        n = c['n']
        if c['target_power'] is not None or (d == 'logrank' and not (type(n) is int and n >= 2)) or (d != 'logrank' and not (type(n) is int and n >= 2)):
            raise ValueError('Solving for power needs an integer n >= 2 (per group; total subjects for logrank) and no target_power')
    a = c['assumptions']
    if not isinstance(a, dict):
        raise ValueError('assumptions must be an object')
    allowed = set(NEEDED[d]) | {k + '_source' for k in NEEDED[d]} | (
        {'allocation_ratio'} if d in ('two_sample_t', 'two_proportions', 'logrank', 'tost_two_means') else set()) | (
        {'method'} if d == 'two_proportions' else set())
    if set(a) - allowed:
        raise ValueError(f'Assumptions not used by {d}: {sorted(set(a) - allowed)}')
    for k in NEEDED[d]:
        if k not in a:
            raise ValueError(f'{d} requires assumptions.{k}')
        if not isinstance(a.get(k + '_source'), str) or not a[k + '_source'].strip():
            raise ValueError(f'Declare assumptions.{k}_source (pilot data, literature, historical runs or a protocol); never the planned study itself')
    a.setdefault('allocation_ratio', 1.) if 'allocation_ratio' in allowed else None
    if 'allocation_ratio' in a:
        _number(a['allocation_ratio'], 'allocation_ratio', 0, 100)
    if d in ('two_sample_t', 'paired_t'):
        _number(a['difference'], 'difference')
        if a['difference'] == 0:
            raise ValueError('A zero difference has no power to detect; declare the smallest relevant difference')
        _number(a['sd' if d == 'two_sample_t' else 'sd_difference'], 'sd', 0)
    if d == 'one_way_anova':
        m = a['group_means']
        if not isinstance(m, list) or len(m) < 3 or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) for v in m):
            raise ValueError('group_means must list at least three finite means (use two_sample_t for two groups)')
        if np.ptp(m) == 0:
            raise ValueError('Declared group means are all equal')
        _number(a['sd'], 'sd', 0)
    if d == 'two_proportions':
        for k in ('p1', 'p2'):
            _number(a[k], k, 0, 1)
        if a['p1'] == a['p2']:
            raise ValueError('p1 equals p2')
        if a.get('method') not in ('arcsine', 'normal_approximation'):
            raise ValueError('Declare assumptions.method: arcsine (Cohen h, pwr) or normal_approximation (pooled, power.prop.test)')
        if a['method'] == 'normal_approximation' and a['allocation_ratio'] != 1:
            raise ValueError('normal_approximation supports equal groups only (as power.prop.test)')
    if d == 'logrank':
        _number(a['hazard_ratio'], 'hazard_ratio', 0)
        if a['hazard_ratio'] == 1:
            raise ValueError('hazard_ratio 1 has no power')
        _number(a['event_probability'], 'event_probability', 0, 1 + 1e-12)
    if d == 'tost_two_means':
        _number(a['difference'], 'difference')
        _number(a['sd'], 'sd', 0)
        mg = a['margin']
        if not isinstance(mg, list) or len(mg) != 2 or any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in mg) or not mg[0] < 0 < mg[1]:
            raise ValueError('margin must be [lower, upper] with lower < 0 < upper, declared before the study')
        if not mg[0] < a['difference'] < mg[1]:
            raise ValueError('The assumed true difference lies outside the margin; equivalence cannot be powered')
    s = c['sensitivity']
    if s is not None:
        if not isinstance(s, dict) or set(s) != {'assumption', 'values'} or s['assumption'] not in NEEDED[d] or s['assumption'] in ('group_means', 'margin'):
            raise ValueError('sensitivity needs {assumption, values} naming one scalar declared assumption')
        if not isinstance(s['values'], list) or not 1 <= len(s['values']) <= 25:
            raise ValueError('sensitivity.values must list 1-25 alternative values')
    if c['report']['plot_style'] not in ('prism_like', 'standard'):
        raise ValueError('Invalid style')
    return c


def load_data(path, cfg):
    return None


# ---------- power functions (n is per group unless stated) ----------

def _t_power(ncp, df, alpha, sides):
    if sides == 'two_sided':
        q = stats.t.ppf(1 - alpha / 2, df)
        # Lower tail by symmetry: SciPy's nct.cdf(-q) returns NaN for large ncp (0.13.1 benchmark).
        power = stats.nct.sf(q, df, ncp) + stats.nct.sf(q, df, -ncp)
    else:
        power = stats.nct.sf(stats.t.ppf(1 - alpha, df), df, abs(ncp))
    if not np.isfinite(power):
        raise ArithmeticError(f'Noncentral t power not finite (df={df}, ncp={ncp})')
    return float(power)


def power_two_sample_t(n1, a, alpha, sides):
    n2 = int(np.ceil(a['allocation_ratio'] * n1))
    ncp = a['difference'] / (a['sd'] * np.sqrt(1 / n1 + 1 / n2))
    return _t_power(ncp, n1 + n2 - 2, alpha, sides), n2


def power_paired_t(n, a, alpha, sides):
    return _t_power(a['difference'] / (a['sd_difference'] / np.sqrt(n)), n - 1, alpha, sides), n


def cohen_f(a):
    m = np.asarray(a['group_means'], float)
    return float(np.sqrt(np.mean((m - m.mean()) ** 2)) / a['sd'])


def power_anova(n, a, alpha, sides):
    k = len(a['group_means']); df1, df2 = k - 1, k * (n - 1)
    q = stats.f.ppf(1 - alpha, df1, df2)
    return float(stats.ncf.sf(q, df1, df2, k * n * cohen_f(a) ** 2)), n


def power_two_proportions(n1, a, alpha, sides):
    n2 = int(np.ceil(a['allocation_ratio'] * n1)); p1, p2 = a['p1'], a['p2']
    if a['method'] == 'arcsine':
        h = abs(2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p2))); scale = h * np.sqrt(n1 * n2 / (n1 + n2))
        if sides == 'two_sided':
            z = stats.norm.ppf(1 - alpha / 2)
            return float(stats.norm.sf(z - scale) + stats.norm.cdf(-z - scale)), n2
        return float(stats.norm.sf(stats.norm.ppf(1 - alpha) - scale)), n2
    # power.prop.test (strict = TRUE): pooled variance under H0, unpooled under H1
    pbar, q = (p1 + p2) / 2, abs(p1 - p2)
    z = stats.norm.ppf(1 - alpha / 2) if sides == 'two_sided' else stats.norm.ppf(1 - alpha)
    se1 = np.sqrt(p1 * (1 - p1) + p2 * (1 - p2)); se0 = np.sqrt(2 * pbar * (1 - pbar))
    upper = stats.norm.cdf((np.sqrt(n1) * q - z * se0) / se1)
    if sides == 'two_sided':
        upper += stats.norm.cdf((-np.sqrt(n1) * q - z * se0) / se1)
    return float(upper), n2


def schoenfeld_events(a, alpha, power, sides):
    r = a['allocation_ratio']; p = 1 / (1 + r)  # proportion in group 1
    z = stats.norm.ppf(1 - alpha / 2) if sides == 'two_sided' else stats.norm.ppf(1 - alpha)
    return (z + stats.norm.ppf(power)) ** 2 / (p * (1 - p) * np.log(a['hazard_ratio']) ** 2)


def power_logrank(n_total, a, alpha, sides):
    r = a['allocation_ratio']; p = 1 / (1 + r); events = n_total * a['event_probability']
    z = stats.norm.ppf(1 - alpha / 2) if sides == 'two_sided' else stats.norm.ppf(1 - alpha)
    scale = np.sqrt(events * p * (1 - p)) * abs(np.log(a['hazard_ratio']))
    pw = stats.norm.sf(z - scale) + (stats.norm.cdf(-z - scale) if sides == 'two_sided' else 0.)
    return float(pw), n_total


def power_tost(n1, a, alpha, sides):
    """Exact TOST power for two independent means with pooled SD (df = n1 + n2 - 2)."""
    n2 = int(np.ceil(a['allocation_ratio'] * n1)); df = n1 + n2 - 2
    k = np.sqrt(1 / n1 + 1 / n2); sigma = a['sd'] * k; low, high = a['margin']; delta = a['difference']
    tc = stats.t.ppf(1 - alpha, df)

    def given_s(u):  # u = chi-square(df) draw; sample SD s = sd * sqrt(u / df)
        half = tc * a['sd'] * np.sqrt(u / df) * k
        lo, hi = low + half, high - half
        return max(0., stats.norm.cdf((hi - delta) / sigma) - stats.norm.cdf((lo - delta) / sigma)) * stats.chi2.pdf(u, df)
    # the integrand vanishes once the acceptance interval closes: lo >= hi
    u_max = df * ((high - low) / (2 * tc * a['sd'] * k)) ** 2
    value, _ = integrate.quad(given_s, 0, u_max, limit=400, epsabs=1e-12, epsrel=1e-10)
    return float(min(1., value)), n2


POWER = {'two_sample_t': power_two_sample_t, 'paired_t': power_paired_t, 'one_way_anova': power_anova,
         'two_proportions': power_two_proportions, 'logrank': power_logrank, 'tost_two_means': power_tost}
MINIMUM_N = {'two_sample_t': 2, 'paired_t': 2, 'one_way_anova': 2, 'two_proportions': 2, 'logrank': 2, 'tost_two_means': 2}


def solve_n(design, a, alpha, sides, target):
    """Smallest integer n (per group; total subjects for logrank) with power >= target."""
    f = lambda n: POWER[design](n, a, alpha, sides)[0]
    lo = MINIMUM_N[design]
    if f(lo) >= target:
        return lo
    hi = lo
    while f(hi) < target:
        hi *= 2
        if hi > N_MAX:
            return None
    while hi - lo > 1:  # power is monotone in n for these designs
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if f(mid) < target else (lo, mid)
    return hi


def _effect(design, a):
    if design == 'two_sample_t':
        return {'cohen_d': a['difference'] / a['sd']}
    if design == 'paired_t':
        return {'cohen_dz': a['difference'] / a['sd_difference']}
    if design == 'one_way_anova':
        return {'cohen_f': cohen_f(a)}
    if design == 'two_proportions':
        return {'cohen_h': float(2 * np.arcsin(np.sqrt(a['p1'])) - 2 * np.arcsin(np.sqrt(a['p2'])))}
    if design == 'logrank':
        return {'log_hazard_ratio': float(np.log(a['hazard_ratio']))}
    return {'standardized_margin_lower': a['margin'][0] / a['sd'], 'standardized_margin_upper': a['margin'][1] / a['sd'],
            'standardized_difference': a['difference'] / a['sd']}


def evaluate(cfg, a):
    d, alpha, sides = cfg['design'], cfg['alpha'], cfg['sidedness']
    if cfg['solve_for'] == 'n':
        n = solve_n(d, a, alpha, sides, cfg['target_power'])
        if n is None:
            return {'status': 'not_attainable', 'n_group_1': None, 'reason': f'power target not reached below n = {N_MAX}'}
    else:
        n = cfg['n']
    power, n2 = POWER[d](n, a, alpha, sides)
    row = {'status': 'computed', 'power': power, **_effect(d, a)}
    if d == 'logrank':
        row.update(total_subjects=n, expected_events=n * a['event_probability'],
                   schoenfeld_events_required=float(schoenfeld_events(a, alpha, cfg['target_power'], sides)) if cfg['solve_for'] == 'n' else None,
                   group_1=int(round(n / (1 + a['allocation_ratio']))), group_2=n - int(round(n / (1 + a['allocation_ratio']))))
    elif d == 'paired_t':
        row.update(pairs=n, total_subjects=n)
    elif d == 'one_way_anova':
        row.update(n_per_group=n, groups=len(a['group_means']), total_subjects=n * len(a['group_means']))
    else:
        row.update(n_group_1=n, n_group_2=n2, total_subjects=n + n2)
    return row


METHOD = {'two_sample_t': 'Exact noncentral t, pooled SD, both rejection tails (pwr.t.test / power.t.test strict).',
          'paired_t': 'Exact noncentral t on paired differences, both rejection tails.',
          'one_way_anova': 'Exact noncentral F, lambda = k n f^2 with Cohen f from the declared means (pwr.anova.test).',
          'logrank': 'Schoenfeld approximation: required events from the log hazard ratio; subjects = events / declared event probability.',
          'tost_two_means': 'Exact TOST power (two one-sided t tests at alpha each, pooled SD), integrated over the sample-SD distribution.'}


def compute(_, cfg):
    d, a = cfg['design'], cfg['assumptions']
    primary = {'design': d, 'solve_for': cfg['solve_for'], 'alpha': cfg['alpha'], 'sidedness': cfg['sidedness'],
               'target_power': cfg['target_power'], 'purpose': cfg['purpose'],
               'method': METHOD.get(d) or ('Cohen arcsine h with normal power (pwr.2p2n.test); an approximation to the chi-square test.'
                                           if a['method'] == 'arcsine' else 'Pooled normal approximation (power.prop.test, strict); an approximation to the chi-square test.'),
               'assumptions': a, **evaluate(cfg, a)}
    curve = []
    if primary['status'] == 'computed':
        n0 = primary.get('n_group_1') or primary.get('n_per_group') or primary.get('pairs') or primary.get('total_subjects')
        grid = sorted({int(v) for v in np.unique(np.round(np.geomspace(max(MINIMUM_N[d], n0 / 4), max(n0 * 2.5, MINIMUM_N[d] + 3), 40)))})
        curve = [{'n': n, 'power': POWER[d](n, a, cfg['alpha'], cfg['sidedness'])[0]} for n in grid]
    sensitivity = []
    if cfg['sensitivity']:
        key = cfg['sensitivity']['assumption']
        for v in cfg['sensitivity']['values']:
            alt = {**a, key: v}
            try:
                probe = deepcopy(cfg); probe['assumptions'] = alt; resolve_config({k: v2 for k, v2 in probe.items()})
                sensitivity.append({'assumption': key, 'value': v, **evaluate(cfg, alt)})
            except ValueError as e:
                sensitivity.append({'assumption': key, 'value': v, 'status': 'invalid_assumption', 'reason': str(e)})
    must = ['The result is only as good as the declared assumptions and their sources; it is a planning calculation, not evidence about the effect.',
            'Report the assumed values and sources together with the sample size; do not quote the sample size alone.']
    if not cfg['sensitivity']:
        must.append('No sensitivity analysis was declared; the required n can change sharply with the SD or effect assumption.')
    if d in ('two_proportions', 'logrank'):
        must.append('This design uses an approximation; check the operating characteristics by simulation when n is small or proportions are extreme.')
    if d == 'logrank':
        must.append('The event probability (accrual, follow-up, dropout) drives the subject count; the number of events is what provides power.')
    if d == 'tost_two_means':
        must.append('Equivalence margins must be justified before the study (for CMC, from reference-lot variability or a protocol), never chosen to make n feasible.')
    if primary['status'] == 'not_attainable':
        must.append(primary['reason'])
    return {'analysis_type': 'sample_size', 'primary': primary, 'sensitivity': sensitivity,
            'power_curve': {'n_definition': 'subjects in total' if d == 'logrank' else 'pairs' if d == 'paired_t' else 'per group; group 1 if unequal', 'points': curve},
            'must_mention': must,
            'failing_items': [primary] if primary['status'] != 'computed' else [],
            'limitations': ['Independent observations and the stated test are assumed; clustering, repeated measures, dropout and multiplicity need separate adjustment.',
                            'Power refers to the test named in the method line, analysed as planned.']}


def render(run, result, style):
    """Power curve from saved values only."""
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .report import font_setup
    from . import plot_style as pstyle
    pts = result['power_curve']['points']
    if not pts:
        return []
    token = pstyle.THEMES[style]; out = Path(run)/'figures'; out.mkdir(exist_ok=True); p = result['primary']
    with plt.rc_context(pstyle.rc(token, *font_setup())):
        fig, ax = plt.subplots(figsize=(5.6, 3.6), layout='constrained')
        ax.plot([q['n'] for q in pts], [q['power'] for q in pts], '-', color=pstyle.color(token, 0), lw=1.6)
        if p.get('target_power'):
            ax.axhline(p['target_power'], color=pstyle.MUTED, lw=.8, ls=(0, (3, 2)))
        n = p.get('n_group_1') or p.get('n_per_group') or p.get('pairs') or p.get('total_subjects')
        ax.plot([n], [p['power']], 'o', color=pstyle.color(token, 0), ms=6)
        ax.set(xlabel=f"n ({result['power_curve']['n_definition']})", ylabel='Power', ylim=(0, 1.02), title=p['design'].replace('_', ' '))
        pstyle.save(fig, out/'power_curve', token, 300)
    return ['power_curve']
