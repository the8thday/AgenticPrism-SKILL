"""Tolerance intervals and process capability for one quality attribute.

A tolerance interval describes where a stated proportion (content) of the
individual units lie, with a stated confidence. It summarizes observed
manufacturing variability; it is not by itself a specification, which also
needs clinical relevance and regulatory input.

Normal intervals (independent units assumed Gaussian):
- one-sided exact: k = t'(conf; n-1, z_p sqrt(n)) / sqrt(n) (noncentral t);
- two-sided exact: the k solving sqrt(2n/pi) Int_0^inf P(chi2_nu > nu R(x)^2 / k^2)
  exp(-n x^2/2) dx = conf, with Phi(x+R) - Phi(x-R) = p (Odeh 1978; Eberhardt,
  Mee and Reeve 1989), and the Howe (1969) approximation in the NIST/SEMATECH
  form (without the second-order correction) as an option.
Nonparametric intervals: order statistics X(r) and X(n-r+1) (two-sided) or X(r)
(one-sided), with the largest r whose exact binomial/beta confidence reaches the
declared level; otherwise the minimum sample size is reported.
Capability: Pp/Ppk from the overall SD; Cp/Cpk from the pooled within-subgroup
SD when subgroups are declared. Cp/Pp intervals are chi-square exact; Cpk/Ppk
use the normal approximation C +/- z sqrt(1/(9n) + C^2/(2(n-1))) quoted by the
NIST/SEMATECH handbook, which it recommends for n >= 25.
"""
import numpy as np
from scipy import stats, integrate, optimize
from .cmc_common import merge, text, number, common, table

DEFAULTS = {'analysis_type': 'specification', 'schema_version': 1, 'input': None, 'source': None,
    'attribute': {'name': None, 'unit': None, 'rationale': None},
    'design': {'independent_units': None, 'unit_description': None},
    'tolerance': {'method': None, 'sides': 'two_sided', 'content': None, 'confidence': None, 'normal_factor': 'exact'},
    'capability': {'lsl': None, 'usl': None, 'limits_source': None, 'confidence': .95, 'subgroup_column': None},
    'report': {'plot_style': 'prism_like'}}


def resolve_config(raw):
    c = merge(raw, DEFAULTS); common(c, 'specification')
    a, d, t, p = c['attribute'], c['design'], c['tolerance'], c['capability']
    for k in ('name', 'unit', 'rationale'):
        text(a[k], 'attribute.' + k)
    if d['independent_units'] is not True:
        raise ValueError('Declare design.independent_units as literal true (for example one value per lot)')
    text(d['unit_description'], 'design.unit_description')
    if t['method'] is None and p['lsl'] is None and p['usl'] is None:
        raise ValueError('Request a tolerance interval, capability limits, or both')
    if t['method'] is not None:
        if t['method'] not in ('normal', 'nonparametric'):
            raise ValueError('tolerance.method must be normal or nonparametric')
        if t['sides'] not in ('two_sided', 'lower', 'upper'):
            raise ValueError('tolerance.sides must be two_sided, lower or upper')
        for k in ('content', 'confidence'):
            number(t[k], 'tolerance.' + k, positive=True)
            if not t[k] < 1:
                raise ValueError(f'tolerance.{k} must lie in (0, 1)')
        if t['normal_factor'] not in ('exact', 'howe'):
            raise ValueError('tolerance.normal_factor must be exact or howe')
    elif t['content'] is not None or t['confidence'] is not None:
        raise ValueError('tolerance content/confidence need tolerance.method')
    if p['lsl'] is not None or p['usl'] is not None:
        for k in ('lsl', 'usl'):
            if p[k] is not None:
                number(p[k], 'capability.' + k)
        if p['lsl'] is not None and p['usl'] is not None and p['lsl'] >= p['usl']:
            raise ValueError('lsl must be below usl')
        text(p['limits_source'], 'capability.limits_source')
        number(p['confidence'], 'capability.confidence', positive=True)
        if not .5 < p['confidence'] < 1:
            raise ValueError('capability.confidence must lie in (.5, 1)')
    elif p['limits_source'] is not None or p['subgroup_column'] is not None:
        raise ValueError('capability settings need lsl and/or usl')
    return c


def load_data(path, cfg):
    need = ['observation_id', 'value', 'exclude'] + ([cfg['capability']['subgroup_column']] if cfg['capability']['subgroup_column'] else [])
    d = table(path, need, ['value'])
    if set(d.exclude) - {'true', 'false'}:
        raise ValueError('exclude must be true/false')
    d['exclude'] = d.exclude == 'true'
    if 'exclusion_reason' not in d or (d.exclude & (d.exclusion_reason.str.strip() == '')).any():
        raise ValueError('Excluded rows need an exclusion_reason')
    return d


def k_one_sided(n, p, conf):
    return float(stats.nct.ppf(conf, n-1, stats.norm.ppf(p)*np.sqrt(n))/np.sqrt(n))


def k_howe(n, p, conf):
    # NIST/SEMATECH 7.2.6.3 form, without Howe's second-order correction factor
    # (R tolerance::K.factor method "HE" includes that factor and differs slightly).
    nu = n-1
    return float(stats.norm.ppf((1+p)/2)*np.sqrt(nu*(1+1/n)/stats.chi2.ppf(1-conf, nu)))


def _r(x, p):
    return optimize.brentq(lambda r: stats.norm.cdf(x+r)-stats.norm.cdf(x-r)-p, 1e-12, x+40., xtol=1e-14)


def k_two_sided_exact(n, p, conf):
    nu = n-1

    def confidence(k):
        f = lambda x: stats.chi2.sf(nu*_r(x, p)**2/k**2, nu)*np.exp(-n*x*x/2)
        val, _ = integrate.quad(f, 0, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200)
        return np.sqrt(2*n/np.pi)*val
    h = k_howe(n, p, conf)
    return float(optimize.brentq(lambda k: confidence(k)-conf, h*.8, h*1.25, xtol=1e-12))


def order_statistic(n, p, conf, sides):
    """Largest symmetric rank r with exact confidence >= conf (continuous data); None if impossible."""
    best = None
    for r in range(1, n//2+1 if sides == 'two_sided' else n+1):
        if sides == 'two_sided':
            achieved = float(stats.beta.sf(p, n-2*r+1, 2*r))
        else:
            achieved = float(stats.beta.sf(p, n-r+1, r))
        if achieved >= conf:
            best = (r, achieved)
        else:
            break
    return best


def minimum_n(p, conf, sides):
    n = 2
    while n < 100000:
        if order_statistic(n, p, conf, sides):
            return n
        n += 1
    return None


def capability(x, cfg, subgroups=None):
    c = cfg['capability']; n = len(x); m = float(x.mean()); s = float(x.std(ddof=1)); level = c['confidence']
    z = stats.norm.ppf((1+level)/2); out = {'n': n, 'mean': m, 'overall_sd': s, 'lsl': c['lsl'], 'usl': c['usl'],
                                            'limits_source': c['limits_source'], 'confidence_level': level}

    def indices(sd, df, prefix):
        res = {}
        if c['lsl'] is not None and c['usl'] is not None:
            cp = (c['usl']-c['lsl'])/(6*sd)
            res[prefix+'p'] = cp
            res[prefix+'p_interval'] = [cp*np.sqrt(stats.chi2.ppf((1-level)/2, df)/df), cp*np.sqrt(stats.chi2.ppf((1+level)/2, df)/df)]
        sides = [v for v in ((c['usl']-m)/(3*sd) if c['usl'] is not None else None, (m-c['lsl'])/(3*sd) if c['lsl'] is not None else None) if v is not None]
        cpk = min(sides)
        half = z*np.sqrt(1/(9*(df+1))+cpk**2/(2*df))
        res[prefix+'pk'] = cpk; res[prefix+'pk_interval'] = [cpk-half, cpk+half]
        return res
    out.update({k.replace('Pp', 'Pp'): v for k, v in indices(s, n-1, 'P').items()})
    if subgroups is not None:
        groups = [x[subgroups == g] for g in np.unique(subgroups)]
        dfw = sum(len(g)-1 for g in groups)
        if dfw < 1:
            raise ValueError('Subgroups need replicate measurements for a within-subgroup SD')
        sw = float(np.sqrt(sum(((g-g.mean())**2).sum() for g in groups)/dfw))
        out['within_sd'] = sw; out['within_df'] = dfw
        out.update(indices(sw, dfw, 'C'))
    lo = stats.norm.cdf((c['lsl']-m)/s) if c['lsl'] is not None else 0.
    hi = stats.norm.sf((c['usl']-m)/s) if c['usl'] is not None else 0.
    out['expected_nonconforming_fraction_normal_plugin'] = float(lo+hi)
    out['interval_methods'] = {'Cp/Pp': 'chi-square exact (NIST/SEMATECH 6.1.6)', 'Cpk/Ppk': 'normal approximation C +/- z sqrt(1/(9n) + C^2/(2(n-1)))'}
    return out


def compute(d, cfg):
    used = d[~d.exclude]; x = used.value.to_numpy(float); n = len(x)
    if n < 3:
        raise ValueError('At least three independent values are required')
    must, failing = [], []
    sw = stats.shapiro(x) if n <= 5000 else None
    diag = {'n': n, 'mean': float(x.mean()), 'sd': float(x.std(ddof=1)), 'min': float(x.min()), 'max': float(x.max()),
            'shapiro_w': float(sw.statistic) if sw else None, 'shapiro_p': float(sw.pvalue) if sw else None,
            'skewness': float(stats.skew(x, bias=False)) if n > 2 else None}
    t = cfg['tolerance']; ti = None; primary = {}
    if t['method'] == 'normal':
        if t['sides'] == 'two_sided':
            k = k_two_sided_exact(n, t['content'], t['confidence']) if t['normal_factor'] == 'exact' else k_howe(n, t['content'], t['confidence'])
            limits = [diag['mean']-k*diag['sd'], diag['mean']+k*diag['sd']]
        else:
            k = k_one_sided(n, t['content'], t['confidence'])
            limits = [diag['mean']-k*diag['sd'], None] if t['sides'] == 'lower' else [None, diag['mean']+k*diag['sd']]
        ti = {'method': 'normal', 'factor_method': t['normal_factor'] if t['sides'] == 'two_sided' else 'exact_noncentral_t',
              'k': k, 'limits': limits, 'content': t['content'], 'confidence': t['confidence'], 'sides': t['sides']}
        if diag['shapiro_p'] is not None and diag['shapiro_p'] < .05:
            must.append(f'Shapiro-Wilk p = {diag["shapiro_p"]:.3g}: the normal tolerance interval relies on a Gaussian model the data may not support. Do not switch methods after seeing results without a protocol reason.')
    elif t['method'] == 'nonparametric':
        best = order_statistic(n, t['content'], t['confidence'], t['sides']); xs = np.sort(x)
        if best is None:
            need = minimum_n(t['content'], t['confidence'], t['sides'])
            ti = {'method': 'nonparametric', 'limits': None, 'status': 'insufficient_sample', 'minimum_n': need,
                  'content': t['content'], 'confidence': t['confidence'], 'sides': t['sides']}
            failing.append({'reason': 'nonparametric_interval_not_attainable', 'n': n, 'minimum_n': need})
            must.append(f'{n} values cannot give a nonparametric interval with content {t["content"]} and confidence {t["confidence"]}; at least {need} are needed.')
        else:
            r, achieved = best
            limits = [float(xs[r-1]), float(xs[n-r])] if t['sides'] == 'two_sided' else [float(xs[r-1]), None] if t['sides'] == 'lower' else [None, float(xs[n-r])]
            ti = {'method': 'nonparametric', 'order_statistic_r': r, 'achieved_confidence': achieved, 'limits': limits,
                  'content': t['content'], 'confidence': t['confidence'], 'sides': t['sides']}
    if ti:
        primary['tolerance_interval'] = ti['limits']
        must.append('A tolerance interval describes observed variability of these units; it is not a specification by itself.')
    cap = None
    if cfg['capability']['lsl'] is not None or cfg['capability']['usl'] is not None:
        sg = used[cfg['capability']['subgroup_column']].to_numpy() if cfg['capability']['subgroup_column'] else None
        cap = capability(x, cfg, sg)
        primary['Ppk'] = cap['Ppk']; primary['Ppk_interval'] = cap['Ppk_interval']
        if 'Cpk' in cap:
            primary['Cpk'] = cap['Cpk']; primary['Cpk_interval'] = cap['Cpk_interval']
            must.append(f'Cpk ({cap["Cpk"]:.3g}) uses only the within-subgroup SD and excludes between-subgroup (for example between-lot) variation; Ppk ({cap["Ppk"]:.3g}) uses the overall SD. Report which one answers the question.')
            must.append(f'The {n} measurements come from {len(np.unique(sg))} subgroups and are not independent units; the Ppk interval treats them as independent and is too narrow. For a lot-level capability, use one value per lot.')
        if n < 25:
            must.append(f'n = {n}: the Cpk/Ppk interval is a normal approximation that the NIST/SEMATECH handbook recommends only for n >= 25.')
        if diag['shapiro_p'] is not None and diag['shapiro_p'] < .05:
            must.append('Capability indices and the nonconforming fraction assume a Gaussian distribution, which the Shapiro-Wilk test questions.')
    return {'analysis_type': 'specification', 'schema_version': 1, 'attribute': cfg['attribute'], 'data_summary': diag,
            'tolerance_interval': ti, 'capability': cap, 'primary': primary, 'must_mention': must, 'failing_items': failing,
            'limitations': ['Units are assumed independent (for example one value per lot); repeated measurements of a lot are not independent units.',
                            'Normal-theory results assume a Gaussian distribution; nonparametric intervals assume continuous data.',
                            'Capability indices describe the process relative to the declared limits; they do not validate those limits.']}
