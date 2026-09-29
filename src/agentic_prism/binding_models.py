"""Mass-action binding primitives; all concentrations use the same molar unit."""
import numpy as np
from scipy.optimize import brentq


def complex_concentration(pt, lt, kd):
    """Cancellation-free P + L <-> PL solution, with scale-safe discriminant."""
    p, l, k = np.broadcast_arrays(np.asarray(pt, float), np.asarray(lt, float), np.asarray(kd, float))
    if not all(np.isfinite(a).all() for a in (p, l, k)) or np.any(p < 0) or np.any(l < 0) or np.any(k <= 0):
        raise ValueError('Finite nonnegative totals and positive KD required')
    scale = np.maximum.reduce([p, l, k])
    a, b, c = p/scale, l/scale, k/scale
    # Algebraically S²-4PL; this form also avoids cancellation at P=L, KD<<P.
    root = np.sqrt((a-b)**2 + c*(2*(a+b)+c))
    return scale * (2*a*b / (a+b+c+root))


def free_fraction(pt, lt, kd):
    """Free constant P fraction; avoids subtracting near-unit occupancy."""
    p, l, k = np.broadcast_arrays(np.asarray(pt, float), np.asarray(lt, float), np.asarray(kd, float))
    if np.any(p <= 0): raise ValueError('Constant active species must be positive')
    bound = complex_concentration(p, l, k)
    free_l = np.maximum(l-bound, 0.)
    return k/(k+free_l)


def competition_fraction(inhibitor_total, receptor_total, tracer_total, tracer_kd, ki):
    """Three mutually exclusive states R, RT, RI. Return fraction tracer bound.

    Solve the strictly increasing receptor mass balance for free R, scaled to
    R_total. This avoids choosing an unphysical cubic root.
    """
    it = np.asarray(inhibitor_total, float)
    if not np.isfinite(it).all() or np.any(it < 0): raise ValueError('Invalid competitor total')
    if any(not np.isfinite(v) or v <= 0 for v in (receptor_total, tracer_total, tracer_kd, ki)):
        raise ValueError('Positive receptor, tracer, tracer KD and Ki required')
    p, t, k, i = receptor_total, tracer_total, tracer_kd, ki
    def one(total):
        u = brentq(lambda u: u + t*u/(k+p*u) + total*u/(i+p*u)-1., 0., 1., xtol=1e-15, rtol=1e-14)
        return p*u/(k+p*u)
    return np.array([one(v) for v in it.ravel()]).reshape(it.shape)


def cheng_prusoff(ic50, tracer_total, tracer_kd, receptor_total, max_depleted_fraction=.05):
    """Labelled approximation, refused if either ligand can be materially depleted."""
    if any(not np.isfinite(v) or v <= 0 for v in (ic50, tracer_total, tracer_kd, receptor_total)):
        raise ValueError('Positive concentrations required')
    bound = float(complex_concentration(receptor_total, tracer_total, tracer_kd))
    if bound/tracer_total > max_depleted_fraction or receptor_total/ic50 > max_depleted_fraction:
        raise ValueError('Cheng-Prusoff invalid under tracer or competitor depletion')
    return ic50/(1+tracer_total/tracer_kd)


def regime(pt, kd):
    ratio = float(pt/kd)
    return {'pt_over_kd': ratio, 'regime': 'binding' if ratio <= .1 else 'intermediate' if ratio < 10 else 'titration',
            'thresholds': {'binding_max': .1, 'titration_min': 10.},
            'source': 'Jarmoskaite et al. 2020, doi:10.7554/eLife.57264; operational thresholds'}
