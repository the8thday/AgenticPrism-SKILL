"""Predeclared 4PL shape checks. These diagnose adequacy; they never refit/select a model."""
from math import comb
import numpy as np
from scipy.stats import f

SHAPE_ALPHA = .005  # each of two checks; calibrated as a combined gate


def _choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def runs_lower_tail(signs):
    """Conditional P(R <= observed) for random order of fixed counts of +/- signs.

    When used on fitted residuals this is a diagnostic, not an exact model test.
    """
    signs = np.asarray(signs, bool)
    a = int(signs.sum()); b = len(signs) - a
    if not a or not b:
        return None, None
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    ways = 0
    for r in range(2, runs + 1):
        k = r // 2
        if r % 2 == 0:
            ways += 2 * _choose(a - 1, k - 1) * _choose(b - 1, k - 1)
        else:
            ways += _choose(a - 1, k) * _choose(b - 1, k - 1) + _choose(b - 1, k) * _choose(a - 1, k - 1)
    return runs, float(ways / comb(a+b, a))


def shape_checks(x, y, prediction, n_parameters, weighting):
    x, y, prediction = (np.asarray(a, float) for a in (x, y, prediction))
    doses, inverse = np.unique(x, return_inverse=True)
    means = np.array([y[inverse == i].mean() for i in range(len(doses))])
    fitted = np.array([prediction[inverse == i][0] for i in range(len(doses))])
    tol = 1e-12 * max(1., float(np.max(np.abs(y))))
    residual = means - fitted
    informative = residual[np.abs(residual) > tol]
    runs, p = runs_lower_tail(informative > 0)
    run_result = {'method': 'lower_tail_conditional_runs_on_dose_mean_residual_signs',
                  'status': 'computed' if p is not None else 'not_estimable',
                  'n_doses': len(doses), 'n_informative_signs': len(informative),
                  'runs': runs, 'p_value': p, 'alpha': SHAPE_ALPHA,
                  'limitation': 'Fitted residuals are not independent random signs; use the registered calibration, not an exact-test claim.'}
    pe = float(np.sum((y - means[inverse]) ** 2))
    df_pe = len(y) - len(doses); df_lof = len(doses) - n_parameters
    lack = {'method': 'replicate_pure_error_F', 'alpha': SHAPE_ALPHA,
            'df_numerator': df_lof, 'df_denominator': df_pe, 'p_value': None, 'F': None,
            'status': 'not_applicable_relative_weighting' if weighting != 'unweighted' else 'insufficient_pure_error_df',
            'limitation': 'Approximate nonlinear lack-of-fit F under independent homoscedastic Gaussian replicate errors.'}
    if weighting == 'unweighted' and df_pe > 0 and df_lof > 0:
        sse = float(np.sum((y-prediction)**2)); ss_lof = max(0., sse-pe)
        lack.update(pure_error_sse=pe, lack_of_fit_sse=ss_lof)
        if pe <= tol**2 * len(y):
            lack.update(status='pure_error_not_estimable')
            # Exact replicates with material model discrepancy still refute adequacy.
            if ss_lof > tol**2 * len(y):
                lack.update(status='lack_of_fit_without_measurable_pure_error', p_value=0.)
        else:
            statistic = (ss_lof / df_lof) / (pe / df_pe)
            lack.update(status='computed', F=float(statistic), p_value=float(f.sf(statistic, df_lof, df_pe)))
    reasons = []
    if lack['p_value'] is not None and lack['p_value'] < SHAPE_ALPHA:
        reasons.append('four_pl_lack_of_fit')
    if p is not None and p < SHAPE_ALPHA:
        reasons.append('four_pl_structured_residuals')
    return {'passed': not reasons, 'lack_of_fit': lack, 'runs': run_result,
            'diagnostics': reasons,
            'limitation': 'A passing or ineligible check does not establish a 4PL mechanism or exclude a hook outside the measured range.'}
