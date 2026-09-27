# Model and scientific interpretation

Version 0.1.0: `Y = B + A*C/(KD+C)`, positive KD and A. Linear-mode fitted B may
be negative; log-mode B is positive. KD is parameterized in log10 M, amplitude in
log10 response-scaled units. Five deterministic initial KD values by default.
Numerical bounds and optimizer success are distinct from assay range and
identifiability. Increasing saturation is supported; decreasing curves require a
separately specified model, not flipping signs silently.

## Assay applicability

Require evidence for equilibrium, signal proportional to binding and a justified
single-site description. Record titrated species, molecular formats/valency,
temperature, presentation/immobilization, controls and signal processing in the
source/rationale. Total concentration requires a supported free-concentration
approximation. Ligand depletion, concentration-dependent nonspecific binding,
multivalency or wash-dependent nonequilibrium responses are not solved by a
constant baseline. Reject an inapplicable model. `apparent_KD` is allowed only
with a specific supported operational interpretation, not as a blanket fallback.

## Baseline, loss and intervals

- `fitted`: fit B jointly. `fixed`: supplied finite B. `zero_control`: arithmetic
  mean of included zero-concentration responses; exclude those controls from the
  optimization objective and residual degrees of freedom. Control-derived B is
  treated as fixed, so label its intervals conditional and do not claim baseline
  uncertainty has been propagated. Multiple controls do not eliminate this issue.
- `linear`: additive normal residuals (scaled internally for numerical stability).
  `log`: additive normal errors in log10 response, implying multiplicative errors;
  all included responses and model baselines must be positive. No pseudocount.
- `profile_f`: approximate individual KD confidence interval using
  `SS(KD)/SS_min = 1 + F(1,n-p;level)/(n-p)`. Other parameters are reoptimized for
  each fixed log KD; residual variance is estimated. For weights this uses
  weighted SS and treats supplied SDs as relative. It is not a simultaneous
  confidence region, nor an exact finite-sample coverage guarantee in a nonlinear
  model. The name deliberately specifies the F convention.
- Endpoint search extends to configured numerical bounds; an unclosed endpoint
  is null with `lower_open`, `upper_open` or `both_open`. These mean unclosed
  within the numerical search domain, not proven infinite intervals. A profile
  computation failure is distinct. Near-zero residual SS yields
  `noise_scale_not_estimable`, not a fabricated zero-width statistical interval.
- Only KD intervals are computed. Baseline/amplitude estimates are saved without
  intervals. Curve confidence bands and prediction intervals are not implemented.

## Repeats and reportability

Fit per curve; no automatic point pooling across curves. Technical curves within
an explicitly identified independent experiment are averaged on log10 KD, then
experiments receive equal weight. Across comparable independent experiments,
report geometric mean and a Student-t interval of the mean log KD, back-transformed.
This assumes independent approximately normal log experimental estimates and
does not perform a meta-analysis that propagates each within-curve interval.
With one experiment, no between-experiment interval. If any curve in a sample
fails or is not reportable, withhold the whole sample estimate rather than
silently excluding it. Dates alone are never experiment IDs.

Point reportability requires successful optimization, local Jacobian full rank,
no numerical boundary hit, KD inside positive concentration range, and a closed
requested interval (or no interval requested / synthetic noiseless limit).
These mechanical criteria do not certify biological suitability. Out-of-range
numeric optima are audit-only; no automatic `KD < min concentration` claim.
Retain open interval endpoints for audit, not as assay limits.

Sensitivity scenarios are explicit (`alternate_loss`, `exclude_highest`), do not
change the primary result and do not supply alternate CIs. Alternate loss is not
supported with weighting. Compare estimates only as sensitivity observations.

Sources: [one-site binding](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_one_site_specific.htm),
[Prism profile interval description](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_how_confidence_intervals_are_c.htm).
This release has no Prism numerical reference project and does not claim equivalence.
