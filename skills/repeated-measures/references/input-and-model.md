# Repeated-measures contract (0.8.0)

## Applicability and estimands

One continuous outcome; independent donors/animals/experiments measured under
three or more declared conditions. Each condition is a within-unit category,
not a separate treatment arm. Estimands are condition population means and
predeclared differences B−A. Input values remain in the supplied physical unit.
If the scientific endpoint requires log transformation, do it upstream and
record its meaning; the runtime does not transform responses or back-transform
contrasts silently.

One between-unit factor (`arm`) crossed with the within-unit condition is supported
by the two-way designs below. No random slope, crossed or nested effect,
covariate, heteroscedastic residual, AR(1), non-Gaussian link or unstructured
repeated-measures covariance is implemented. The random-intercept
model is not a general MMRM. Changing to it does not by itself solve violation
of sphericity or informative dropout.

## CSV

Required: `observation_id,independent_unit_id,group,outcome,value,unit`.
`group` is the within-unit condition. Optional `exclude=true|false` plus
`exclusion_reason`; exclusions require a reason. Values must be finite, IDs
nonempty and unique at observation level. Included rows are unique by unit and
condition. Technical replicates must be aggregated explicitly upstream.

Missing measurements are absent rows, not zeros or strings such as `NA` in the
numeric value column. Literal string identifiers such as `NA` are preserved.
No mean imputation, interpolation or silent complete-case deletion occurs.
Exclusions are retained in the input snapshot and preprocessing log. Conditions
in the file must exactly equal the declared conditions. All included values
must have the declared outcome and unit. Entirely absent subjects cannot be
inferred from this table; provide a separate study disposition record when
accounting for all enrolled subjects.

## JSON

```json
{
  "analysis_type": "repeated_measures",
  "input": "measurements.csv",
  "source": "Experiment identifier and upstream preparation history",
  "comparison": {
    "design": "random_intercept",
    "groups": ["baseline", "condition-A", "condition-B"],
    "outcome": "response",
    "unit": "pg/mL",
    "rationale": "One unit-level value per donor per condition; differences against baseline were specified in advance",
    "post_hoc": "vs_control",
    "control": "baseline",
    "confidence_level": 0.95,
    "missing_policy": "available_case",
    "missingness_rationale": "Replace with experimental evidence supporting MAR given this model",
    "covariance_rationale": "Replace with why donor intercepts and homogeneous independent residual errors fit this experiment",
    "inference": "parametric_bootstrap",
    "bootstrap_reps": 999,
    "seed": 20260927
  }
}
```

For RM ANOVA set `design=one_way_rm_anova`, `missing_policy=require_complete`,
and `inference=greenhouse_geisser` (also the default for this design).
Mixed-model inference must be explicitly specified. `post_hoc` is `none`,
`vs_control` or `all_pairs`; `control` is allowed only for `vs_control`.
Confidence level lies strictly between 0.5 and 1. Bootstrap repetitions are an
integer 199–10000, seed an unsigned 32-bit integer. More repetitions reduce
Monte Carlo noise, not model misspecification. Near a decision threshold,
consider the recorded Monte Carlo SE and a prespecified higher repetition count;
do not rerun different seeds until significant.

RM requires ≥3 complete independent units. Mixed models require ≥6 observed
units per condition and ≥6 units with repeated values. Included unit/condition
overlaps must form a connected design. These are computability guards, not
sample-size sufficiency guarantees. Entirely excluded units do not enter the fit.

## Repeated-measures ANOVA

With n subjects and k conditions, decompose SS into subject, condition and
subject-by-condition error. F = MS_condition / MS_error with uncorrected df
(k−1, (n−1)(k−1)). For S, the sample covariance across conditions, and
H=I−11'/k, epsilon_GG = tr(HSH)^2 / [(k−1) tr((HSH)^2)], clipped to
[1/(k−1),1]. Multiply both df by epsilon; retain F unchanged. Always report the
GG-corrected p rather than switching according to a sphericity pretest.
`partial_eta_squared` = SS_condition/(SS_condition+SS_error); it is descriptive,
with no interval. It is not generalized eta squared.

Contrast SE comes from each paired difference and uses n−1 df. Holm adjusts
p values across the declared family; Bonferroni t intervals give simultaneous
coverage under the paired normal-difference assumption without requiring
sphericity. Zero-variance paired differences have no reported test/interval;
the original family size is preserved for the other comparisons.

## Linear mixed model

`y_ij = mu_j + b_i + e_ij`, with b_i ~ N(0,tau²), e_ij ~ N(0,sigma²), mutually
independent. The implied within-unit covariance is compound symmetry with
nonnegative correlation. REML variance components come from statsmodels
`MixedLM` (BFGS; Powell fallback on failed convergence), after centering/scaling
for numerical stability. The result is restored to the original units.

For V_i = sigma² I + tau² 11', beta_hat is GLS and its plug-in covariance is
C=(X'V^-1 X)^-1. Cell-means coding gives one coefficient per condition. SEs do
not themselves correct for uncertainty in variance components. BLUP donor
intercepts, marginal/conditional predictions and residuals are saved.

- `satterthwaite`: for contrast l, t = l'beta / sqrt(l'Cl) with
  df = 2 (l'Cl)^2 / (g'Ag), g the gradient of l'C(theta)l in the log variance
  parameters and A = 2 H^-1, H the Hessian of the REML deviance. Multi-df tests
  (omnibus, type III effects) eigen-decompose L C L' and combine per-component df
  (Fai and Cornelius 1996) exactly as lmerTest does. Derivatives are central finite
  differences on an REML optimum polished on the exact deviance. Withheld at a
  zero random-intercept variance. Contrasts: Holm p and Bonferroni t intervals.
- `wald_asymptotic`: overall W = (L beta)'(L C L')^-1(L beta), compared with
  chi-square(k−1). Pairwise z uses c'beta/sqrt(c'Cc); Holm p and Bonferroni normal
  intervals. No denominator df, Satterthwaite or KR adjustment is supplied.
- `parametric_bootstrap`: draw independent Gaussian random intercepts and errors
  using fitted REML variances, on the same observed design/missingness pattern.
  Refit REML each time. Center beta* on beta_hat; compute W* and the t* vector
  with each draw's refitted GLS covariance. Overall p is
  (1 + count(W* ≥ W_observed))/(B_success+1). Contrast p uses |t*|; family p uses
  max |t*| across all declared contrasts. Symmetric simultaneous intervals use
  the empirical confidence-level quantile (`higher`) of max |t*| times each
  observed SE. This is inversion of approximate centered bootstrap pivots,
  **not** a null-model ML likelihood-ratio bootstrap. At finite B, quantile-based
  interval endpoints and plus-one p thresholds may differ by a Monte Carlo rank.
  Bootstrap p resolution and omnibus Monte Carlo SE are recorded. This remains
  a model-based approximation, including at small n or a variance boundary.

Each bootstrap record retains success, boundary, optimizer attempts, W* and
contrast t* values. Boundary estimates are retained, not discarded. If fewer
than 199 or 95% of requested fits succeed, p values and contrast intervals are
withheld while point estimates remain for audit. The bootstrap is conditional
on the observed missingness pattern, not a sensitivity analysis for MNAR.

## Between-unit arm × within-unit condition (0.8.0)

CSV adds a required `arm` column (each unit in exactly one arm); `group` is the
within-unit condition. JSON: `design` `two_way_rm_anova` or `two_way_mixed`,
`arms` (≥2, declared order), `groups` (≥2 conditions), `post_hoc` `none`,
`arm_vs_control_each_condition` (needs `control_arm`) or
`condition_vs_first_each_arm`. See the [synthetic body-weight
example](../../../fixtures/repeated_two_way_synthetic/config_two_way_mixed.json).
Each arm needs ≥3 units and every arm × condition cell ≥3 observed units.

`two_way_rm_anova` (complete data): arm = one-way ANOVA on unit means
(df a−1, N−a). Condition and arm × condition use the multivariate linear model on
orthonormal within-unit contrasts Z = Y C (C: k × (k−1)), with type III
hypotheses (arms weighted equally): F = [tr(H)/df_h] / [tr(E)/df_e],
df_e = (N−a)(k−1), both df multiplied by GG epsilon from the pooled within-arm
covariance E/(N−a). This equals car::Anova type III univariate tests with GG
(afex default). Contrasts: Welch t between arms at a condition, paired t within
an arm between conditions; Holm p, Bonferroni intervals.

`two_way_mixed`: one fixed mean per arm × condition cell plus a unit random
intercept (compound symmetry), REML. Type III tests are equal-weight cell-mean
hypotheses (arm, condition, arm × condition) with Satterthwaite F (identical to
lmerTest `anova(type=3)` with sum-to-zero coding) or Wald chi-square. Contrasts
are cell-mean differences with Satterthwaite t (identical to emmeans with
`lmer.df = "satterthwaite"`) or Wald z. `tests.csv` holds the three tests; the
interaction p is recorded as the primary test. Parametric bootstrap is not yet
available for two-way designs.

## Sources and validation scope

- [statsmodels MixedLM](https://www.statsmodels.org/stable/mixed_linear.html):
  ML/REML implementation and its inference interfaces.
- [statsmodels AnovaRM](https://www.statsmodels.org/stable/generated/statsmodels.stats.anova.AnovaRM.html):
  balanced within-subject F benchmark; this API does not supply GG corrections.
- [Pingouin epsilon source](https://github.com/raphaelvallat/pingouin/blob/main/src/pingouin/distribution.py):
  GG covariance formula reference; not a runtime dependency.
- [R nlme](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/lme.html):
  independent numerical benchmark, used only in validation.
- [Halekoh and Højsgaard, 2014](https://www.jstatsoft.org/article/view/v059i09):
  small-sample mixed-model inference context. Their KR and likelihood-ratio
  bootstrap methods are not implemented or claimed equivalent here.
- [Kuznetsova, Brockhoff and Christensen, 2017](https://www.jstatsoft.org/article/view/v082i13)
  (lmerTest): Satterthwaite df method reproduced here; R lmerTest, afex, car and
  emmeans are validation-only references.

See [0.8.0 validation](../../../validation/RELEASE_0.8.0.md) for actual observed
agreement and simulations. Synthetic data do not establish suitability for an
assay, informative dropout, non-Gaussian responses, random slopes or serial
correlation.
