---
name: multivariable-regression
description: Fit one outcome against multiple prespecified continuous or categorical predictors with Gaussian linear regression, binomial Logistic regression, or Poisson count/rate regression. Supports explicit reference levels, pairwise interactions, exposure, HC3 or quasi-Poisson uncertainty, diagnostics and conditional mean predictions. Use for independent units; repeated or clustered outcomes require another model.
---

# Multivariable regression and GLM

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input/model contract](references/input-and-model.md) before constructing
a config. Use `analysis_type: generalized_regression`. Locate the installed
runtime with the [ordered procedure](../agentic-prism/references/runtime.md).

## Match the question and design

Establish the outcome, independent unit, predictor meanings/units, intended
conditional effect and prior model specification. Read actual data before
claiming distribution, missingness or independence checks. Multiple predictors
do not make this a multi-outcome or repeated-measures analysis.

- **Gaussian identity:** continuous outcome and additive conditional mean;
  classical OLS inference or explicitly chosen HC3 with approximate t reference.
- **Binomial logit:** one declared positive event, encoded 0/1, or successes and
  integer trial counts. Clarify which event is positive. A percentage alone does
  not supply its denominator. Grouped rows assume a binomial trial mechanism;
  donor-level repeated cells/wells cannot manufacture independent evidence.
- **Poisson log:** integer event counts; known follow-up exposure enters as its
  natural logarithm with coefficient one. An exposure is not an arbitrary
  normalization divisor or a predictor with an estimated coefficient.

Declare predictors and reference categories for the scientific question, not
to search for the smallest p value. Main effects participating in interactions
are conditional on the partner's zero/reference value. Record any preprocessing
and its scientific source; do not silently standardize raw data, log outcomes,
impute, add pseudocounts, delete outliers or turn BLQ into zero.

Repeated donor observations, clustered animals/cells, complete batch-treatment
confounding or a post-treatment adjustment variable can change the estimand or
invalidate inference. Ask only the missing material question. HC3 and Pearson
scaling do not correct dependence, causal confounding or a wrong mean model.
Explain when GLMM/GEE, Firth/penalization or negative binomial is scientifically
appropriate but unavailable here; do not force data into this module.

## Execute and interpret

Use a new output directory, then `verify` and read `interpretation_facts.json`.
Inspect coefficient coding, model rank, convergence, separation, dispersion,
influence and prediction-specific support. A successful process is not a passed
scientific gate. Parent-withheld coefficients remain under `audit`; never use
them as reportable ORs, p values, significant findings or decision thresholds.

Report the conditional coefficients, pointwise uncertainty and independent
sample count. Logistic exponentiated slopes are **odds ratios**, not relative
risks or probability differences. Poisson exponentiated slopes are expected-count
ratios, or rates when an exposure/log offset is declared. An exponentiated
intercept is a baseline odds/count, and an interaction exponentiates to a ratio
of conditional ratios, not a population-average effect.

The Gaussian heteroscedasticity and fixed Poisson/grouped-binomial dispersion
gates are approximate diagnostics. Explicit HC3 and Pearson quasi-Poisson are
separate choices, not automatic repairs. Preserve registered calibration misses,
sparse-data limitations, warnings and failed predictions from the saved facts.
Do not describe Wald intervals or post-selection diagnostic gates as exact.

Requested predictions are conditional **mean** responses with pointwise
intervals. They do not validate predictive performance, classify patients or
provide future-individual prediction intervals. Preserve support withholding;
do not calculate an extrapolated answer elsewhere to bypass it.

Save/report the input, resolved config, coefficient/prediction CSVs, diagnostics,
design matrix, preprocessing log, results, facts and report. Rendering changes
style only. See [0.13.5 evidence](../../validation/RELEASE_0.13.5.md) for method-
specific numerical, calibration, worked-example and behavioral scope.
