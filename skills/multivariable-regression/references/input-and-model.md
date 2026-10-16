# Independent-unit multivariable regression contract (0.13.5)

Search terms: generalized_regression, gaussian, binomial, poisson, HC3, pearson,
reference, trials, exposure, offset, missingness, interactions, separation.

## Inputs and configuration

CSV: exactly one row per independent unit, `unit_id` (nonempty, unique string),
numeric `outcome`, declared predictor columns, and optional declared trials or
exposure/log-offset column. `exclude` is optional true/false with a required
`exclusion_reason` for true. Raw bytes are copied before fitting. Units and the
meaning of the positive event belong in the source/design declaration.

Example for a continuous response; binomial and count fixtures are also in
[`fixtures/regression_0135`](../../../fixtures/regression_0135).

```json
{
  "analysis_type": "generalized_regression",
  "input": "study.csv",
  "source": "Study/protocol/instrument export and known prior processing",
  "design": {
    "unit": "donor",
    "rationale": "One measurement per independent donor; no repeated samples",
    "outcome": "Functional response",
    "outcome_unit": "AU",
    "estimand": "Conditional treatment difference adjusted for pretreatment baseline"
  },
  "model": {
    "family": "gaussian",
    "rationale": "Prespecified mean model; baseline measured before treatment",
    "predictors": {
      "baseline": {"kind": "continuous", "unit": "AU"},
      "arm": {"kind": "categorical", "levels": ["control", "treated"], "reference": "control"}
    },
    "interactions": [],
    "covariance": "model_based",
    "confidence_level": 0.95
  },
  "missing": {"policy": "error", "rationale": null},
  "prediction": {"rows": []},
  "report": {"plot_style": "prism_like"}
}
```

- `schema_version` defaults to 1. Unknown settings are refused. Required design
  fields and model rationale must be grounded in the experiment, not a copied fixture.
- Predictors use safe column identifiers (letters, digits, underscores; start
  with a letter); `unit_id`, `outcome`, `exclude`, `exclusion_reason` and
  `prediction_id` and `Intercept` are reserved. Formula strings/executable expressions
  are not accepted. Encoded coefficient labels must remain unique; ambiguous
  categorical labels that collide with interaction names are refused.
- Each continuous predictor requires its unit. Each categorical predictor
  requires distinct string levels and an explicit reference; every declared
  level must occur among included complete cases. Labels such as `NA` remain
  categorical labels when declared. Numeric blanks, `NA` and `NaN` are missing;
  `<LLOQ`, infinity and other nonnumeric qualifiers are errors, not zeros.
- Always includes an intercept. Indicator columns omit the reference. Optional
  `interactions: [["baseline", "arm"]]` adds all cross-products of these two
  declared main effects; duplicate pairs and undeclared partners are refused.
  No three-way terms, splines, automatic transformations or variable selection.
- Missing policy: `error` by default; `complete_case` requires a rationale.
  Save every dropped ID separately from explicit exclusions; no imputation or
  inference of MAR/MCAR from the table. Inference describes the retained sample
  under the declared selection assumptions.
- One outcome only. No frequency/variance/survey weights or robust clustered SE.

## Families and uncertainty

| Family | Outcome and link | Covariance and reference |
|---|---|---|
| `gaussian` | Real-valued outcome; identity | `model_based`: OLS residual variance, Student t; `HC3`: leverage-adjusted sandwich, residual-df t (approximate) |
| `binomial` | 0/1 without trials, or integer successes `0 <= outcome <= trials`; logit | `model_based` only; binomial variance, asymptotic normal Wald |
| `poisson` | Nonnegative integer counts; log | `model_based`: Poisson variance, normal Wald; `pearson`: quasi-Poisson Pearson/residual-df scale and t reference |

`trials` is a column name, only for binomial; every included denominator must be
a positive integer. It is not the number of technical replicates, a percentage
or an inferred effective sample size. All-zero/all-one binary data cannot fit a
finite unpenalized Logistic model.

For Poisson, optionally declare **one** of `exposure` or `offset` (column name),
plus `offset_rationale` including the meaning, unit and source. Exposure must be
finite and positive and is logged; offset is already on the natural-log scale.
They are fixed additions to the linear predictor, not fitted predictors. A rate
as outcome must not be substituted for integer counts. A zero-count observation
is valid; zero exposure is not. Pearson scaling does not fit a negative-binomial
distribution and does not supply a likelihood or AIC.

All coefficient tests are two-sided versus zero on the linear-predictor scale.
Intervals are pointwise at the declared confidence level; no adjustment across
predictors/terms. Grouped coefficient Wald tests use F for Gaussian/quasi-Poisson,
chi-square for fixed-dispersion GLMs. They test the coded coefficients jointly;
with interactions these main-effect hypotheses remain reference-conditional,
not Type II/III tests or marginal population effects.

## Diagnostics and reportability

The implementation uses scaled design columns for numerical fitting and returns
coefficients/covariance in original input units. Saves the original design matrix.
Rank deficiency, residual df < 2, scaled condition number > 1e10, likelihood
separation/recession, failed convergence, invalid uncertainty or leverage one
withhold reportable inference. The linear-program separation check includes
quasi-separation and Poisson zero-count strata with no finite MLE.

Classical Gaussian inference is withheld when the robust Breusch-Pagan diagnostic
p < .01. Fixed Poisson or grouped-binomial inference is withheld when the upper-
tail Pearson chi-square dispersion diagnostic p < .01. These are approximate
diagnostics, especially with sparse counts, not formal certificates or exact
post-selection coverage guarantees. No meaningful grouped-dispersion test is
claimed for Bernoulli rows. No automatic HC3/quasi substitution follows a failure.

High leverage (>2p/n), approximate Cook distance (>4/n), and fewer than ten
successes/failures per coefficient are **review heuristics**, not automatic
exclusions or sample-size prescriptions. Diagnostics retain unit IDs. Passing
checks cannot establish independence, causal identification or correct functional
form. Extra-binomial data may require a model not available in this release.

## Prediction and saved artifacts

`prediction.rows` contains explicit dictionaries with a unique `prediction_id`,
all declared predictors and any exposure/log-offset column, and no outcome/trials.
Unknown categories or invalid values fail validation. Mean-response intervals
are obtained on the link scale then transformed; binomial means are probabilities,
Poisson means are counts for the supplied exposure, Gaussian means use outcome units.
No individual prediction interval or prediction-performance metric is produced.

Only rows within the convex hull of the expanded observed design, including any
offset, receive reportable predictions. This conservative interpolation check
may withhold an otherwise scientifically justified prediction; it does not prove
transportability or validate an unobserved joint predictor distribution. A failed
parent also withholds predictions; each requested row keeps its own reason.

`results.json` stores reportable `primary`, separate `audit`, matrix/covariance,
diagnostic observations and notes. `coefficients.csv` contains reportable rows
only; withheld estimates remain explicitly audit-only in results. Also save
`predictions.csv`, `design_matrix.csv`, `diagnostics.json`, `preprocessing_log.json`,
raw input, resolved config, extraction-only facts and hash manifest. Re-render
uses saved results and never refits or changes scientific artifacts.

## Primary references

- [statsmodels GLM API](https://www.statsmodels.org/stable/generated/statsmodels.genmod.generalized_linear_model.GLM.html)
- [R glm](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/glm.html)
- [R summary.glm dispersion and inference](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/summary.glm.html)
- [Published Star98 worked example](https://www.statsmodels.org/stable/examples/notebooks/generated/glm.html)

Actual validation and retained misses are in the [release record](../../../validation/RELEASE_0.13.5.md),
not established by these documentation links.
