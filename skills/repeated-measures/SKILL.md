---
name: repeated-measures
description: Analyze a continuous outcome measured repeatedly on the same donors, animals, or independent experiments - one within-unit factor (three or more conditions), or a between-unit arm by within-unit condition design such as treatment arm x day (body weight, tumor volume at scheduled days). Supports complete RM ANOVA and split-plot ANOVA with Greenhouse-Geisser correction, and Gaussian random-intercept REML models with available-case fitting and Satterthwaite (lmerTest-equivalent), parametric bootstrap or asymptotic Wald inference. Also supports marginal MMRM with unstructured or ordered-visit AR(1) covariance, numeric baseline covariates, Python Satterthwaite or optional R Kenward-Roger. Not for random slopes, continuous-time growth rates, nested/crossed effects or non-Gaussian outcomes.
---

# Repeated measurements

Read the [input and model contract](references/input-and-model.md). Establish the
estimand (condition mean differences), independent unit, within-unit conditions,
response scale and missingness mechanism before selecting a model. Technical
wells and cells do not become independent donors. One observation per unit and
condition is required; record any upstream aggregation or transformation.

For a marginal model with unstructured or AR(1) covariance and optional
Kenward–Roger inference, follow the separate [MMRM contract](references/mmrm.md).
It uses `analysis_type=mmrm`; do not add MMRM options to legacy configs.

## Choose the model from the experiment

- **Complete one-factor within-unit design:** `one_way_rm_anova`. Always use the
  Greenhouse–Geisser corrected omnibus p; the uncorrected p is audit information.
  Predeclare no comparisons, a control family, or all pairs. Contrasts use paired
  t tests with Holm p values and Bonferroni simultaneous intervals.
- **Gaussian random-intercept design:** `random_intercept`. Condition means are
  fixed effects; one intercept per independent unit is random. Declare why equal
  residual variance and no residual time correlation are reasonable. Missing
  measurements may use `available_case` with an explicit MAR rationale; MAR is
  conditional on this model and cannot be verified from observed data alone.
  Inference: `satterthwaite` (small-sample t/F with Satterthwaite df; numerically
  equal to R lmerTest; deterministic) or `parametric_bootstrap` (999 refits by
  default, model-based, not distribution-free) are the small-sample choices. In the
  null simulations of the validation record (8 and 12 units, four conditions, 15%
  MCAR) Satterthwaite stayed within the prespecified bounds (omnibus 4.7%/5.0%,
  simultaneous coverage 96.0%/94.9%); the bootstrap, run there with 199 refits,
  held the error bounds but missed the coverage bound at 12 units (93.4% against a
  93.6% minimum). Prefer Satterthwaite; keep the default 999 refits if the
  bootstrap is used. The 0.8.2 B=999 results are given below. Satterthwaite is withheld
  when the random-intercept variance is at its zero boundary (a third of the
  simulated datasets with weak clustering, subject SD 0.3 vs residual SD 1);
  bootstrap still runs there and held its bounds. `wald_asymptotic` has no small-sample correction and was liberal
  (about 7–9% at nominal 5% with 8–24 units); use it only when explicitly
  justified. Kenward–Roger is available only through the separate MMRM contract.
- **Between-unit arm × within-unit condition** (for example treatment arm × scheduled
  day, each animal in one arm): `two_way_rm_anova` for complete data (split-plot
  ANOVA, GG-corrected condition and interaction tests), or `two_way_mixed` (cell
  means + unit random intercept, `satterthwaite` or `wald_asymptotic`) when
  measurements are missing. Declare `arms`, and the family: `arm_vs_control_each_condition`
  (with `control_arm`) or `condition_vs_first_each_arm`. The CSV needs an `arm` column;
  `group` stays the within-unit condition.
- Days are treated as categories. A growth-rate question (slope over continuous
  time, tumor-growth modelling, TGI), random slopes, irregular-time
  autocorrelation, batch plus donor effects, crossed/nested random effects,
  covariates, counts or binary outcomes need another model. Do not flatten these
  designs, discard the grouping factor, or use independent-group ANOVA. Two
  complete conditions can use the paired t module.

Never select a covariance structure, contrast family, missing-data policy or
inference method after comparing p values. Assess the scientific design, not
just whether the software converges. A missingness rationale must come from the
experiment; never copy the synthetic fixture's MCAR declaration into real data.

## Run and inspect

Start from [complete RM](../../fixtures/repeated_synthetic/config_rm.json) or
[mixed model](../../fixtures/repeated_synthetic/config_mixed.json). Resolve the
runtime using [shared instructions](../agentic-prism/references/runtime.md).
Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`, then
`agentic-prism verify --run DIRECTORY`.

Inspect `omnibus.csv`, `contrasts.csv`, `missingness.csv`, `diagnostics.json`,
`results.json` and `report.html`. Mixed models also save `fixed_effects.csv`,
`residuals.csv` and `bootstrap_samples.json`. All raw inputs, exclusions,
configuration, dependency versions and scientific hashes are retained. Render
only restyles saved results and does not refit or resample.

## Report

For new 0.8.1 runs, read `interpretation_facts.json` after verification and base
the narrative on its primary results, reportability and required disclosures
([shared contract](../agentic-prism/references/interpretation-facts.md)). Lead
with the corrected omnibus/interaction and the declared difference family with
intervals. If inference is withheld, say why; retained fixed-effect means or
audit SEs do not restore a p value or interval. If the experimental unit,
covariance or missingness rationale is unclear or contradicts the data, ask
the user before interpreting: a successful run cannot validate those assumptions.

Name the independent unit and its count, observations and missing cells, design,
covariance assumption, estimation/inference methods, contrast direction B−A and
family correction. For arm × condition designs report the interaction test
first: if it is significant, the arm difference changes across conditions, so
interpret arm differences through the per-condition contrasts rather than the
arm main effect; a non-significant interaction is not evidence of parallel
profiles. Say whether a day-0 difference exists before treatment effects are
discussed. Report corrected omnibus p and simultaneous contrast CIs
without screening contrasts by omnibus significance. Mixed-model estimated
condition means differ from observed means when data are incomplete.

Surface near-zero random-intercept variance, optimizer failures, fewer than 30
units, failed bootstrap refits and Monte Carlo uncertainty. The 30-unit warning
is a heuristic, not a reliability threshold. At least 95% of bootstrap refits
and at least 199 must succeed; otherwise inferential p values/CIs are withheld.
Audit point estimates and SEs may remain. Never promote withheld inference.
A boundary estimate does not prove that clustering is absent. Residual and QQ
plots help inspect assumptions but do not validate them.
When discussing a bootstrap fallback after Satterthwaite failure, retain its
evidence limits: B=199 missed one coverage bound; the B=999 results below
apply only to the stated Gaussian compound-symmetry/MCAR simulation designs. Do not call it validated for the current experiment or
promise valid inference merely because it runs. A revised analysis needs a
scientific rationale; keep the original withheld result and identify the
revision. Do not infer that paired observations became independent at a fitted
variance boundary.

[Validation](../../validation/RELEASE_0.8.0.md) describes the R numerical
benchmark, simulated designs and remaining limits. Do not extrapolate its
coverage results to other covariance structures or missingness mechanisms.

### 0.8.2 bootstrap B=999 observations

Each fixed row generated 1,000 datasets (four conditions, residual SD 1, 15% MCAR, comparisons versus the first; seed 20260929). Bounds are rejection ≤ 6.3784% and simultaneous coverage ≥ 93.6216% when all 1,000 fits are reportable.

| Units / subject SD | Reportable | Omnibus / family rejection | Simultaneous coverage | Outcome |
|---|---:|---|---:|---|
| 8 / 2 | 1000 | 4.4% / 4.4% | 95.5% | passed |
| 12 / 2 | 1000 | 4.8% / 5.9% | 94.0% | passed |
| 12 / 0.3 | 1000 | 4.5% / 4.4% | 95.6% | passed |

The 12-unit/SD=2 row now has 94.0% coverage; its B=199 value was 93.4%. Preserve the historical miss. These simulation results do not validate the current experiment, another covariance or MNAR dropout. Full denominators and bounds are in [0.8.2 evidence](../../validation/RELEASE_0.8.2.md).

### Explain variance-boundary withholding precisely

A zero random-intercept estimate triggers this implementation's Satterthwaite
withholding rule. It does not prove that fixed-effect SEs are universally invalid,
that every alternative method fails, or that a hypothetical Wald calculation
necessarily treats observations as independent. Contrasts require their full
covariance and a justified inference method; marginal SE columns alone do not
supply them. Describe the actual saved withholding rule. Distinguish failed
optimizer attempts from the final converged fit. Clustering is within donors,
not between donors.

中文表述中，“聚集”指**同一供体内**重复观测的相关性。不要写成“供体之间的
聚集性”：这里供体是声明的独立单位。方差边界只触发本实现的推断扣留规则。

Do not list `fewer_than_30_units_small_sample_inference_requires_caution` as a
cause of withholding. It is a caution only, even when it appears beside genuine
blocking diagnostics in the same facts/result object. In a "why withheld"
table, list the boundary/singular-inference rule or actual fit failure; put the
sample-count warning separately. 中文：少于30个独立单位只是警告，不能写成扣留
p值/区间的原因；把它与真正触发扣留的方差边界或拟合失败分开说明。
