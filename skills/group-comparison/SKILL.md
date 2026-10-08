---
name: group-comparison
description: Compare independent-unit outcomes with one-sample, Welch, paired or ratio t tests, one/two-way ANOVA and declared multiple-comparison families, scoped rank tests, nested comparisons and baseline ANCOVA. Use for antibody candidates vs controls, recovery vs a reference and geometric-mean ratios. Preserves biological units, pairing, effect scales and suitability gates; general repeated measures use their specialist.
---

# Group comparison

For a single sample versus a declared reference, or an explicitly multiplicative
effect (geometric means / treatment-to-control ratio), use the
[single-sample and ratio-test contract](references/location-tests-0134.md),
`analysis_type=location_test`. Preserve positive-value and complete-pair gates.

Read the [input and design contract](references/input-and-model.md). Establish the true independent unit (animal, donor, independent experiment, etc.), outcome and units before choosing a design. Technical wells or cells are not independent units. When the data hold several technical replicates per independent unit, use the nested contract below instead of averaging by hand or treating wells as the sample size. Preserve any upstream aggregation and exclusion decisions in the source record; the other designs require one value per unit per condition.

For rank-based inference, read the [nonparametric contract](references/nonparametric.md):
MW/signed-rank with HL, KW with Dunn, or complete-block Friedman. Establish the
shift/symmetry assumptions and pairing before running. Results include verified
interpretation facts.

## Choose the design before looking at results

- **Two groups** (`analysis_type=group_comparison`): `independent_welch` for disjoint units, `paired_t` when every unit is measured in both conditions.
- **Three or more independent groups** (`analysis_type=multi_group_comparison`): choose `design` and `post_hoc` from the scientific question, and record it in `rationale`:
  - many candidates against one control → `one_way_anova` + `dunnett_vs_control`, or `welch_anova` + `holm_welch_vs_control` when group variances may differ;
  - all pairwise differences → `one_way_anova` + `tukey_all_pairs` (Tukey-Kramer), or `welch_anova` + `games_howell_all_pairs`;
  - omnibus only → `post_hoc: none`.
  Do not pick the family, the control or the design after seeing which gives smaller p values. Never run several unadjusted t tests as a substitute for a family.
- The same unit measured under three or more conditions routes to [repeated-measures](../repeated-measures/SKILL.md), which supports one-factor RM ANOVA and scoped random-intercept mixed models. Covariates and nested/crossed effects remain outside these modules; do not force such data into independent-group designs.

## Run and inspect

Create the canonical CSV and JSON config ([two-group example](../../fixtures/groups_synthetic/paired_config.json), [multi-group example](../../fixtures/multigroup_synthetic/config_dunnett.json)). Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY` with the runtime located and checked per the [runtime instructions](../agentic-prism/references/runtime.md), then `agentic-prism verify --run DIRECTORY`. Inspect `comparison.csv` (two groups) or `omnibus.csv`, `group_summary.csv` and `contrasts.csv` (multi-group), plus `diagnostics.json` and the offline `report.html`. `render` changes style without recomputing statistics.

## Report

- Two groups: the difference **B minus A**, two-sided p value, confidence interval, sample sizes and design.
- Multi-group: omnibus F, its degrees of freedom and p value; then each contrast in the declared family as **later group minus earlier group** (group minus control for vs-control families) with its simultaneous interval and family-adjusted p value. Contrasts are reported regardless of the omnibus result, because the family was predeclared. Name the family and adjustment method.
- Surface diagnostics: `unequal_group_sd_consider_welch_design` (classic ANOVA with SD ratio above `sd_ratio_warning`, default 3) and `small_group_welch_procedures_may_be_liberal` (Welch design with a group of fewer than 6 units, where simulation showed error slightly above nominal). Do not switch design after the fact to remove a diagnostic; state it.

At least three units per group (or three complete pairs) are required; zero-variance t/ANOVA designs are withheld; rank-test handling is described in the separate contract. None of these tests establish normality, biological importance or equivalence. Error-rate evidence is in [validation](../../validation/README.md) and [the 0.7.0 record](../../validation/RELEASE_0.7.0.md).

## 0.11.1 independent factorial and categorical designs

Use [routine contract](references/routine-0111.md) for `factorial_anova` or
`contingency`. Require one row per independent subject. Refuse repeated units
in two-way independent ANOVA and route to repeated-measures. Declare SS II/III
and its reason; lead with interaction, and withhold marginal contrasts when it
is material. Contrast families and thresholds must be declared before analysis.
For categorical outcomes, technical wells cannot supply independent subject
counts. Disclose expected counts below five and exact versus Monte Carlo
methods. All results, including older two/multi-group runs, have saved facts.

## Nested replicates (0.13.2, `analysis_type=nested_comparison`)

Use when each independent unit (mouse, donor, independent culture) contributes several technical
replicates (wells, fields, organoids) and each unit belongs to one group. Declare the unit and the
replicate with rationales; a unit appearing in two groups is a repeated-measures design and is refused.
The run fits units as a random intercept (REML, Satterthwaite t/F as lmerTest) and always reports the
unit-means analysis; at a boundary fit (between-unit variance estimated as zero) the unit-means analysis
is primary. Report the number of units as n, the ICC and the design effect, never the number of wells.
In calibration a pseudo-replicated ANOVA on the same null data rejected 26% of the time.
Contract: [nested and ANCOVA](references/nested-ancova-0132.md).

## Baseline-adjusted comparison (0.13.2, `analysis_type=ancova`)

Use when one value per independent unit is compared between groups after adjusting for 1-3 covariates
measured **before** treatment (randomization tumor volume, baseline body weight). Refuse covariates
measured after dosing: they can absorb the treatment effect. The slope-homogeneity alpha is declared
before analysis; if slopes differ, adjusted comparisons are withheld (the group difference then depends on
the covariate). Report adjusted means at the overall covariate mean, adjusted differences with the
declared family, the covariate slope, and the unadjusted comparison beside them; mention covariate
imbalance or extrapolation flags. Repeated outcomes with a baseline belong to MMRM.

Evidence (0.13.2): nested agrees with lmerTest within 1.7e-7 and ANCOVA with lm/emmeans/car within 4.2e-12;
all registered calibration rows passed (nested coverage 0.946-0.963, type I 0.041; ANCOVA coverage
0.947/0.954, type I 0.048, slope-test withholding 0.044, power 0.38 adjusted vs 0.21 unadjusted).
Published worked example and live-agent gates: see [the 0.13.2 record](../../validation/RELEASE_0.13.2.md).
