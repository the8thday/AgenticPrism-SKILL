---
name: group-comparison
description: Compare one outcome between independent experimental units - two groups (Welch or paired t test) or three or more independent groups (classic or Welch one-way ANOVA with a predeclared Dunnett, Tukey, Games-Howell or Holm family). Reports differences, simultaneous intervals, adjusted p values and design checks. Use for unit-level antibody-project measurements such as candidates vs isotype control; not repeated-measures, mixed-effects, covariate or nonparametric designs.
---

# Group comparison

Read the [input and design contract](references/input-and-model.md). Establish the true independent unit (animal, donor, independent experiment, etc.), outcome and units before choosing a design. Technical wells or cells are not independent units. Preserve any upstream aggregation and exclusion decisions in the source record; the runtime requires one value per unit per condition.

## Choose the design before looking at results

- **Two groups** (`analysis_type=group_comparison`): `independent_welch` for disjoint units, `paired_t` when every unit is measured in both conditions.
- **Three or more independent groups** (`analysis_type=multi_group_comparison`): choose `design` and `post_hoc` from the scientific question, and record it in `rationale`:
  - many candidates against one control → `one_way_anova` + `dunnett_vs_control`, or `welch_anova` + `holm_welch_vs_control` when group variances may differ;
  - all pairwise differences → `one_way_anova` + `tukey_all_pairs` (Tukey-Kramer), or `welch_anova` + `games_howell_all_pairs`;
  - omnibus only → `post_hoc: none`.
  Do not pick the family, the control or the design after seeing which gives smaller p values. Never run several unadjusted t tests as a substitute for a family.
- The same unit measured under three or more conditions (repeated measures), covariates, nested/clustered units, or rank-based tests are **not implemented**: explain the gap and request a suitable method; do not force such data into these designs.

## Run and inspect

Create the canonical CSV and JSON config ([two-group example](../../fixtures/groups_synthetic/paired_config.json), [multi-group example](../../fixtures/multigroup_synthetic/config_dunnett.json)). Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY` with the runtime located and checked per the [runtime instructions](../agentic-prism/references/runtime.md), then `agentic-prism verify --run DIRECTORY`. Inspect `comparison.csv` (two groups) or `omnibus.csv`, `group_summary.csv` and `contrasts.csv` (multi-group), plus `diagnostics.json` and the offline `report.html`. `render` changes style without recomputing statistics.

## Report

- Two groups: the difference **B minus A**, two-sided p value, confidence interval, sample sizes and design.
- Multi-group: omnibus F, its degrees of freedom and p value; then each contrast in the declared family as **later group minus earlier group** (group minus control for vs-control families) with its simultaneous interval and family-adjusted p value. Contrasts are reported regardless of the omnibus result, because the family was predeclared. Name the family and adjustment method.
- Surface diagnostics: `unequal_group_sd_consider_welch_design` (classic ANOVA with SD ratio above `sd_ratio_warning`, default 3) and `small_group_welch_procedures_may_be_liberal` (Welch design with a group of fewer than 6 units, where simulation showed error slightly above nominal). Do not switch design after the fact to remove a diagnostic; state it.

At least three units per group (or three complete pairs) are required; zero-variance designs are withheld. None of these tests establish normality, biological importance or equivalence. Error-rate evidence is in [validation](../../validation/README.md) and [the 0.7.0 record](../../validation/RELEASE_0.7.0.md).
