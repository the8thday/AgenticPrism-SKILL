# Nonparametric unit-level comparisons (0.8.2)

Use `analysis_type: nonparametric` and the existing canonical unit CSV:
`observation_id, independent_unit_id, group, outcome, value, unit`, optional
`exclude` and `exclusion_reason`. Exactly one observation per unit/group;
no silent aggregation or removal of incomplete pairs/blocks. Exclusions need
reasons. The same units must occur in every paired/Friedman condition;
independent groups must have disjoint IDs. At least three units per group.

Start from `fixtures/nonparametric/config_fixture_*.json`. These are synthetic
examples; establish applicability from the actual experiment.

`comparison` fields (unknown keys are refused):

- `design`: `mann_whitney`, `wilcoxon_signed_rank`, `kruskal_wallis`, `friedman`.
- `groups`: ordered names, two for MW/signed-rank, ≥3 for KW/Friedman.
- `outcome`, `unit`, `rationale`: nonempty strings. Identify units, pairing,
  planned hypothesis and distribution assumptions in the rationale.
- `independent_units`, `rank_model_appropriate`: must be literal `true`.
  For paired data independence is between pairs/blocks, not within them.
- `location_shift_or_symmetry`: literal `true` for MW/signed-rank, otherwise null.
  MW location inference assumes comparable distribution shapes and a common
  shift; signed ranks require symmetric paired differences. A non-significant
  normality test does not establish these assumptions.
- `inference`: `auto` (default), `exact`, or `asymptotic` for MW/signed-rank.
  Auto uses exact inference when each sample has <50 observations and has no
  pooled ties (MW), or no absolute-difference ties/zero differences (paired).
  Otherwise it explicitly uses the normal approximation, tie correction and
  continuity correction. Exact mode refuses ties/zeros or larger samples.
  This scope differs from R 4.6's conditional exact option with ties. KW and
  Friedman require `auto` and use chi-square approximation, with tie correction.
- `confidence_level`: default .95, between .5 and 1 for MW/signed-rank.
  KW/Friedman do not produce such intervals; keep the default .95.
- `post_hoc`: `none` or `dunn_all_pairs` for KW only. The family is predeclared,
  reported regardless of the omnibus result. No Friedman post-hoc inference.
- `adjustment`: `holm` (default) or `bonferroni` for Dunn. Without Dunn, keep
  the default. Canonical monotone Holm adjusted p values match `stats::p.adjust`;
  the `dunn.test` package's raw Holm step values require its stopping rule.

## Interpretation

Verify the run, then read `interpretation_facts.json`, `omnibus.csv`,
`contrasts.csv`, `group_summary.csv` and diagnostics. Lead with the estimand and
independent sample/block count. MW/signed-rank direction is **B minus A**.
`hodges_lehmann` is the median of pairwise B−A differences (MW) or Walsh averages
of paired B−A differences. It is not generally a difference of group medians.
`estimate` matches the R location-estimate convention: HL for exact inference,
normal rank-inversion approximation otherwise; retain this distinction.

Report the interval method and saved `achieved_confidence_level` with its
`confidence_level_basis`: attained confidence for exact inference, but nominal
normal-approximation confidence for asymptotic inference, not measured coverage. Small samples
may yield an unbounded exact interval or an approximation below the requested
confidence. Null endpoints never authorize invented finite limits. All-tied or
all-zero designs may withhold inference. Normal/chi-square approximations are
flagged for small samples; see the release's calibration scope and misses.
Dunn contrasts are **pooled mean-rank differences**, with adjusted p values;
no raw-unit shift intervals are provided. Friedman tests a complete block
condition effect, not a treatment-by-time interaction.

Stop and ask when the independent unit, matching, exclusions, distribution
assumptions or intended estimand cannot be established. Do not switch tests,
round data, remove zeros/ties, choose comparisons or relabel technical wells as
independent units to obtain significance. Non-significance does not establish
equivalence. Rank tests do not correct informative dropout or confounding.

For death-related or outcome-related missingness, stop and establish the target
population before excluding incomplete blocks, even if the user asks to delete
them. A documented deletion does not by itself justify a complete-case estimand.
Do not suggest separate pairwise signed-rank runs as a multiplicity-controlled
Friedman post-hoc family: that family is not implemented here.

The interval targets the population location parameter, not a particular point
estimator. Reporting HL with a clearly named approximate rank-inversion interval
is appropriate under the declared model. Do not say that HL has no interval
merely because `estimate` and `hodges_lehmann` differ numerically. Endpoints within
1e-4 of zero are within solver precision; do not claim exclusion of zero from
such numerical tails.
The normal rank-inversion point estimate is not necessarily the interval's
midpoint: the interval may be asymmetric. Do not call `estimate` its center.
If the saved interval spans zero, it is numerically inclusive of zero. Solver
precision limits an exclusion claim; it does not justify saying the displayed
interval cannot be described as including zero. Prefer "compatible with zero"
for an endpoint within solver tolerance (中文：与零相容，不能据此排除零).
