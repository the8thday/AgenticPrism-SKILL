# Independent factorial and categorical contracts

See complete runnable configurations under `fixtures/routine_0111/`.
Both types require schema_version=1, input, source, confidence_level,
assay.independent_units=true, assay.predeclared=true, unit_definition, rationale.
No automatic exclusion; original inputs and upstream exclusions are retained.

`factorial_anova`, method=two_way: CSV observation_id, independent_unit_id,
factor_a, factor_b, value. All crossed cells populated with >=2 unique units;
no repeated IDs anywhere. Declare factorial.ss_type II/III, ss_reason,
interaction_alpha, factor (factor_a/factor_b), scope (simple/marginal),
adjustment (tukey/dunnett/sidak/holm), and control for Dunnett.
EMMs weight the other factor's levels equally. Separate families per stratum;
Holm p values have conservative Bonferroni simultaneous intervals. Shared
Gaussian residual variance is assumed; heteroscedastic stress misses remain.
Interaction is reported first. Material interaction withholds marginal contrasts.

`contingency`: one row per subject, outcome; group for independent tables;
McNemar has before plus outcome in the same row. Binary outcomes use strings
0/1; Fisher and chi-square also accept R×C categories. Methods fisher,
chi_square, trend, mcnemar, mcnemar_exact, wilson, clopper_pearson.
Categorical config declares correction, seed, resamples (>=999),
null_proportion, and scores={group:number} for trend. Fisher 2×2 uses conditional
MLE OR with exact limits; larger tables use fixed-margin Monte Carlo with seed
and MCSE. Binary tables include Newcombe difference and Katz log-RR intervals;
zero-event RR intervals are withheld without an automatic correction.
McNemar units are independent pairs, not independent pre/post measurements.

Read facts and diagnostics before writing a conclusion. A nonsignificant test
is not equivalence; no causal claim follows without a causal design.
