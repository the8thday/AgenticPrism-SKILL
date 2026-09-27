# Precision contract — 0.9.1

`analysis_type=variance_components`. Start from
`fixtures/variance_reml/config_unbalanced.json`; do not reuse its applicability
assertions as real evidence. All defaults and validation live in `precision.py`.

CSV: unique nonblank `observation_id`, finite numeric `value`, and each declared
factor/covariate column. Factor labels are strings (including literal NA).
Missing values and duplicate IDs refuse. Preserve units and upstream provenance
in the required `source` and design rationale. No automatic exclusions.

Required config: `input`, `source`, `response_unit`; `design.random_terms` (1–8 lists of column
names), `design.rationale`; literal true `independent_random_effects`,
`gaussian_homoscedastic_errors`, `representative_levels`. Optional
`fixed_numeric=[]`, `fixed_categorical=[]`; categorical treatment coding uses
sorted first level as reference. No formula strings or automatic interactions.
Nested effects use explicit composite terms, crossed effects separate terms.
Each random term needs >=3 levels with replication, full-rank fixed design and
identifiable variance kernels. Maximum 1000 observations for this dense solver.
These are computation gates, not power or adequacy guarantees.

`inference.confidence_level=.95`; `cv_reference=null` disables CV. Positive
`cv_reference` requires `cv_rationale`; its uncertainty is not propagated.
`report.plot_style` is standard or prism_like; report rendering never refits.
Unknown keys refuse. No silent regulatory acceptance limits.

Model: y = X beta + sum Z_j u_j + epsilon, independent scalar Gaussian
u_j ~ N(0, variance_j I), epsilon ~ N(0, residual I). Nonnegative REML optimized
in Python. Expected information is trace(P K_j P K_k)/2. VCA-compatible
Giesbrecht–Burns covariance is the inverse on active components. For variance
v or its declared total, df=2*v^2/Var(v); chi-square inversion gives the
Satterthwaite interval. Total variance includes all covariance terms in Var(v).
Zero-boundary individual intervals are unavailable; total inference conditions
on the active set. Infinite upper endpoints are JSON null, explicitly unbounded. An unrepresentable
lower endpoint instead makes the entire interval unavailable, with
`interval_status=unavailable_numeric_range`.
No random slopes, correlated effects, serial/heterogeneous residuals, residual
normality assessment, fixed-effect hypothesis tests or profile likelihood CIs.

Artifacts: raw input, resolved config, results, source-hashed interpretation facts,
manifest, offline HTML. Components include variance, SD, CV%, df and interval;
repeatability is residual, intermediate precision is all components combined.
Missing individual interval is unavailable; null upper endpoint within an existing
interval is unbounded. These two cases must not be conflated.

## 0.9.2 optional intervals

`inference.sum_interval` is `mls` (default from 0.9.3) or `satterthwaite`.
`mean_square_intervals_applicable` and `interval_rationale` are optional; declaring `false` with `mls`, or `true` with `satterthwaite`, is refused as contradictory.
Use the latter to justify independent Gaussian effects, the declared ordering
of random terms and the scope of the mean-square approximation. Without that
declaration the runtime refuses; it never chooses the interval after comparing widths.
`positive_sums` maps names to positive finite component weights, for example
`{"within_lot":{"lot:run":1,"residual":1}}`; it requires MLS and cannot replace
`total`. Unlisted components have weight zero. CV uses the same declared reference.

Orthogonality is checked against every covariance kernel. Exact orthogonal
strata use their integer df; otherwise df = 2 E(MS)^2 / Var(MS) and correlations
come from 2 tr(A_i V A_j V), with fitted REML V. The mean-square coefficients
are weights times the inverse expected-MS matrix. MLS/MOVER uses signed
lower/upper recovered errors with their correlation matrix. Output preserves
mean squares, mapping, df, correlations, interval center and old Satterthwaite
interval. Point estimates remain constrained REML. All component one-sided
upper bounds are computed, avoiding selective reporting only at zero estimates.
No simultaneous interval guarantee. The 0.9.2 nested-crossed lot upper-bound
miss must be stated when relevant. Unbalanced MOVER is an approximate application,
checked by 1000 fitted-model parametric bootstrap samples on each of two fixtures.
