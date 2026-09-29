# Contract (0.11.1)

CSV: observation_id, independent_unit_id, x, y. All pairs complete and finite;
no repeated subject IDs. Linear inverse-variance WLS also requires `weight>0`.
Missing pairs must be removed upstream under a documented rule with originals
retained. Examples: `fixtures/routine_0111/config_*.json`.

Top-level keys: schema_version=1, analysis_type (`correlation_regression` or
`method_comparison`), input, source, method, confidence_level, assay, regression,
comparison, report. Defaults for other optional sections are in `routine.DEFAULTS`;
unknown keys fail. Assay: independent_units=true, predeclared=true,
unit_definition and rationale. Method comparison additionally requires
same_scale=true; Deming/Passing-Bablok require linear_relation_supported=true.

Correlation methods: pearson (Fisher-z interval), spearman (exact untied n<=9,
otherwise R exact=FALSE t approximation with average tied ranks), linear
(slope/intercept t intervals, mean and individual prediction bands).
`regression.weighting`: unweighted or inverse_variance; latter requires a
nonempty weight_basis and is only allowed for linear regression.

Method comparison: deming, passing_bablok, bland_altman.
`comparison.error_variance_ratio`: Var(error Y)/Var(error X)>0, ratio_source
required for Deming. This is the reciprocal of mcr's error.ratio convention.
Jackknife uses t(n-2). Passing-Bablok uses mcr's tangent median/rank interval,
relative ties and positive-association convention; CUSUM 5% threshold 1.36.
Bland-Altman reports 95% population limits, each with a confidence_level exact
noncentral-t confidence interval, not a simultaneous interval or a tolerance
interval. Acceptance_limits=[lower,upper] and acceptance_source are optional
user declarations; without them no acceptance verdict is available.

The reported curve range is the observed range. No extrapolated agreement,
latent true-value inference from OLS, or automatic outlier exclusion.
