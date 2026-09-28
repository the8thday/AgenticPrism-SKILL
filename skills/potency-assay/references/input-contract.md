# Potency-assay contract

Declared estimates CSV:
`observation_id,run_id,nominal_rp,relative_potency,reference_quality,parallelism,fit_source`.
RP and nominal RP are finite and positive. Statuses are `pass`, `fail`, or
`not_evaluable`; `fit_source` identifies the actual source record. No exclusions.
Each run includes every nominal level with at least two independent RP estimates.
`combine` uses one level; `validate` requires at least three levels.

For `input_mode=dose_runs`, use
`observation_id,run_id,nominal_rp,dose_run,comparison_id`. `dose_run` is a path
relative to the config to a verified dose-response output directory using
prespecified equivalence parallelism. The importer reads the comparison's RP,
reportability, reference fit and equivalence status; copies hashed source
artifacts under `upstream/`; rejects repeated comparisons/shared reference fits;
and never substitutes an RP for a failed fit. Reruns use the portable copied
paths. Use dose-response to create each independent curve comparison first.

See [combination example](../../../fixtures/potency_assay/config_combine.json),
[validation example](../../../fixtures/potency_assay/config_validate.json), and
[prespecified screening example](../../../fixtures/potency_assay/config_outlier.json).
Synthetic declarations do not establish a real assay's applicability.

Required config: `analysis_type: potency_assay`, `input`, `source`, `mode`
(`combine` default or `validate`), `input_mode` (`declared_estimates` default or
`dose_runs`), and `assay` containing literal true `independent_runs`,
`independent_determinations`, `common_log_variance`, `same_material_per_level`,
plus `rationale`. `suitability.source` is mandatory.

For validation, `criteria.source`, `max_abs_log_bias`, `max_gcv_percent` are
mandatory. Optional `linearity_equivalence: true` also requires
`margins_prespecified: true`, `slope_limits: [lo,hi]`, and
`intercept_limits: [lo,hi]`. Criteria use natural-log RP, not log10.
`statistics.confidence_level` defaults .95; use .90 for TOST alpha .05 per side.
The software does not select a confidence level to improve acceptance.

`outliers.method` defaults `none`. `grubbs_run_means` requires declared `alpha`,
`prespecified: true`, and `source`, equal run sizes and `mode=combine`.
There is no exclusion option. Unknown settings are rejected.

Saved `results.json`, `interpretation_facts.json`, config, input, source hashes
and offline HTML are produced. `render` cannot refit or alter scientific files.
