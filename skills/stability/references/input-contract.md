# Stability contract

CSV: `observation_id,batch,time,value`. IDs are nonempty and unique; preserve
batch labels (including `NA`). Each batch needs baseline 0, at least three
distinct times and four observations. Time is months and values are finite.
No exclusions are supported inside this module. Additional columns are kept
in the copied input; `exclude` may only be `false`.

Use the complete [synthetic config](../../../fixtures/stability/config_pooled.json)
and [declared extrapolation example](../../../fixtures/stability/config_extrapolated.json)
as schema examples, never as evidence for the user's applicability flags.

Required config:

- `analysis_type: stability`, `input`, `source`.
- `assay`: `attribute`, `unit`, `direction` (`decreasing`, `increasing`,
  `two_sided`), `direction_source`, `storage_condition`, `storage_source`,
  `storage_class` (`room`, `refrigerated`, `frozen`), literal true
  `linear_model_appropriate`, `independent_errors`, `common_variance`, `rationale`.
- `specification`: `source` plus `lower`, `upper`, or both according to direction.
- Optional `extrapolation`: `requested` defaults false. If requested, supply
  `source`; `accelerated_change` is `none_6_months`, `within_3_months`, or
  `between_3_6_months`; `intermediate_no_significant_change` is boolean.
  `supporting_data`, `change_pattern_continues`, `commitment_study` must all be
  explicitly true for any extension. Missing conditions cap the result at X.

Unknown config keys are rejected. Alpha 0.25 and confidence 0.95 are fixed by
this Q1E contract. Artifact `results.json` contains poolability, selected and
individual fits, bound crossings, decision-tree declarations, curvature,
failures and the capped primary result. `interpretation_facts.json` copies the
saved result and calibration disclosures before the manifest is written.

```sh
<collection>/.venv/bin/agentic-prism analyze --config config.json --output NEW_DIR
<collection>/.venv/bin/agentic-prism verify --run NEW_DIR
```
