# Combination and HTS contracts prepared for 0.13.0

Combination input: one complete rectangular dose matrix per declared experiment,
including zero-dose axes and at least four positive doses per agent. One condition
and one response scale; canonical percent inhibition with sourced normalization,
no implicit viability conversion. Technical wells must be summarized explicitly.
Repeated public rows without independent-experiment identifiers support point
comparisons only, not uncertainty or synergy claims.

Bliss uses independent fractional inhibition, HSA the larger observed single-agent
effect, Loewe equal-effect dose additivity over compatible response ranges, ZIP
conditional4PL fits along both directions minus fitted single-agent Bliss. Reuse
the existing separable4PL solver. Single-agent or conditional fits outside their
supported range withhold the affected model. No silent choice of favorable model.
Uncertainty uses equal-weight independent matrix scores and Student t intervals;
unreplicated matrices have descriptive scores only. All requested models and
incomplete states remain visible.

HTS input: full rectangular canonical plate with well/compound/unit identities,
explicit positive/negative/sample roles, >=4 distributed controls of each kind.
Declare independent wells, randomized layout and majority-inactive assumption.
Z-prime and plug-in SSMD describe control separation. B-score uses base-R-style
median polish and an explicitly declared MAD scale (default1.4826). Median polish
may be inappropriate for structured biological layouts or dense activity.
Predeclared QC failures withhold hit claims. One-sided predictive t tests against
negative controls, followed by BH within plate, are exploratory well-level calls;
plate correction changes dependence, so FDR calibration limits must be disclosed.
No hit is automatically a confirmed biological compound effect.

## Canonical columns and configuration

Combination CSV: `experiment_id`, `condition_id`, `dose_a`, `dose_b`, `response`.
Config `analysis_type=drug_combination` declares `agent_a`, `agent_b`,
`dose_unit_a`, `dose_unit_b` (M/mM/uM/nM/pM), one `condition_id`, `source`,
`readout_description` (what the detector/endpoint measures, or an explicit unknown),
`normalization_source`, `same_response_scale=true`, `independent_experiments`,
`independence_supported` and `independence_source`. `models` explicitly names
Bliss/Loewe/HSA/ZIP. Selecting Loewe requires
`loewe.compatible_maximal_effects=true` and a source; fitted plateau differences
must also be within its declared tolerance (default10 percentage points).
4PL fits are unweighted on the stated inhibition scale. Optional fixed bottom/top
must be declared; midpoint must lie in the positive tested range and boundary,
rank and flatness gates apply. Default mean-score interval level0.95.
See [synthetic config](../../../fixtures/pharmacology_0130/config_combination.json).
Do not copy its synthetic applicability declarations into measured data.

HTS CSV: `plate_id`, `well_id`, `compound_id`, `independent_unit_id`, `row`,
`column`, `role`, `value`. Row/column are positive integer grid coordinates.
Roles are `positive`, `negative`, `sample`. Config `analysis_type=hts_qc`
declares `source`, `readout`, `response_unit`, `direction`, independence,
randomization and majority-inactive flags with sources. Correction is
`median_polish`, with epsilon0.01 and max_iterations10 by default; MAD scale1.4826
requires a source. `criteria.predeclared=true` and `criteria.source` are mandatory;
default criteria Z-prime>=0.5, |SSMD|>=3 and FDR0.05 must be justified for the assay.
See [synthetic config](../../../fixtures/pharmacology_0130/config_hts.json).

For a grid import, select `target=hts_qc` in the existing `plate_grid_csv`
manifest. Provide `role`, `compound_id` and `independent_unit_id` layers; each
well's row/column and unmodified readout are copied automatically. The emitted
config is deliberately incomplete and needs actual assay declarations before
analysis. No skipped or unassigned well is imputed to complete an HTS plate.

HTS `plate_shape` must be declared as `[8,12]` or `[16,24]`. Every declared
coordinate must be present, including entire rows/columns; a smaller observed
rectangle does not silently redefine the plate edge. The import manifest records
source plate sizes; fill this required field in its configuration template.
