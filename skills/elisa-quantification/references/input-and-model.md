# ELISA input and method contract

CSV columns: `plate_id,well_id,role,sample_id,concentration,concentration_unit,dilution_factor,response,response_unit`. `role` is `standard` or `unknown`. Standards need a nonnegative physical concentration and a blank `dilution_factor`; unknowns need a blank `concentration` and a positive `dilution_factor` (1 for neat samples). Optional `exclude=true|false` and `exclusion_reason` retain predeclared exclusions. Well IDs are unique within a plate. Each plate has one concentration unit and one response unit. Supported concentration units are those of the 4PL specialist; no automatic conversion between molar and mass units. Enter measured response values without silent blank subtraction or normalization.

The JSON config declares `analysis_type=elisa_quantification`, the CSV `input`, `direction=increasing|decreasing`, and nonempty `assay.analyte`, `response_definition`, `matrix` and `rationale`. Optional `fit.weighting=unweighted|relative`, `fixed_bottom`, `fixed_top`; fixed plateaus need independent control evidence. Optional `qc` (defaults shown) sets the standard acceptance criteria:

```json
"qc": {"recovery_limits_percent": [80, 120], "edge_recovery_limits_percent": [75, 125],
       "max_cv_percent": 20, "min_passing_fraction": 0.75, "min_passing_levels": 6}
```

`max_cv_percent: null` disables the replicate CV criterion. The edge limits must contain the interior limits. See [the synthetic example](../../../fixtures/elisa_synthetic/config.json).

The standard curve uses the shared symmetric 4PL fitter. It needs at least six distinct positive standard concentrations and ten included standard observations, plus a reportable in-range fit. `calibration_C50` is a model midpoint, not antibody binding EC50 or KD. For increasing response, an unknown's fitted-curve fraction is `(Y−Bottom)/(Top−Bottom)` and its log concentration is `log10(C50)+logit(fraction)/(ln(10)·Hill)`; decreasing response uses the signed negative Hill. The inverse is multiplied by the supplied dilution factor. **Standard QC.** Every included positive standard well is back-calculated with the same inverse. Per level: mean back-calculated concentration, recovery = mean/nominal × 100, and CV% of the back-calculated replicates. A well whose response lies outside the fitted plateaus cannot be back-calculated, so its level fails. The LLOQ is the lowest and the ULOQ the highest level passing the edge limits (recovery and CV); levels between them must pass the interior limits; levels outside them count as failing. The plate is accepted when LLOQ < ULOQ and the number of passing levels is at least `max(min_passing_levels, ceil(min_passing_fraction × number of levels))`. Zero-concentration standards (blanks) are fitted but not back-calculated. Failing levels stay in the fit; excluding one requires an explicit, reasoned exclusion and a new run.

An unknown is withheld when the calibration failed QC, the standard fit is not reportable, its response falls outside observed standard responses, the fitted fraction is outside (0,1), or its in-well concentration is below the LLOQ or above the ULOQ. Withheld wells carry the reason (`calibration_qc_failed`, `concentration_below_lloq`, `concentration_above_uloq`, …). No numeric limit-of-detection or upper/lower statistical bound is manufactured.

## Independent QC extension (0.6.0)

`role` also accepts `qc`. `concentration` is its positive nominal concentration
in the original QC material; `dilution_factor` is required and must be positive.
Use separate QC preparation from calibration standards. `sample_id` records its
identity. Controls are grouped by nominal concentration within each plate;
replicates are technical, not independent experiments. Controls outside the
accepted in-well standard range fail, and exclusions remain recorded.

```json
"independent_qc": {
  "required": true,
  "preparation_independent": true,
  "rationale": "Describe the actual independently prepared QC material",
  "recovery_limits_percent": [80, 120],
  "max_cv_percent": 20,
  "min_levels": 3,
  "min_replicates": 2
}
```

The example is a JSON config fragment. See `fixtures/elisa_independent_qc` for a
complete synthetic example. Standard-derived `lloq_input_unit`/`uloq_input_unit`
are per-plate screening bounds, not experimentally validated assay limits.

## 0.7.0 extensions

### Calibration model, unknown intervals and dilution linearity

```json
"fit": {"model": "5pl", "weighting": "relative"},
"uncertainty": {"unknown_interval": "delta_method", "level": 0.95},
"dilution_linearity": {"required": true, "summary_policy": "in_range_dilutions_with_linearity",
                       "recovery_limits_percent": [80, 120], "max_cv_percent": 20}
```

Defaults: `model=4pl`, `unknown_interval=delta_method`, `required=false`, `summary_policy=all_wells_required`, 80–120% and 20%. See the [5PL dilution example](../../../fixtures/elisa_dilution_5pl/config.json).

**5PL.** `Y = Bottom + (Top−Bottom)·[expit(ln10·s·h·(log10 x − log10 C))]^g`; g=1 is the 4PL. Only log10 C, log10 h and log10 g (bounded 0.05–20) are optimized; plateaus come from the same bounded linear solve. `c_parameter_canonical` is C; `half_response_input_unit` is the concentration at the mid-plateau response. Inverse: `x = C·10^(logit(f^(1/g))/(s·h))` with `f=(Y−Bottom)/(Top−Bottom)`. No profile interval is computed for the 5PL midpoint (`ci_status=not_computed_for_5pl`). `model_comparison` reports the extra-sum-of-squares F test of g=1 and AICc; it is informational.

**Intervals.** Free parameters θ (plateaus, log10 C, log10 h[, log10 g]); `Cov(θ) = s²(JᵀJ)⁻¹` from the weighted residual Jacobian at the fitted solution with `s² = SSE/df`. For a sample, the estimate is the mean over its wells of dilution × inverse(Y). Its variance is `∇θᵀ Cov ∇θ + Σ (∂/∂Yᵢ)² σᵢ²`, where σᵢ² = s² (unweighted) or s²·Yᵢ² (relative). The interval is `log10(estimate) ± t(df)·SE_log10`. Withheld wells and incomplete samples get no interval (`ci_status` states why). The interval treats standard-curve error as the only source of uncertainty.

**Dilution linearity.** Per sample and plate, a dilution is *in range* when every well at that dilution is quantified. With ≥2 in-range dilutions: each dilution mean's recovery relative to the equal-weight mean of dilution means, their CV, and with ≥3 dilutions a regression of log10(mean) on log10(dilution factor). Pass = all recoveries inside the limits and CV ≤ limit. A significant positive slope, or a failed series whose means rise monotonically with dilution, is flagged as a possible hook/matrix pattern; a falling pattern is flagged separately.

### Plate-reader grids (0.7.0)

`agentic-prism import-plate --manifest plate_manifest.json --output NEW_DIRECTORY` turns grid CSVs into this contract. A grid has a header row (free first cell, then `1..12` or `1..24`) and rows `A..H` or `A..P`. See the [example manifest](../../../fixtures/plate_import_example/plate_manifest.json):

```json
{"format": "plate_grid_csv", "target": "elisa_quantification", "source": "...",
 "constants": {"concentration_unit": "ng/mL", "response_unit": "OD450"},
 "plates": [{"plate_id": "plate-1", "readout": "reader_od450.csv",
             "layers": {"role": "layout_role.csv", "sample_id": "layout_sample.csv",
                        "concentration": "layout_concentration.csv", "dilution_factor": "layout_dilution.csv"},
             "constants": {}}]}
```

Every output column comes from exactly one layer grid, plate constant or global constant; duplicates or gaps fail. The `role` layer decides which wells are used; wells with a blank role are skipped and counted, including wells that have a reading. Readout text is copied unchanged; non-numeric reads fail. Source files are hashed and snapshotted, and `config.template.json` leaves every assay fact unset. The importer only understands this generic grid layout, not vendor-specific exports (SoftMax Pro, Gen5 and so on); reshape those into grids without changing values, and keep the original export.
