---
name: elisa-quantification
description: Quantify unknown ELISA sample concentrations from plate-specific 4PL or predeclared 5PL standard curves, with standard back-calculation and independent-control QC, dilution correction, dilution-linearity checks, calibration-conditional unknown intervals and explicit interpolation-range checks. Imports plate-reader grids with a plate map. Use when standards with known concentrations and unknown wells are supplied; not for antibody binding EC50 or equilibrium KD.
---

# ELISA quantification

Read the [input and method contract](references/input-and-model.md) before analyzing a plate. Confirm analyte, standard matrix, response scale, standards and unknowns, dilution factors, blank/reference processing, and whether a monotone 4PL is appropriate. An ELISA binding titration without a known-concentration calibration standard belongs in dose-response or equilibrium binding, depending on its scientific question.

If the data are plate-reader grids (rows A–H/A–P, columns 1–12/1–24) plus a plate map, convert them with `agentic-prism import-plate --manifest MANIFEST --output NEW_DIRECTORY` (see the [plate import section](references/input-and-model.md#plate-reader-grids-070)); the importer copies values exactly, skips wells with a blank role and refuses non-numeric reads such as `OVRFLW`. Resolve saturated or missing reads in the source record yourself; never type a replacement number. Complete the generated `config.template.json` from the actual assay.

Create a JSON config and canonical CSV, then run the version-matched `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`. Locate and check the executable per the [runtime instructions](../agentic-prism/references/runtime.md). Run `agentic-prism verify --run DIRECTORY`, inspect `unknown_wells.csv`, `sample_summary.csv`, the standard-fit diagnostics and `report.html`. For appearance-only changes, `render --run DIRECTORY --style standard|prism_like` does not refit.

Check `standards_qc.csv` (or the report's standard-QC table) before anything else. Each positive standard level is back-calculated through the fitted curve. By default a level passes at 80–120% mean recovery (75–125% at the LLOQ and ULOQ) with replicate CV ≤ 20%. The plate is accepted only if at least 75% of levels, and at least six, pass (software screening rules, not full ICH M10 validation; set `qc` in the config to apply your SOP's limits). The legacy JSON keys LLOQ/ULOQ denote the lowest and highest passing levels on this plate; they do not establish assay-validated quantification limits. Failing standards are reported, never dropped automatically. If your SOP excludes a failing standard, mark it in the input with `exclude=true` and an `exclusion_reason`, and rerun into a new directory.

Report only wells marked `quantified`: their in-well concentration lies within LLOQ–ULOQ and the observed standard responses, and the value includes the declared dilution factor. Wells below the LLOQ or above the ULOQ are withheld with a reason; a failed calibration withholds every unknown on that plate. The C50 profile interval concerns the curve midpoint, **not** an unknown's uncertainty.

**Calibration model.** `fit.model` is `4pl` (default) or `5pl`. Choose 5PL only when asymmetry is expected from method development or prior plates, before seeing this plate's fit; the report adds an informational 5PL-vs-4PL F test and AICc (`asymmetry_not_supported_over_4pl` when the extra parameter is not supported) but never switches the model. 5PL needs ≥7 distinct positive standard levels.

**Unknown intervals.** By default (`uncertainty.unknown_interval=delta_method`) each quantified well and each complete sample summary gets a delta-method interval on log concentration that combines the plate's calibration-parameter covariance with response noise implied by the standard residuals. Quote it as a *calibration-conditional* interval: it excludes plate-to-plate, matrix, pipetting and dilution error, and it is not a validated assay precision. Simulated coverage is in the validation record.

**Dilution linearity.** When a sample has several dilution factors on one plate, `sample_summary.csv` and `dilution_linearity.csv` compare dilution-corrected means of in-range dilutions (recovery vs their mean, CV, and a log-log trend when ≥3 dilutions). `dilution_linearity.summary_policy`: `all_wells_required` (default; any withheld well withholds the sample) or `in_range_dilutions_with_linearity` (estimate from dilutions whose wells all quantify, only if linearity passes). `required=true` withholds samples whose linearity is not demonstrated. `corrected_concentration_increases_with_dilution_possible_hook_or_matrix_interference` is a pattern flag, not a diagnosis; report it and recommend further dilution or spike-recovery work. Never drop a dilution to make linearity pass.

Current scope: one 4PL or 5PL calibration per plate, raw OD or other finite response, unweighted or 1/Ŷ² relative fit, optional independently justified fixed plateaus, standard back-calculation QC, independent controls and descriptive cross-plate QC, analytical inverse, dilution correction, within-plate dilution linearity, calibration-conditional intervals and strict range checks. No automatic blank subtraction or re-fit after excluding standards, cross-plate parallelism, matrix spike-recovery model or mixed-effects intermediate precision. See [validation](../../validation/README.md) and [the 0.7.0 record](../../validation/RELEASE_0.7.0.md).

Independent controls use `role=qc`, positive nominal concentration before dilution,
and positive `dilution_factor` (1 for neat). They never enter the calibration fit.
For new assay workflows set `independent_qc.required=true`, document separately
prepared controls via `preparation_independent=true` and experiment-specific
`rationale`, and prespecify thresholds. Defaults require at least three positive
QC levels, two replicates each, recovery 80–120% and CV ≤20%; these are configurable
software screening rules, not a full ICH M10 implementation. Supplied failing or
unconfirmed controls withhold unknowns even if `required=false`.

Inspect `independent_qc.csv`, `cross_plate_qc.csv`, per-well `independent_qc_status`
and `quantification_reportable`. Historical configs without independent QC may
still return exploratory interpolation, explicitly marked `not_assessed`.
Three or more plates allow descriptive CV of plate means; this is not total
intermediate precision. Matrix effects, hook effect, stability, selectivity and
independent real-assay validation still need experimental evidence; the dilution
linearity check covers only the dilutions supplied on one plate.
