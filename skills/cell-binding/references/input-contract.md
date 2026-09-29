# Cell-binding contract

`analysis_type=cell_binding`, `model=cell_binding_apparent`. Defaults and strict
validation: `src/agentic_prism/cell_binding.py`. Unknown keys fail.

CSV columns: sample_id, experiment_id, curve_id, concentration,
concentration_unit (M/mM/uM/nM/pM), response, response_unit, cell_line, density_id,
series (`total` or `nonspecific`). Optional observation_id, exclude and
exclusion_reason. Each curve contains the paired total/control titrations from
one experiment; control identity/source is configured in `nonspecific`.
At least three control concentrations span the complete total-series range.
Excluded observations remain in the snapshot; no automatic exclusion or smoothing.

`cell_assay` declares cell_line, equilibrium_incubation_supported, incubation_s,
temperature_C, internalization_controlled, wash_protocol,
wash_dissociation_supported, detection (`direct_label`/`secondary_antibody`),
valency (`monovalent`/`bivalent_IgG`/`multivalent`) and rationale. Flags must be
literal true with evidence. Common assay flags require equilibrium, proportional
signal and a justified operational single-site response, total concentration,
interpretation `apparent_KD`, and rationale. `model_rationale` is required.

`receptors`: depletion (`negligible`/`quadratic`), cells_per_mL,
receptors_per_cell, source and density_id. Quadratic receptor concentration is
cells/mL × 1000 × receptors/cell / 6.02214076e23 M. This is active accessible
site concentration conditional on the supplied quantitation, not estimated from
binding MFI. Negligible mode requires free_approximation_supported=true and source.

`background`: mode `joint_baseline` (fit a common baseline using both series) or
`subtract_declared` (subtract explicit value from both series, then fit residual
baseline), value, source. Declared subtraction is conditioned on; uncertainty in
that external background is not propagated. Total = baseline + amplitude ×
occupancy + slope × total ligand concentration. Control = baseline + slope ×
concentration (+ control offset). Separate control affinity is not fitted.

`nonspecific.baseline` must be declared: `shared` (control measured on the same
cells, e.g. isotype or excess competitor) or `separate` (one extra control
intercept per curve; the nonspecific slope stays shared). Antigen-negative cells
can differ in autofluorescence and therefore require `separate`; `shared` is
refused for them. `separate` has an R numerical comparison but no coverage study.

Fit: fitted baseline, linear residuals, weighting explicitly `unweighted` or
`relative`. Relative minimizes sum((predicted-observed)/predicted)^2, as in
 dose-response; positive responses required after background handling.
95% profile-F reoptimizes all nuisance coefficients. Intervals are approximate,
not a calibrated guarantee for arbitrary noise, valency or assay design.

`replicates`: independent_unit `none` or `experiment_id`, conditions_comparable
boolean. Technical curves first averaged on log KD, then experiments equally
weighted; Student-t interval across experiments, none with one experiment.
Any missing/limited estimate withholds the sample summary.
