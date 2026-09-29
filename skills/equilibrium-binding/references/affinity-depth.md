# Opt-in affinity depth (0.11.0)

All models require model_rationale before fitting. Defaults and validation live
in `affinity_schema.py`. Retain old `one_site_with_baseline` for justified free
concentrations. Total ligand with material depletion cannot be forced through
that model without a supported free approximation; use mass balance instead.

New models use total concentrations, linear residuals, fitted baselines,
unweighted or inverse_sd weighting, and profile_f intervals. Unknown settings
fail. Source concentration/response columns are the equilibrium input contract.
No model selection by improved residual SSE alone.

## Direct depletion

`model=one_site_depletion`; `active_sites` includes mode `declared` or `fitted`,
concentration_M, basis `active_sites`, source, activity_basis, bounds_M.
A nominal protein concentration without an activity basis is refused. In fitted
mode the concentration is an initial declaration only, not fixed truth.
The exact quadratic complex / Pt determines occupancy. Pt/KD is always saved:
operational binding ≤0.1, intermediate >0.1 to <10, titration ≥10, following the
regime distinction of Jarmoskaite et al. 2020 (doi:10.7554/eLife.57264).
These cutoffs are design diagnostics, not universal calibration guarantees.
When Pt is fitted, profile KD while reoptimizing Pt and all signal coefficients.
A titration regime or lower-open profile withholds the KD point. Only a closed
upper profile endpoint explicitly saved as supported_upper_bound_M can be cited;
otherwise no numerical affinity bound is supported. Numerical search boundaries
are never confidence limits. Audit estimates remain available but unreportable.

`design_warnings` (non-blocking, all mass-balance models except competition):
`titrant_max_below_twice_constant_species` when the highest titrant concentration
is below twice Pt (or the SET constant species). The equivalence region was then
not spanned and KD rests on sub-saturation curvature; the 0.11.0 known-Pt
Pt/KD=100 calibration miss came from exactly this design. Mention it and
recommend extending the titration; it does not change reportability.

## Solution equilibrium titration (SET/KinExA style)

`model=solution_equilibrium_titration`; declared active sites and at least two
curves with different constant_species_M. CSV also requires fit_group_id and
constant_species_source per curve. Curves in a group share one KD, with separate
positive signal amplitudes and baselines; one sample, independent experiment and
response unit per group. Constant species are independently sourced active
concentrations, not fitted nominal total protein. Model signal is proportional
to free constant-species fraction. No pooling across independent experiments.

`valency` (SET only; other models refuse it): constant_species `monovalent` or
`bivalent`; detected_signal `free_sites` or `molecules_with_any_free_site`;
constant_concentration_basis `binding_sites` or `molecules`;
titrant_concentration_basis `binding_sites` (e.g. 2 × IgG) or
`monovalent_molecules`; source. A bivalent constant species (IgG) whose capture
step detects every molecule with at least one free site is refused: that signal
is not proportional to the free-site fraction and no bivalent SET model is
implemented. Use a Fab/monovalent format, or a readout justified as free sites
with site concentrations; the latter is fitted as independent sites and the
assumption is added to must_mention. Ask the user which species the plate or
bead captures; do not choose `free_sites` to make the run proceed.

`equilibration`: incubation_s, method `koff` or `measured_time_course`,
koff_per_s or time_to_95_s, and source. Incubation must reach conservative
-ln(0.05)/koff or the experimentally established time_to_95_s. The bound from
koff assumes a reversible 1:1 process and appropriate initial conditions; it
cannot establish equilibration in aggregation/internalization models.

## Competition

`model=competition_exact`: mutually exclusive R, RT and RI states, no ternary
complex. Active receptor Pt must be declared and sourced. `competition` requires
tracer_total_M, tracer_kd_M, source. Primary parameter is **Ki**, not direct KD.
Four-state/allosteric competition is not implemented or claimed benchmarked.
Signal is proportional to bound tracer fraction; fluorescence anisotropy with
binding-dependent quantum yield needs prior justified correction, not this raw
linear signal model. Baseline/amplitude are fitted, tracer KD and concentration
uncertainty are conditioned on, not propagated.

Optional cheng_prusoff=true requires an externally supplied ic50_M. This is a
labelled diagnostic approximation, never the default model. Refuse it when
baseline tracer-bound fraction >5% or receptor_total/ic50 >5%; these are
conservative applicability checks, not proof of negligible all assay effects.

## Evidence and interpretation

Read [release evidence](../../../validation/RELEASE_0.11.0.md). Published
worked-example gates are UNMET for these new models. Do not call synthetic
fixtures a published example or claim native KinExA/Prism equivalence.
Use facts for reportability, upper bounds and calibration misses. Report the
fitted range, weighting, profile-F convention and source-dependent limitations.

Observed registered calibration misses: known Pt profile coverage 93.1% at
Pt/KD=1 and 80.0% at 100 (bound 93.62%). Fitted-Pt titration/flat-profile
diagnostic withholding 48.5% at ratio 10 and 91.5% at 100, below the same bound.
At ratio 100 any point withholding was 95.7%; that broader count does not erase
the original diagnostic-specific miss. The exact ratio=10 decision boundary is
unstable under noise. Do not promise reliable titration classification. SET
coverage 95.2% and three-state Ki coverage 95.8% passed in their registered designs.
