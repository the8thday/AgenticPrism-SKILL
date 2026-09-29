# Advanced surface models (0.12.0)

Select a model in config; the analysis type remains `binding_kinetics`.

| Model | Meaning and limits |
|---|---|
| `one_to_one_drift` | Linear instrument drift per curve; sourced instrument rationale required. Does not repair an invalid reference. |
| `heterogeneous_ligand` | Two independent sites with ordered KD1 < KD2 and shared capacity/fraction. These describe surface heterogeneity, not two epitopes. |
| `bivalent_analyte` | First binding step includes factor 2 for free analyte arms. Second association uses response^-1 s^-1; state signal is X1+X2 and conserved sites are free+X1+2X2. Declare second rebinding during dissociation. No single KD. |
| `mass_transport` | Binding plus finite surface delivery. Declare conversion between response and surface-compartment nM, with source. Identified km alone does not establish transport limitation. |
| `off_rate_screening` | Per-curve apparent koff ranking from dissociation; no kon or KD. |

Multi-cycle, known concentrations, regenerated surfaces and unweighted response
are required. Global models require shared Rmax and a justified common response
scale. At least three concentrations and eight points per phase are required;
these minimum counts do not establish an informative design. Titrate and observe
long enough to resolve curvature and dissociation. Wells, curves and time points
are not independent experiments.

`advanced` declarations: `predeclared`, `mechanistic_rationale`, `reference_valid`,
`response_is_bound_analyte`, `response_scales_comparable`, `valency` (monovalent or
bivalent). Bivalent requires `second_rebinding`; drift requires `drift_source`;
transport requires `surface_conversion_nM_per_response` and its source.

Reliability gates: numerical boundaries, parameter correlation, profile endpoint
support (default one log unit each way), tolerance sensitivity, residual serial
correlation, irregular sampling, little dissociation, 75%/50% dissociation-window
sensitivity and >=90% accepted segment-wise block bootstrap refits. These
profiles check whether a declared span is supported; they are not full CIs.
AICc comparison uses identical windows. `comparison_rule=aicc_and_identifiability`
requires the prespecified improvement (default10); `diagnostic_only` does not
select a mechanism. A better AICc does not prove a biological mechanism.

Examples: `fixtures/kinetics_0120/*_config.json`. Real export originals,
licences and SHA256 are in its `public/` folder. Source software may fit different
capacity-sharing or baseline conventions; export compatibility does not establish
numerical equivalence. Read `validation/RELEASE_0.12.0.md` for observed gates.
