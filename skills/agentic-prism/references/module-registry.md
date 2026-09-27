# Module registry — development 0.7.1

| Purpose | Skill | Capability state |
|---|---|---|
| Single-site equilibrium binding, KD and diagnostics | [equilibrium-binding](../../equilibrium-binding/SKILL.md) | Implemented; see release validation record |
| Re-render an existing equilibrium run | Same specialist, render command | Implemented; no refitting |
| BLI/SPR 1:1 kon/koff and kinetic KD | [binding-kinetics](../../binding-kinetics/SKILL.md) | Implemented for independent cycles and (0.7.0) single-cycle series in canonical CSV, explicit reference-channel or double referencing, and the scoped Octet Results.txt layout; reliability gates with audit-only intervals for blocked fits. No Biacore/native parser, mass-transport or heterogeneous models |
| ELISA 4PL/5PL calibration and unknown concentration interpolation | [elisa-quantification](../../elisa-quantification/SKILL.md) | Implemented for per-plate standards with back-calculation QC, independent-control gates, descriptive cross-plate QC, range-checked inverse and dilution, (0.7.0) predeclared 5PL, calibration-conditional delta-method unknown intervals and within-plate dilution linearity. No mixed-model intermediate precision or matrix model |
| Relative EC50/IC50 from symmetric 4PL curves; replicate summaries; reference-vs-test relative potency | [dose-response](../../dose-response/SKILL.md) | Implemented for long-form CSV with explicit endpoint/direction/units; optional fixed plateaus and relative weighting; parallel-line RP with F-test or (0.7.0) predeclared-margin equivalence parallelism and optional RP acceptance limits |
| Two-group independent or paired comparisons | [group-comparison](../../group-comparison/SKILL.md) | Implemented: Welch or paired two-sided t test on independent-unit data |
| One-way comparison of ≥3 independent groups | [group-comparison](../../group-comparison/SKILL.md) | Implemented (0.7.0): classic or Welch ANOVA with predeclared Dunnett, Tukey-Kramer, Games-Howell or Holm-Welch family. No repeated-measures, mixed-effects, covariate or nonparametric model |
| Install / locate / check the runtime | `install.py`, `agentic-prism doctor` ([runtime instructions](runtime.md)) | Implemented (0.7.1): per-clone `.venv` with lock-constrained dependencies and self-check; macOS evidence only |
| Plate-reader grid + plate map to long table | `agentic-prism import-plate` (used by the ELISA and dose-response Skills) | Implemented (0.7.0) for a generic 96/384-well grid CSV layout; no vendor-specific export parser |

Published validation evidence: [validation/README.md](../../../validation/README.md).
“Implemented” is not Prism numerical equivalence. No Prism project benchmark has
been run. New modules must supply their own data contract, estimand, assumptions,
tests and release status. A requested mixture of available and unavailable tasks
may complete the available work while explicitly identifying the remaining gap.

Route sensorgrams to the kinetic specialist; its Skill records current limits. Route
functional EC50/IC50 curves to the 4PL specialist, not the equilibrium KD fitter.

Routing examples:

- “Here is a titration table, calculate KD”: specialist checks equilibrium,
  concentration basis, signal meaning and valency before choosing the model.
- “Here is an ELISA plate, infer unknown concentrations”: calibration specialist;
  do not send its OD values through the KD workflow.
- “Compare three independent KD experiments”: equilibrium specialist summarizes
  log KD descriptively. For a predeclared two-group contrast, use the group
  specialist on one independent-unit value per group/experiment.
- “Compare four candidate antibodies with the isotype control”: group specialist,
  multi-group design with a vs-control family chosen before testing.
- “Same donors under three treatments”: repeated measures, not available; do not
  run a one-way ANOVA or unadjusted paired t tests instead.
- “Compute kon/koff from an endpoint table”: insufficient time-course data;
  request sensorgrams rather than manufacturing rate constants.
- “Single-cycle kinetics from Biacore”: kinetic specialist with
  `injection_design=single_cycle` after mapping the export to the canonical CSV.
- “Estimate IC50 from viability versus inhibitor concentration”: dose-response
  specialist; determine whether relative plateau midpoint is the intended IC50.
- “Is our candidate more potent than the reference antibody?” / relative potency:
  dose-response comparisons, pairing reference and test curves from the same
  plate, one comparison per independent experiment. Equivalence margins and RP
  acceptance limits must come from the user's historical data or protocol.
- “Make existing plots Prism-like”: render-only command on a saved run.

Current reliability scope: [0.6.0 evidence](../../../validation/RELIABILITY_0.6.0.md)
and [0.7.0 additions](../../../validation/RELEASE_0.7.0.md).
