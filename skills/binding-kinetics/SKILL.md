---
name: binding-kinetics
description: Analyze BLI or SPR time-response sensorgrams with explicit phases - independent regenerated cycles or single-cycle (sequential injection) series - using a global 1:1 kon/koff model and kinetic KD, with explicit reference-channel or double referencing. Supports canonical CSV and a verified Octet Results.txt export layout.
---

# Binding kinetics

Use this specialist for time-resolved association and dissociation data, not endpoint titrations. Inspect the actual export and experimental design before choosing the model. Read the [input and model contract](references/input-and-model.md) for new data, and the [Octet adapter guide](references/octet-results-txt.md) for the supported text export.

1. Establish analyte concentration and units, phase boundaries, response units, reference/baseline processing, the injection design (`assay.injection_design`: `multi_cycle` with `independent_cycles=true`, or `single_cycle` with `independent_cycles=false` for sequential injections without regeneration), sample/experiment identity, and evidence for an ideal 1:1 interpretation. Choose `preprocessing.reference_mode` from what was actually measured: `reference_column` (reference channel/sensor, same cycle) or `double_reference` (also a matched buffer blank cycle on both channels). Do not construct a blank or reference trace that was not recorded. The `assay` flags and rationale require experiment-specific justification; do not copy them from a demonstration.
2. Locate and check the runtime per the [runtime instructions](../agentic-prism/references/runtime.md) (collection root two folders above this Skill after resolving symlinks; `<root>/.venv`; `doctor` once; ask before running `install.py`). Commands below use `agentic-prism` for that resolved executable.
3. For the supported Octet `Time1/Data1` text layout, create an explicit import manifest and run `agentic-prism import-octet --manifest ... --output ...`. Inspect `import_manifest.json` and `config.template.json`; complete the assay gate and fitting choices before analysis. For other exports, map measured values to the canonical CSV contract without claiming the instrument adapter supports them.
4. Run `agentic-prism analyze --config ... --output NEW_DIRECTORY`, then `agentic-prism verify --run ...`. Inspect `results.json`, `diagnostics.json`, `preprocessing_log.json`, source snapshots, residuals, and the offline HTML report. `render` changes appearance without refitting.
5. Report kon (M⁻¹ s⁻¹), koff (s⁻¹), and kinetic KD = koff/kon (M) only for reportable fits. State the global sharing, selected time windows, reference/baseline handling, interval method, fit diagnostics, and whether the model is biologically plausible. A fit is not proof of 1:1 binding. When reliability blockers are present, the implementation marks the fit limited and withholds reporting intervals. `audit_intervals` and raw point estimates remain for audit only; do not quote them as scientifically reportable results, and do not suggest an approximate value for a report. For a limited fit, the reportable statement is that no reliable kinetic estimate was obtained, with the reasons. Inspect `window_sensitivity.csv` and residuals. Prefer appropriately designed independent experiments for experimental reproducibility.

Current model: independent association/dissociation cycles or single-cycle series (each row carries its own injection's start, end and concentration), at least two distinct positive concentrations, one shared kon/koff per declared fit group, configurable shared or per-curve Rmax and fixed-zero or fitted offsets, unweighted least squares, and segment-wise block bootstrap. For single-cycle data the analytic 1:1 state is carried across injections and the window-sensitivity check shortens only the final dissociation. No mass-transport, avidity, heterogeneous-ligand or drift model, no native `.frd` or Biacore export parser: Biacore/other SPR data must be mapped to the canonical CSV by hand, keeping the original export, and that mapping has no instrument-specific validation. Do not infer phase boundaries from trace appearance. The public Octet validation is a scoped comparison with published TitrationAnalysis results, not a native Octet or Prism equivalence claim. See [validation](../../validation/README.md).

The router is [agentic-prism](../agentic-prism/SKILL.md). Use this specialist directly when the task is clearly kinetic; it shares the versioned Python implementation with the router and equilibrium skill.

Reliability checks include absolute phase-wise residual lag-1 correlation >0.5,
model-implied dissociation decay <2%, irregular sampling, and deterministic refits
using the first 75% and 50% of the declared dissociation window. More than twofold
rate/KD change, boundary hits or failed/insufficient sensitivity refits block reporting.
These are conservative engineering heuristics, not a validated classifier for mass
transport, avidity or model truth. Passing them does not establish nominal CI coverage.
See [0.6.0 evidence](../../validation/RELIABILITY_0.6.0.md).

Single-cycle (0.7.0): in one simulated iid design, 95% bootstrap intervals covered
kon, koff and KD in 92.7%, 91.3% and 94.0% of 150 fits (koff just below the
prespecified bound), similar to the slight undercoverage seen for multi-cycle data;
under strongly correlated noise every interval was withheld by the reliability gates.
Short intermediate dissociations can make the default block length infeasible;
choose `block_length` from the sampling design before fitting, not after seeing
intervals. See [0.7.0 evidence](../../validation/RELEASE_0.7.0.md).
