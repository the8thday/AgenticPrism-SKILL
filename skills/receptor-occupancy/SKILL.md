---
name: receptor-occupancy
description: Calculate percent receptor occupancy (RO) from flow-cytometry MFI - free receptor against the subject's predose sample, free over total receptor, or bound drug over an ex vivo saturated aliquot - with declared background, no clipping to 0-100%, subject-level summaries by group and timepoint, and an optional exposure-occupancy EC50 with a cluster-robust (subject-level) small-sample interval. Use for target-engagement PD biomarkers of antibodies; not for in vitro binding KD (equilibrium-binding, cell-binding) or PK models.
---

# Receptor occupancy

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: receptor_occupancy`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **The format comes from the assay design, before the data.** Ask which reagents
   were used and what evidence shows their behaviour:
   - `free_baseline`: a reagent that competes with the drug measures unoccupied
     receptor; needs a predose sample per subject; assumes receptor levels per cell
     do not change.
   - `free_over_total`: adds a non-competing total-receptor reagent; use it when the
     drug may internalize or down-regulate the receptor, which free/baseline would
     count as occupancy.
   - `bound_saturation`: an anti-drug reagent measures bound drug, relative to the
     same sample spiked ex vivo with saturating drug.
   Do not compute several formats and keep the one that looks best.
2. **Background.** Each signal needs its declared background (isotype, FMO or
   unstained) MFI from the same sample and run. Gating and compensation happen before
   this analysis; the table holds population MFIs.
3. **Do not clip.** Values below 0% or above 100% reflect assay noise and are kept;
   clipping biases group means near 0% or 100%.
4. **Units.** The subject is the unit; one sample per subject and timepoint
   (average technical tubes first by a declared rule). Do not pool timepoints as if
   they were independent subjects.
5. **Exposure-occupancy (optional).** Needs drug concentrations from the same
   samples, at least 6 subjects, and data on both sides of 50% occupancy. The EC50 is
   an in vivo apparent value for that matrix and population, not an intrinsic KD,
   and must not be extrapolated beyond the observed concentrations. Set
   `"interval": "cluster_robust_t"` (0.14.3): it covered 0.95-0.96 with 8 subjects
   and 0.948 with 6 in calibration. The default `cluster_bootstrap` remains only so
   older configurations reproduce; it was too narrow with 6-8 subjects (0.86-0.89),
   so if a run used it, call that interval approximate and optimistic.

## Run and report

Report the format and reagent evidence, per-group and per-timepoint mean RO with
intervals and the number of subjects, values outside 0-100%, samples that could not
be computed (signal not above background), and, if requested, the EC50 with its
interval (state which method and its degrees of freedom) and the observed concentration range. Respect
`must_mention` and failing items in `interpretation_facts.json`.

## Evidence

See [the 0.14.2 record](../../validation/RELEASE_0.14.2.md): base-R oracle, the
registered calibration (including the internalization and clipping comparisons) and
the unmet published-example gate. The cluster-robust EC50 interval: [the 0.14.3
record](../../validation/RELEASE_0.14.3.md) (R clubSandwich oracle, registered
calibration after a disclosed pilot).
