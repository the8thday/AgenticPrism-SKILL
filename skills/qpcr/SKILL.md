---
name: qpcr
description: Relative quantification of qPCR data - efficiency-corrected delta-delta Cq (normalised to the geometric mean of declared reference genes), technical replicates averaged within biological replicates, Welch or paired t on log2 relative quantities with fold-change intervals, plus NTC, technical-SD, undetermined-Cq and reference-stability checks. Use for gene-expression changes after treatment; not absolute copy number (use standard-curve).
---

# qPCR relative quantification

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[contract](references/input-and-model.md). Implemented in 0.13.2.

## Before analysis

1. **Biological replicates are the sample size.** Technical wells are averaged within each biological
   replicate; never compute statistics on wells. Declare pairing: independent (different donors or
   cultures per condition) or paired (the same source in every condition, with `pair_id`).
2. **Reference genes** need stability evidence for these conditions (`reference_rationale`, e.g. geNorm
   or NormFinder results). Two or more references are recommended (MIQE). The run reports how far the
   references shift between conditions and flags shifts above the declared cycles.
3. **Efficiencies**: declare per-gene amplification factors with their source (standard-curve slope) or
   state `assume_100_percent`. Declared values are treated as known.
4. **QC limits** (`cq_max`, technical SD flag) come from the lab SOP. Undetermined wells are never
   replaced by 40 or the cycle limit.

## Run and report

Report fold change versus the control with its interval and the method (Welch or paired), Holm p within
each target, the number of biological replicates, and every QC flag (NTC amplification, partial
undetermined replicates, technical SD, reference shift). Multiplicity across many genes is not adjusted.

## Evidence (0.13.2)

An independent base-R delta-delta Cq with `t.test` agrees within 1.9e-13 (the R pcr package is not
available). Published worked example (e.g. Livak and Schmittgen 2001): UNMET, full text not retrievable.
Calibration results are in [the 0.13.2 record](../../validation/RELEASE_0.13.2.md).
