---
name: standard-curve
description: Fit linear or quadratic standard curves per plate (BCA, Bradford, absorbance or fluorescence assays, qPCR copy-number standards) with declared weighting and recovery acceptance, lack-of-fit and monotonicity checks, and interpolate unknown concentrations with inversion intervals and dilution correction, never extrapolating. Not for sigmoid immunoassay curves (use elisa-quantification).
---

# Linear and quadratic standard curves

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[contract](references/input-and-model.md). Implemented in 0.13.2.

## Before fitting

1. **Form** (linear or quadratic) and **weighting** (none, 1/x, 1/x²) come from the assay's validated range
   or SOP, recorded in `form_rationale`; never switch forms because one fits this plate better. Sigmoid
   immunoassays belong to elisa-quantification.
2. **Acceptance**: back-calculated recovery limits (general and lowest standard), the minimum passing
   fraction and their source are declared. A plate also fails on significant lack of fit (replicated
   standards) or a curve that is not monotone over the standard range; its unknowns are withheld.
3. **Dilutions**: one dilution per unknown `sample_id`; give each dilution its own ID.

## Run and report

Report the plate verdict and why, standard recoveries, and each unknown's concentration (after dilution)
with its inversion interval. Unknowns outside the standard range are reported as below/above range,
never extrapolated: suggest re-assay at another dilution. Intervals cover calibration and replicate noise
on this plate only, not between-plate or matrix effects. Weighted fits report point values only.

## Evidence (0.13.2)

Coefficients and inversion intervals agree with R lm and investr 1.4.2 `invest` within 5.1e-8 (root-finding
tolerance), including the replicate-pooled variance convention. Calibration results are in
[the 0.13.2 record](../../validation/RELEASE_0.13.2.md).
