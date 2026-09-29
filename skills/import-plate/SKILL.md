---
name: import-plate
description: Convert explicitly mapped plate grids into canonical data and analyze HTS control separation, plate effects, median-polish B scores and exploratory FDR hits with declared control roles, layout and independence. Keeps raw files and refuses hit calling without justified correction and passing plate QC.
---

# Plate import and HTS QC

Read [input contract](references/input-and-model.md) and shared
[runtime](../agentic-prism/references/runtime.md). Keep raw inputs and explicit
well maps. Import establishes a data layout, not assay applicability; fill the
emitted configuration template before analysis.

For `analysis_type=hts_qc`, declare readout, direction, independent wells,
randomized layout, a majority of inactive wells and the sources for QC/FDR
criteria. Repeated measurements of one unit are not independent wells. Require
complete rectangular data and dispersed positive/negative controls; never
impute missing wells to complete a plate.

Z-prime and SSMD describe controls. Median polish estimates row/column effects;
B scores scale residuals by declared MAD. Biological layouts confounded with
position invalidate this correction. Inspect row/column and edge diagnostics.
Failed QC, nonconvergence or zero residual variation withhold hit tests.

BH-adjusted one-sided predictive t tests are exploratory, assumption-dependent
well screens. They do not establish validated compound activity. Carry observed
calibration misses into reports and confirm hits in independent experiments.
Use `interpretation_facts.json` and [release evidence](../../validation/RELEASE_0.13.0.md).

0.13.0 calibration observed a global-null false-hit probability of0.100
(100/1000; MCSE0.009487), above the registered0.063784 bound. Lead with this
limitation: nominal BH adjustment has not demonstrated FDR control here.
The saved exploratory hit indicators must not be called validated FDR-controlled
hits. No threshold was retuned after these results.
