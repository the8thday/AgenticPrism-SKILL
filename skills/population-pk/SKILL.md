---
name: population-pk
description: Population PK with R nlmixr2 (FOCEi) for two-compartment linear, Michaelis-Menten (nonlinear) or quasi-steady-state TMDD antibody models after IV bolus, IV infusion or SC dosing, with declared structure, random effects, residual model and fixed parameters, identifiability gates (covariance, correlation, saturation inside the observed range, RSE, shrinkage) and descriptive model comparison. Use for animal or human concentration-time data across dose levels; NCA belongs to pk-nca, PD response models are not covered.
---

# Population PK with nlmixr2

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Use `analysis_type: population_pk` (0.15.1). R nlmixr2 estimates every value
through a fixed packaged script with model templates; it needs R and the pinned
packages (`python3 install.py --with-r-pop`; doctor shows `optional_r_pop`). If they
are missing, say so; do not fit the model another way. Read the
[input and model contract](references/contract.md). Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **One population across dose levels.** All animals and doses are fitted
   together. Never fit each dose separately and compare to claim nonlinearity.
2. **The structural model is declared before fitting, with a reason**
   (`linear_2cmt`, `linear_mm_2cmt`, `qss_tmdd_2cmt`). Comparison models are
   optional, declared in advance and descriptive; never pick the model with the
   best AIC and report it as the finding.
3. **What the assay measures.** For QSS-TMDD declare whether concentrations are
   total or free drug. KSS and R0 are then in the drug's concentration unit as
   drug-binding capacity (1:1), not measured target levels.
4. **Identifiability.** PK data alone rarely identify every TMDD parameter. Fix
   parameters that the design cannot inform (often kdeg) at declared values with a
   source. SC-only data cannot identify F: fix it with a source, or set 1 and call
   the clearances CL/F and V/F.
5. **Starting values** for every estimated parameter, with a source (typical
   species values, NCA results). Declared random effects (diagonal), residual
   model (proportional or combined), seed. BLQ only as `drop` (M1) in this release.
6. **A nonlinear claim needs data at saturating and non-saturating
   concentrations.** One dose level, or a Km outside the observed range, cannot
   establish TMDD.

## Run and report

Report the declared model and why, the typical values with 95% intervals and
units, IIV (CV%) with shrinkage, residual error, objective and covariance status.
Withheld parameters are not results: say why (no covariance, |correlation| > 0.95,
saturation not inside the observed range, RSE > 100%). Individual (EBE)
parameters are withheld when shrinkage exceeds 30%. Comparisons: quote ΔAIC as
descriptive support only. Show the goodness-of-fit and individual-fit figures.
Respect `must_mention` and failing items in `interpretation_facts.json`.

Covariates, interoccasion variability, PD/effect models, full kon/koff TMDD, M3
BLQ, VPC and bootstrap are not part of this analysis type. For other structures,
write a saved R script following [the R workflow](references/r-workflow.md) and
label it as a custom analysis outside the validated runtime.

## Evidence

See [the 0.15.1 record](../../validation/RELEASE_0.15.1.md): independent IPRED
recheck (with its missed tolerance), the published Bolus_2CPTMM truth, the
registered simulation-re-estimation calibration (including fit identifiability
rates and the first attempt's defect) and the live-agent scenarios.
