---
name: population-pk
description: Guide population pharmacokinetic modeling in R with nlmixr2 and rxode2, including dosing-event preparation, basic compartmental models, between-subject variability and fit diagnostics. Use for population PK or PopPK requests, not NCA or pharmacodynamic response modeling.
---

# Population PK with nlmixr2

Status: local guidance prototype, 2026-10-13. Use **nlmixr2** (not `nlmix2`) for
estimation and its rxode2 infrastructure for model/event handling. Do not write
a substitute NLME estimator or repurpose independent-unit regression. This is
not an `agentic-prism run` analysis type or a validated production workflow.

Preserve the study's [discussion/direct-analysis choice](../agentic-prism/references/analysis-entry.md).
Read the [data/model contract](references/contract.md) before fitting and the
[R workflow](references/r-workflow.md) for an educational starting point.

## Scope and design

First-version guidance is for one analyte, linear one- or two-compartment PK,
declared IV bolus/infusion or first-order extravascular input, a justified small
random-effects structure and additive/proportional/combined residual error.
Only a one-compartment oral example is supplied; other structures/routes need
their own verified event mapping and example before they are claimed tested.

The estimands are population typical parameters, between-subject variability,
residual variability and explicitly requested predictions. Multiple records
from one subject are not independent subjects. Sparse sampling may be useful
across a population, but does not guarantee identification of individual Ka,
multiple compartments or every random effect.

PD/Emax, effect compartments, indirect response, TMDD, nonlinear elimination,
automatic covariate searches, interoccasion variability and censored-likelihood
BLQ models are outside this initial guide. A later extension can use mature R
packages; do not force a complex antibody PK pattern into linear PK to obtain
an answer. For NCA use [pk-nca](../pk-nca/SKILL.md).

## Workflow

1. Inspect source data and dosing history. Identify subjects, analyte/matrix,
   amount/concentration/time units, covariates and BLQ. Construct a documented
   event table and inspect per-subject timelines before fitting.
2. State the structural, variability and residual models, parameter meanings,
   initial estimates, estimation method and rationale in an analysis plan.
   Keep the model fixed for the initial fit; comparisons are separate saved
   runs with explicit scientific reasons, not automated model shopping.
3. Check R, packages and compiler support. Save model code, data mapping,
   seeds, controls, package versions and provenance. Use installed-version help;
   do not guess slot names or silently change algorithms after an error.
4. Fit with nlmixr2 in a fresh output directory. Retain warnings, failed attempts,
   model/input and fit RDS. A successful call does not establish convergence,
   identifiability or usable uncertainty.
5. Inspect optimizer/SAEM behavior, covariance/SE availability, parameter
   plausibility and correlations, random-effect boundaries, shrinkage and
   subject-level fits. Review observed versus PRED/IPRED and residuals versus
   time/prediction; CWRES must actually be computed if it is reported.
6. For predictive claims, specify and run an appropriate simulation/VPC design
   including sampling/dosing, strata, binning, replicates and seed. Record any
   unperformed VPC/bootstrap as not run. Confidence in a typical parameter is
   different from uncertainty for a new individual's concentration.
7. Extract reportability and source-linked facts from saved outputs before
   writing the report. Preserve unsupported parameters and failed diagnostics;
   do not present conditional individual estimates as precise measurements.

## Interpretation

For extravascular-only data with F fixed to one, clearance and volume estimates
are apparent CL/F and V/F. Label population typical values, individual empirical
Bayes estimates, variance parameters and residual SD separately. High shrinkage
limits individual/covariate interpretations. A VPC is a diagnostic, not proof
that the model is true; package use does not establish NONMEM equivalence.

See the [prototype evidence record](../../validation/PK_R_SKILLS_PROTOTYPE_2026-10-13.md).
Unmet evidence gates and missing diagnostics must remain visible in the report.
