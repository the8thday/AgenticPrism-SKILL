# Population PK data, model and output contract

This is an agent guidance contract, not an automatic validator. Keep the exact
source-to-analysis mapping and do not execute model code taken from data cells.

## Input and events

| Item | Declaration/check |
|---|---|
| Study | Species/population, analyte/free-total definition, matrix, route/formulation, concentration/amount/time units, subject and period/occasion definitions |
| Subject | Stable `ID`; preserve repeated samples and dose history within subject; covariate units, timing and missingness |
| Time | `TIME` on one continuous documented clock per subject; keep actual dose/sample times and same-time ordering; resetting time per dose must not erase history |
| Observations | `DV`, observation `EVID`, source row, missing/BLQ/LLOQ fields; do not use missing/BLQ as measured zero |
| Doses | `AMT`, dosing `EVID`, model-specific `CMT`; route, infusion rate/duration where relevant; explicit records or documented `ADDL`/`II` expansion |
| Regimen | Relevant previous doses, actual versus prescribed administration, missed doses, and explicit steady-state initialization if used |

Use the installed rxode2/nlmixr2 event convention. Standard NONMEM-style
observation `EVID=0` and dose `EVID=1` are not interchangeable with legacy
rxode2 codes (the upstream theophylline example uses codes such as 101).
Do not change legacy codes by a blanket replacement. Compartment numbers depend
on the model: oral depot, IV central input and observation compartments differ.
Infusion rate and duration must agree with amount and units. Dose rows are not
concentration observations, even when a source uses `DV=0` as their placeholder.
Retain documented same-time dose/sample ordering; do not just sort on TIME and
assume a predose sample became postdose.

Prefer explicit dose records for the first adaptation. For complex event codes,
infusions, steady state or expanded regimens, verify the event trajectory against
the source schedule before fitting. Missing dosing history blocks the affected
model; covariate values or records cannot be invented to get a fit.

If BLQ exists, stop the basic-example path and establish a justified treatment
and sensitivity plan. Never silently drop BLQ or use LLOQ/2. A censored-likelihood
extension needs version-specific coding and its own verification; this prototype
does not demonstrate it.

## Model plan and diagnostics

Save the target population and requested predictions; structural model and
route; parameter scales/units and F convention; fixed/random effects including
omega covariance; residual error; starting values and any bounds; estimator and
controls; random seeds; exclusion policy; and prespecified comparisons. Covariate
extensions require explicit hypotheses and assessment of available information.
Starting values from the educational example are not study-derived priors.

Record package diagnostic outputs and human/agent review separately. At minimum
review convergence/SAEM traces, covariance and SE, correlation/conditioning when
available, boundary estimates, shrinkage, individual trajectories and residual
patterns. Do not invent universal cutoff values. Missing SE/covariance withholds
precision claims; unstable or implausible fits remain exploratory/audit only.
Save reasons at parameter and prediction level, not just a single pass flag.

Model selection, sensitivity analysis, bootstrap and VPC each get a declared
purpose and budget. Retain failed replicates and alternative runs. No VPC is
produced merely by plotting observed versus fitted concentrations. Do not compare
objective values across different likelihood approximations without justification.

## Artifacts and narrative

Each run saves source files and hashes, mapped event data, mapping log, study
plan, `model.R`/`analysis.R`, dependency versions (prefer a project lockfile),
`sessionInfo.txt`, seed/controls, console/warnings, fit RDS, parameter and
prediction tables, diagnostics and figures. A failed fit still retains its inputs
and error. Never overwrite a previous fit directory.

Before report and manifest, create extraction-only `interpretation_facts.json`:
engine/version, `validation_status: prototype_unvalidated`, subject/observation
counts, parameter scales/units, reportable/withheld/audit states with reasons,
diagnostics, missing checks, and `must_mention`. Each value or diagnostic refers
to an actual saved table row or fit-object component. Preserve CL/F labeling,
shrinkage, missing uncertainty and unmet evidence gates in `must_mention`.
This is agent-authored extraction code, not an existing nlmixr2/AgenticPrism
automatic feature. Do not refit or introduce new inference during extraction.

Then save a manifest with relative artifact paths, hashes, invocation, timestamp
and environment. Base the narrative on saved facts; do not declare a fit validated
because the package is widely used. Re-rendering reads saved outputs only.
