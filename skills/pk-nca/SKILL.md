---
name: pk-nca
description: Guide pharmacokinetic concentration-time review and noncompartmental analysis in R with PKNCA, including Cmax, Tmax, exposure and terminal-phase assessment. Use for individual PK profiles and NCA study summaries; population compartmental modeling belongs to population-pk.
---

# PK NCA with PKNCA

Status: local guidance prototype, 2026-10-13. PKNCA performs the calculations;
do not implement an alternative NCA engine in Python or R. This is not a new
`agentic-prism run` analysis type. Package-backed numerical execution and the
collection's five evidence gates are not yet validated for this prototype.

For a new study, follow the collection's [analysis entry](../agentic-prism/references/analysis-entry.md)
and preserve any existing discussion/direct-analysis choice. Read the
[input and reporting contract](references/contract.md) before preparing a run.
Use the [R workflow](references/r-workflow.md) for environment checks and a
small, explicitly educational template.

## Applicability

Start with individual, adequately sampled concentration-time profiles and an
identified analyte, matrix, route, dose history and time origin. Separate subjects,
periods and analytes. Summarize PK parameters calculated per subject, not a
pooled mean concentration curve. Destructive/sparse sampling cannot be turned
into individual profiles; sparse NCA requires a separately specified extension.

Initial guidance covers single-dose Cmax/Tmax/AUClast, and terminal-dependent
parameters only when their assumptions and diagnostics are documented.
Repeated-dose/steady-state, partial-window imputation, renal PK, bioequivalence
and individual compartmental fitting are follow-ups, not demonstrated here.
For population PK, use [population-pk](../population-pk/SKILL.md). Ordinary
curve-auc does not supply PK terminal extrapolation or dose-dependent parameters.

## Workflow

1. Inspect actual concentration and dosing files. Establish the requested
   endpoints and units. Ask only for missing facts that affect the analysis;
   unknown dose need not block dose-independent observed-profile summaries.
2. Save original files and a mapping/preprocessing log. Resolve duplicates,
   missingness, BLQ and time-zero handling explicitly before calculation.
   Never set all BLQ values to zero or subtract predose concentration silently.
3. Declare interval boundaries, AUC integration method, terminal-selection
   policy and any reportability criteria with sources. Use explicit intervals
   in `PKNCAdata`; inspect the installed package help and resolved options.
   A package default is a software choice, not a study acceptance criterion.
4. Run PKNCA in a new output directory using a saved R script. Retain the
   package result object, full result table, warnings and exclusions. Do not
   replace excluded package estimates with hand calculations.
5. Review linear/log concentration plots and the terminal fit. Preserve selected
   points, selection rule, lambda-z diagnostics and extrapolated fraction. A
   numeric half-life alone does not establish a credible terminal phase.
6. Classify each subject/endpoint as reportable, withheld or not requested.
   Unsupported lambda-z withholds dependent half-life/AUCinf/CL/Vz, while
   Cmax/Tmax/AUClast may remain usable if their own inputs are adequate.
7. Save extraction-only interpretation facts and provenance before writing the
   narrative. State observed versus extrapolated exposure and per-endpoint n.

## Interpretation boundaries

After extravascular dosing with unknown bioavailability, label dose-based
clearance and volume as **CL/F** and **Vz/F**. Units must follow the dose basis:
mg/kg and mg are not interchangeable. An absorption-limited terminal phase can
reflect absorption rather than elimination. Do not infer linear PK, steady state,
bioequivalence or population variability from this workflow alone.

The [prototype evidence record](../../validation/PK_R_SKILLS_PROTOTYPE_2026-10-13.md)
states what was checked. Upstream credibility does not establish correctness of
our input mapping, agent behavior or study-specific scientific conclusions.
