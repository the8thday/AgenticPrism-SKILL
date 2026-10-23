---
name: pk-nca
description: Single-dose noncompartmental PK analysis (NCA) computed by R PKNCA - Cmax, Tmax, AUClast, lambda-z and half-life, AUCinf, CL or CL/F, Vz, MRT and Vss for IV bolus (back-extrapolated C0), IV infusion or extravascular dosing - with declared BLQ, AUC-method, terminal-phase and reportability rules, subject-level summaries and an optional Smith power-model dose-proportionality assessment. Use for individual animal or human PK profiles; population compartmental models belong to population-pk.
---

# PK NCA with PKNCA

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Use `analysis_type: pk_nca` (0.15.0). R PKNCA computes every NCA value through a
fixed packaged script; the runtime checks the contract, applies your declared
reportability, summarises and reports. It needs R and the pinned PKNCA packages
(`python3 install.py --with-r-pk`; `agentic-prism doctor` shows `optional_r_pk`). If
they are missing, say so; do not compute NCA by hand or in another tool. Read the
[input and output contract](references/contract.md). Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Design facts come from the study, not the curve.** Analyte (free/total),
   matrix, species, route (`intravascular_bolus`, `intravascular_infusion`,
   `extravascular`), dose and its unit (mg or mg/kg), time origin, actual versus
   nominal times, infusion duration. Single dose at time 0 per subject only.
2. **Serial sampling per subject.** Each subject must contribute its own profile.
   Destructive or sparse sampling (one sample per animal) cannot be assembled into
   pseudo-subjects; say it is unsupported here.
3. **BLQ.** Mark BLQ rows `blq = true` with an empty concentration. When any
   post-dose sample is BLQ, the user declares first/middle/last rules
   (`keep_as_zero` or `drop`) with a rationale. Never substitute LLOQ/2, delete
   values or set them to zero without that declaration.
4. **AUC method** (`linear` or `lin_up_log_down`) is declared with a rationale.
   WinNonlin's default is linear trapezoidal; results differ, so state it.
5. **Terminal phase.** Default: PKNCA's best-adjusted-R² rule over points after
   Tmax (and after the end of infusion). `declared_points` uses a `lambda_z_point`
   column and still needs at least `min_points`. Never choose points to obtain a
   preferred half-life; report the rule actually used.
6. **Reportability.** Either criteria with a source (minimum adjusted R², span
   ratio, maximum % AUCinf extrapolated; any may be null) or `descriptive_only`
   with a rationale. Criteria are the study's choice, not a universal rule.
7. **Dose proportionality (optional)** needs at least 2 dose levels (3 is better)
   and an acceptance ratio interval with a source (Smith et al. 2000 commonly use
   0.80-1.25). It describes exposure scaling across the tested range only.

## Run and report

Report per subject and per group: n reportable and n withheld with reasons;
geometric mean and CV% for exposure parameters; Tmax as median and range. For
extravascular dosing say CL/F and Vz/F (bioavailability unknown) and that MRT
includes absorption. For IV bolus say C0 was back-extrapolated and check the
`c0_not_log_back_extrapolated` flag (the first two samples were not declining, so
C0 is the first observed value). Withheld terminal parameters are not results;
do not quote them in summaries. Dose proportionality: quote b, its interval, the
acceptance interval for b and the three-way conclusion. `inconclusive` is the
usual result of a small study and is not evidence of nonlinearity or of
proportionality; a mechanism (TMDD, saturable absorption) needs a model.
Respect `must_mention` and failing items in `interpretation_facts.json`.

Multiple-dose/steady-state, partial AUCs, urine, sparse NCA and bioequivalence
are not supported by this analysis type. For them, write a saved R script with
PKNCA following [the R workflow](references/r-workflow.md) and label it as a
custom analysis outside the validated runtime.

## Evidence

See [the 0.15.0 record](../../validation/RELEASE_0.15.0.md): the published
WinNonlin worked example reproduced, an independent formula recheck, the
registered dose-proportionality calibration (including its low power) and the
live-agent scenarios.
