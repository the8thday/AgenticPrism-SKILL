---
name: potency-assay
description: Combine replicated relative-potency determinations across independent runs and validate relative accuracy, intermediate precision, log linearity and tested range. Uses a log-RP random-run model with REML and MLS/MOVER, declared acceptance criteria and per-run suitability. Does not identify separate run and residual variance from one RP per run.
---

# Relative potency across runs

Follow [runtime discovery](../agentic-prism/references/runtime.md), then the
[input contract](references/input-contract.md). AgenticPrism 0.10.0 uses
`analysis_type=potency_assay`. Write each analysis to a new directory.

## Applicability and inputs

Establish the estimand: geometric mean RP across independent runs, or validation
at nominal RP levels. Require at least three runs and two independently fitted
RP determinations per run per nominal level. Technical wells and RP estimates
sharing a reference fit are not independent determinations. One RP per run
cannot separate run and residual variance in this model; stop and explain the
required design. Do not duplicate a run to make the model fit.

Require the user's rationale for independent runs/determinations, common log
residual variance and the same material at each nominal level. Heterogeneous
known-SE meta-analysis and correlated/shared-reference estimates are unavailable.

For raw curves, use [dose-response](../dose-response/SKILL.md) first, retaining
its `dose_fit.py` parallel model and `equivalence.py` parameter-equivalence
checks. `input_mode=dose_runs` imports verified saved comparisons and preserves
their scientific artifacts. Do not refit with a new curve model or reinterpret
an F-test's nonsignificance as equivalence. Alternatively, use an explicitly
sourced table of RP estimates with reference-quality and parallelism statuses.

## Suitability and criteria

Every determination needs a source and `pass`, `fail` or `not_evaluable` for
reference quality and parallelism. The user supplies the suitability protocol.
Any failing or unevaluable supplied run withholds the combined result. Never
silently drop a failing run or turn an audit RP into an accepted result.
Prespecified single Grubbs screening of equal-sized run means is available;
report every flag, its threshold and source, and retain all runs. Do not screen
post hoc or iteratively remove flagged runs. This is scoped screening, not a
claim to implement every USP <1010> procedure.

Validation criteria always require a source: absolute natural-log bias and
intermediate geometric CV limits. They are never defaulted from USP chapters.
For linearity equivalence, require slope and intercept margins that bracket
1 and 0 respectively, prespecification and justification in `criteria.source`.
Do not choose margins after seeing the current observations.

## What to report

Read verified `interpretation_facts.json`; lead with status, combined RP and
interval, or the reason it is withheld. Include all must-mention disclosures
and every failing item. Report per-run diagnostics and every outlier flag.

The model is natural-log RP = fixed mean (or intercept + slope × log nominal RP)
+ random run intercept + independent residual. Report log bias and transformed
relative bias, run/residual components, intermediate precision and MLS interval
(MOVER with fitted correlation when unbalanced). Geometric CV is
100 × sqrt(exp(total log variance) − 1), for one determination in a future run.
Fixed-effect intervals use GLS covariance and t with runs−1 df; balanced mean
inference is conservative at a boundary and other designs are approximate.

Bias/precision verdicts use user-declared point criteria. Intervals are labelled
statistical supplements. Linearity reports slope/intercept intervals; the
optional equivalence decision requires the entire intervals inside the declared
limits (alpha per side = (1−confidence)/2). Without margins there is no validated
range. Range requires at least three contiguous passing tested levels plus
global linearity equivalence; no interpolation across failing levels or extension
beyond tested levels. Passing this scope does not establish full USP validation.

## Evidence limits

All 12 prespecified potency coverage/type-I rows passed the fixed bounds in
1,000-dataset simulations. The intercept-margin rejection rate was 5.7%
(MCSE 0.73 percentage points); the slope-margin rate was 2.0% (0.44).
No published worked example matching this replicated random-run model was
obtained: that evidence gate is **not met**. Do not call synthetic fixtures a
published reproduction. See [0.10.0 evidence](../../validation/RELEASE_0.10.0.md).
The numerical runtime is Python; R lme4 is a validation oracle only.
