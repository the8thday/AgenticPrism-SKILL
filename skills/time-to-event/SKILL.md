---
name: time-to-event
description: Analyze time-to-event outcomes such as survival or time to a humane endpoint in in-vivo efficacy studies (xenograft, syngeneic, infection models) - Kaplan-Meier curves with confidence limits, median survival, survival at landmark times, log-rank tests (asymptotic or exact-by-permutation for small arms) with predeclared pairwise families, Cox proportional-hazards hazard ratios with optional covariates, and a proportional-hazards test. Numerically equal to R survival. Not for recurrent events, competing risks, time-varying covariates, clustered subjects or interval-censored data.
---

# Time-to-event analysis

Read the [input and model contract](references/input-and-model.md). The endpoint
definition is a scientific decision, not a data column: establish it before any
analysis.

## Establish the design first

- **Event:** what exactly counts (death, humane endpoint such as tumor volume or
  body-weight threshold, relapse). Humane-endpoint euthanasia is usually the
  event itself; write the rule in `study.endpoint`.
- **Time zero:** randomization, first dose or inoculation. It must be the same
  for every subject (`study.time_origin`).
- **Censoring:** who is censored and why (alive at study end, death judged
  unrelated, technical loss). Kaplan–Meier, log-rank and Cox all assume
  censoring is unrelated to prognosis. If animals were removed because they
  were doing badly, that is an event, not censoring; ask the user rather than
  choosing. State the reason in `study.censoring_rationale`.
- **Unit:** one row per animal/subject. Cage or litter clustering, repeated
  events and competing causes of death are not handled; say so if present.

## Choose the analysis

- Comparisons: `arm_vs_control` (declare `control_arm`) or `all_pairs`, chosen
  before looking at curves. Pairwise tests are log-rank with Holm-adjusted p;
  hazard ratios come from one Cox model with Bonferroni simultaneous limits.
- `logrank_inference`: `permutation` (exact up to Monte Carlo error when arm
  labels are exchangeable; recommended for small animal arms) or `asymptotic`
  (chi-square, equal to R `survdiff`). In null simulations with 6 animals per arm
  the asymptotic test rejected 6.6% at nominal 5% while permutation held 4.9%.
  Do not switch after seeing p values.
- `conf_type`: `log-log` (SAS/Prism asymmetrical, default) or `log` (R survfit
  default). With 10 animals per arm, log-log limits at the median covered 96.4%
  and log limits only 91.9% in simulation; keep log-log unless matching R output. `landmarks`: prespecified times for survival probabilities.
- `covariates`: numeric baseline variables (for example baseline tumor volume)
  for an adjusted Cox model, only if prespecified; log-rank stays unadjusted.

Locate the runtime with the [shared runtime instructions](../agentic-prism/references/runtime.md),
then run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY` and
`agentic-prism verify --run DIRECTORY`. Inspect `arm_summary.csv`, `km_curves.csv`,
`logrank.csv`, `contrasts.csv`, `cox_terms.csv`, `ph_test.csv`, `diagnostics.json`
and `report.html`. Start from the [synthetic example](../../fixtures/survival_synthetic/config.json);
never copy its endpoint or censoring statements.

## Report and interpret

- For new 0.8.1 runs, read `interpretation_facts.json` after verification
  ([shared contract](../agentic-prism/references/interpretation-facts.md)). Lead
  with arm-specific KM summaries and the prespecified log-rank comparison.
  Keep Cox reportability separate. With zero informative log-rank df, the saved
  audit p=1 supplies no test of arm equality. With reduced rank, name the
  limitation; pairwise tests lacking saved rank information are withheld in
  the facts. Do not recover them from CSV to bypass that state.
  Ask the user when endpoint or censoring definitions contradict removal notes.
  "Not reached" stays unavailable; a lower median confidence limit is not a
  point median. Include all material disclosures, even when another result
  in the same run is reportable.
- Per arm: n, events, censored, median with its limits, or "not reached" (then
  give the lower limit and landmark survival instead of inventing a median).
- Overall and pairwise log-rank p (state asymptotic or permutation), and the
  hazard ratio with limits. A hazard ratio summarizes the whole curve only under
  proportional hazards: if `ph_test` global p < 0.05 or the curves cross or
  separate late, call the HR a time-averaged effect and describe the curves and
  landmark survival instead of leaning on the HR.
- Surface diagnostics: few events per arm, fewer than 10 events per Cox
  parameter, an arm without events (HR withheld), landmarks beyond follow-up.
  A saved Greenwood interval of 1–1 before any event is degenerate, not evidence
  of certain survival; retain that limitation alongside the saved estimate.
- Log-rank tests equality of the whole curves; "not significant" does not show
  equivalence. With small arms, confidence limits are wide: report them.
- Do not describe animals as "cured", and do not extrapolate beyond follow-up.

Validation against R survival and small-sample simulations: [0.8.0 evidence](../../validation/RELEASE_0.8.0.md).
