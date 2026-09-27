---
name: tumor-growth
description: Analyze longitudinal tumor volumes from in-vivo efficacy studies (xenograft, syngeneic) - arm-specific exponential growth rates from a log-volume random-slope mixed model (Satterthwaite, numerically equal to R lmerTest), rate differences, doubling times and model-based T/C versus control, plus observed TGI% and T/C% at a prespecified day with Fieller limits and explicit handling of humane-endpoint dropout. Use with time-to-event for the survival endpoint. Not for tumor regression kinetics models, bioluminescence counts, or non-exponential growth fitted as curves.
---

# Tumor growth

Read the [input and model contract](references/input-and-model.md). Tumor-growth
readouts are easy to overstate: decide the estimand and how dropout is handled
before running anything.

## Establish the design first

- One row per animal and measurement day (`animal_id, arm, day, volume`); each
  animal in one arm; caliper replicates averaged upstream with provenance.
- Time zero (`study.time_origin`, usually first dose), the humane-endpoint /
  removal rule (`study.removal_rule`), and the prespecified readout day
  (`analysis.analysis_day`) and baseline day. Every animal needs a baseline value.
- `log_offset`: volumes are analysed as log(V + offset). With complete
  regressions (V = 0) the offset must be positive; state its value and treat
  results as conditional on it. Never pick the offset that gives the best p.
- `study.dropout_rationale`: why removal can be treated as depending only on
  observed volumes (MAR). If animals left for reasons tied to unobserved
  worsening, the growth model is also biased; say so.

## What the Skill reports and how to read it

For new 0.8.1 runs, read `interpretation_facts.json` after verification
([shared contract](../agentic-prism/references/interpretation-facts.md)). Use its
primary growth-model results and their intervals first; preserve withheld
model inference even if observed ratios are available. Report observed ratios
only for the animals measured on the analysis day, including their counts and
dropout caveats. Ask the user if removal reasons or the MAR rationale are
unknown or contradictory. The facts record declarations; they cannot establish
MAR. Missing Fieller bounds are not a narrow or zero-width interval.

1. **Growth model (primary):** per-arm log growth rate with Satterthwaite limits,
   doubling time, the equal-rates F test, rate difference vs control (Holm p,
   Bonferroni limits) and the model-estimated geometric-mean T/C at the analysis
   day. This uses every measurement and remains valid under MAR dropout, which is
   why it leads: with about two animals per study removed before the readout day,
   the rate-difference family still covered 95.1% in simulation.
2. **Observed TGI% and T/C% (secondary):** at the analysis day, using animals
   still measured that day; TGI = 100 × (1 − mean ΔT / mean ΔC) with ΔV from
   baseline; Fieller limits with a Bonferroni family. If any animal was removed
   before that day (`dropout.csv`, diagnostic flags), these are biased toward the
   survivors (family coverage fell to 92.2% in that simulation, versus 95.4%
   without dropout); report them only with that caveat.
3. **Survival:** time to the humane endpoint belongs in the
   [time-to-event](../time-to-event/SKILL.md) Skill; recommend it whenever
   removal is common.

Other rules: a growth-rate difference is on the log scale (per day); T/C is a
geometric-mean ratio from the model or an arithmetic-mean ratio when observed:
name which. With a positive offset, model T/C is the geometric-mean ratio of
V + offset, and doubling time concerns V + offset, not raw V. Do not attach
the rate-difference p value to model T/C; those are different contrasts.
Do not describe TGI > 100% as "cure" (it means regression). Flag
fewer than four measurement days, random-slope boundary (Satterthwaite withheld),
or a control arm that did not grow (TGI undefined).

Locate the runtime with the [shared runtime instructions](../agentic-prism/references/runtime.md),
run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`, then
`agentic-prism verify --run DIRECTORY`. Inspect `growth_rates.csv`,
`growth_tests.csv`, `model_contrasts.csv`, `observed_tgi.csv`, `dropout.csv`,
`arm_day_summary.csv`, `diagnostics.json` and `report.html`. Start from the
[synthetic example](../../fixtures/tumor_growth_synthetic/config.json); never
copy its removal or dropout statements. Evidence: [0.8.0 validation](../../validation/RELEASE_0.8.0.md).
