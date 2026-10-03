---
name: sample-size
description: Plan prospective sample size or power from declared, sourced assumptions for two-sample or paired t tests, one-way ANOVA, two proportions, log-rank (events) and TOST equivalence of two means (e.g. CMC comparability lot counts), with exact noncentral t/F and TOST power, a saved power curve and a sensitivity grid. Use when designing an in vivo, in vitro or comparability study; not for post-hoc ("observed") power of a completed study.
---

# Sample size and power

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[input and method contract](references/input-and-model.md). Implemented in 0.13.1.

## Before calculating

1. **Name the planned analysis.** The power belongs to one test analysed as
   planned. Choose the design that matches it: independent groups
   (`two_sample_t`, `one_way_anova`, `two_proportions`), the same subjects twice
   (`paired_t`), time to event (`logrank`), or showing equivalence
   (`tost_two_means`). Repeated measures over time, animals clustered in cages,
   multiple endpoints or interim looks are not covered: say so and stop, or ask
   for a statistician; never feed per-timepoint or per-well n into these formulas.
2. **Every assumption needs a source** (`<key>_source`): pilot or historical data,
   literature, or a protocol. The difference to detect is the smallest difference
   that would change a decision, not the effect the team hopes for. Never take
   the difference or SD from the study being planned, and never use "observed
   power" of a finished study.
3. **Declare a sensitivity grid** for the least certain assumption (usually the
   SD or the event probability). Ask the user for plausible values if they have
   none; the report lists how n moves.
4. **TOST margins** come from a predeclared protocol or reference-lot data (see
   the comparability Skill). Never pick the margin that makes n feasible.
5. For `logrank`, the declared event probability (follow-up, accrual, dropout)
   converts required events into subjects. Power comes from events.

## Run and report

Run `analyze --config ... --output NEW`, `verify`, then read
`interpretation_facts.json`. Report: the design and test, alpha and sidedness,
each assumption with its source, the n per group (and total), the power reached
at that integer n, the sensitivity rows, and every `must_mention` line. Add
expected attrition yourself only if the user gives a rate, and show it as a
separate step. Proportions and log-rank are approximations; for small n or
extreme proportions, suggest a simulation check.

## Evidence (0.13.1)

Power functions agree with pwr 1.3.0, `power.prop.test` and PowerTOST 1.5.7
(exact) in 306 cases, maximum absolute difference 9.5e-10; the integer n matches
`ceiling(pwr.t.test(...)$n)`. Simulating 4000 studies at the solved n gave
rejection rates within two binomial SE of the computed power for all ten
registered rows, including the approximate proportion and Schoenfeld log-rank
methods against Pearson chi-square and the package log-rank. The first
benchmark run hid three NaN powers (fixed and retained). The published
worked-example gate is unmet; live-agent misuse scenarios passed ([live-agent results](../../validation/agent-scenarios/0.13.1/RESULTS.md)). See
[the 0.13.1 record](../../validation/RELEASE_0.13.1.md).
