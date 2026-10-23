---
name: quantal-dose-response
description: Fit quantal (responder, mortality, yes/no) dose-response data by probit, logit or cloglog maximum likelihood and report ED50/LD50 and other EDp with Fieller limits, with declared natural (background) response (none, Abbott known rate, or estimated from controls), Finney heterogeneity scaling and parallel-line relative potency. Use for counts of responding animals or subjects per dose; continuous responses belong to dose-response (4PL), ordered scores to ordinal-scores.
---

# Quantal dose response

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Use `analysis_type: quantal_dose_response` (0.15.3). Read the
[input contract](references/input-and-model.md). Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Counts, not percentages.** Each dose needs n tested and the number responding
   (or one 0/1 row per independent subject). Percentages without n, or responses
   averaged over replicate wells, are not quantal data.
2. **Independent units.** Animals, not cages or litters. Clustering inflates the
   χ²; the heterogeneity policy (declared before fitting) scales the intervals but
   does not repair a design problem.
3. **Link and metameter are declared** (probit on log10 dose is classical) with a
   reason; do not try links and keep the one that fits best.
4. **Natural response.** If controls respond, declare `abbott_known` (a rate from
   historical controls, with source) or `estimate` (concurrent control groups);
   `none` is refused when controls show responders. EDp is then the dose giving p%
   additional response beyond background.
5. **EDp within the tested range.** Values outside the tested doses and Fieller
   limits that are unbounded (g ≥ 1) are withheld: say so; never extrapolate.

## Run and report

Report the slope, χ² with df and whether the covariance was scaled, each EDp with
its Fieller limits in dose units, withheld EDp with the reason, and the natural
response used. Relative potency needs parallel lines: report the parallelism test
or equivalence result first; a non-significant LR test does not prove
parallelism, and a failed one withholds ρ. Respect `must_mention` and failing
items in `interpretation_facts.json`.

## Evidence

See [the 0.15.3 record](../../validation/RELEASE_0.15.3.md): Finney's printed
rotenone example, the R glm / MASS oracle, the registered calibration and its
first attempt, and the live-agent scenarios.
