---
name: potency-ranking
description: Rank many antibody or compound candidates against a benchmark tested in the same independent experiments - paired log EC50/IC50 differences, relative potency with Holm-adjusted p values and Bonferroni simultaneous intervals, optional prespecified fold-threshold tiers, and withholding of candidates with non-estimable runs. Use for discovery triage of clones, variants or hits from saved dose-response estimates; not for fitting the curves themselves (dose-response) or for a single test-versus-reference potency assay (potency-assay).
---

# Multi-candidate potency ranking

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: potency_ranking`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Curves first.** EC50/IC50 values must come from curve fits that passed
   their own checks (run dose-response first). Every candidate-experiment
   combination is entered, including runs whose EC50 was not estimable
   (`above_range`, `below_range`, `failed`). Never delete them: a candidate whose
   EC50 is above the top dose in some runs would otherwise look more potent than it
   is. Such candidates are withheld, not ranked.
2. **Same-experiment benchmark.** Each candidate is compared with the benchmark
   tested in the same independent experiment (separate cell thaw, plate day).
   Technical replicate wells or curves within one plate are summarized first by
   the dose-response step; they are not experiments. At least 3 shared
   experiments are needed.
3. **One family.** All candidates versus the benchmark are one family: report
   Holm-adjusted p values and simultaneous intervals. Do not pick the candidate
   with the smallest unadjusted p value, and do not drop "boring" candidates
   before the analysis to make the family smaller.
4. **Thresholds before results.** A fold threshold (for example "at least 3-fold
   more potent") is optional and must come from a triage rule set before the
   data; it is applied to simultaneous intervals only.

## Run and report

Report each candidate's relative potency (10^−mean log10 difference; > 1 = more
potent), its simultaneous interval, Holm p value, number of experiments and any
withholding reason; the reference's failed experiments; and the threshold tier
if declared. Point ranks are descriptive: candidates are not compared with each
other, and overlapping intervals do not order them. With few experiments power is
low; "not significant" does not mean "no better than the benchmark". Respect
`must_mention` and failing items in `interpretation_facts.json`.

## Evidence (0.14.0)

See [the 0.14.0 record](../../validation/RELEASE_0.14.0.md): base-R oracle
agreement, the registered family-wise error and coverage calibration with its
low-power rows, and the unmet published-example gate.
