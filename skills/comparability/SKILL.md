---
name: comparability
description: Compare one quality attribute between reference and test lots for process comparability or biosimilarity - equivalence of lot means by TOST against a declared margin, a quality range (reference mean plus or minus k SD) or a descriptive summary. Tier, method, margin and k are user declarations with a source; the lot is the unit. Not a totality-of-evidence judgement.
---

# Comparability and biosimilarity

Follow [runtime discovery](../agentic-prism/references/runtime.md), then read the
[input contract](references/input-contract.md). AgenticPrism 0.10.1;
`analysis_type=comparability`, one attribute per run, a new output directory.

## Gate before analysis

- **The lot is the unit.** Replicate measurements of a lot are averaged first; they
  are not extra lots. Ask how lots were selected (all available, or a subset chosen
  how and when).
- **Tier, method, margin and k come from the user's protocol or risk assessment,**
  with a source. Never assign a tier yourself, never infer the method from the tier,
  and never choose a margin or k after seeing the data. If they are missing, stop
  and ask. The FDA's 2017 draft (withdrawn in 2018) used 1.5 × reference SD margins
  for tier 1 and mean ± k·SD ranges for tier 2; the user must say which convention
  their protocol follows.
- An SD-multiple margin is estimated from the reference lots, so it varies with
  them; say so.

## Run and interpret

Run `analyze`, `verify`, then read `interpretation_facts.json`. Lead with the
attribute, tier and its source, the method, then:

- **TOST:** difference of lot means (test − reference) with its 100(1 − 2α)%
  interval, the margin and its source, and whether the interval lies inside.
  Equivalent means are not identical lots: individual test lots can still differ.
  Non-equivalence is not proof of a difference when lots are few.
- **Quality range:** the range, how many test lots fall inside, which ones do not,
  and the required fraction if declared. A quality range has no nominal error rate:
  with identical products, 10 reference and 6 test lots, all test lots fell inside
  a k = 3 range in 88% of simulations and inside k = 2 in only 61%. A lot outside
  the range needs investigation, not an automatic conclusion.
- **Descriptive:** summaries only; no verdict.

Mention every `must_mention` item. Never state that products are "biosimilar" or
"comparable" overall from one attribute: that is a totality-of-evidence judgement
across attributes, functional and clinical data.

## Numerical evidence (0.10.1)

TOST intervals, degrees of freedom and p values equal R `t.test` (Welch) to 1e-10;
see [the 0.10.1 record](../../validation/RELEASE_0.10.1.md). Calibration: TOST type I
error at a fixed margin 4.5% (equal SD) and 3.1% (test SD twice the reference),
bound 6.38%; with an estimated 1.5 × SD margin 5.8% (descriptive).
