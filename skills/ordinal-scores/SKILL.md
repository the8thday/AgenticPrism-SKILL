---
name: ordinal-scores
description: Compare ordered categorical scores (histopathology grades, clinical or arthritis scores, severity ratings) between groups of independent subjects with a cumulative-logit proportional-odds model - odds ratios with profile-likelihood intervals, a group likelihood-ratio test, declared pre-treatment covariates and a proportional-odds check. Use instead of t tests or ANOVA on score means; repeated scores over time are not covered.
---

# Ordinal scores

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Use `analysis_type: ordinal_regression` (0.15.3). Read the
[input contract](references/input-and-model.md). Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Scores are ordered categories, not numbers.** Do not average grades or run t
   tests on them; the distance between grade 1 and 2 is not that between 3 and 4.
2. **One score per independent subject.** Weekly scores from the same animal are
   not independent; summarise per animal by a rule declared in advance (for
   example the score at a fixed day or the maximum) or use a mixed model, which is
   not available here. Cage or litter clustering is not modelled.
3. **Declare the category order and every possible category** before seeing the
   data. A declared category nobody reached is refused: merge categories only by
   a pre-declared rule.
4. **Covariates** must be pre-treatment (for example baseline weight) and declared.

## Run and report

Report the odds ratio for each group against the reference with its profile
interval (higher scores when OR > 1), the group LR test, the score distribution
per group and the proportional-odds check. With few subjects that check is often
unavailable or has little power: not rejecting it does not establish
proportional odds, and a violation means a single OR does not summarise the
effect. Respect `must_mention` and failing items in `interpretation_facts.json`.

## Evidence

See [the 0.15.3 record](../../validation/RELEASE_0.15.3.md).
