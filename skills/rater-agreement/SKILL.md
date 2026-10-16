---
name: rater-agreement
description: Inter-rater or inter-operator agreement for independently rated units - Cohen's kappa (unweighted, linear/quadratic or sourced custom weights) for two raters and categorical scores, and intraclass correlation (ICC) in the declared McGraw-Wong/Shrout-Fleiss form for numeric scores, with confidence intervals. Use for pathology scores, IHC grades, image reads or operator repeatability; not for diagnostic accuracy against a reference (diagnostic-accuracy) or comparing two measurement methods (method comparison / Bland-Altman).
---

# Rater agreement (kappa, ICC)

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: rater_agreement` with `method: cohen_kappa` or `icc`. Locate the
runtime with the [ordered procedure](../agentic-prism/references/runtime.md).

## Choose by the question, before seeing results

- **Categorical scores, two raters → Cohen's kappa.** Declare the category list
  (all possible categories, in order if ordinal). Unordered categories use no
  weights. Ordinal categories may use linear or quadratic agreement weights, or a
  custom weight matrix, only with a source (protocol, guideline). Never try
  several weightings and keep the highest kappa.
- **Numeric scores → ICC.** Declare the form with a rationale:
  - `one_way_random`: each unit rated by a different set of raters.
  - `two_way_random`: the same raters rate every unit and represent a larger
    pool you want to generalize to.
  - `two_way_mixed`: these specific raters are the only ones of interest.
  - `agreement` (absolute values must match) or `consistency` (only ranking or
    additive shifts matter; systematic rater offsets are ignored).
  - `single` (one rater's score will be used in practice) or `average` (the mean
    of k raters will be used).
  Forms answer different questions and give different values; never pick the
  form by its size. Run only the declared form; do not tabulate the other forms
  "for context" unless the user asks a separately declared question. Two measurement methods on a continuous scale are a method
  comparison (Bland-Altman/Deming in correlation-regression), not ICC.
- **Units.** Each rated unit is independent (slide, animal, patient). Technical
  re-reads by the same rater are not extra raters; average them by a declared
  rule first or treat them as a separate repeatability question.
- **Complete data.** Every unit is rated once by every rater (two-way) or by the
  same number of raters (one-way). Missing ratings are refused; exclusions are
  whole units with reasons. No imputation.

## Run and report

Kappa: report kappa with its interval, observed and expected agreement, the
rating table, the weighting and its source. Benchmark labels ("moderate",
"substantial") are conventions, not inference. With skewed prevalence, kappa can
be low despite high observed agreement; report prevalence/bias indices, PABAK and
maximum kappa as context, not as substitutes. The large-sample interval is
approximate; it undercovered for binary kappa in calibration (see evidence).

ICC: report the form label (for example ICC(A,1)), estimate, interval, F test,
mean squares and rater means. Consistency ignores rater offsets; average-measure
ICC does not describe a single rater. Negative estimates are reported as
observed. Agreement is reliability, not accuracy: raters can agree and both be
wrong.

Respect `must_mention`, withheld states and failing items in
`interpretation_facts.json`.

## Evidence (0.13.6)

See [the 0.13.6 record](../../validation/RELEASE_0.13.6.md): base-R formula oracle
and psych agreement, Shrout & Fleiss (1979) ICCs and Cohen (1968) kappas
reproduced, and the registered calibration including the binary-kappa miss.
