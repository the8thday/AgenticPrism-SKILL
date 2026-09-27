---
name: method-validation
description: Summarize bioanalytical method-validation experiments for ligand binding assays (ICH M10-style) from back-calculated concentrations - accuracy and precision per QC level with total error, dilution linearity and hook effect, parallelism, selectivity, specificity and stability - against acceptance criteria the user declares with their source. Adds labelled statistical supplements (bias and precision intervals, beta-expectation accuracy profile, parallelism trend, exact pass-rate intervals). Not a calibration-curve fit and not a regulatory verdict.
---

# Bioanalytical method validation

Read the [runtime procedure](../agentic-prism/references/runtime.md), then the
[input and numerical contract](references/contract.md). Implemented in 0.9.3.

## Before analysis

1. **Establish which experiment the data are.** One run of the tool analyzes one
   experiment: `accuracy_precision`, `dilution_linearity`, `parallelism`,
   `selectivity`, `specificity` or `stability`. Do not merge experiments to get
   one "pass".
2. **Inputs are back-calculated concentrations** from each run's own calibration
   curve. Raw signals need calibration first (elisa-quantification or the lab's
   system). Calibration-curve acceptance and run acceptance are upstream and
   outside this tool.
3. **Acceptance criteria come from the user**, with a source: their SOP, or the
   guideline edition they follow. Never fill them in because the user said
   "ICH M10". For reference only, ICH M10 (2022, Step 4) section 4.2 states for
   LBA: accuracy within ±20% (±25% at LLOQ and ULOQ), within- and between-run
   precision ≤20% (≤25%), total error ≤30% (≤40%); selectivity in at least 10
   individual matrices with at least 80% of sources meeting ±25% at LLOQ and ±20%
   at high QC, blanks below LLOQ in at least 80%; dilution QC mean within ±20%
   and precision ≤20% for at least 3 dilution factors; stability mean within
   ±20% of nominal. Section 7.2: parallelism CV ≤30%. Ask the user to confirm the
   values against their document before running; say that you quoted them.
4. **Independent runs and exclusions.** A run is an independent analytical run
   with its own calibration curve. Exclusions need a documented, obvious error
   (M10 includes all results otherwise). Do not exclude results to make a level
   pass.
5. **CV denominator** (`observed_mean` or `nominal`) must be declared for
   accuracy/precision. Ask if unknown.

## Run and interpret

Run `analyze --config ... --output NEW`, then `verify --run ...`, then read
`interpretation_facts.json`. Lead with the experiment, the declared criteria
and their source, and the per-level (or per-group) result. Then report every
item in `must_mention` and `failing_items`.

- **Accuracy and precision:** report bias %, within-run CV, between-run
  (intermediate) CV and total error per level, and which check failed. A level
  that fails is a failed result, not a reason to drop a run. The supplements
  (90% bias interval, MLS interval for between-run CV, beta-expectation
  tolerance interval) show uncertainty; the acceptance rule uses point estimates.
  When the accuracy profile narrows the range at an end (facts say so), report
  both: the level-based rule passes, the profile suggests the end is marginal.
  The profile is a supplement, not an M10 requirement.
- **Dilution linearity:** the headline is false when a hook effect is suspected,
  even if every in-range dilution passes. Report the dilution-factor range that
  was actually validated. Diluting study samples outside that range is not
  supported by this experiment.
- **Parallelism:** the CV rule can pass while concentrations drift with dilution.
  Report per-sample log slopes and any sample flagged against the declared slope
  margin, as M10 asks to watch for trends. Non-parallelism may come from matrix
  components or binding proteins; do not attribute a cause without evidence.
- **Selectivity / specificity:** report passing/total per role (and interferent),
  the required fraction and the exact interval. With 10 sources, 8/10 meets an
  80% rule but the interval is wide; say so when the user wants assurance about
  future matrices.
- **Stability:** report the mean accuracy per condition and level. A condition
  that passes on the mean with an interval extending beyond the limit is
  marginal; say so, and do not extend the claim to longer storage.

A single A&P pass or fail is itself uncertain. In the 0.9.3 simulation (6 runs ×
3 replicates, Gaussian), a level with no bias and a true intermediate CV of 15%
met the ±20%/20%/30% rules 91.7% of the time, and at 18% only 72.8%; with a true
bias of 15% and CV of 15%, 66.0%. Use the intervals to say how close a level is
to its limits instead of treating the verdict as exact.

Never state that the method is "validated" or "M10 compliant" from one
experiment. Passing is evidence for that experiment under the declared criteria.
Keep synthetic-fixture results separate from real data.

## Numerical evidence (0.9.3)

Every computed quantity is checked against an independent R implementation
(VCA::anovaVCA, lm, binom.test, t.test, approx): 188 fields, maximum relative
difference 3.4e-12. Calibration (1000 datasets per row, seed 20261005): the
beta-expectation interval (β = .80) had mean content 0.785–0.821 across five
designs, the 90% bias interval covered 90.8–95.2% and the 90% MLS interval for
between-run CV 89.5–92.0%; all rows met their prespecified criteria, including
3 runs × 3 replicates. Calibration of the statistical supplements is in
[the 0.9.3 record](../../validation/RELEASE_0.9.3.md); read the rows before
relying on an interval. No printed worked example of the Mee (1984) interval
was available; its formula was derived and checked independently.
