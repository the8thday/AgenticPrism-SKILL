---
name: diagnostic-accuracy
description: Empirical ROC analysis of a continuous or ordinal marker against a declared binary reference standard - AUC with DeLong logit interval, Mann-Whitney test, paired DeLong comparison of two markers on the same units, and sensitivity/specificity at prespecified cut-offs with exact intervals. Use for biomarker, assay or rating discrimination between condition-present and condition-absent independent units; not for agreement between raters or methods.
---

# Diagnostic accuracy (ROC)

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: roc_analysis`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Reference standard.** Each independent unit (patient, animal, sample) has
   one condition classification (1/0) from a named reference standard. ROC
   accuracy is relative to that reference; an imperfect reference biases it.
   Units must be independent: several lesions, wells or timepoints from one
   subject are not separate units. Clustered designs are outside this module:
   offer a declared one-value-per-subject analysis through it (for example the
   most suspicious lesion, chosen by a rule that does not use the reference
   result), or state that a clustered ROC method is needed. Do not write your
   own clustered-ROC or bootstrap computation and present it as the result; any
   such calculation is an unvalidated, clearly labelled exploration at most.
2. **Direction before data.** Declare whether higher or lower marker values
   indicate the condition, with its source (biology, assay design, prior study).
   Never flip the direction because AUC came out below 0.5; report it as
   observed.
3. **Cut-offs.** Sensitivity and specificity are reported only at cut-offs
   declared before analysis with a source. A Youden "optimal" cut-off may be
   requested as `exploratory_youden`; it is selected from these data, optimistic,
   has no interval and is never a validated decision threshold.
4. **Two markers.** Compare AUCs with `comparison: paired_delong` only when both
   markers were measured on the same units. Different cohorts are not paired.
   A non-significant difference does not show the markers are equivalent.
5. **Missing values.** Missing or non-numeric marker values, BLQ text and
   duplicated units are refused. Do not impute, set BLQ to a number or drop
   units silently; excluded units need a reason.

## Run and report

Report per marker: n with and without the condition, AUC with its DeLong logit
interval, the Mann-Whitney p value, the declared direction, and flags (fewer
than 10 per class, many ties, AUC below 0.5). For cut-offs report sensitivity and
specificity with Clopper-Pearson intervals and the 2x2 counts. Predictive values
depend on prevalence and are not produced; case-control sampling makes them
meaningless. AUC is discrimination, not calibration or clinical utility.

If the interval is withheld (complete separation, zero DeLong variance), report
the observed AUC as descriptive only. Respect `must_mention` and failing items in
`interpretation_facts.json`.

## Evidence (0.13.6)

See [the 0.13.6 record](../../validation/RELEASE_0.13.6.md): base-R formula oracle
agreement on every fixture, the Hanley & McNeil (1982) AUC reproduced, and the
registered calibration including its stress rows.
