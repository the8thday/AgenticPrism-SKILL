---
name: panel-comparison
description: Compare two groups across a panel of many outcomes measured on the same units - cytokine/chemokine panels, flow-cytometry marker panels, multiplex biomarkers - with every measured outcome assigned to a declared family and Holm (family-wise error, simultaneous intervals) or Benjamini-Hochberg (false discovery rate) adjustment within each family. Use to avoid reporting only the outcomes that happened to reach significance; not for a single prespecified primary outcome (group-comparison). Outcomes with values below or above quantification limits stay censored through an opt-in censoring block.
---

# Panel comparisons with multiplicity control

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: panel_comparison`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Every measured outcome counts.** All outcomes in the file go into a declared
   family (for example "cytokines", "chemokines"); the runtime refuses data with
   undeclared outcomes. Do not remove outcomes that "did not work" or look flat;
   if an analyte failed QC, record that as a data decision with its reason before
   seeing group results.
2. **Families and procedure from the plan.** Choose families by biology or the
   analysis plan, not by which grouping makes more results significant. Choose
   Holm when any false positive matters (claims about specific outcomes) and
   Benjamini-Hochberg for screening where a controlled proportion of false
   discoveries is acceptable. Do not run both and keep the more favorable one.
   A single prespecified primary outcome may be analysed on its own in
   group-comparison; that does not remove it from the panel family.
3. **Design.** Declare pairing (the same donors/animals in both groups versus
   independent groups) and the scale (log for concentration-type outcomes, with a
   reason, before results). Technical replicates are averaged first.
4. **Values beyond quantification limits stay censored.** Enter the
   sample-specific limit as `value` with `censored: left` (or `right` above the
   ULOQ) and add the `censoring` block with the limit's source; the scale must be
   log. Never substitute 0, LLOQ or LLOQ/2 and never drop an analyte because many
   values are below the limit. Outcomes with censored values use censored lognormal
   maximum likelihood with a small-sample t correction (independent groups assume
   a common SD; paired designs accept left censoring only); the others keep the t
   test. An outcome with too few detected values is withheld but still counted as
   measured.

## Run and report

Report every outcome: estimate (difference or ratio), unadjusted p, adjusted p,
the decision at the declared alpha, and interval (simultaneous for Holm;
unadjusted, labelled, for BH). Report how many outcomes were tested in each
family. Under BH some discoveries may be false. Non-significant outcomes are not
shown to be unaffected. Respect `must_mention` and failing items in
`interpretation_facts.json`.

## Evidence

See [the 0.14.0 record](../../validation/RELEASE_0.14.0.md): base-R oracle agreement,
Benjamini & Hochberg (1995) reproduced, and the registered family-wise error and FDR
calibration. Censored outcomes: [the 0.14.1 record](../../validation/RELEASE_0.14.1.md)
(R survival oracle and a registered calibration; published-example gate unmet).
