---
name: outlier-flagging
description: Prespecified outlier flagging within declared groups of replicate values - Grubbs single-outlier test, generalized ESD for up to a declared number of outliers, or ROUT (robust fit with false discovery control) - reporting flags and results with and without flagged values side by side, never removing data. Use when a protocol requires an outlier screen of replicate determinations, lots or donors; not for choosing which points to drop to reach significance, and not for curve-fit residuals.
---

# Outlier flagging (Grubbs, generalized ESD, ROUT)

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: outlier_screen`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Rules that apply before any screen

1. **Flag, never delete.** The module removes nothing, and neither do you. A
   flagged value may be excluded from a later analysis only with a documented
   assignable cause (recorded pipetting, instrument or sample error). A test
   result alone is not a cause. Report the main result with all values, and any
   sensitivity analysis without flagged values beside it.
2. **Prespecified.** The method, alpha (or Q), sidedness, maximum number of
   outliers and the scale must come from a protocol, SOP or analysis plan named in
   `design.prespecified_source`. If the user has none, say the screen is post hoc
   and exploratory; never pick a method or alpha after seeing which one flags the
   "inconvenient" point, and never screen only the group that spoils a result.
3. **Right sample.** Each declared group must be one homogeneous sample:
   replicate determinations of the same thing. Values from different doses,
   timepoints or conditions are not one sample; residuals of a fitted curve are
   outside this module. Technical wells inside one donor are not donors.
4. **Right method.**
   - Grubbs tests exactly one outlier. Do not repeat it after removing a point:
     that inflates the error rate and suffers masking (in calibration Grubbs never
     detected both of two planted outliers).
   - Generalized ESD handles up to a declared r suspects.
   - ROUT uses a robust fit and a false discovery level Q (1% is the authors'
     default); it follows the published description and is not Prism-equivalent.
5. **Distribution.** Grubbs and ESD assume normal data apart from outliers;
   skewed or heavy-tailed data give false flags (clean t(3) data were flagged in
   about 40–48% of samples). Use the log scale only when declared with a reason
   before screening (for example log-normal concentrations); do not switch scales
   to make a flag appear or disappear.

## Run and report

Report per group: n, the test statistic and critical value (or the ROUT table),
the flagged units and values, and summaries with and without flagged values.
Small groups (fewer than 10) have little power: no flag does not mean no outlier.
Respect `must_mention` and withheld groups in `interpretation_facts.json`.

## Evidence (0.13.6)

See [the 0.13.6 record](../../validation/RELEASE_0.13.6.md): NIST/SEMATECH
Grubbs and generalized ESD examples reproduced, base-R formula oracle agreement,
and the registered false-flag calibration with masking and heavy-tail results.
ROUT has no published numeric worked example (gate UNMET).
