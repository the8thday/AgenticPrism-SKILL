---
name: censored-data
description: Compare groups when some values are below a lower limit of quantification or detection (or above an upper limit) - cytokines, biomarkers, ADA titers, drug levels - keeping those values as censored at their limit instead of substituting zero, LLOQ or LLOQ/2. Censored lognormal maximum likelihood with geometric-mean ratios, or the Peto-Peto rank test for left-censored data. Use for independent units with one value each.
---

# Values below (or above) quantification limits

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: censored_comparison`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Keep censored values censored.** A "<LLOQ" result means the true value is
   somewhere below that sample's reporting limit. Enter the limit as `value` with
   `censored: left`. Never replace it with 0, the LLOQ, LLOQ/2, LLOQ/√2 or a
   random draw, and never delete those units; each choice invents data or biases
   the comparison. Above-ULOQ results are `censored: right` at the ULOQ unless the
   sample was re-assayed at a higher dilution (then use the measured value).
2. **Use the right limit.** The limit is sample-specific: LLOQ × dilution factor
   for that sample. Different batches or dilutions can have different limits; keep
   each one. Name the limit's source (assay validation report, kit insert).
3. **Choose the method from the question, before results.**
   - `censored_lognormal`: when a ratio of geometric means is the quantity of
     interest and log values are plausibly normal with similar spread across groups.
     Gives geometric means, ratios versus the reference with intervals, and a
     likelihood-ratio test.
   - `peto_peto`: a rank test that needs no distribution, for left-censored data
     only; it says whether distributions differ but gives no effect size.
   Do not run both and report the smaller p value.
4. **Independent units.** One value per donor, animal or sample. Technical
   replicate wells are averaged by a declared rule first (if any replicate is
   below the limit, say how the assay SOP reports the sample).

## Run and report

Report per group: n, detected, censored fraction and limits, and Kaplan-Meier
quartiles where defined (a median is undefined when more than half the values are
censored; say so rather than quoting the limit). Report the ratio or test with
its method. With more than half the values censored in a group, results depend
heavily on the distributional assumption; with fewer than 3 detected values in a
group the parametric fit is withheld. Respect `must_mention`, withheld states and
failing items in `interpretation_facts.json`.

## Evidence (0.14.0)

See [the 0.14.0 record](../../validation/RELEASE_0.14.0.md): R survival oracle
agreement, the registered calibration (including a retained convergence-check
defect, fixed and rerun post hoc) and the unmet published-example gate.
