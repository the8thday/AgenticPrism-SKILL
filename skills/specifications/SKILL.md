---
name: specifications
description: Summarize the variability of a quality attribute across independent lots for specification work - normal (exact or Howe) and nonparametric tolerance intervals with declared content and confidence, and process capability (Pp/Ppk overall, Cp/Cpk within subgroups) with confidence intervals against declared limits. Tolerance intervals describe data; they do not set specifications.
---

# Tolerance intervals and process capability

Follow [runtime discovery](../agentic-prism/references/runtime.md), then read the
[input contract](references/input-contract.md). Introduced in AgenticPrism 0.10.1;
`analysis_type=specification`, one attribute per run, a new output directory.

## Gate before analysis

- **Independent units:** normally one release value per lot. Repeated measurements
  of the same lot are not independent units; use subgroups only for within-lot
  capability (Cpk) and say what they are.
- **Content and confidence** (for example 99% of lots with 95% confidence) and the
  choice of normal versus nonparametric come from the user's protocol. Do not switch
  methods after seeing the result. A failed normality check is a reason to discuss
  the protocol, not to pick whichever interval looks better.
- **Specification limits** for capability come from the user with a source.

## Run and interpret

Run `analyze`, `verify`, read `interpretation_facts.json`. Report:

- The tolerance interval with content, confidence, sides, method and n. It says where
  the stated proportion of lots from this process lies, given these data. It is not
  a specification: specifications also need clinical relevance, stability, assay
  variability and regulatory input.
- Nonparametric intervals need many lots: 93 for 95% content with 95% confidence,
  46 for 90%/95% using the extremes. When too few, the tool reports the minimum n
  instead of an interval; relay it.
- Capability: Ppk uses the overall SD; Cpk uses the within-subgroup SD and excludes
  between-lot variation, so Cpk can look much better than Ppk. Report which one
  answers the question, with its interval. Below n = 25 the Cpk/Ppk interval is an
  approximation that NIST recommends against; say so.
- The expected nonconforming fraction is a normal plug-in estimate, not an observed rate.

Mention every `must_mention` item.

## Numerical evidence (0.10.1)

Exact two-sided factors equal R `tolerance` (method EXACT) over 45 combinations of
n, content and confidence (max relative difference 3.3e-9); one-sided factors,
order-statistic ranks, TOST and capability formulas match R; the NIST/SEMATECH
handbook's printed values (k = 2.217 and 1.8740 for n = 43; 46 and 473 lots for
nonparametric intervals; Cp 1.0, Cpk 0.6667) are reproduced. Calibration (1000
datasets per row): tolerance intervals attained 91.7–98.5% confidence against
90–95% targets. The n = 10 Pp and Ppk rows missed their bound (93.6% vs 93.62%); a
later 100,000-dataset check gave 94.9% and 95.4%, so the misses reflect Monte Carlo
error, but they remain recorded as misses. See [the 0.10.1 record](../../validation/RELEASE_0.10.1.md).
