# 0.14.1 — closing the 0.14.0 follow-ups

Baseline `52cfef0` (0.14.0 as published). Plan: [plans/RELEASE_0.14.1_PLAN.md](../plans/RELEASE_0.14.1_PLAN.md).
All additions are opt-in; existing default outputs are checked byte for byte at the end.

## Item 1 — paired censored designs (`censored_comparison` with `paired`)

[Skill](../skills/censored-data/SKILL.md), [contract](../skills/censored-data/references/input-and-model.md).
Helsel's censored paired difference (NADA2 `cen_paired`: interval-censored normal ML of
paired differences; original or log scale) and censored sign test (NADA2 `cen_signtest`,
with Fong's tie correction). Prespecified in `a6b9bbe`.

| Gate | Observed scope |
|---|---|
| Established implementation | NADA2 2.0.2 (CRAN 2026-03-23, MIT; tarball SHA-256-pinned, sources sourced at validation time, not installed or redistributed) and R survival 3.8.6 `survreg` refitted at full precision: 10/10 cases (7 synthetic fixtures and NADA2's `atrazine` data on the original scale, log scale and sign test); full-precision differences ≤ 2.6e-7; NADA2's printed Z, p, mean and sign-test p values reproduced exactly at their printed precision. R `survreg` with default control did not converge on the atrazine original-scale fit (warning) but reached the same optimum (log-likelihood −146.178) as this implementation. [Record](censored_paired_benchmarks_0.14.1.json) |
| Published worked example | UNMET: no printed numeric result for a paired censored comparison was located (NADA2's atrazine example prints none in its documentation). |
| Calibration, first registration | Wald z (the NADA2 inference), seed 141101, 2000 draws per scenario, all kept: **2/7 bounded metrics passed (5 misses).** Misses retained: log-scale coverage 0.926 (null) and 0.925 (shift), type I 0.074; original-scale coverage 0.936 and 0.927 with 20 pairs (bounds about 0.940). Original-scale type I 0.064 and sign-test type I 0.030 passed. The size matches the ML SD with z quantiles (about 0.93 without censoring at 20 pairs). [Record](censored_paired_calibration_0.14.1.json) |
| Calibration, second registration | After the first, a small-sample t default was added (ML SE × √(k/(k−1)), t with k − 1 df over k informative pairs; identical to the paired t test without censoring) and a second calibration was committed before running (`c72f235`; seed 141102; same scenarios and bounds): **7/7 bounded metrics passed** — coverage 0.948–0.958, type I 0.042 (log) and 0.049 (original), sign test 0.028. Descriptive: Fong's corrected p was very conservative (type I 0.002; power 0.195 versus 0.418 for the exact binomial); with 10 pairs and 50% censoring 0.801 of datasets were reportable (coverage 0.960); the original-scale method on lognormal data had type I 0.033; script-side L/2 substitution covered 0.945 in this equal-limit design, which does not show substitution bias. [Record](censored_paired_t_calibration_0.14.1.json) |
| Live-agent misuse | 6/6 for the 0.14.1 set (criteria `d37bf7f`, not blinded); paired items: `paired_blq_as_independent` passed, `drop_blq_pairs` passed but weaker than designed (complete-case result shown first, labelled biased, censored analysis recommended as primary). [Record](agent-scenarios/0.14.1/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results; paired runs carry the 0.14.1 evidence block |

The first registration's code is `a6b9bbe`; the inference change is `c72f235`. The
Wald z option remains (`paired.inference: wald_z`) for matching NADA2, with its
calibration miss stated in every such run.

## Item 2 — censored outcomes inside panels (`panel_comparison` with `censoring`)

[Skill](../skills/panel-comparison/SKILL.md), [contract](../skills/panel-comparison/references/input-and-model.md).
Outcomes with any censored value use censored lognormal ML with a small-sample t
correction (independent: two-group Tobit, common SD, N − 2 df; paired: interval-censored
paired log differences, k − 1 df); other outcomes keep Welch or paired t; Holm/BH as in
0.14.0. Prespecified in `52639b3`.

| Gate | Observed scope |
|---|---|
| Established implementation | R survival 3.8.6 `survreg` (interval2, gaussian on ln values) with the same correction, `t.test` and `p.adjust`: 24/24 outcomes across 4 fixtures (12 censored, including right-censored values and a withheld outcome), largest difference 6.1e-11. [Record](panel_censored_benchmarks_0.14.1.json) |
| Published worked example | UNMET: none located for censored outcomes inside a multiplicity-controlled panel. |
| Calibration | Seed 141201, 2000 draws per scenario, all kept: **9/9 bounded metrics passed** — Holm FWER 0.043 (independent) and 0.034 (paired), BH FDR 0.028, censored-outcome coverage 0.948–0.969, simultaneous coverage 0.957–0.964 (15 units per group, 30–40% censoring). Descriptive: with 8 per group simultaneous coverage 0.916 (withheld outcomes count as not covered; one censored outcome reportable in 0.977); with 60% censoring and 10 per group a censored outcome was reportable in 0.697; script-side L/2 substitution covered 0.953 in this equal-limit design, which does not show substitution bias. [Record](panel_censored_calibration_0.14.1.json) |
| Live-agent misuse | `panel_drop_blq_analytes` and `paired_panel_half_lloq` passed (all analytes kept, censored outcomes through the new block). [Record](agent-scenarios/0.14.1/RESULTS.md) |
| Interpretation facts | Extraction-only; panels with censoring carry `validation_evidence_0141` in addition to the 0.14.0 panel evidence |

## Item 3 — fitted GRmax definition (`growth_rate_inhibition`, `fit.grmax_definition`)

[Skill](../skills/growth-rate-inhibition/SKILL.md), [contract](../skills/growth-rate-inhibition/references/input-and-model.md).
Opt-in `fitted_minimum_top_two`: GRmax from the selected fit (sigmoid or flat level) at
the two highest tested concentrations; the authors' observed value is kept as
`GRmax_observed_audit`. Default output unchanged. Prespecified in `b424cf9`.

| Gate | Observed scope |
|---|---|
| Established implementation | Identity tests: the fitted value equals the logistic (or flat level) at the top two concentrations; every other per-curve field equals the default analysis; observed value unchanged. No external implementation of this variant. |
| Published worked example | UNMET for this variant (the authors' toy example applies to their observed definition, reproduced in 0.14.0). |
| Calibration | Seed 141301, 1000 draws per scenario, 0.14.0 data model, both definitions on the same draws: **3/3 bounded rows passed** — fitted GRmax summary coverage 0.958 (cytostatic), 0.940 (cytotoxic; bound 0.936) and 0.956 (flat); the observed definition covered 0.924, 0.943 and 0.931. Mean per-curve error, fitted versus observed: −0.001 vs −0.027, +0.000 vs −0.007, +0.008 vs −0.037. Descriptive: with 0.7 control doublings the fitted value ran high (+0.065; observed −0.064), coverage 0.960 versus 0.921 over 692 estimable summaries; GRinf 0.6 coverage 0.960 versus 0.918; GR50 near the top dose 0.945 versus 0.960. [Record](gr_grmax_calibration_0.14.1.json) |
| Live-agent misuse | `grmax_switch_for_one_drug` and `grmax_extrapolate_top` passed (one definition for both drugs; GRmax only at tested doses). [Record](agent-scenarios/0.14.1/RESULTS.md) |
| Interpretation facts | Extraction-only; runs with the fitted definition carry `validation_evidence_0141` |

The 0.14.0 GRmax miss (0.913) remains the record for the authors' definition; this
calibration's observed-definition rows (0.924, 0.943, 0.931) are a different seed and do
not replace it.

## Corrections before release

The first draft of this record (commit `1099722`) and the CHANGELOG counted the paired
calibrations as 2/8 and 8/8 bounded metrics; the registered scripts have 7 bounded metrics
(the shift scenarios bound coverage only). Correct counts: first registration 2/7 (5 misses),
second 7/7. The evidence module always computed 7. Found when all six live agents quoted 7/7.

## Checks across 0.14.1

- Full suite: `.venv/bin/python -m pytest -q` 937 passed (2026-10-17).
- Unchanged methods: 235 configurations, 1,689 byte-identical scientific artifacts against
  the published 0.14.0 (`52cfef0`), including every 0.14.0 fixture ([record](legacy_fixture_checks_0.14.1.json)).
  No earlier benchmark or calibration was rerun.
