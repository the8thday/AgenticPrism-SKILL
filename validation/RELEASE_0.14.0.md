# 0.14.0 — discovery-stage batch

Development baseline: `16613d8` (0.13.6 as published). Plan `plans/RELEASE_0.14.0_PLAN.md`
(development repository only). Items are recorded in the order they were built.

## Item 1 — values below/above quantification limits: `censored_comparison`

[Skill](../skills/censored-data/SKILL.md), [contract](../skills/censored-data/references/input-and-model.md).
Censored rows are kept at their sample-specific limit; no substitution or deletion.
Censored lognormal maximum likelihood (geometric means, ratios versus the reference,
likelihood-ratio test) or the Peto-Peto k-sample rank test for left-censored data.
Scope and calibration were prespecified in `1b051a1`, with `censored.py` frozen there.

| Gate | Observed scope |
|---|---|
| Established implementation | R survival 3.8.6 oracle (`scripts/censored_0140_oracle.R`: `survreg` gaussian on log values with `interval2` coding; `survdiff(rho = 1)` on flipped values): 6/6 fixtures, largest difference 6.3e-11, including left-plus-right censoring, two limits, three groups and a withheld case. No NADA or Prism comparison. [Record](censored_benchmarks_0.14.0.json) |
| Published worked example | **UNMET**: Helsel's CuZn example has code but no printed numeric result in reachable sources; no other printed censored comparison was located |
| Calibration | Seed 140101, 1000 draws per each of 12 scenarios, every draw kept. Registered run: 10/10 bounded criteria passed among reportable fits, but a too-strict convergence check wrongly withheld converged fits (n=30 30% censoring 0.138; n=30 60% 0.067; alternative 0.034; two limits 0.023; gamma stress 0.122). [Record](censored_calibration_0.14.0.json) |
| Calibration, post hoc | After the fix (analytic gradient, gradient-based convergence; estimates unchanged, still equal to survreg) the same script and seed reported all of those fits: coverage 0.940–0.955, LR type I 0.051–0.056, Peto-Peto type I 0.038 and 0.056; 10/10 bounded criteria. Supplements, does not replace, the registered run. [Record](censored_calibration_0.14.0_posthoc_convergence_fix.json) |
| Live-agent misuse | 6/6 passed (criteria `811a677`, not blinded). The generator mislabelled three files' column names (my error); agents flagged it. [Record](agent-scenarios/0.14.0-censored/RESULTS.md) |
| Interpretation facts | Extraction-only through the extension workflow; facts equal saved results; rendering neither recomputes nor changes hashed artifacts |

Descriptive rows retained: 80% censoring with n=12 left only 0.187 of datasets with at
least 3 detected values per group (the rest withheld); a non-lognormal (gamma) stress
kept coverage 0.944 and type I 0.054; script-side LLOQ/2 substitution with equal
censoring in both groups kept ratio coverage (0.949, 0.936), so this design does not
show substitution bias and is not evidence for it. Follow-up: paired censored designs.

## Item 2 — multi-candidate potency ranking: `potency_ranking`

[Skill](../skills/potency-ranking/SKILL.md), [contract](../skills/potency-ranking/references/input-and-model.md).
Paired log10 EC50 differences against a benchmark from the same independent
experiments; relative potency with Holm p values and Bonferroni simultaneous
intervals; optional prespecified fold-threshold tiers; non-estimable runs withhold the
candidate. Prespecified in `0bdb35c` with `potency_ranking.py` frozen there.

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R oracle (`scripts/potency_ranking_0140_oracle.R`: paired `t.test`, `p.adjust(holm)`, Bonferroni-level `t.test` intervals): 4/4 fixtures, largest difference 3.2e-14, including identical withholding. [Record](potency_ranking_benchmarks_0.14.0.json) |
| Published worked example | **UNMET**: no published multi-candidate ranking with printed numbers was located |
| Calibration | Seed 140201, 1000 draws per each of 8 scenarios, all kept: 11/11 bounded criteria passed (Holm family-wise error 0.041–0.057; simultaneous coverage 0.943–0.959; no null candidate classified beyond a 3-fold threshold). [Record](potency_ranking_calibration_0.14.0.json) |
| Live-agent misuse | 6/6 passed (criteria `8a62f3a`, not blinded). [Record](agent-scenarios/0.14.0-potency/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results; rendering neither recomputes nor changes hashed artifacts |

Descriptive rows retained: with 10 candidates and 4 experiments a true 3-fold candidate
reached Holm significance in 0.243 of cases, and a true 10-fold candidate was classified
at least 3-fold in 0.275; unadjusted p < .05 across 10 null candidates gave family-wise
error 0.312; heavy-tailed residuals kept FWER 0.036 and coverage 0.964; candidates whose
EC50 exceeded the top dose in some runs were withheld in 0.932 of cases.

## Item 3 — growth-rate inhibition metrics: `growth_rate_inhibition`

[Skill](../skills/growth-rate-inhibition/SKILL.md), [contract](../skills/growth-rate-inhibition/references/input-and-model.md).
GR values (Hafner et al., Nat Methods 2016) and per-curve GR50/GEC50/GRinf/GRmax/GR_AOC
following the authors' reference code (bounded sigmoid, flat-fit F test), no
extrapolation of GR50, summaries across declared independent experiments. Scope and
calibration prespecified in `2d0727f`, with `growth_rate.py` frozen after a pre-freeze
check against the authors' example (two pre-freeze fixes: flat level = mean as in R
GRmetrics; optimizer status flags not used).

| Gate | Observed scope |
|---|---|
| Established implementation / published worked example | The authors' toy example (`datarail/gr_metrics` commit 9c18132, files SHA-256-pinned and downloaded at validation time because the repository has no license) against their R GRmetrics output: 112/112 curves with the same fit type; GRmax and GR_AOC within 4.2e-15; GEC50, GRinf, h and GR50 within 6.2e-4 relative; this implementation's residual sum of squares was no larger than R's in every sigmoid fit. [Record](gr_published_0.14.0.json) |
| Calibration | Seed 140301, 1000 draws per each of 7 scenarios, all kept: 6/7 bounded criteria passed. **Miss retained:** cytostatic GRmax summary coverage 0.913 below 0.936 — GRmax is the minimum observed GR at the two highest doses (reference definition) and is biased low under noise. False sigmoid-fit rate with no drug effect 0.018; GR50 summary coverage 0.953 and 0.940; GR_AOC 0.947 and 0.943. [Record](gr_calibration_0.14.0.json) |
| Calibration, post hoc | The registered division-rate row averaged infinite audit values (mean −inf, kept). A supplement on the same seed using only estimated values: with identical GR truth, a 3-doubling versus a 1-doubling line gave GR50 ratio 1.07 but relative-count IC50 ratio 0.43. [Record](gr_division_rate_supplement_0.14.0.json) |
| Live-agent misuse | 6/6 passed (criteria `8f475cb`, not blinded). Finding: a control that grew 1.03-fold passed the original gate and produced GR up to +29; the agent overrode it. [Record](agent-scenarios/0.14.0-gr/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results; `gr_values.csv` saved before the manifest |

Post-freeze changes (neither changes a fitted value): non-finite GR50 audit values are
stored as text for JSON; the control-growth gate was raised from "did not grow" to
"less than half a doubling" after the live-agent finding. The published-example check
was rerun after the gate change and is identical; the calibration was not rerun (its
lowest-growth scenario, 0.7 doublings, is above the new threshold). Other descriptive
rows: GRinf 0.6 withheld GR50 in 0.914 of curves; with 0.7 control doublings GR50 was
withheld in 0.230 and GRmax coverage was 0.926.

## Item 4 — multi-outcome panels with multiplicity control: `panel_comparison`

[Skill](../skills/panel-comparison/SKILL.md), [contract](../skills/panel-comparison/references/input-and-model.md).
Two groups across a panel; every measured outcome in exactly one declared family;
Welch or paired t on a declared scale; Holm (with Bonferroni simultaneous intervals)
or Benjamini-Hochberg within each family. Prespecified in `b159608` with `panel.py`
frozen there.

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R oracle (`scripts/panel_0140_oracle.R`: Welch/paired `t.test`, `p.adjust`): 3/3 fixtures, largest difference 6.7e-15. [Record](panel_benchmarks_0.14.0.json) |
| Published worked example | Benjamini & Hochberg (1995) 15 p values: BH rejects 4 as printed (Holm 3); adjusted values equal to R. This checks the procedure; per-outcome panel tests have no published example (UNMET) |
| Calibration | Seed 140401, 1000 draws per each of 8 scenarios, all kept: 9/9 bounded criteria passed (Holm family-wise error 0.029–0.036, conservative under positive correlation; simultaneous coverage 0.964–0.971; BH FDR 0.030). [Record](panel_calibration_0.14.0.json) |
| Live-agent misuse | 6/6 passed (criteria `369ae9a`, not blinded). [Record](agent-scenarios/0.14.0-panel/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results |

Descriptive rows retained: BH had at least one false discovery in 0.095 of mixed panels
(it controls the proportion, not any error); per-effect power BH 0.631 versus Holm 0.460;
unadjusted p < .05 across 20 null outcomes gave family-wise error 0.442; under strong
dependence (rho 0.9) BH FDR was 0.020. Follow-up: censored outcomes inside the panel.

## Checks across 0.14.0

- Full suite: `.venv/bin/python -m pytest -q` 911 passed (2026-10-17, after all four items).
- Unchanged methods: 220 configurations, 1,612 byte-identical scientific artifacts against
  the published 0.13.6 (`16613d8`) ([record](legacy_fixture_checks_0.14.0.json)). No earlier
  benchmark or calibration was rerun.
- Clean build from `git archive` of `dadb8a5`: 113 Python source files in the wheel are
  byte-identical to `src/` (wheel SHA-256 `992b135c95d2007c…`; [record](archive_build_0.14.0.json)).
- Isolated wheel-only install with locked dependencies: 25 configurations and 135
  scientific artifacts identical to the source (the four 0.14.0 fixture folders plus
  `outliers_0136` and `groups_synthetic`); six new-type runs re-rendered and verified
  with unchanged scientific artifacts ([record](clean_install_0.14.0.json)).
  Same macOS host only; not a cross-platform claim.

Follow-ups recorded by this batch: paired censored comparisons; censored outcomes inside
panel comparisons; a GRmax definition or interval that is not biased low under noise.
