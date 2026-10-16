# 0.13.6 — ROC, rater agreement and outlier flagging

Development baseline: `c48570f` (0.13.5 plus the 0.13.3 live-agent review).
Scope and calibration were prespecified in `ddee73d` (plan
`plans/RELEASE_0.13.6_PLAN.md`, development repository only). Outlier flagging
is recorded in its own section below.

## ROC diagnostic accuracy — `roc_analysis`

[Skill](../skills/diagnostic-accuracy/SKILL.md), [contract](../skills/diagnostic-accuracy/references/input-and-model.md).
Empirical AUC in a declared direction against a declared reference standard,
DeLong logit interval, Mann-Whitney test, paired DeLong comparison, Clopper-Pearson
sensitivity/specificity at declared cut-offs, exploratory Youden only.

## Rater agreement — `rater_agreement`

[Skill](../skills/rater-agreement/SKILL.md), [contract](../skills/rater-agreement/references/input-and-model.md).
Cohen's kappa (unweighted, linear/quadratic or sourced weights) with Fleiss-Cohen-Everitt
variances; McGraw-Wong ICC forms declared before analysis; complete balanced designs.

## Five evidence gates (ROC and agreement)

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R formula oracle (`scripts/agreement_roc_0136_oracle.R`: wilcox.test, explicit DeLong, binom.test, explicit FCE kappa, anova(lm) mean squares and F limits): 19/19 fixtures within 1e-9 (largest difference 2e-14). psych 2.6.5: six Shrout-Fleiss ICCs with limits and three kappa variances agree. [Record](agreement_roc_benchmarks_0.13.6.json). The [initial attempt](agreement_roc_benchmarks_0.13.6_initial.json) is kept: its psych checks failed because the oracle passed proportion tables plus `n.obs`, for which psych returns a wrong variance; base-R checks were already 19/19 |
| Published worked example | Hanley & McNeil (1982) CT ratings AUC 0.893 (their own SE formula gives 0.032, checked script-side; runtime reports DeLong); Shrout & Fleiss (1979) ICC .17/.29/.71/.44/.62/.91; Cohen (1968) kappa .492 and weighted kappa .348 as quoted in the psych documentation (primary article not inspected). 10/10. [Record](agreement_roc_published_0.13.6.json) |
| Calibration | Prespecified in `ddee73d`: seed 136001, 1000 draws per each of 19 scenarios through the package wrappers, all draws kept; 18/19 bounded criteria passed. [Record](agreement_roc_calibration_0.13.6.json) |
| Live-agent misuse | 6/6 passed (criteria `719ec54`, not blinded). Two findings led to Skill wording changes; post-hoc reruns of both passed (one first aborted by an API usage limit, kept). [Record](agent-scenarios/0.13.6/RESULTS.md) |
| Interpretation facts | Extraction-only through the extension workflow; tests check facts equal saved results and that rendering neither recomputes nor changes hashed artifacts |

## Calibration observations retained

- ROC: AUC coverage 0.961 (m=n=50, AUC .75), 0.948 (30/90, .85), 0.958 (100/100, .95),
  0.952 at AUC .5 with Mann-Whitney type I 0.053; paired DeLong difference coverage
  0.966, type I 0.034; Clopper-Pearson sensitivity 0.979 and specificity 0.961.
  Stress m=n=10, AUC .9: 953/1000 reportable (complete separation withholds the
  interval), conditional coverage 0.977.
- Kappa: **miss** — binary kappa n=100 coverage 0.924 below the bound 0.936
  (MCSE 0.0084). Stress rows: binary n=40 0.914, skewed prevalence n=100 0.913.
  Three-category unweighted 0.949 and quadratic 0.940 passed; null kappa coverage
  0.937 and type I 0.052 passed. The large-sample interval is too narrow for binary
  kappa; runtime results say so.
- ICC: A1 n=15 0.952, A1 n=40 0.951, C1 n=15 0.952, one-way n=20 k=2 0.951, A,k
  n=20 0.948; stress n=6 0.962.
- The first execution finished every row but crashed while writing JSON (NumPy
  booleans). Its log is kept as `agreement_roc_calibration_0.13.6_attempt1.log`; after
  changing only the serializer, the same seed reproduced all 19 rows exactly.

## Unchanged methods

193 configurations from earlier fixtures produced 1,473 byte-identical scientific
artifacts against `c48570f` ([record](legacy_fixture_checks_0.13.6.json)). No
earlier benchmark or calibration was rerun.

No Prism, pROC or software-equivalence claim beyond these checks.

## Outlier flagging — `outlier_screen`

[Skill](../skills/outlier-flagging/SKILL.md), [contract](../skills/outlier-flagging/references/input-and-model.md).
Grubbs (single outlier), generalized ESD (declared r) and ROUT (Q, Motulsky & Brown
2006 for a location model) within declared groups. Flags only: nothing is removed,
and results with and without flagged values are reported side by side. Scope and
calibration prespecified in `25414d2`, with `outliers.py` frozen in that commit.

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R formula oracle (`scripts/outliers_0136_oracle.R`; qt, pt, optimize): 14/14 fixture groups give identical flags; Grubbs/ESD statistics within 1e-9; ROUT robust location and RSDR within 1e-7 relative (largest 3e-9, two independent 1-D optimizers). No contributed R outlier package (outliers 0.15 dates from 2011). [Record](outliers_benchmarks_0.13.6.json) |
| Published worked example | NIST/SEMATECH e-Handbook: Grubbs G 2.4687 and one-sided critical value 2.032 (Tietjen & Moore data); generalized ESD ten R_i and λ_i (printed values are truncated, so agreement is within 0.001) and three outliers (Rosner data): 24/24. **ROUT: UNMET** — no published numeric example; Motulsky & Brown's simulation statements are context only. [Record](outliers_published_0.13.6.json) |
| Calibration | Seed 136301, 1000 draws per each of 11 scenarios, all kept: 5/5 bounded false-flag criteria passed (Grubbs n=10 0.050, n=30 0.058; ESD n=30 r=3 0.049; ROUT Q=1% n=20 0.019, n=50 0.024). [Record](outliers_calibration_0.13.6.json) |
| Live-agent misuse | 6/6 passed (criteria `143c088`, not blinded); two scenarios were weaker than designed (no flag occurred in the iterated-Grubbs case; the "skewed" sample was only mildly skewed). [Record](agent-scenarios/0.13.6-outliers/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results; rendering does not recompute or change hashed artifacts |

Descriptive calibration rows retained:
- One +5 SD contaminant, n=20: Grubbs detected it in 0.891 with no other flags.
- Two contaminants (masking): Grubbs detected one in 0.586 and **both in 0.000**;
  generalized ESD (r=3) both in 0.822 (exactly the two in 0.779, any false flag
  0.044); ROUT both in 0.508 (exactly the two 0.497, false flag 0.011).
- Clean heavy-tailed t(3) data, n=30: flagged in 0.482 (Grubbs) and 0.400 (ROUT)
  of datasets. A flag in non-normal data is not evidence of an error.

## Checks across 0.13.6

- Full suite: `.venv/bin/python -m pytest -q` 868 passed (2026-10-16, after all three items).
- Unchanged methods rechecked after the outlier addition: 193 configurations, 1,473
  byte-identical scientific artifacts against `c48570f`
  ([record](legacy_fixture_checks_0.13.6_final.json)); the earlier record taken
  after item 2 is kept.
- Clean build from `git archive` of `36385f7`: 108 Python source files in the wheel are
  byte-identical to `src/` (wheel SHA-256 `39d01f4789360dfa…`; [record](archive_build_0.13.6.json)).
  The first attempt failed offline because the dedicated cache used for 0.13.5 no longer
  existed; the build then used the default local uv cache, still offline.
- Isolated wheel-only install with locked dependencies: 42 configurations and 287
  scientific artifacts identical to the source; six new-type runs re-rendered and
  verified with unchanged scientific artifacts ([record](clean_install_0.13.6.json)).
  Same macOS host only; not a cross-platform claim.
