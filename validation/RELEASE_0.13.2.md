# Release 0.13.2 — nested comparisons, ANCOVA, curve AUC, standard curves and qPCR

Local development release on `dev/0.13.2` from `dev/0.13.1`. Five new analysis types through the
extension workflow; earlier scientific outputs are byte-identical (`legacy_0.13.2.json`).

| Gate | nested_comparison | ancova | curve_auc | standard_curve | qpcr_relative |
|---|---|---|---|---|---|
| Established implementation | lmerTest (lmer REML, contest Satterthwaite) | R lm, emmeans, car | independent base-R trapezoid, t.test | R lm, investr invest | independent base-R ΔΔCq, t.test (pcr package unavailable) |
| Published worked example | UNMET | UNMET | UNMET | UNMET | UNMET (Livak & Schmittgen 2001 full text not retrievable from PMC) |
| Calibration | 5/5 bounded rows passed | 4/4 | 3/3 | 2/2 | 3/3 |
| Live-agent misuse | Passed (2 scenarios) | Passed (2) | Passed (2) | Passed (2) | Passed (2) |
| Interpretation facts | `evidence_0132.py` | same | same | same | same |

Every calibration script was committed before it ran (nested e32587c, ANCOVA 3423d22, AUC 48fc2fe,
standard curve and qPCR de18c79).

## Development findings retained

- The nested omnibus F matches lmerTest only with the same hypothesis basis (`contest`); `anova(type = 3)`
  uses another basis and its Fai–Cornelius denominator df differ by about 2e-4 relative.
- ANCOVA fixture seed 13202 rejected slope homogeneity by chance (p = 0.047); it is kept as
  `config_slope_flag.json` and the main fixture uses seed 13203.
- A test caught that AUC withholding fractions were NaN for groups without withheld units, silently
  disabling the unequal-withholding flag; fixed before any benchmark or calibration record.
- Standard-curve inverse intervals first disagreed with investr because the unknown's replicate variance
  was not pooled; the investr (Graybill) convention was adopted before any record.

## Cross-implementation

**Nested** (`nested_benchmarks_0.13.2.json`)

| Dataset | Max relative difference | Worst field | Passed |
|---|---|---|---|
| balanced_2groups | 2.3e-08 | sigma2 | True |
| unbalanced_2groups | 8.6e-08 | trt - ctrl p | True |
| unbalanced_3groups_vs_control | 4.2e-09 | sigma2 | True |
| unbalanced_4groups_all_pairs | 2.8e-08 | tau | True |
| fixture | 1.7e-07 | mAb-high - vehicle p | True |

**ANCOVA** (`ancova_benchmarks_0.13.2.json`)

| Dataset | Max relative difference | Worst field | Passed |
|---|---|---|---|
| two_groups_one_covariate | 4.2e-12 | homog_F | True |
| three_groups_unequal_n | 9.1e-13 | group_p | True |
| four_groups_two_covariates | 1.2e-13 | b - c p | True |
| heterogeneous_slopes | 4.2e-14 | homog_p | True |
| fixture | 1.9e-14 | mAb-low - vehicle p | True |

**Curve AUC** (`curve_auc_benchmarks_0.13.2.json`)

| Dataset | Max relative difference | Worst field | Passed |
|---|---|---|---|
| regular_no_baseline | 3.6e-14 | t - c p | True |
| irregular_first_value | 1.5e-14 | t - c p | True |
| declared_constant | 2.9e-15 | auc | True |
| fixture_with_dropout | 1.8e-13 | TCE-low - isotype p | True |

**Standard curve and qPCR** (`standard_curve_qpcr_benchmarks_0.13.2.json`)

| Dataset | Max relative difference | Worst field | Passed |
|---|---|---|---|
| standard_curve linear_m3 | 5.1e-08 | U0 upper | True |
| standard_curve linear_m1 | 2.3e-08 | U1 lower | True |
| standard_curve quadratic_m2 | 1.5e-08 | U2 est | True |
| standard_curve fixture_plate_P1 | 2.7e-08 | P1-lysate1 upper | True |
| qpcr_relative efficiency_100_one_ref | 4.9e-14 | p | True |
| qpcr_relative declared_efficiencies_two_refs | 1.9e-13 | p | True |
| qpcr_relative fixture (partial undetermined, NTC) | 5.2e-14 | IL6 p | True |

## Calibration

**Nested** (seed 131262001)

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| two_balanced | coverage | 1000 | 0.9480 | 0.0070 | 0.9362 | True |
| two_balanced | fraction_mixed_model_primary | 1000 | 0.9570 | — | — | descriptive |
| two_unbalanced | coverage | 1000 | 0.9460 | 0.0071 | 0.9362 | True |
| two_unbalanced | fraction_mixed_model_primary | 1000 | 0.9360 | — | — | descriptive |
| low_icc | coverage | 1000 | 0.9560 | 0.0065 | 0.9362 | True |
| low_icc | fraction_mixed_model_primary | 1000 | 0.5660 | — | — | descriptive |
| three_null | type_I_error_primary_omnibus | 1000 | 0.0410 | 0.0063 | 0.0638 | True |
| three_null | type_I_error_pseudoreplicated_anova | 1000 | 0.2580 | 0.0138 | — | descriptive |
| three_null | fraction_mixed_model_primary | 1000 | 0.9550 | — | — | descriptive |
| three_vs_control | simultaneous_coverage | 1000 | 0.9630 | 0.0060 | 0.9362 | True |
| three_vs_control | fraction_mixed_model_primary | 1000 | 0.9840 | — | — | descriptive |

**ANCOVA** (seed 131262002)

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| two_groups | coverage | 1889 | 0.9471 | 0.0052 | 0.9400 | True |
| two_groups | power_adjusted | 1889 | 0.3774 | 0.0112 | — | descriptive |
| two_groups | power_unadjusted | 2000 | 0.2085 | 0.0091 | — | descriptive |
| three_vs_control | simultaneous_coverage | 1888 | 0.9544 | 0.0048 | 0.9400 | True |
| null_three | type_I_error_adjusted_F | 1880 | 0.0484 | 0.0049 | 0.0601 | True |
| slope_test_null | slope_test_withholding_rate | 2000 | 0.0440 | 0.0046 | 0.0597 | True |

**Curve AUC** (seed 131262003)

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| two_hetero | coverage | 2000 | 0.9490 | 0.0049 | 0.9403 | True |
| three_null | familywise_error_holm | 2000 | 0.0390 | 0.0043 | 0.0597 | True |
| three_vs_control | simultaneous_coverage | 2000 | 0.9595 | 0.0044 | 0.9403 | True |
| informative_dropout | coverage_under_informative_dropout | 1956 | 0.8211 | 0.0087 | — | descriptive |
| informative_dropout | mean_bias | 1956 | -343.4818 | — | — | descriptive |

**Standard curve and qPCR** (seed 131262004)

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| sc_linear_m3 | coverage_per_unknown | 6000 | 0.9482 | 0.0029 | 0.9444 | True |
| sc_quadratic_m2 | coverage_per_unknown | 6000 | 0.9468 | 0.0029 | 0.9444 | True |
| qpcr_independent | coverage | 2000 | 0.9615 | 0.0043 | 0.9403 | True |
| qpcr_paired | coverage | 2000 | 0.9505 | 0.0049 | 0.9403 | True |
| qpcr_null | type_I_error | 2000 | 0.0470 | 0.0047 | 0.0597 | True |
| qpcr_efficiency_misdeclared | coverage | 2000 | 0.8675 | 0.0076 | — | descriptive |
| qpcr_efficiency_misdeclared | mean_bias_log2 | 2000 | 0.1896 | — | — | descriptive |

Descriptive rows show what each module guards against: a pseudo-replicated ANOVA rejected a true null
26% of the time; adjusting for baseline nearly doubled power (0.377 vs 0.208); outcome-related dropout
dropped AUC interval coverage to 0.821; assuming 100% efficiency when the target's is 85% dropped qPCR
coverage to 0.868 with a +0.19 log2 bias.

## Not done

Published worked examples, wheel build and clean-install checks. Live-agent review: 10/10 passed
([results](agent-scenarios/0.13.2/RESULTS.md), implementer-scored); it found that curve AUC lacks a paired design.
