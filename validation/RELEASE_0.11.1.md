# Release evidence — local development 0.11.1

Base: dev/0.11.0 at 1a12574. Local development; no push/publication.

## Phase 0: sources and licences (before numerical implementation)

- R datasets ToothGrowth (60 guinea pigs, Crampton 1947 / Bliss 1952) and
  warpbreaks (54 looms, Tippett 1950 p106) have numeric tables in the installed
  datasets package. Base R distributes datasets under GPL-2/GPL-3; data will
  retain attribution and the distribution licence, with original table hashes.
  https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/ToothGrowth.html
  https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/warpbreaks.html
  https://www.r-project.org/doc/R-FDA.pdf
- mcr's creatinine table and documented analysis are available from CRAN;
  package licence GPL >=3. Package download and actual numeric data verification
  are in progress. https://search.r-project.org/CRAN/refmans/mcr/html/mcreg.html
- Bland & Altman 1986 PEF and Agresti tables: no separate usable data licence
  verified at this stage. Published worked-example gates remain UNMET until
  a matched usable numeric example has actually been reproduced.
- Synthetic fixtures are for implementation/calibration, never published gates.

## Implemented scope and evidence

Python runtime: independent SS II/III two-way ANOVA, equal-weight EMMs with
Tukey, Dunnett, Sidak or Holm families; conditional Fisher OR/exact interval,
RxC fixed-margin Monte Carlo, chi-square, Newcombe RD/Katz RR, trend, McNemar,
Wilson/Clopper-Pearson; Pearson/Spearman and OLS/WLS; Deming jackknife,
Passing-Bablok with mcr tangent/tie convention and CUSUM, Bland-Altman exact
normal-quantile intervals. R is validation-only for new methods.

Cross-implementation: 66 checks passed against car 3.1.5, emmeans 2.0.4,
mcr 1.3.3.1 and base R. Includes ToothGrowth, warpbreaks, complete-pair
creatinine and synthetic fixtures. mcr creatinine has 110 original rows, with
missing pairs explicitly removed in a derived table; original rows, hashes,
exclusion IDs and GPL licences preserved. Read the source manifest for counts.
Dunnett stochastic integration differed by 2.80e-8 absolute (4.47e-5 relative);
the first comparison using 2e-5 relative alone failed. Final comparison records
a 1e-7 absolute integration tolerance separately. RxC Monte Carlo differs from
R's exact p by 11.43% relative, within four MCSE plus Monte Carlo resolution;
this is a statistical check, not deterministic byte equality. See benchmark JSON.

Published gates: public factorial worked datasets and the documented mcr
creatinine comparison reproduced. The same creatinine data exercise correlation,
linear regression and normal LoA numerics, but no separate published matched
correlation/LoA result is claimed. Bland-Altman 1986 PEF and Agresti table gates
remain UNMET (licence/data not verified). No synthetic example substitutes.

## Registered calibration

Preregistered commit efe03df; 32 scenarios × 1000 draws = 32000 simulations.
Original run completed numerics but failed JSON serialization of NumPy scalar
counts. The complete stdout/error log is preserved. c675a26 fixes serialization
only; same seeds, scenarios, bounds and algorithms were rerun. No retuning.
Both all-evaluable and reportable-only denominators are registered; first-run
printed outcomes match the final JSON. All draws are retained.

| Scenario | Metric | Population | Rate | MCSE | Bound |
|---|---|---|---:|---:|---:|
| anova_2_heteroscedastic | C(a, Sum) | all_evaluable | 0.2160 | 0.0130 | 0.0638 |
| anova_2_heteroscedastic | C(a, Sum) | reportable_only | 0.2160 | 0.0130 | 0.0638 |
| anova_2_heteroscedastic | C(a, Sum):C(b, Sum) | all_evaluable | 0.2840 | 0.0143 | 0.0638 |
| anova_2_heteroscedastic | C(a, Sum):C(b, Sum) | reportable_only | 0.2840 | 0.0143 | 0.0638 |
| anova_2_heteroscedastic | C(b, Sum) | all_evaluable | 0.1880 | 0.0124 | 0.0638 |
| anova_2_heteroscedastic | C(b, Sum) | reportable_only | 0.1880 | 0.0124 | 0.0638 |
| anova_3_equal | C(a, Sum):C(b, Sum) | all_evaluable | 0.0730 | 0.0082 | 0.0638 |
| anova_3_equal | C(a, Sum):C(b, Sum) | reportable_only | 0.0730 | 0.0082 | 0.0638 |
| anova_3_heteroscedastic | C(a, Sum) | all_evaluable | 0.2190 | 0.0131 | 0.0638 |
| anova_3_heteroscedastic | C(a, Sum) | reportable_only | 0.2190 | 0.0131 | 0.0638 |
| anova_3_heteroscedastic | C(a, Sum):C(b, Sum) | all_evaluable | 0.2650 | 0.0140 | 0.0638 |
| anova_3_heteroscedastic | C(a, Sum):C(b, Sum) | reportable_only | 0.2650 | 0.0140 | 0.0638 |
| anova_3_heteroscedastic | C(b, Sum) | all_evaluable | 0.2590 | 0.0139 | 0.0638 |
| anova_3_heteroscedastic | C(b, Sum) | reportable_only | 0.2590 | 0.0139 | 0.0638 |
| wilson_n20_p0.05 | coverage | all_evaluable | 0.9020 | 0.0094 | 0.9362 |
| wilson_n20_p0.05 | coverage | reportable_only | 0.9020 | 0.0094 | 0.9362 |
| wilson_n100_p0.5 | coverage | all_evaluable | 0.9340 | 0.0079 | 0.9362 |
| wilson_n100_p0.5 | coverage | reportable_only | 0.9340 | 0.0079 | 0.9362 |
| loa_n30_skew | limit_0_coverage | all_evaluable | 0.0150 | 0.0038 | 0.9362 |
| loa_n30_skew | limit_0_coverage | reportable_only | 0.0150 | 0.0038 | 0.9362 |
| loa_n30_skew | limit_1_coverage | all_evaluable | 0.4300 | 0.0157 | 0.9362 |
| loa_n30_skew | limit_1_coverage | reportable_only | 0.4300 | 0.0157 | 0.9362 |
| loa_n100_skew | limit_0_coverage | all_evaluable | 0.0000 | 0.0000 | 0.9362 |
| loa_n100_skew | limit_0_coverage | reportable_only | 0.0000 | 0.0000 | 0.9362 |
| loa_n100_skew | limit_1_coverage | all_evaluable | 0.1780 | 0.0121 | 0.9362 |
| loa_n100_skew | limit_1_coverage | reportable_only | 0.1780 | 0.0121 | 0.9362 |

Every failed criterion above is copied into evidence_0111.py and relevant
results/facts. Heteroscedastic ANOVA and skewed normal-LoA rows are stress
failures. Wilson at sparse p and one nominal ANOVA interaction row also missed.
Fisher observed sizes are conservative; no promise of exact 5% size.
Correlation/regression, RR, trend and McNemar do not have standalone calibration
rows; their numerical checks do not establish coverage across designs.

## Live-agent and facts gates

Six live-agent scenarios PENDING, per user instruction. No Claude CLI or
substitute agents executed. Exact command in
[RESULTS](agent-scenarios/0.11.1/RESULTS.md).

All four new analysis types save source-linked facts. Older two/multi-group
runs add only facts; no refit or changed scientific tables. State tests cover
withheld nonlinearity, failed estimates, sparse counts and interaction.

## Self-review

- **Valency/avidity:** these are statistical measurements, not affinity models.
  Correlation/regression never re-labels a measured readout as intrinsic KD.
- **Readout:** x/y must be the declared measurements; association is not
  agreement. Bland-Altman assumes normal, constant-spread differences;
  proportional trend and normality assumptions are disclosed, not repaired.
- **Common scale:** method-comparison gate requires same_scale=true; Deming
  explicitly records Var(error Y)/Var(error X), reciprocal to mcr's convention.
- **Independence:** duplicate subject IDs fail; paired methods have one row
  per independent subject. Technical wells need upstream justified aggregation.
- **Design range:** every factorial cell needs >=2 units; empty/rank-deficient
  designs fail. No extrapolated agreement. Passing-Bablok requires positive
  association, nonnegative data and finite rank interval; CUSUM failure withholds.
- **Interactions/multiplicity:** interaction first; marginal contrasts withheld
  when material. Families are per stratum, not an undeclared global family;
  Holm intervals use conservative Bonferroni bounds and say so.
- **Sparse categorical data:** expected counts and Monte Carlo SE retained;
  zero-event RR intervals withheld with no automatic correction.

## Deviations and limits

- Larger Fisher tables use seeded Monte Carlo, an explicitly allowed option.
- Spearman is exact only for untied n<=9; larger/tied samples use base R's
  exact=FALSE convention, not its AS89 expansion.
- Known assumptions have disclosed stress misses rather than altered estimators.
- Manifest version/code/time necessarily changes; prior scientific files stay
  byte-identical. HTML includes current metadata, not a claimed legacy byte match.
- New methods are Python; no R runtime dependency was added.
- R compilation required MACOSX_DEPLOYMENT_TARGET=14.0 to work around the host
  Fortran toolchain's unsupported default. Packages stay in ignored .r-lib.
- Initial calibration serialization failure and Dunnett integration tolerance
  adjustment are disclosed above; no calibration decision rule changed.

## Release checks

Final test, legacy, clean archive, wheel-only smoke and exact-install results
are recorded in release_checks_0.11.1.json after execution. No unchanged-method
benchmarks/calibrations rerun. No runs/, build/ or historical evidence modified.

Final observed checks: **585 passed**, eight existing warnings (223.42 s);
three changed Skills validate. Legacy **110 configs / 922 prior artifacts**
byte-identical to 1a12574. Clean-environment passed (26 new configurations
exact). Separate wheel-only smoke **136 configs** passed; its final stdout JSON
is saved unchanged. Broader install-exact check **30 configs / 170 artifacts**
includes older two/multi-group facts. All 71 Python files match wheel bytes.
Wheel built from clean archive 1d7cec8. A later clean-archive sdist from 06d2fe2
corrects README Skill count to 18; original sdist preserved in dist/superseded,
runtime unchanged. Hashes are in release_checks_0.11.1.json. No native Prism
comparison or cross-platform claim. No live-agent execution by implementer.
