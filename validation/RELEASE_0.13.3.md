# Release 0.13.3 — R-package audit and common nonlinear models

Local development release on `dev/0.13.3` from `dev/0.13.2`.

## R-package audit

Prompted by the user's guidance that long-unmaintained or short-lived R packages should not be trusted
by default, every R package used as an oracle was checked for its latest CRAN publication and
maintainer. These installed versions are maintained with recent releases: survival 3.8-6, lme4 2.0-6,
lmerTest 3.2-1, pbkrtest 0.5.5, emmeans 2.0.4, car 3.1-5, nlme 3.1-169 (R Core), mmrm 0.3.18,
afex 1.5-1, VCA 1.5.2, tolerance 3.0.0, PowerTOST 1.5-7, cmprsk 2.2-12 (2024) and synergyfinder 3.20.0
(Bioconductor). Three stand out:

| Package | Latest CRAN release | Earlier use | Re-check with base R only |
|---|---|---|---|
| drc 3.0-1 | 2016-08-30 | 0.13.1 5PL benchmark and ELISA 5PL check | `nls` (port) on the same bounded problem |
| pwr 1.3-0 | 2020-03-17 | sample-size ANOVA and arcsine power | `pf`/`qf` noncentral F, `pnorm` |
| investr 1.4.2 | 2022-03-31 | standard-curve inversion intervals | independent `uniroot` implementation |

Results (`base_r_rechecks_0.13.3.json`):

| Item | Dataset | Max difference | Passed |
|---|---|---|---|
| one-way ANOVA power (previously pwr) | 18 cases | 9.5e-10 | True |
| arcsine two-proportion power (previously pwr) | 72 cases | 7.8e-16 | True |
| inversion interval (previously investr) | linear_m3 | 2.7e-15 | True |
| inversion interval (previously investr) | linear_m1 | 2.9e-15 | True |
| inversion interval (previously investr) | quadratic_m2 | 9.8e-13 | True |
| inversion interval (previously investr) | fixture_plate_P1 | 6.8e-13 | True |

5PL: 37/43 curves agree with bounded `nls`; the remaining
6 are `nls` convergence failures at the asymmetry bound on fits the package already withholds, with no
numerical disagreement. The initial re-check (`base_r_rechecks_0.13.3_initial.json`) refitted without
bounds: 22/43 agreed; of the 21 that did not, 6 were `nls` failures and the rest were fits where the
package estimate sits at a parameter bound and unconstrained `nls` moved past it. The comparison was
changed to the same bounded problem after that run.
The nonlinear library below uses base R only.

## Nonlinear model library (`nonlinear_fit`, Skill `nonlinear-models`)

| Gate | Status |
|---|---|
| Established implementation | base R `nls` (R-side starts) and `confint`: 18/18 datasets |
| Published worked example | UNMET |
| Calibration | 11/11 bounded rows passed (script committed in bb469b5 before running) |
| Live-agent misuse | Not run |
| Interpretation facts | `evidence_0133.py` |

Benchmark (`nonlinear_benchmarks_0.13.3.json`; the initial run, which scaled interval differences by
the end value and failed a dataset when R confint failed, is kept as `_initial`):

| Dataset | SSE rel. diff | Max estimate rel. diff | Max interval-end diff / width | Passed |
|---|---|---|---|---|
| one_phase_decay_seed1 | 3.0e-16 | 1.2e-08 | 1.2e-05 | True |
| one_phase_decay_seed2 | -7.3e-16 | 1.3e-09 | 1.8e-05 | True |
| one_phase_association_seed1 | -2.0e-15 | 6.7e-08 | 2.1e-05 | True |
| one_phase_association_seed2 | 3.8e-16 | 3.7e-08 | 3.0e-05 | True |
| two_phase_decay_seed1 | -4.9e-16 | 1.6e-08 | 1.5e-04 | True |
| two_phase_decay_seed2 | -1.6e-16 | 2.6e-08 | 1.2e-04 | True |
| exponential_growth_seed1 | 3.5e-16 | 1.4e-09 | 2.0e-05 | True |
| exponential_growth_seed2 | -1.9e-16 | 1.4e-09 | 1.8e-05 | True |
| michaelis_menten_seed1 | 1.8e-15 | 6.8e-10 | 3.3e-05 | True |
| michaelis_menten_seed2 | -4.8e-15 | 2.3e-09 | 3.8e-05 | True |
| logistic_growth_seed1 | -4.8e-15 | 2.7e-10 | 1.3e-05 | True |
| logistic_growth_seed2 | -1.2e-15 | 2.4e-09 | 1.5e-05 | True |
| fixture mAb-X-E1-r1 | -1.5e-16 | 2.4e-08 | 9.5e-05 | True |
| fixture mAb-X-E1-r2 | -4.4e-14 | 1.3e-07 | 1.5e-04 | True |
| fixture mAb-X-E2-r1 | 3.3e-15 | 4.1e-08 | 1.1e-04 | True |
| fixture mAb-X-E2-r2 | -1.4e-14 | 2.4e-08 | 1.4e-04 | True |
| fixture mAb-X-E3-r1 | -3.3e-15 | 1.7e-08 | R confint failed | True |
| fixture mAb-X-E3-r2 | -8.5e-15 | 2.4e-08 | 9.4e-05 | True |

Calibration (`nonlinear_calibration_0.13.3.json`, seed 131330001):

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| decay | coverage_k | 1000 | 0.9670 | 0.0056 | 0.9362 | True |
| decay | reportable_fraction_k | 1000 | 1.0000 | — | — | descriptive |
| decay | coverage_plateau | 1000 | 0.9560 | 0.0065 | 0.9362 | True |
| decay | reportable_fraction_plateau | 1000 | 1.0000 | — | — | descriptive |
| association | coverage_k | 1000 | 0.9600 | 0.0062 | 0.9362 | True |
| association | reportable_fraction_k | 1000 | 1.0000 | — | — | descriptive |
| association | coverage_plateau | 1000 | 0.9520 | 0.0068 | 0.9362 | True |
| association | reportable_fraction_plateau | 1000 | 1.0000 | — | — | descriptive |
| two_phase | coverage_k_fast | 1000 | 0.9420 | 0.0074 | 0.9362 | True |
| two_phase | reportable_fraction_k_fast | 1000 | 1.0000 | — | — | descriptive |
| two_phase | coverage_k_slow | 1000 | 0.9590 | 0.0063 | 0.9362 | True |
| two_phase | reportable_fraction_k_slow | 1000 | 1.0000 | — | — | descriptive |
| exp_growth_weighted | coverage_k | 1000 | 0.9500 | 0.0069 | 0.9362 | True |
| exp_growth_weighted | reportable_fraction_k | 1000 | 1.0000 | — | — | descriptive |
| michaelis_menten | coverage_vmax | 1000 | 0.9480 | 0.0070 | 0.9362 | True |
| michaelis_menten | reportable_fraction_vmax | 1000 | 1.0000 | — | — | descriptive |
| michaelis_menten | coverage_km | 1000 | 0.9530 | 0.0067 | 0.9362 | True |
| michaelis_menten | reportable_fraction_km | 1000 | 1.0000 | — | — | descriptive |
| logistic | coverage_xmid | 1000 | 0.9420 | 0.0074 | 0.9362 | True |
| logistic | reportable_fraction_xmid | 1000 | 1.0000 | — | — | descriptive |
| logistic | coverage_scal | 1000 | 0.9600 | 0.0062 | 0.9362 | True |
| logistic | reportable_fraction_scal | 1000 | 1.0000 | — | — | descriptive |
| decay_short | coverage_k | 946 | 0.9514 | 0.0070 | — | descriptive |
| decay_short | reportable_fraction_k | 1000 | 0.9460 | — | — | descriptive |
| decay_short | coverage_plateau | 0 | — | — | — | descriptive |
| decay_short | reportable_fraction_plateau | 1000 | 0.0000 | — | — | descriptive |
| mm_unsaturated | coverage_vmax | 2 | 0.0000 | 0.0000 | — | descriptive |
| mm_unsaturated | reportable_fraction_vmax | 1000 | 0.0020 | — | — | descriptive |
| mm_unsaturated | coverage_km | 2 | 0.0000 | 0.0000 | — | descriptive |
| mm_unsaturated | reportable_fraction_km | 1000 | 0.0020 | — | — | descriptive |

After the calibration, an overflow warning in the rate transform was silenced (np.errstate); 83 recorded
draws recomputed identically. In the unsaturated Michaelis-Menten stress row, 2 of 1000 fits passed the
saturation gate with Km underestimated; this is disclosed in the facts.
