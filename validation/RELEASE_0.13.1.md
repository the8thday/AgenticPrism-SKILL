# Release 0.13.1 — HTS null reference, ISR and carry-over, 5PL and bell-shaped dose response

Local development release on `dev/0.13.1` from `dev/0.13.0` (77cdd17). Three
opt-in additions; every default output is unchanged. No PK/PD added.

## Scope

1. **HTS layout-conditional hit reference** (`hts.py`, `hit_reference.method=layout_simulation`):
   simulated null plates keep the fitted additive row/column pattern, the positive-control
   offset and the exact layout; each is median-polished and scored by the same
   predictive statistic, standardized per well, and observed wells are referred to the
   pooled null. Addresses deferred follow-up 3 of the 0.11.1–0.13.0 review.
2. **Method validation: incurred sample reanalysis and carry-over** (`method_validation.py`):
   ISR percent difference against declared limits with a declared policy for pairs that
   are not both quantified, Clopper–Pearson and mean-difference t supplements and the M10
   extent diagnostic; carry-over as blank response after ULOQ relative to same-sequence
   LLOQ response from raw responses.
3. **Dose response: asymmetric 5PL and Prism bell-shaped models** (`dose_models.py`):
   5PL optimized on log EC50 directly so the profile-F interval is on EC50; bell-shaped
   with linear plateaus, profile-F intervals for both EC50s and a middle-plateau-reached
   gate. Both need a declared `model_rationale` and refuse RP comparisons.

## Evidence gates

| Gate | HTS simulated reference | ISR | Carry-over | 5PL | Bell-shaped |
|---|---|---|---|---|---|
| Established implementation | Same median polish already checked against base R `medpolish` (0.13.0); new code is the simulation reference | Arithmetic plus Clopper–Pearson (checked vs `binom.test` in 0.9.3) | Arithmetic ratio | drc 3.0.1 LL.5: passed | base-R `nls` (port) on Prism's equation: passed |
| Published worked example | UNMET | UNMET | UNMET | UNMET (S.alba is published data, not a printed 5PL result) | UNMET |
| Calibration | Passed all 4 registered rows | Passed all 8 rows | Not applicable (no inference) | Passed all 3 rows | Passed both separated rows; stress row descriptive (see below) |
| Live-agent misuse | Passed (hts_reseed) | Passed (isr_policy_choice) | Not scenario-tested | Not scenario-tested | Passed (bell_after_bad_4pl, overlapping_bell_ec50) |
| Interpretation facts | `must_mention` from `evidence_0131.py` | facts `must_mention` | facts `must_mention` | `calibration_evidence` items | `calibration_evidence` items |

Calibration scripts were committed before execution: HTS and the HTS implementation in
63a9c3b, ISR in d89e46d, dose models in b9648e3. Pre-execution design changes are written
into each script's docstring (HTS edge/gradient positive controls raised to +30 so
raw-scale Z′ QC can pass; u_separated EC50₂ widened from 1e-7 to 3e-7 because the true
curve reached only 93% of the middle plateau, next to the 90% gate). Both changes used
smoke-run withholding reasons only, never coverage.

## Observed calibration

### HTS hit reference (`validation/hts_calibration_0.13.1.json`, seed 131261003)

| Scenario | Reference | Metric | Population | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|---|---|
| hts_null_plain | layout_simulation | false_hit | evaluable | 1000 | 0.0410 | 0.0063 | 0.0638 | True |
| hts_null_plain | layout_simulation | false_hit | reportable | 999 | 0.0410 | 0.0063 | 0.0638 | True |
| hts_null_plain | predictive_t | false_hit | evaluable | 1000 | 0.0950 | 0.0093 | — | descriptive |
| hts_null_plain | predictive_t | false_hit | reportable | 999 | 0.0951 | 0.0093 | — | descriptive |
| hts_null_edge | layout_simulation | false_hit | evaluable | 1000 | 0.0440 | 0.0065 | 0.0638 | True |
| hts_null_edge | layout_simulation | false_hit | reportable | 1000 | 0.0440 | 0.0065 | 0.0638 | True |
| hts_null_edge | predictive_t | false_hit | evaluable | 1000 | 0.0930 | 0.0092 | — | descriptive |
| hts_null_edge | predictive_t | false_hit | reportable | 1000 | 0.0930 | 0.0092 | — | descriptive |
| hts_null_gradient | layout_simulation | false_hit | evaluable | 1000 | 0.0270 | 0.0051 | 0.0638 | True |
| hts_null_gradient | layout_simulation | false_hit | reportable | 738 | 0.0366 | 0.0069 | 0.0660 | True |
| hts_null_gradient | predictive_t | false_hit | evaluable | 1000 | 0.0720 | 0.0082 | — | descriptive |
| hts_null_gradient | predictive_t | false_hit | reportable | 738 | 0.0976 | 0.0109 | — | descriptive |
| hts_four_actives | layout_simulation | mean_fdp | evaluable | 1000 | 0.0404 | 0.0035 | 0.0571 | True |
| hts_four_actives | layout_simulation | active_detection | evaluable | 1000 | 0.6710 | 0.0117 | — | descriptive |
| hts_four_actives | layout_simulation | mean_fdp | reportable | 1000 | 0.0404 | 0.0035 | 0.0571 | True |
| hts_four_actives | layout_simulation | active_detection | reportable | 1000 | 0.6710 | 0.0117 | — | descriptive |
| hts_four_actives | predictive_t | mean_fdp | evaluable | 1000 | 0.0685 | 0.0047 | — | descriptive |
| hts_four_actives | predictive_t | active_detection | evaluable | 1000 | 0.7795 | 0.0097 | — | descriptive |
| hts_four_actives | predictive_t | mean_fdp | reportable | 1000 | 0.0685 | 0.0047 | — | descriptive |
| hts_four_actives | predictive_t | active_detection | reportable | 1000 | 0.7795 | 0.0097 | — | descriptive |

### ISR supplements (`validation/isr_calibration_0.13.1.json`, seed 131261004, level 0.90)

| Scenario (n, s, delta) | Truth: mean % diff / pass prob. | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|---|
| isr_small (12, 0.08, 0.0) | -0.000 / 0.9925 | mean_coverage | 1000 | 0.9140 | 0.0089 | 0.8810 | True |
| isr_small (12, 0.08, 0.0) | -0.000 / 0.9925 | pass_coverage | 1000 | 0.9960 | 0.0020 | 0.8810 | True |
| isr_shift (40, 0.15, 0.1) | 9.882 / 0.8009 | mean_coverage | 1000 | 0.9060 | 0.0092 | 0.8810 | True |
| isr_shift (40, 0.15, 0.1) | 9.882 / 0.8009 | pass_coverage | 1000 | 0.9440 | 0.0073 | 0.8810 | True |
| isr_large (100, 0.2, 0.0) | 0.000 / 0.7148 | mean_coverage | 1000 | 0.9010 | 0.0094 | 0.8810 | True |
| isr_large (100, 0.2, 0.0) | 0.000 / 0.7148 | pass_coverage | 1000 | 0.9240 | 0.0084 | 0.8810 | True |
| isr_noisy (40, 0.25, -0.05) | -4.852 / 0.6027 | mean_coverage | 1000 | 0.9050 | 0.0093 | 0.8810 | True |
| isr_noisy (40, 0.25, -0.05) | -4.852 / 0.6027 | pass_coverage | 1000 | 0.9150 | 0.0088 | 0.8810 | True |

### 5PL and bell-shaped profile-F intervals (`validation/dose_models_calibration_0.13.1.json`, seed 131261005, level 0.95)

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| 5pl_asym_low | ec50_1_coverage_reportable | 960 | 0.9500 | 0.0070 | 0.9359 | True |
| 5pl_asym_low | ec50_1_coverage_all_two_sided | 1000 | 0.9510 | 0.0068 | — | descriptive |
| 5pl_asym_low | reportable_fraction | 1000 | 0.9600 | — | — | descriptive |
| 5pl_asym_high | ec50_1_coverage_reportable | 797 | 0.9573 | 0.0072 | 0.9346 | True |
| 5pl_asym_high | ec50_1_coverage_all_two_sided | 1000 | 0.9590 | 0.0063 | — | descriptive |
| 5pl_asym_high | reportable_fraction | 1000 | 0.7970 | — | — | descriptive |
| 5pl_symmetric | ec50_1_coverage_reportable | 982 | 0.9440 | 0.0073 | 0.9361 | True |
| 5pl_symmetric | ec50_1_coverage_all_two_sided | 1000 | 0.9440 | 0.0073 | — | descriptive |
| 5pl_symmetric | reportable_fraction | 1000 | 0.9820 | — | — | descriptive |
| bell_separated | ec50_1_coverage_reportable | 1000 | 0.9560 | 0.0065 | 0.9362 | True |
| bell_separated | ec50_1_coverage_all_two_sided | 1000 | 0.9560 | 0.0065 | — | descriptive |
| bell_separated | ec50_2_coverage_reportable | 1000 | 0.9480 | 0.0070 | 0.9362 | True |
| bell_separated | ec50_2_coverage_all_two_sided | 1000 | 0.9480 | 0.0070 | — | descriptive |
| bell_separated | reportable_fraction | 1000 | 1.0000 | — | — | descriptive |
| u_separated | ec50_1_coverage_reportable | 990 | 0.9515 | 0.0068 | 0.9361 | True |
| u_separated | ec50_1_coverage_all_two_sided | 1000 | 0.9500 | 0.0069 | — | descriptive |
| u_separated | ec50_2_coverage_reportable | 990 | 0.9455 | 0.0072 | 0.9361 | True |
| u_separated | ec50_2_coverage_all_two_sided | 999 | 0.9449 | 0.0072 | — | descriptive |
| u_separated | reportable_fraction | 1000 | 0.9900 | — | — | descriptive |
| bell_overlap_stress | ec50_1_coverage_reportable | 33 | 0.3939 | 0.0851 | — | descriptive |
| bell_overlap_stress | ec50_1_coverage_all_two_sided | 1000 | 0.9780 | 0.0046 | — | descriptive |
| bell_overlap_stress | ec50_2_coverage_reportable | 33 | 0.3030 | 0.0800 | — | descriptive |
| bell_overlap_stress | ec50_2_coverage_all_two_sided | 1000 | 0.9700 | 0.0054 | — | descriptive |
| bell_overlap_stress | reportable_fraction | 1000 | 0.0330 | — | — | descriptive |

### Cross-implementation (`validation/dose_models_benchmarks_0.13.1.json`)

| Dataset | Oracle | SSE rel. diff (Python − oracle) | EC50 rel. diff | EC50₂ rel. diff | Python lower objective | Passed | Initial run passed |
|---|---|---|---|---|---|---|---|
| 5pl:S1-E1 | drc 3.0.1 LL.5 | -4.56e-07 | 6.01e-05 | — | False | True | True |
| 5pl:S1-E2 | drc 3.0.1 LL.5 | -5.84e-07 | 7.54e-05 | — | False | True | True |
| 5pl:S1-E3 | drc 3.0.1 LL.5 | -1.18e-07 | 1.55e-04 | — | False | True | True |
| bell:S1-E1 | R 4.6.0 stats::nls port, Prism bell-shaped equation | -1.68e-14 | 9.69e-09 | 7.12e-10 | False | True | True |
| bell:S1-E2 | R 4.6.0 stats::nls port, Prism bell-shaped equation | -5.67e-15 | 7.02e-09 | 5.08e-09 | False | True | True |
| bell:S1-E3 | R 4.6.0 stats::nls port, Prism bell-shaped equation | 1.42e-15 | 6.55e-09 | 6.54e-09 | False | True | True |
| bell_overlap:S1-E1 | R 4.6.0 stats::nls port, Prism bell-shaped equation | -2.44e-05 | 4.83e-02 | 4.77e-02 | True | True | False |
| drc::S.alba bentazone (Christensen et al. 2003) | drc 3.0.1 LL.5 | -1.67e-04 | 2.16e-02 | — | True | True | False |

Initial-run datasets: 5pl:S1-E1=True, 5pl:S1-E2=True, 5pl:S1-E3=True, bell:S1-E1=True, bell:S1-E2=True, bell:S1-E3=True, bell_overlap:S1-E1=False, drc::S.alba bentazone (Christensen et al. 2003)=False

Key observations:

- HTS: the simulated reference keeps the global-null false-hit probability at or below
  0.044 (plain, edge, gradient) where the 0.13.0 predictive-t reference gives 0.093–0.098
  on the same plates. It costs power: 0.671 versus 0.779 of +5 SD actives detected.
  Calibrated only for independent Gaussian wells with additive plate effects.
- ISR: Clopper–Pearson is conservative when the pass probability is near 1 (0.996 at
  n = 12).
- Bell-shaped stress: with phases one decade apart, 96.7% of fits were withheld, but the
  33 that passed every gate covered the true EC50s only 0.394 and 0.303. The gate selects
  the worst draws among overlapping fits. Saved facts tell users to treat EC50s from
  phases less than about two decades apart as unsupported even when reported. This is
  not hidden by the passing separated rows.

## Cross-implementation notes

- The initial benchmark run (`dose_models_benchmarks_0.13.1_initial.json`) failed S.alba:
  the 5PL multistart had no steep Hill starts and stopped at the asymmetry upper bound
  (g = 20, SSE 3.8096 versus drc 3.8075). h = 4 starts were added and the run repeated;
  Python then reached a lower objective than drc (3.80688) with an interior solution.
  The fit is still withheld by the workflow (sparse transition doses, nuisance at a
  bound). ELISA `calibration.fit_5pl` had the same start grid; see the ELISA section below.
- `python_lower_objective` (Python SSE below the oracle by more than 1e-6 relative) was
  added as a pass condition after the initial run; it applies to bell_overlap and S.alba.
- A first draft of the ryegrass benchmark was replaced by S.alba before any record was
  written: ryegrass has six positive doses, below the 5PL gate of seven.

## Byte identity and tests

- `scripts/validate_legacy_0131.py` against `dev/0.13.0` (77cdd17): 154 configurations,
  1167 scientific artifacts byte-identical (`legacy_0.13.1.json`). The new options are
  popped from resolved configs when unused.
- Full suite: 691 passed (macOS, Python 3.13). New tests cover config gates, legacy
  config bytes, determinism and row-order invariance (HTS), ISR policies and shift
  flag, carry-over ordering, 5PL/bell recovery, Prism-equation identity, U-shape
  direction, overlap withholding, facts withholding and render invariance.
- Not run for this local release: wheel build, clean-environment and isolated install
  checks. Live-agent scenarios: 10/10 passed (see the live-agent section).

## Self-review

- HTS: the reference is a parametric model of the plate; nonadditive or biological
  layout effects, non-Gaussian tails and correlated wells are outside its calibration.
  Seed and simulation count are declared before analysis; changing them to move hits
  would be p-hacking and the Skill says so.
- ISR: the unquantified-pair policy changes the denominator; it must be declared, and
  the facts state how many pairs it touched.
- Carry-over: raw responses only; no concentration claim.
- Dose models: shape is declared before fitting with a rationale; 4PL remains the
  default and the only model for relative potency. Bell EC50s are model parameters
  relative to the fitted middle plateau.

## ELISA 5PL start grid (changed method)

A targeted scan (300 steep synthetic standard curves, h 2.5–8, development check, not a
registered record) found the 0.13.0 ELISA 5PL grid missed a lower objective in 5 curves:
3 stopped at the asymmetry bound (withheld by the boundary gate) and 2 reportable fits
sat at a worse point on a flat surface. `calibration.fit_5pl` now appends h = 4 starts;
the new grid is a superset of the old one (commit f289d87, before any registered run).

- **Registered drc check** (`elisa_5pl_benchmarks_0.13.1.json`, seed 131261006, S.alba plus
  40 steep curves): **31/41 passed; the registered gate is
  not met.** In every miss drc's optimum lies outside the declared bounds (9 with Hill
  slope 8.8–10.7 above the default `hill_bounds` upper limit 8, one with g = 21.6 above
  20). All 10 production fits were withheld by `numerical_boundary_hit`; no miss produced
  a reportable value. The 0.13.0 grid was no worse than the new grid on any of these 41
  curves (descriptive). The pass rule did not anticipate declared bounds; it is not
  changed after the fact.
- **Post-hoc supplement** (`elisa_5pl_bounds_supplement_0.13.1.json`): refitting the 10
  misses with `hill_bounds` [0.05, 20] matched or beat drc in
  9/10; the remaining curve (steep-02) is held by the g
  bound, as stated before the supplement ran. This supplements, and does not replace, the
  registered miss.
- **Interval calibration rerun** (`elisa_interval_coverage_0.13.1.json`, same seed and
  design as 0.7.0): all 12 rows identical to `elisa_interval_coverage.json` in reported
  count, coverage and pass status (5PL rows 0.9367 / 0.9467 / 0.9606).
- **Fixture change**: `elisa_dilution_5pl` artifacts move within optimizer tolerance only
  (maximum relative difference 5.0e-7, all statuses and flags identical);
  `legacy_0.13.1.json` records them as numerically changed and the other 1161 artifacts
  as byte-identical. `release_validation.legacy` gained a `changed` tolerance map for
  such intended changes; embedded SHA-256 hashes are masked in that comparison.
- The ELISA Skill now says that `hill_bounds` may be widened only for a reason declared
  before seeing the plate.

## Sample size and power (new specialist)

`sample_size.py`, `skills/sample-size/`; design-only runs through the extension workflow,
which now accepts configs without an input file.

| Gate | Status |
|---|---|
| Established implementation | pwr 1.3.0, `power.prop.test` (strict), PowerTOST 1.5.7 exact: 306 cases, max absolute difference 9.5e-10; integer n equals `ceiling(pwr.t.test()$n)` in 6/6 (`sample_size_benchmarks_0.13.1.json`). Log-rank has no package oracle. |
| Published worked example | UNMET |
| Calibration | 10/10 rows passed (`sample_size_calibration_0.13.1.json`, seed 131261007, script committed in 197cbc2 before running) |
| Live-agent misuse | Passed (see the live-agent section) |
| Interpretation facts | `must_mention` from `evidence_0131.py['sample_size']` |

The initial benchmark (`sample_size_benchmarks_0.13.1_initial.json`) reported `passed`
while three Python powers were NaN: SciPy's `nct.cdf(-q)` fails for large noncentrality
and pandas skipped the rows in its summary. The lower tail now uses symmetry
(`nct.sf(q, df, -ncp)`), non-finite power raises, and the benchmark fails on NaN rows.

| Scenario | Design | Solved n | Computed power | Simulated rejection | MCSE | Bound (abs) | Passed |
|---|---|---|---|---|---|---|---|
| t_equal | two_sample_t | 26 | 0.8075 | 0.8120 | 0.0062 | 0.0125 | True |
| t_ratio2 | two_sample_t | 45 | 0.9036 | 0.9083 | 0.0046 | 0.0093 | True |
| paired | paired_t | 34 | 0.8078 | 0.8125 | 0.0062 | 0.0125 | True |
| anova4 | one_way_anova | 17 | 0.8036 | 0.8090 | 0.0062 | 0.0126 | True |
| prop_arcsine | two_proportions | 93 | 0.8013 | 0.7980 | 0.0063 | 0.0126 | True |
| prop_normal | two_proportions | 93 | 0.8000 | 0.8037 | 0.0063 | 0.0126 | True |
| logrank_complete | logrank | 121 | 0.8022 | 0.7935 | 0.0064 | 0.0126 | True |
| logrank_censored | logrank | 199 | 0.8009 | 0.7927 | 0.0064 | 0.0126 | True |
| tost_equal | tost_two_means | 11 | 0.8344 | 0.8390 | 0.0058 | 0.0118 | True |
| tost_ratio15 | tost_two_means | 9 | 0.8343 | 0.8400 | 0.0058 | 0.0118 | True |

The proportion rows check the arcsine and pooled-normal approximations against the
Pearson chi-square test without continuity correction; the log-rank rows check the
Schoenfeld formula against the package log-rank, with and without exponential censoring.

## Competing risks (time-to-event extension)

`competing.py` (`analysis_type=competing_risks`): ports of cmprsk 2.2-12 `cinc`, `crst` and
`crr` (Newton with backtracking; `crrvv` sandwich variance, whose triple sum over censoring
times is evaluated through the separable censoring weights) plus cause-specific Efron Cox.

| Gate | Status |
|---|---|
| Established implementation | cmprsk 2.2-12 and survival 3.8-6 on 6 datasets (`competing_benchmarks_0.13.1.json`): |
| Published worked example | UNMET (mgus2 is public data; no printed analysis reproduced) |
| Calibration | 6/6 metrics passed (`competing_calibration_0.13.1.json`, seed 131261009, script committed in ef05429 before running) |
| Live-agent misuse | Passed (see the live-agent section) |
| Interpretation facts | `must_mention` from `evidence_0131.py['competing_risks']` |

| Dataset | n | CIF/var max abs | Gray max abs | crr coef max abs | crr var max rel | Cox coef max abs | Passed |
|---|---|---|---|---|---|---|---|
| synthetic seed 3 | 150 | 5.0e-16 | 8.9e-15 | 6.7e-16 | 6.5e-14 | 1.6e-15 | True |
| synthetic seed 11 | 150 | 4.7e-16 | 1.8e-15 | 3.3e-16 | 3.7e-14 | 9.4e-16 | True |
| synthetic seed 29 | 150 | 5.0e-16 | 6.9e-14 | 3.3e-16 | 1.1e-14 | 2.1e-15 | True |
| synthetic seed 47 | 150 | 5.0e-16 | 1.2e-14 | 3.1e-16 | 2.5e-14 | 7.7e-16 | True |
| synthetic seed 83 | 150 | 5.0e-16 | 1.8e-15 | 3.3e-16 | 3.0e-14 | 1.8e-15 | True |
| survival::mgus2 (Kyle et al.; progression vs death) | 1384 | 5.6e-16 | 1.1e-14 | 2.6e-15 | 1.8e-14 | 8.4e-16 | True |

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---|---|---|---|---|
| cif_landmark | coverage_1 | 1000 | 0.9540 | 0.0066 | 0.9362 | True |
| cif_landmark | coverage_2 | 1000 | 0.9430 | 0.0073 | 0.9362 | True |
| fine_gray_200 | coverage_1 | 1000 | 0.9520 | 0.0068 | 0.9362 | True |
| fine_gray_60 | coverage_1 | 1000 | 0.9570 | 0.0064 | 0.9362 | True |
| gray_null | type_I_error | 1000 | 0.0510 | 0.0070 | 0.0638 | True |
| cause_specific | coverage_1 | 1000 | 0.9640 | 0.0059 | 0.9362 | True |

The Fine–Gray generator was checked against the closed-form F1(t|x) before the run
(400 000 uncensored draws, differences within Monte Carlo error).

## Thermal unfolding (new specialist)

`thermal.py`, `skills/thermal-stability/` (`analysis_type=thermal_unfolding`).

| Gate | Status |
|---|---|
| Established implementation | base-R `nls(algorithm = "plinear")` on the same model, 7 fixture curves, SSE within 2.4e-13 relative and Tm within 1.4e-7 °C (`thermal_benchmarks_0.13.1.json`) |
| Published worked example | UNMET |
| Calibration | Run 1: 3 misses (retained). Run 2 after gate changes: 1 miss (retained). See below |
| Live-agent misuse | Passed (see the live-agent section) |
| Interpretation facts | `must_mention` from `evidence_0131.py['thermal_unfolding']`, both runs |

Development checks (not registered): the first start grid stopped at a wrong optimum for
a three-transition curve; derivative-seeded starts and local restarts were added, after
which a 120-curve scan found no fit worse than the generating parameters.

**Run 1** (`thermal_calibration_0.13.1.json`, seed 131261008, script committed in 3b2086b):

| Row | Transition | n reportable | Coverage | Reportable fraction | Bound | Passed |
|---|---|---|---|---|---|---|
| single | 1 | 1000 | 0.955 | 1.000 | 0.9362 | True |
| mab3 | 1 | 1000 | 0.928 | 1.000 | 0.9362 | False |
| mab3 | 2 | 944 | 0.944 | 0.944 | 0.9358 | True |
| mab3 | 3 | 971 | 0.953 | 0.971 | 0.9360 | True |
| mab3_noisy | 1 | 922 | 0.914 | 0.922 | 0.9356 | False |
| mab3_noisy | 2 | 765 | 0.919 | 0.765 | 0.9342 | False |
| mab3_noisy | 3 | 767 | 0.948 | 0.767 | 0.9343 | True |
| overlap_stress | 1 | 637 | 0.890 | 0.671 | — | descriptive |
| overlap_stress | 2 | 55 | 0.055 | 0.058 | — | descriptive |
| delta_tm | 1 | 497 | 0.970 | 0.994 | 0.9304 | True |
| delta_tm | 2 | 367 | 0.959 | 0.734 | 0.9272 | True |
| delta_tm | 3 | 370 | 0.968 | 0.740 | 0.9273 | True |

In the overlap stress row 51 draws raised an error in the calibration script's coverage
scoring (a half-open interval compared with `None`); those transitions were withheld
by `profile_interval_open`, but the other transition of each draw was not scored.

**Post-hoc supplement** (`thermal_profile_supplement_0.13.1.json`): for every run-1 miss
in the multi-transition rows, the profile SSE at the true Tm was recomputed from 8 starts.

| Row and transition | Misses | Truth inside after re-optimization | Median relative excess over threshold |
|---|---|---|---|
| mab3_transition_1 | 72 | 2 | 0.0087 |
| mab3_transition_2 | 53 | 0 | 0.0121 |
| mab3_transition_3 | 46 | 0 | 0.0101 |
| mab3_noisy_transition_1 | 79 | 16 | 0.0061 |
| mab3_noisy_transition_2 | 62 | 1 | 0.0114 |
| mab3_noisy_transition_3 | 40 | 3 | 0.0094 |

Most misses are genuine profile-F approximation limits (the truth stays just outside);
16/79 noisy transition-1 misses came from under-optimized profiles. Inspection of the
overlap row showed the leaked transition-2 estimates came from a spurious low-amplitude
component next to a merged transition.

**Changes before run 2** (commit e9f981a): refit with k − 1 transitions and withhold every Tm
unless k beats k − 1 (extra-sum-of-squares F, p < 0.01); withhold transitions adjacent to an
unresolved one; five nuisance starts per profile point; half-open intervals score as not
covering.

**Run 2** (`thermal_calibration_run2_0.13.1.json`, seed 131261010, same rows and bounds):

| Row | Transition | n reportable | Coverage | Reportable fraction | Bound | Passed |
|---|---|---|---|---|---|---|
| single | 1 | 1000 | 0.950 | 1.000 | 0.9362 | True |
| mab3 | 1 | 941 | 0.949 | 0.941 | 0.9358 | True |
| mab3 | 2 | 941 | 0.949 | 0.941 | 0.9358 | True |
| mab3 | 3 | 941 | 0.963 | 0.941 | 0.9358 | True |
| mab3_noisy | 1 | 721 | 0.943 | 0.721 | 0.9338 | True |
| mab3_noisy | 2 | 703 | 0.947 | 0.703 | 0.9336 | True |
| mab3_noisy | 3 | 701 | 0.932 | 0.701 | 0.9335 | False |
| overlap_stress | 1 | 4 | 0.750 | 0.004 | — | descriptive |
| overlap_stress | 2 | 9 | 0.111 | 0.009 | — | descriptive |
| delta_tm | 1 | 298 | 0.973 | 0.596 | 0.9247 | True |
| delta_tm | 2 | 298 | 0.980 | 0.596 | 0.9247 | True |
| delta_tm | 3 | 298 | 0.966 | 0.596 | 0.9247 | True |

Run 2 still misses the noisy third transition (0.932 against 0.9335), still lets 0.9% of
overlapped second transitions through (covering 0.11), and withholds far more: about 72%
of noisy three-transition curves are reportable and delta-Tm summaries exist in 60% of
simulated comparisons, because one replicate failing the structure test withholds the
sample. No further changes were made after run 2.

## Live-agent misuse review

Ten scenarios with fresh inputs and criteria committed before the run (518ccfb):
**10/10 passed**, scored by the implementer (not blinded). Details and every response:
[agent-scenarios/0.13.1/RESULTS.md](agent-scenarios/0.13.1/RESULTS.md). Findings: the default 4PL
does not gate lack of fit on hook data (follow-up); one overlapped bell curve passed the gates
in practice and was withheld by the agent from the facts disclosure; one scenario dataset was
borderline significant by design error; the harness writes traces inside the agent workspace.
Carry-over and the ELISA start-grid change were not scenario-tested.

