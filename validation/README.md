# Validation — AgenticPrism development 0.13.1

## 0.13.1

[Release evidence](RELEASE_0.13.1.md): opt-in layout-conditional HTS hit reference
(registered global-null false-hit 0.041 / 0.044 / 0.037 versus 0.093–0.098 for the default;
detection 0.671 versus 0.779), ISR and carry-over in method validation (ISR supplements
covered 0.901–0.996 at 90%), and asymmetric 5PL / bell-shaped dose response (drc LL.5 and
base-R nls 8/8 after a retained failed initial run; 95% coverage 0.944–0.957; the
overlapping-phase stress row covered 0.39 / 0.30 among the 3.3% of fits that passed the
gates). Unchanged scientific artifacts: 154 configurations, 1167 files byte-identical to
dev/0.13.0. Published worked examples UNMET; live-agent review PENDING; no build or
clean-install check for this local release. ELISA 5PL start grid widened (changed method):
registered drc check 31/41 with every miss outside the declared Hill/asymmetry bounds and
withheld (post-hoc wider Hill bound matched 9/10); ELISA interval calibration rerun
identical to 0.7.0 in all 12 rows; ELISA 5PL fixture changes at most 5.0e-7 relative.
New sample-size specialist: 306 power cases agree with pwr/power.prop.test/PowerTOST
(9.5e-10; initial run hid three NaN powers, retained) and 10/10 simulated-power rows passed.
Competing risks: cmprsk/survival agreement within 7e-14 on six datasets including mgus2;
6/6 calibration metrics passed. Thermal unfolding: nls agreement 2.4e-13; registered run 1
missed multi-transition coverage (0.914–0.928) and leaked overlapped fits; after a k vs k−1
structure test, neighbour gate and multistart profiles, run 2 passed all rows but the noisy
third transition (0.932) with a 0.9% overlap leak. Both runs and a post-hoc supplement kept.

## 0.13.0

[Release evidence](RELEASE_0.13.0.md): Named Bliss, Loewe, HSA and ZIP combination references with independent-matrix uncertainty; declared HTS plate QC, median-polish B scores and exploratory FDR hits. Dose-response and ELISA facts complete the oldest-six retrofit.
All observed misses and unmet gates retained; live-agent review PENDING.

[Four-release review, all calibration misses and pending commands](RELEASE_SERIES_0.11.1-0.13.0.md).

## 0.12.1

[Release evidence](RELEASE_0.12.1.md): Directed epitope binning with declared controls and thresholds, asymmetric-pair diagnostics, conditional bootstrap clustering stability and reciprocal-block communities.
All observed misses and unmet gates retained; live-agent review PENDING.

## 0.12.0 — Surface kinetic mechanisms

[Release evidence](RELEASE_0.12.0.md): declared complex surface models, deSolve
checks, retained calibration misses, scoped real-export imports and kinetics facts.
Public worked-example fit reproduction remains limited; live-agent review PENDING.

## 0.11.1 — Routine statistics

[Release evidence](RELEASE_0.11.1.md): independent factorial and categorical
methods, correlation/regression and method comparison; numerical R comparisons,
registered calibration with retained misses, facts and exact legacy checks.
Live-agent gate PENDING.


## 0.11.0 — equilibrium affinity depth and separate cell-binding specialist

[Release evidence](RELEASE_0.11.0.md): opt-in mass-balance affinity models,
three-state Ki, cell apparent KD and plateau steady state. BindCurve and R
comparisons cover the new implementations. 14,000 simulations preserve five
failed criteria; published worked-example gates remain UNMET. Unchanged
scientific artifacts: 100 configurations, 800 byte-identical files against
`dev/0.10.1`; equilibrium adds saved-result facts. Full suite 554 passed; six changed Skills valid; six live-agent scenarios passed
(after preserving Claude authentication failures). Clean install has 98 exact
facts matches; isolated smoke passes 109 configurations. All 109 new scientific
artifacts across nine configs match source. See release_checks_0.11.0.json.
Review fixes (SET valency/readout gate, declared cell control baseline, shared
Rmax for steady state, non-blocking titration-span warning) left every registered
calibration result unchanged (56 draws recomputed); a post-hoc supplement shows
the known-Pt Pt/KD=100 row titrated only up to Pt (reportable fits 95.0%).
Benchmarks 36/36; independent Claude rerun of the six scenarios 6/6; after the
fixes 558 tests pass, clean install has 99 exact facts matches, isolated smoke
110 configurations and 121 new artifacts across ten configs match source.

## 0.10.1 — comparability, tolerance intervals and process capability

[Release evidence](RELEASE_0.10.1.md): new comparability (TOST, quality range) and specifications
(exact/Howe/one-sided normal and order-statistic tolerance intervals, Pp/Ppk, Cp/Cpk) Skills. 115 fields
agree with R `tolerance` and `t.test` (max relative 3.3e-9); six NIST/SEMATECH printed values reproduced.
Calibration (script committed before running): 13 of 15 rows pass; the n = 10 Pp and Ppk rows miss (93.6%
vs 93.62%), which a post-hoc 100,000-dataset check attributes to Monte Carlo error (94.9%, 95.4%). Five
live-agent scenarios pass.


## 0.10.0 — CMC stability and replicated potency across runs

[Release evidence](RELEASE_0.10.0.md): new stability and potency-assay Skills,
Python Q1E regression/ANCOVA and log-RP random-run REML/MLS/MOVER. All 318
benchmark fields pass (308 R and 10 printed regression/ANCOVA/SE fields).
Calibration: 18/22 rows pass, four stability coverage misses disclosed in Skills
and saved facts. Printed stability shelf life and a matching published potency
example are unmet gates. Five live-agent scenarios are pending reviewer run
because Claude CLI authentication was unavailable. Legacy: 76 configurations,
680 byte-identical artifacts against dev/0.9.3. Comparability, specifications,
Arrhenius and oldest-module facts are not started. See
[release checks](release_checks_0.10.0.json) for tests and clean installation.

## 0.9.3 — method validation, ADA sensitivity/drug tolerance, two-way bootstrap bound

[Release evidence](RELEASE_0.9.3.md): new method-validation Skill (accuracy/precision with total
error, dilution linearity and hook, parallelism with trend, selectivity, specificity, stability),
ADA sensitivity and drug tolerance, a two-way bootstrap nonparametric lower bound, and MLS as the
default precision interval. 188 fields agree with an independent R oracle (max relative 3.4e-12).
Calibration (script committed before running): 17 of 18 rows with criteria pass; the sensitivity
prediction limit with 3-fold dilution spacing misses (92.1% vs 92.2%). Five live-agent scenarios
pass. No printed worked example of the Mee interval was available.

## 0.9.2 — precision intervals and ADA lower confidence bounds

[Release evidence](RELEASE_0.9.2.md): orthogonal MLS, unbalanced correlated
MOVER, boundary upper bounds, and parametric/independent-pair ADA lower limits.
Reviewer 84.6% / 95.6% coverage reproduced. Nested-crossed lot upper-bound
miss (92.8%) and conservative nonparametric titer FPR are disclosed.
All calibration rows include MCSE; historical records are preserved.
See [release checks](release_checks_0.9.2.json) and
[live-agent attempts](agent-scenarios/0.9.2/RESULTS.md).


## 0.9.1 — ADA recalibration and general precision REML

[Release record](RELEASE_0.9.1.md): reviewer ADA numbers reproduced using exact
conditional FPR; old floating/IQR failure labels superseded. New Python nested/
crossed REML with unbalanced data and Satterthwaite intervals, regenerated VCA/
lme4 comparisons and EP05 worked example. Small three-lot total coverage and
true-zero component inference miss their bounds; all rows and MCSE are recorded.
Live-agent gate is pending explicit external-Claude authorization following an
automatic review rejection. See [checks](release_checks_0.9.1.json).


> The published distribution contains these records but not `runs/`, `tests/`,
> `scripts/` or the agent-session transcripts. Links to those items below resolve
> only in the development repository; the results they document are summarized
> in the text.

## Scoped ADA cut points and shared variance components — 2026-09-27 (0.9.0)

Adds a Python ADA specialist for complete balanced drug-naive negative panels:
point screening/confirmatory/titer percentiles, declared normalization and
transformations, audited outlier policies, fixed/floating gates and a reusable
subject/run crossed ANOVA primitive. Dynamic deployment, lower cut-point bounds,
separate analyst/day/plate REML, sensitivity/drug tolerance and full method
validation are not implemented. New ADA runs save interpretation facts.

Nine regenerated R comparison rows pass; positive component estimates agree
with lme4 within 8.072245e-7 relative. Published rADA scalar percentiles and the
Bates Penicillin variance example are reproduced at printed precision. A second
lecture's printed sample-variance discrepancy is retained as a miss. Nine
cut-point and three component-interval calibration rows each use 1,000 datasets:
floating parametric screening (6.8349%), biological-IQR exclusion (7.1579%) and
unanticipated future shift (61.5628%) miss their fixed FPR bounds. All three
component simultaneous-coverage rows pass (95.9%, 96.1%, 96.9%).

Full pytest: 349 passed. Both changed Skills validate. Legacy regression retains
428 byte-identical artifacts across 34 configurations; installed wheel/source
facts match in 31 configurations. All four final live-agent stopping guards
pass, but A1/A3 narrative overgeneralizations remain. See
[0.9.0 evidence](RELEASE_0.9.0.md), [release checks](release_checks_0.9.0.json)
and [scenario assessment](agent-scenarios/0.9.0/RESULTS.md) for exact scope,
misses, defects and installation results. Work remains local on dev/0.9.0.

## Rank tests, MMRM and runtime discovery — 2026-09-27 (0.8.2)

Adds MW/signed-rank with HL intervals, KW/Dunn and Friedman, plus marginal
US/AR(1) MMRM with Python Satterthwaite or optional pinned R Kenward–Roger.
The runtime guide now checks the collection `.venv` before PATH. Both new
analysis types save interpretation facts. New R benchmarks include public
worked examples and fixtures; the KR/direct-mmrm maximum relative difference
is 6.0e-15, and Python MMRM/direct-mmrm is 1.810916e-6. An exact-MW tail cancellation
was fixed; the tiny-tail R comparison is within 9.86e-15 relative. Approximate
rank endpoints have limited near-zero precision, detailed in the release record.

All 20 rank-test calibration rows met their fixed evaluable-run bounds; two
tied signed-rank rows had 384 and 31 reduced-confidence intervals that were
excluded from the 95% coverage denominator. US/Satterthwaite missed interaction
rejection (7.0%); US/KR missed interaction/family rejection (7.1%/6.7%) and
coverage (93.3%). Both AR(1) rows passed. B=999 bootstrap coverage was
95.5%/94.0%/95.6% across the three fixed rows (3/3 rows passed all bounds).
Each row generated 1,000 datasets; no bounds, seeds or scenarios were retuned.

The full suite passed 306 tests. The final wheel's facts match source exactly
for 23 configurations; isolated install smoke passed 33 configurations.
Optional-R checks covered absent R and an existing pinned R library. Legacy
fixtures retain 295 byte-identical artifacts across 22 configs; one facts
limitation was updated without changing numerical results. See
[0.8.2 release evidence](RELEASE_0.8.2.md) for every row, R scope, preserved
failures, live-agent scores and remaining limits. Historical evidence below is
unchanged; no 0.9 work was started.

## Interpretation facts — 2026-09-27 (0.8.1)

Adds a source-linked facts artifact to repeated-measures, time-to-event and
tumor-growth, with no new statistical method. Twelve configurations reproduce
170 legacy scientific artifacts byte for byte against dev/0.8.0; 487 facts
values match their saved sources exactly. The published R ovarian log-rank
example matches within 9.18e-15 relative. Existing R benchmarks were rerun.
The full suite passed 250 tests; four changed Skills passed quick validation.
All three new RM/RI simulation rows used 1,000 draws and passed their fixed
inference bounds; all 3,000 facts/state checks passed. Survival and tumor
calibrations reproduce the earlier limits and comparator shortfalls.
Live-agent core misuse checks passed in the final three scenarios, with one
partial workflow result (runtime location); there was no independent human
re-score. See [0.8.1 release evidence](RELEASE_0.8.1.md) for every simulation
row, omissions and installation evidence. Other 0.8 follow-ups remain open.

## Satterthwaite, arm × time, time to event and tumor growth — 2026-09-27 (0.8.0)

Adds Satterthwaite inference and arm × condition designs to repeated measures,
and two new Skills: time-to-event and tumor-growth (nine Skills in total). All
new estimators agree with R (lmerTest, afex, emmeans, survival, nlme) on public
datasets and fixtures to within 1e-6 relative (most to 1e-8 or better). See
[RELEASE_0.8.0.md](RELEASE_0.8.0.md) for the observed numbers.
Calibration simulations with bounds fixed before running:
repeated-measures Satterthwaite held 4.7–5.5% type I error and 94.9–96.0%
simultaneous coverage (8–12 units, including weak clustering, where it was
withheld in 33% of datasets at the variance boundary); split-plot GG and the
two-way Satterthwaite tests held 4.5–5.2%. One miss: the one-factor parametric
bootstrap (199 refits) covered 93.4% against a 93.6% bound at 12 units.
Survival: permutation log-rank 4.3–4.9% at n = 6–10 per arm, compared with 6.6%
for the asymptotic test at n = 6. Log-log KM limits covered 96.4% and log
limits 91.9%. Tumor growth: rate-difference family coverage 95.1% under MAR
dropout, where survivor-based observed TGI fell to 92.2%. Twelve live-agent
scenarios passed. The appendix of that record keeps the first 0.8.0 commit's
evidence; the entry below describes that commit only, and its 100-experiment
smoke check is superseded.

## Repeated measurements — 2026-09-27 (0.8.0)

One-factor RM ANOVA with GG correction, and Gaussian random-intercept REML
with explicit missingness policy and parametric bootstrap or asymptotic Wald
inference. Independent R/nlme fitting comparisons, prespecified simulations,
input/refusal tests and installed-wheel checks are in
[RELEASE_0.8.0.md](RELEASE_0.8.0.md). The 100-experiment mixed-model simulation
is a coarse smoke check (93% simultaneous coverage); it does not establish
nominal coverage in other designs. No KR/Satterthwaite, random slopes,
treatment-by-time, AR(1) or general MMRM implementation is claimed.

## Installable distribution — 2026-09-27 (0.7.1)

Packaging only: `install.py`, `agentic-prism doctor`, shared runtime instructions for the
Skills, and `requires-python >=3.12`. Fresh-copy installs passed on macOS with uv (Python
3.13) and with pip (Python 3.14), and 3.14 reproduced the reference results exactly. Two
live-agent installation scenarios passed. See [RELEASE_0.7.1.md](RELEASE_0.7.1.md).

## New modules — 2026-09-26 (0.7.0)

See [RELEASE_0.7.0.md](RELEASE_0.7.0.md) for the checks and observed numbers:
one-way multi-group ANOVA with predeclared families, plate-reader grid import,
ELISA 5PL, unknown intervals and dilution linearity, equivalence-margin
parallelism with RP acceptance limits, single-cycle kinetics and double
referencing, plus recorded live-agent sessions. Four simulations used bounds
fixed before running. Three results missed their bound and are reported as
misses, with no retuning: Welch ANOVA and Games-Howell with an n=4
high-variance group (0.069, 0.077 vs 0.064, 0.075), and single-cycle koff
bootstrap coverage (0.913 vs 0.914). Multi-cycle kinetic artifacts from the
0.6.0 reference library reproduce byte for byte ([regression record](regression_0.7.0.json)).

## Reliability extension — 2026-09-26 (0.6.0)

See [RELIABILITY_0.6.0.md](RELIABILITY_0.6.0.md) for current gates, executable
reference-library results, independent QC and installation evidence. Sections below
are historical release evidence; their earlier reportability decisions and test
counts do not describe 0.6.0. In particular the earlier public kinetic intervals
are now audit-only, and ELISA standard-derived LLOQ/ULOQ keys are per-plate
screening bounds, not experimentally established assay limits.

## ELISA standard QC and review fixes — 2026-09-24 (0.5.1)

Follow-up to an independent review of the 0.5.0 ELISA and two-group modules.
That review found both numerically correct: the ELISA standard fit and unknown
inverses matched an independent full-parameter fit with root-finding
(≤3.9e-10), and the Welch and paired t statistics, SE and intervals matched
first principles, with 95.3% / 95.8% simulated coverage of 95% intervals
(4000 simulations each).

| Change | Evidence | Result / boundary |
|---|---|---|
| Standard back-calculation QC (ICH M10 ligand-binding convention, configurable) | `tests/test_elisa_qc.py` (6 tests) | Back-calculated standards, recovery and CV match an independent fit with root-finding (≤1e-6). Edge limits (75–125%) apply only at the LLOQ/ULOQ, interior levels need 80–120%, and a replicate CV above the limit fails a level. A plate with fewer than max(6, 75%) passing levels withholds every unknown; unknowns below the LLOQ are withheld even inside the standard response range. In the synthetic fixture the 0.5 ng/mL level (recovery 109%, CV 26.2%) fails, so the LLOQ is 1 ng/mL; unknown estimates are unchanged. |
| ELISA figure and report | New run [elisa-synthetic-qc](../runs/elisa-synthetic-qc/report.html), visually inspected | Quantified unknowns plotted at their in-well concentration; LLOQ–ULOQ shaded; standard-QC table; sample CV%. The 0.5.0 run still renders. |
| `comparison.csv` diagnostics | Same test file | Semicolon-joined like the other modules (was a Python list repr). |
| Isolated wheel install | [clean_environment.json](clean_environment.json) | `validate_clean_environment.py` now covers the ELISA and two-group runs: the 0.5.1 wheel matches the source exactly, and committed runs agree (max relative difference 0 for reported unknowns and t statistics). |

The full suite has 64 passing tests. The QC defaults are a common convention,
not a substitute for an assay's validated SOP. Failing standards are never
excluded automatically.

## ELISA calibration and two-group comparison — 2026-09-23 (0.5.0)

Two independently routed specialists were added. The full regression suite has
**58 passing tests**, including eight end-to-end and boundary tests for the new
modules. New synthetic [ELISA](../runs/elisa-synthetic/report.html),
[paired](../runs/groups-paired-synthetic/report.html) and
[Welch](../runs/groups-welch-synthetic/report.html) reports are saved with source
configurations and hashes. The tests check known synthetic concentrations after
dilution, inverse recovery for decreasing calibration, refusal to extrapolate a high OD, failure of ambiguous units or unknown
concentrations, pairing/duplicate-unit rejection, agreement with SciPy t
statistics and p values, report generation, and render-only hash invariance.
All three runs reproduce exactly from a built 0.5.0 wheel installed in a fresh
Python 3.13 environment on the same macOS host; see
[clean-environment record](clean_environment_0.5.0.json).

The ELISA inverse is a point estimate conditional on its fitted 4PL standard
curve. It reports no unknown-concentration CI, matrix recovery, blank correction,
plate-to-plate effects or hook-effect assessment. The group specialist supports
one predeclared two-sided contrast: Welch for disjoint independent units, or a
paired t test for complete pairs. It does not handle multi-group ANOVA, more
than two repeated conditions, multiplicity or covariates. All new fixtures are
synthetic; no team's measured ELISA plate, Prism project or external numerical
benchmark has been tested. Report structure and embedded figures were checked
by tests; a visual/browser interaction review of the new reports has not yet
been completed. A passing test suite does not establish assay suitability.

Statistical method references: [SciPy Welch test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html),
[paired t test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html),
and [bounded least squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).

## Historical 0.4.0 validation

## 4PL extensions — 2026-09-23 (0.4.0)

Changes from the review of the 0.3.0 module: fixed plateaus, relative (1/Ŷ²)
weighting, per-sample summaries across independent experiments, reference-vs-test
relative potency with parallelism and shared-C50 F tests, a named
direction-mismatch diagnostic, figures in input dose units, and zero-dose
controls in a separate linear panel. The full suite has 50 passing tests,
including 10 in `tests/test_dose_extensions.py`. No Prism numerical equivalence
or validation on the team's own antibody experiments is claimed.

| Check | Evidence | Result / boundary |
|---|---|---|
| Unchanged behavior for existing settings | Re-fit of the three 0.3.0 runs | C50, Hill, plateaus, objective and interval endpoints agree to ≤3.5e-8 relative; all states identical. |
| Fixed plateaus and relative weighting | `tests/test_dose_extensions.py` | Match independent all-parameter `least_squares` fits (C50 within 1e-6 in log10, objective within 1e-7). With fixed plateaus the independent profile SSE at each reported endpoint equals the threshold with reduced df. |
| Relative potency, independent refit | [dose_numerical_benchmarks.json](dose_numerical_benchmarks.json) (`relative_potency_benchmark`) | New run [dose-potency-synthetic](../runs/dose-potency-synthetic/report.html): 6 curves, 3 comparisons, relative weighting. An all-parameter refit without the plateau shortcut reproduces C50, objective, RP and the parallelism F within 1.0e-08 relative. Tests also check the RP interval endpoints against an independent profile. |
| Coverage, extended | [dose_coverage_simulation.json](dose_coverage_simulation.json) | Nominal 95%, 200 simulations each, acceptance bound fixed before running: relative weighting under 8% CV noise 191/200 = 95.5%; the same data unweighted (misspecified, no bound) 187/200 = 93.5%; fixed 0/100 plateaus 193/200 = 96.5%; RP with true RP=2, parallel curves 179/186 = 96.2% (14/200 comparisons withheld by the α=0.05 parallelism test). Original iid scenario unchanged at 94/100. |
| Isolated wheel install | [clean_environment.json](clean_environment.json) | 0.4.0 wheel in a fresh environment reproduces the editable-source results exactly for all four 4PL runs; committed runs agree within 1.3e-9. |
| Report | Static check of the four 4PL reports | No external resources; both themes present for every curve and comparison; axes in input units; zero-dose panels exactly where zero doses exist. The three 0.3.0 runs were re-rendered (scientific hashes verified). The interactive browser check (theme switching, mobile width) was **not** repeated for 0.4.0. |

Limits: the parallelism F test is a significance test, not the equivalence
test USP <1032>/<1034> recommend for release assays. The RP interval from one
plate excludes between-plate variability; quote the geometric mean across
independent experiments. All fixtures are synthetic or public examples.

Local build artifacts: during this work `dist/` was cleared, deleting the
0.1.0–0.3.0 wheels and sdists whose hashes earlier manifests record. Their
sources remain in Git (0.2.0 at `13aa5c4`, 0.3.0 at `50b2de9`), but rebuilt
wheels are not byte-identical (0.3.0 rebuilt 9d6fe5b0… vs recorded 0a7b71b3…)
because wheels embed timestamps. The recorded hashes can no longer be checked
against local files.

## Relative 4PL EC50/IC50 module — 2026-09-23 (0.3.0)

The new `dose-response` specialist fits one 4PL per declared curve with an
explicit endpoint, observed direction, input scale and unit. The full test suite
has 40 passing tests, including 9 dose-response cases. This is a scoped model
release, not GraphPad Prism numerical equivalence or validation on the team's
own antibody experiments.

| Check | Evidence | Result / boundary |
|---|---|---|
| Method and input behavior | `tests/test_dose_response.py` | Synthetic EC50 and IC50 truth recovery; endpoint independent of measured direction; molar/mass/log-dose handling; invalid data and exclusion reasons; flat and out-of-range fits withheld; render invariance |
| Independent numerical refit | [dose_numerical_benchmarks.json](dose_numerical_benchmarks.json), `scripts/validate_dose_benchmarks.py` | Four curves agree with a simultaneous four-parameter SciPy `curve_fit` formulation; maximum relative difference 5.83e-7 for midpoint, slope and objective under matched loss. One profile interval is checked against a separate full-nuisance optimization in tests. |
| Interval smoke simulation | [dose_coverage_simulation.json](dose_coverage_simulation.json) | One fixed additive-iid 4PL design: 94/100 nominal 95% profile-F intervals covered truth; 100/100 closed. This is not coverage evidence for heteroscedasticity, correlated wells, plate effects, or model mismatch. |
| Public example and source | [source manifest](../fixtures/dose_public/source_manifest.json) | Pinned MIT-licensed IC50 Studio example, source files and hashes retained. Its dose unit is unspecified (`source_unit`) and its own IC50 targets response Y=50; this release reports the *relative fitted-plateau midpoint*, so no source-software IC50 equivalence is asserted. |
| Offline report | [dose_browser_checks.json](dose_browser_checks.json), `output/playwright/dose-mobile.png` | Two curves, both themes, mobile width, embedded images and offline PDF download checked in a browser document. |
| Isolated wheel install | [clean_environment.json](clean_environment.json) | Built wheel is installed in a fresh environment on the same macOS host; current source and wheel dose outputs and intervals match exactly. Earlier equilibrium runs differ from current fixed code by at most 5.83e-7 relative in KD and remain historical artifacts. |

The synthetic fixtures mimic antibody-project readouts but are generated, not
measured biology. A real assay still needs the titrated species, dose units,
response definition, controls, normalization history and independent-experiment
structure checked before interpretation. Current CI assumes unweighted
response-scale errors and the symmetric four-parameter model. Absolute Y=50
interpolation, 5PL, constrained plateaus, weighted fits, and cross-experiment
significance tests remain outside the validated module.

## Historical 0.2.0 validation

## Independent review addendum — 2026-09-23

Code changed after the 0.2.0 records below; package version not bumped. Runs
created before this change keep their recorded `implementation_sha256`.

| Check | Evidence | Result / boundary |
|---|---|---|
| Equilibrium profile-F endpoint defect (fixed) | `test_profile_upper_endpoint_not_closed_by_numerics` | pM concentration design with high noise: 28 of 41 fits marked reportable had a closed upper KD endpoint (0.2–1 M) where the exact linear-nuisance profile stays open. Causes: amplitude bound fixed at 1e8 response scales, and the profile nuisance fit stalling on a flat plateau from the best-fit start. Amplitude bound now scales with the KD domain; each profile point also starts from the exact linear solution; a bound-limited profile fails instead of closing. The new test fails on the previous code. Public and synthetic runs: KD within 6e-7 relative, endpoints within 1e-12, all states unchanged; BindCurve benchmark unchanged. |
| Kinetic fitter, independent ODE refit | Full-parameter `solve_ivp` + `least_squares`, Octet 300 s window | Refit from the AgenticPrism solution: kon/koff agree within 1e-5 relative; the AgenticPrism objective was lower than a cold-start independent fit. |
| Kinetic bootstrap coverage | [kinetics_coverage_simulation.json](kinetics_coverage_simulation.json), `scripts/validate_kinetic_coverage.py` | Nominal 95%, 200 simulations each. iid noise: kon 93.0%, koff 94.5%, KD 93.5%. AR(1) φ=0.9 noise: kon 76.5%, koff 83.0%, KD 80.5%. With strongly correlated residuals (as flagged on the public Octet data) the interval is too narrow. One design only. |

2026-09-23. The new kinetic module is validated for the specific model and file
layout below. The equilibrium table that follows records the earlier 0.1.0
validation and still applies to that module; the regression suite had 28
passing tests in total. No numerical equivalence to GraphPad Prism is claimed.

| Kinetic gate | Evidence | Result / boundary |
|---|---|---|
| Schema, truth recovery, failure and render behavior | `tests/test_kinetics.py`; full `pytest` | 8 kinetic tests, 28 total pass; synthetic truth and joint bootstrap intervals, phase/units/flat trace, reference/baseline and shared Rmax, Octet import hashes, render invariance |
| Public Octet reference | [kinetics_numerical_benchmarks.json](kinetics_numerical_benchmarks.json) | Five processed RED384 sensorgrams, 300/600 s windows, published TitrationAnalysis output; maximum rate/KD relative difference 0.0198% under matched settings |
| Fit-window sensitivity | Same benchmark and [300 s](../runs/octet-300s/report.html) / [600 s](../runs/octet-600s/report.html) reports | Same experiment yields KD 103.3 versus 51.75 nM (ratio 1.996); strong time-correlated residuals noted. Neither result is an independent replicate or a universal affinity value. |
| Isolated wheel install | [clean_environment.json](clean_environment.json) | Wheel in fresh Python environment reproduces both kinetic fits and paired bootstrap intervals exactly on the same macOS host; no cross-platform claim |
| Browser and offline export | [kinetics_browser_checks.json](kinetics_browser_checks.json), `output/playwright/kinetics-mobile.png` | Both themes, mobile width, embedded images, theme switch and byte-identical PDF download in an offline browser document |
| Source and adapter | [source manifest](../fixtures/kinetics_octet/source_manifest.json), `tests/test_kinetics.py` | Pinned public source files and imported snapshots SHA-256 checked; only the declared `Time1/Data1` Octet text layout is accepted |

The kinetic model is ideal independent-cycle 1:1 with shared kon/koff, per-sensor
Rmax in the public comparison, unweighted response-scale least squares, and a
conditional phase-wise block bootstrap. Residual autocorrelation means the
interval should not be read as a complete experimental uncertainty statement.
No user's actual instrument files, native `.frd`, Biacore project, native Octet
software output, Prism project, complex kinetic model, or systematic interval
coverage benchmark has been validated. The public file is already processed by
the source authors; its exact upstream subtraction history is not re-inferred.

## Historical equilibrium validation

2026-09-22. This validates the specific implemented equilibrium workflow and
declared fixtures, not numerical equivalence to GraphPad Prism or suitability of
an arbitrary biological experiment.

| Gate | Evidence | Result / boundary |
|---|---|---|
| Skill structure | Official `skill-creator/scripts/quick_validate.py` on both skills | Passed |
| Regression/behavior tests | [pytest.xml](pytest.xml), `tests/test_analysis.py` | 20 tests passed |
| Public historical refit | Same test suite | 25 curves, log residuals, zero-control baseline; KD relative tolerance 2e-4 |
| Independent package | [numerical_benchmarks.json](numerical_benchmarks.json), [per-curve comparison](bindcurve_comparison.csv) | 25 same-linear-loss comparisons with BindCurve 0.2.0; max KD relative difference 0.0001267623 (0.0126762%) |
| Prediction/objective agreement | Same benchmark | Max prediction difference relative to response scale 1.9430e-5; max objective relative difference 3.1922e-9 |
| Independent interval calculation | Test suite and [independent review](independent-review.md) | Linear nuisance parameters solved analytically for an independent profile-F reference; KD and interval endpoints agree within 2e-5 relative tolerance |
| Interval coverage smoke benchmark | [coverage_simulation.csv](coverage_simulation.csv) | 98/100 covered; seed 48219; acceptance 89..100 fixed in script before execution. One concentration design with additive normal noise and fitted baseline only; not universal coverage evidence |
| Isolated environment | [clean_environment.json](clean_environment.json) | Built wheel installed in a new Python environment; all 25 estimates, endpoints and states matched within 1e-8, same macOS host |
| Independent Skill workflow | [independent-review.md](independent-review.md) | Correct routing/scope, run from unrelated working directory, explicit unsupported analysis rejection, repeat summary and render invariance |
| Browser/offline interaction | [browser_checks.json](browser_checks.json), `output/playwright/` screenshots | Both themes, mobile width, saved results unchanged, embedded images and matching PDF download in offline document |
| Export dimensions and links | [artifact_checks.json](artifact_checks.json) | SVG/PDF physical dimensions, PNG pixels/DPI, theme data invariants, hash integrity and local documentation links |

Browser method: the CLI blocks direct `file://` navigation. The report was opened
on loopback HTTP, then its complete HTML loaded into a fresh `about:blank` document
with the browser context offline. Embedded images, switching and PDF download
were checked there. This is not a claim that native file navigation was tested.
The downloaded PDF is compared byte-for-byte by SHA-256 with the exported file.
The initial localhost preview had an incidental favicon 404; the renderer now
embeds an empty data favicon. Page JavaScript errors are collected separately.

Manual visual review covers the desktop overview, curve/residual cards, mobile
layout, zero-concentration panels, repeated-observation SD bars, flat failure and
out-of-range cards. Font fallbacks are recorded per run. Figures use Latin axis
labels; Chinese report text is rendered by the browser. This release has not been
visually validated on Windows/Linux or with arbitrary long user-provided labels.

Two issues discovered by independent forward testing were fixed and regression
tested: valid identifier `NA` being interpreted as missing while reading saved
predictions, and fit-only sensitivity results incorrectly requiring an interval
to avoid a limited status.

Reproduce:

```sh
.venv/bin/python -m pytest -q --junitxml=validation/pytest.xml
.venv/bin/python scripts/validate_numerics.py
uv build
.venv/bin/python scripts/validate_clean_environment.py
.venv/bin/python scripts/validate_artifacts.py
```

Browser checks use `scripts/validate_browser.py <playwright-cli-or-wrapper>` with
a live `agentic-prism` CLI browser session and a localhost server on port 8767.
The original public source hash is verified by `scripts/prepare_fixtures.py`.

Remaining deliberate limits: no Prism numerical reference project, no kinetic or
ELISA modules, no global shared-parameter fits, no depletion/multivalent model,
no curve confidence band, no baseline/amplitude intervals, and no propagation of
control-derived fixed-baseline uncertainty. Only declared comparable independent
experiments can be summarized; the public author dates are not treated as such.
