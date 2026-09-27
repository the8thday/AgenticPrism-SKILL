# Validation — AgenticPrism development 0.7.1

> The published distribution contains these records but not `runs/`, `tests/`,
> `scripts/` or the agent-session transcripts. Links to those items below resolve
> only in the development repository; the results they document are summarized
> in the text.

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
