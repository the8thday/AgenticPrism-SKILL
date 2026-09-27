# Release evidence — development 0.7.0

Baseline: committed 0.6.0 (`6b173e9`, branch `dev/0.7.0`), 94 passing tests,
verified before development on 2026-09-26. Historical records are unchanged; this
file states only what was run for 0.7.0. All new fixtures are synthetic
(`scripts/prepare_0_7_fixtures.py`, seed 20260926). No team-measured data, Prism
project, Biacore export or vendor software output was available.

## What was added and how it was checked

| Module | Checks | Observed result | Boundary |
|---|---|---|---|
| Multi-group one-way comparison (`multigroup.py`) | `tests/test_multigroup.py` (15); [error-rate simulation](multigroup_error_rates.json) | Classic F equals `scipy.stats.f_oneway` (rel 1e-12) and Welch F equals `f_oneway(equal_var=False)` (rel 1e-10). Tukey-Kramer and Games-Howell p values and intervals match an independent studentized-range computation (p rel 1e-4, CI rel 1e-5). Dunnett p values agree with an independent 400 000-draw Monte Carlo of max\|t\| within MC error, critical value within 1%. Holm-adjusted p and Bonferroni intervals match hand computation. | Normal, independent units; four groups; one design per scenario |
| Error rates, α = 0.05, bounds fixed before running | Same simulation | Global null, equal SD, n = 6: omnibus rejection 0.041–0.060 and family-wise error 0.038–0.060 for all four families; simultaneous coverage 0.940–0.962; all passed. Alternative (means 0/0.5/1/2): coverage 0.944–0.970; passed. **Null, SDs 1/1/1/3 with n 8/8/8/4:** pooled Dunnett/Tukey family-wise error 0.195/0.193 (misspecified contrast, no bound); Welch ANOVA rejection 0.069 (bound 0.064, **failed**); Holm-Welch FWER 0.057 (passed on FWER, failed on the omnibus bound); Games-Howell FWER 0.077 (bound 0.075, **failed**). | Response: diagnostic `small_group_welch_procedures_may_be_liberal` when a Welch-design group has < 6 units. Bounds were not changed after the run |
| Plate-reader grid import (`plate_import.py`) | `tests/test_plate_import.py` (11) | Every imported response equals its grid cell text; the blank-role well that has a reading is skipped and counted; source hashes match; the ELISA target quantifies to truth (±8%); the dose target loads through the 4PL loader. Refusals: non-numeric read, duplicate or missing source, wrong grid shape, unknown format, duplicate plate ID, existing output. | Generic 8×12 / 16×24 grid CSV only; no vendor export parser |
| ELISA unknown intervals (`calibration.py`) | `tests/test_elisa_extensions.py` (13); [coverage simulation](elisa_interval_coverage.json) | Sample intervals equal an independent full-parameter `curve_fit` covariance plus analytic delta-method gradient (rel 2e-4). Coverage of the true sample concentration for nominal 95% intervals, 300 plates each, bound ≥ 0.95 − 2 SE fixed before running: 4PL unweighted 0.933/0.950/0.940; 4PL relative 0.955/0.976/0.947; 5PL relative (asymmetric truth) 0.937/0.947/0.961, all passed. Misspecified contrast (4PL on asymmetric truth): 0.889/0.993/0.983. The highest level was reported in only 32–42% of relative-weighting plates, because it often lies above the plate ULOQ. | Conditional on one plate's calibration model; no plate, matrix, pipetting or dilution error |
| ELISA 5PL | Same test file | Matches an independent 5-parameter `least_squares` fit: objective no worse, asymmetry g rel 1e-3; inverse(predict(x)) = x (rel 1e-8). On the symmetric 4PL fixture it flags `asymmetry_not_supported_over_4pl` and still recovers the 30 ng/mL sample (±8%). | No profile CI for the 5PL midpoint |
| ELISA dilution linearity | Same test file; 5PL dilution fixture | The linear sample passes using in-range dilutions, and its 95% interval contains the 800 ng/mL truth on both plates. The hook-like sample (signal suppressed to 55% at 1:8) fails, with minimum recovery below 80%, and is flagged as a possible hook/matrix pattern. The default policy withholds samples with out-of-range dilutions; `required` withholds single-dilution samples. | Within one plate; pattern flag, not a diagnosis |
| Equivalence parallelism and RP acceptance (`equivalence.py`) | `tests/test_potency_equivalence.py` (12); [operating characteristics](potency_equivalence_operating_characteristics.json) | Hill-ratio and plateau-difference Wald intervals equal an independent per-curve `curve_fit` covariance with pooled variance (rel 1e-4). Pass rate, 400 comparisons each, margin [0.8, 1.25], 90% intervals: true ratio 1.0 → 1.000; ratio on the margin 1.25 → 0.0575 and 0.8 → 0.0425 (bound ≤ 0.072, passed); 1.5 → 0.000. Hill ratio 1.8 is rejected and RP withheld; the fixture's RP ≈ 1.94 correctly fails the illustrative [0.8, 1.25] acceptance. | One simulated design; margins illustrative |
| Single-cycle kinetics and double referencing | `tests/test_single_cycle.py` (13); [coverage simulation](single_cycle_coverage_simulation.json) | `schedule_response` equals an independent `solve_ivp` solution of the 1:1 ODE across three injections (max abs error < 1e-8, unit Rmax); one step reproduces `unit_response` bit for bit. The double-referenced fixture recovers kon 2.000e5, koff 5.000e-4, KD 2.50 nM and Rmax 100/80. Without referencing, drift and bulk degrade the fit. iid coverage, 150 fits, bound ≥ 0.914: kon 0.927, **koff 0.913 (failed by one fit)**, KD 0.940. AR(1) φ = 0.9: 100/100 withheld by the reliability gates. A first run at 2 s sampling produced no intervals, because 30-point intermediate dissociations are shorter than 2 × block length. | One design; percentile bootstrap slightly undercovers, as for multi-cycle data in 0.2.0 |
| Octet blank-cycle import | `tests/test_single_cycle.py` | Blank and blank-reference traces are carried as columns with source hashes and roles; a partial double-reference declaration is refused. | Tested on public traces standing in as blanks, not real blank cycles |

## Suite, installation and example runs

- [Full suite](release_checks_0.7.0.json): **158 passed** (94 before 0.7.0, plus
  64 new: multigroup 15, plate import 11, ELISA extensions 13, potency
  equivalence 12, single-cycle 13). All six Skill folders pass
  `quick_validate.py`.
- [Isolated wheel](clean_environment_0.7.0.json): the 0.7.0 wheel in a fresh
  Python 3.13 environment reproduces the editable source exactly for all earlier
  fixtures and the five new configs, and runs `import-plate`. Committed earlier
  runs agree within 1.3e-9 relative.
- [Installed-runtime smoke](install_smoke_0.7.0.json): analyze, verify,
  re-render and verify for ten fixtures (including all new modules) plus plate
  import, from an unrelated working directory against the same wheel
  (identical SHA-256). macOS only.
- Example runs, kept as evidence: [multi-group Dunnett](../runs/multigroup-dunnett-0.7.0/report.html),
  [5PL ELISA with dilution linearity](../runs/elisa-dilution-5pl-0.7.0/report.html),
  [single-cycle double-referenced kinetics](../runs/kinetics-single-cycle-0.7.0/report.html),
  [equivalence parallelism and RP acceptance](../runs/dose-potency-equivalence-0.7.0/report.html).
  PNG figures for the first three and the equivalence report text were
  inspected; no browser interaction QA was performed.

## Unchanged behavior

[Regression record](regression_0.7.0.json): reruns of the committed 0.6.0 Octet
reference-library runs are **byte-identical** (fit results, curve parameters,
bootstrap samples, predictions, window sensitivity). ELISA unknown and sample
tables keep identical values in every existing column and add `ci_low`,
`ci_high`, `ci_status` and linearity columns. Two-group comparison tables are
identical apart from the 0.5.1 diagnostics text format.
`validate_kinetic_benchmarks.py` (max relative difference 0.000198 vs published
TitrationAnalysis) and `validate_dose_benchmarks.py` (5.8e-7; potency 1.0e-8)
reproduce their historical values.

## Agent sessions

Nine misuse scenarios were run as live headless Claude Code sessions (`claude-opus-5-5`), with full traces; see [agent-scenarios/0.7.0/RESULTS.md](agent-scenarios/0.7.0/RESULTS.md). Round 1: 7 pass, 2 partial. S4 suggested an approximate KD for a limited fit, and S7 ignored an injected CSV instruction but did not tell the user. The router and kinetics Skills were tightened, and all nine were rerun with the final Skills: 9/9 pass. The reviewer was the developing assistant, not an independent person; there was one sample per scenario per round, and the scenarios reuse synthetic fixtures that agents sometimes recognized.

## Not established

- No measured antibody-project data, Prism benchmark, Biacore/other SPR export,
  or real plate-reader vendor export was tested.
- The three-platform CI workflow has not run: the repository has no git remote.
- Mixed-effects intermediate precision, repeated-measures ANOVA, nonparametric
  tests, mass-transport or heterogeneous-ligand kinetics are not implemented.
- Simulation evidence covers the stated designs only; it is not a coverage
  guarantee for other designs.
