# Release evidence — local development 0.9.3

Branch `dev/0.9.3` from `dev/0.9.2`. Implemented by the reviewing assistant (Claude) at the user's
request, in place of a Codex round. Scope: the remaining 0.9 roadmap items. No push or publication.

## Scope

| Item | Status |
|---|---|
| Method validation (ICH M10-style LBA experiments) | **Done**: new `method-validation` Skill, `analysis_type=method_validation` with six experiments |
| ADA sensitivity and drug tolerance | **Done**: `ada_sensitivity`, `ada_drug_tolerance` in the ADA Skill |
| Efficient nonparametric ADA lower bound | **Done**: two-way (subject × run) bootstrap; the 0.9.2 six-pair bound averaged 28.6% FPR for a 5% target |
| MLS as the default interval for intermediate precision | **Done** (Satterthwaite covered 83–85% in 0.9.1; MLS 95.1–95.9% in 0.9.2) |
| Dynamic cut point; interpretation facts for the six oldest modules | Not started |

Acceptance numbers quoted in the method-validation Skill were checked against the ICH M10 Step 4 text
(2022-05-24, sections 4.2.1–4.2.7 and 7.2; PDF sha256 e306f3b6…). The tool never defaults them: users
declare criteria and their source.

## Methods

- **Accuracy/precision** per QC level: one-way random-run ANOVA with the classical unbalanced n0
  (as `VCA::anovaVCA`), bias %, within-run CV, between-run CV, total error = |bias| + between-run CV,
  per-run results, M10 minimum-design diagnostics. Supplements: t interval on the mean (df = runs − 1),
  MLS interval for MS_B/n0 + (1 − 1/n0)MS_W, Mee (1984) beta-expectation tolerance interval (the method
  adopted by SFSTP accuracy profiles; formula derived from the prediction variance and Satterthwaite df
  and checked in `tests/test_method_validation.py`), accuracy-profile range interpolated on log
  concentration without extrapolation.
- **Dilution linearity**: dilution-corrected accuracy and CV per factor; hook effect when an above-ULOQ
  QC does not read above the ULOQ; the headline fails on a suspected hook.
- **Parallelism**: CV of corrected concentrations per sample; per-sample log slopes against a declared
  margin; common slope with sample intercepts. M10 7.2 notes that trends can pass the 30% CV rule.
- **Selectivity/specificity**: per-source pass/fail, pass fraction per role with Clopper–Pearson interval.
- **Stability**: mean accuracy per condition and level, t interval as a supplement.
- **ADA sensitivity**: per-run stable crossing of the cut point, linear on log PC concentration;
  log-normal upper prediction limit exp(m + t(conf, n−1)·s·√(1 + 1/n)) for a future run; censored runs
  withhold it. **Drug tolerance**: per PC level and run, highest drug concentration still positive,
  linear on log drug concentration; edge results censored.
- **Two-way bootstrap bound**: pigeonhole bootstrap (Owen 2007) of the pooled type-7 percentile,
  declared seed and replicate count (default 2000); refused when subjects × FPR < 3 (fixed before
  calibration).

All new numerics are Python. R is a validation oracle only (no new `r_bridge` method).

## Numerical agreement with R

[method_validation_benchmarks_0.9.3.json](method_validation_benchmarks_0.9.3.json): R 4.6.0, VCA 1.5.2.
188 fields (ANOVA mean squares and components, bias, CVs, bias interval, MLS CV interval, Mee interval
and df for both A&P fixtures; parallelism slope and interval via `lm`; selectivity intervals via
`binom.test`; stability intervals via `t.test`; per-run sensitivities via `approx` and the prediction
limit; the bootstrap's point percentile via `quantile(type = 7)`). Maximum relative difference
3.39e-12, maximum absolute 4.42e-12; tolerance (1e-8 relative) fixed before running. No printed worked
example of the Mee interval could be obtained (Mee's 1984 cement-briquette data were not accessible), so
gate 2 for that supplement is not met; the precision components themselves were checked on the
published EP05 example in 0.9.1.

## Calibration (seed 20261005, 1000 datasets per row)

The script was committed (`2f8292c`) before it was run. Record: [calibration_0.9.3.json](calibration_0.9.3.json).
Targets are exact functions of the generating model; MCSE is reported for each.

| Check | Designs | Result | Criterion | Outcome |
|---|---|---|---|---|
| A. Mee interval mean content (β = .80) | 6×3 with ratio 0 / .5 / 2; 6 runs unbalanced; 3×3 | 0.821, 0.808, 0.805, 0.805, 0.785 (MCSE ≤ 0.006) | within .80 ± (.02 + 2 MCSE) | 5/5 pass |
| B. 90% bias interval coverage | same | 95.2, 92.4, 90.9, 90.8, 92.6% | ≥ 88.1% | 5/5 pass |
| C. 90% MLS between-run CV coverage | same | 90.2, 91.1, 90.0, 89.5, 92.0% | ≥ 88.1% | 5/5 pass |
| D. 95% sensitivity prediction limit, mean conditional coverage | 6 runs 2-fold; 6 runs 3-fold; 3 runs 2-fold | 94.0% (MCSE 0.28); **92.1% (0.40)**; 94.0% (0.44) | ≥ .93 − 2 MCSE | 2/3; **3-fold spacing misses (bound 92.2%)** |
| E. Two-way bootstrap, P(FPR ≥ target), 90% | 5% target, 80×6, run var 0 / .04; 1% target, 300×6, run var 0 / .04 | 96.1, 94.9, 96.2, 96.7% (mean FPR 9.0, 12.3, 1.97, 3.55%) | ≥ 88.1% | 4/4 pass |
| E. Tail-support gate | 5% target, 50 subjects | 1000/1000 withheld | all withheld | pass |

One miss: with 3-fold PC dilutions, log-linear interpolation of a sigmoid is biased enough to lower
coverage to 92.1%. The ADA Skill recommends 2-fold spacing around the cut point. Row A's 3×3 design is
the lowest (0.785); B at ratio 0 is conservative (95.2%) because df = runs − 1 ignores within-run df.

**Descriptive (no criterion)**: probability that a mid QC meets ±20% / 20% / 30% in a 6×3 design —
bias 0: 100% (CV 10%), 91.7% (15%), 72.8% (18%); bias 10%: 99.9%, 86.1%, 66.0%; bias 15%: 92.0%,
66.0%, 42.0%. A single verdict is uncertain; the Skill says so.

## Live-agent scenarios

[agent-scenarios/0.9.3/RESULTS.md](agent-scenarios/0.9.3/RESULTS.md): MV1 (undeclared M10 verdict), MV2
(hook omission), MV3 (parallel trend), AP1 (PC sensitivity as patient LOD) and a V3 rerun after the MLS
default change. 5/5 pass, one session each, scored by the implementing assistant.

## Preservation, tests and installation

- Full suite 456 passed (33 new method-validation/ADA-performance tests plus MLS-default test).
- [Legacy check](legacy_fixture_checks_0.9.3.json) against `dev/0.9.2`: 65 configurations, 573
  byte-identical artifacts; 48 declared differences (three inert ADA defaults, explicit
  `sum_interval=satterthwaite` in the `variance_reml` fixtures, dependent source hashes), other JSON equal.
- Behaviour change: precision configs without `sum_interval` now get MLS. No committed fixture relied on
  the old default; the VCA benchmark fixtures now state Satterthwaite explicitly.
- Wheel sha256 faeffaf3…, built from `git archive` of the committed tree:
  [clean environment](clean_environment_0.9.3.json) passed (11 new fixtures equal to source; 66 facts
  files byte-identical); [install smoke](install_smoke_0.9.3.json) 74 configurations passed.
  [Release checks](release_checks_0.9.3.json). macOS only.

## Boundaries

Inputs are back-calculated concentrations; calibration-curve fitting, run acceptance, carry-over and
incurred-sample reanalysis are not implemented. Calibrations used Gaussian data and fixed curve shapes.
PC sensitivity is not patient-ADA sensitivity. Passing an experiment is not a complete validation or a
regulatory decision. No Prism, Watson, WinNonlin or SAS equivalence is claimed.
