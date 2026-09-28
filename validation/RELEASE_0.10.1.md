# Release evidence — local development 0.10.1

Branch `dev/0.10.1` from `dev/0.10.0` (`f8c83d9`). Implemented by the reviewing assistant (Claude) at the
user's request. Scope: the two CMC items 0.10.0 did not start. Local commits only; no push or publication.

## Scope

| Item | Status |
|---|---|
| Comparability / biosimilarity (`comparability` Skill) | **Done**: TOST, quality range, descriptive; one attribute, lot as unit |
| Tolerance intervals and capability (`specifications` Skill) | **Done**: normal exact/Howe/one-sided, order-statistic, Pp/Ppk, Cp/Cpk |
| Arrhenius; interpretation facts for the six oldest modules | Not started |

## Methods

- **Comparability**: replicate measurements averaged per lot. TOST on the difference of lot means (test −
  reference), Welch (default) or pooled, 100(1 − 2α)% interval inside ±margin; the margin is absolute or a
  declared multiple of the reference-lot SD (the convention of the FDA 2017 draft, withdrawn 2018), with its
  sampling uncertainty disclosed. Quality range: reference mean ± k·SD with k declared, fraction of test lots
  inside and an optional required fraction. Tier, method, margin and k must be declared with a source.
- **Tolerance intervals**: two-sided exact factor from the Odeh integral (Eberhardt, Mee & Reeve 1989), or
  Howe in the NIST/SEMATECH 7.2.6.3 form (without the second-order correction of R `tolerance` method "HE");
  one-sided exact noncentral-t factor; order-statistic intervals with the largest exact-confidence rank, or
  the minimum n when none exists.
- **Capability**: Pp/Ppk from the overall SD, Cp/Cpk from the pooled within-subgroup SD; Cp/Pp chi-square
  intervals and the Cpk/Ppk normal approximation quoted by NIST/SEMATECH 6.1.6 (recommended for n ≥ 25).
  Subgroup results disclose within- versus overall variation and non-independent measurements.

All numerics are Python; R is a validation oracle only.

## Agreement with R and published values

[cmc_0101_benchmarks_0.10.1.json](cmc_0101_benchmarks_0.10.1.json): R 4.6.0, tolerance 3.0.0. 115 fields:
two-sided exact and one-sided factors over 45 (n, content, confidence) combinations each, order-statistic
ranks, TOST intervals and p values (`t.test`), capability indices and intervals, and fixture tolerance
limits (`normtol.int`). Maximum relative difference 3.3e-9 (R's EXACT factor integrates with 200 points);
tolerances fixed before running. Published values reproduced (NIST/SEMATECH e-Handbook): k = 2.217 (Howe,
n = 43, 90%/99%) and 1.8740 (one-sided), minimum n = 46 and 473 for extreme-value intervals (7.2.6.4), and
Cp = 1.0, Cpk = 0.6667 for the 6.1.6 example (point estimates only; the page prints no interval example).
No published TOST worked example was reproduced.

## Calibration (seed 20261016, 1000 datasets per row)

Script committed before running (`5256e4a`); record [cmc_0101_calibration_0.10.1.json](cmc_0101_calibration_0.10.1.json).
Tolerance content is computed exactly from the true distribution.

| Check | Result (MCSE) | Criterion | Outcome |
|---|---|---|---|
| Normal two-sided exact: n=10 90%/95%; n=30 99%/95%; n=5 95%/90% | 94.9 (0.7), 93.9 (0.8), 91.7 (0.9)% | ≥ conf − 2 SE | pass ×3 |
| Normal two-sided Howe n=10 90%/95% | 95.3 (0.7)% | same | pass |
| Normal one-sided exact n=10 95%/95% | 94.4 (0.7)% | same | pass |
| Nonparametric: lognormal n=120, exponential n=60 (90%/95%) | 98.3, 98.5% | same | pass ×2 (conservative, discrete ranks) |
| Ppk 95% interval: n=10, 25, 50 (true 1.0); n=25 (true 1.33) | **93.6**, 94.3, 95.4, 94.8% | ≥ 93.62% | **n=10 miss**; others pass |
| Pp 95% chi-square interval n=10 | **93.6%** | ≥ 93.62% | **miss** |
| TOST type I at a fixed margin: equal SD; test SD ×2 | 4.5 (0.7), 3.1 (0.5)% | ≤ 6.38% | pass ×2 |

Two misses, both n = 10 and from the same datasets. The Pp interval is exact by construction, so a post-hoc
check was run afterwards and is labelled as such ([supplement](cmc_0101_capability_supplement_0.10.1.json),
seed 20261017, 100,000 datasets): Pp 94.9% and Ppk 95.4% at n = 10 (MCSE 0.07). The registered misses are
Monte Carlo error; they stay recorded as misses and are disclosed in the Skill and facts.

Descriptive rows: TOST with an estimated 1.5 × SD margin at the true 1.5σ difference, 5.8%; quality range with
identical products (10/6 lots), all test lots inside 88% (k = 3) and 61% (k = 2); with a 1-SD shift, 72% and 35%.

## Live-agent scenarios

[agent-scenarios/0.10.1/RESULTS.md](agent-scenarios/0.10.1/RESULTS.md): SP1 (tolerance interval as a
specification), SP2 (Cpk only), CP1 (failed quality range as "not biosimilar"), CP2 (replicates as lots) and
a CMC4 rerun — 5/5 pass. SP2 also caught a wrong design description in the subgroup fixture; corrected
(data unchanged).

## Preservation, tests and installation

- Full suite 523 passed (24 new). Three changed Skills pass quick validation.
- [Legacy check](legacy_fixture_checks_0.10.1.json) against `dev/0.10.0`: 89 configurations, 745
  byte-identical artifacts, no differences.
- Wheel built from `git archive` of `fd0e2c9`: [clean environment](clean_environment_0.10.1.json) passed
  (90 facts files byte-identical to source); [install smoke](install_smoke_0.10.1.json) 98 configurations.
  [Release checks](release_checks_0.10.1.json). macOS only; no local paths in the new records.

## Boundaries

One attribute per comparability run; no totality-of-evidence judgement. Tolerance intervals describe data and
do not set specifications. Normal methods assume Gaussian independent units; capability intervals below
n = 25 are approximations. No JMP, Minitab or SAS equivalence is claimed.
