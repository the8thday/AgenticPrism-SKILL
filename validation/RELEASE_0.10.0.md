# Release evidence — local development 0.10.0

Branch `dev/0.10.0` from `dev/0.9.3` at `c7d97df`. Local commits only; no push,
publication, main merge, public-tree build or 0.11 work. Historical records and
runs are unchanged. Detailed numbers live in linked JSON records.

## Scope and five gates

| Item | Implementation | Evidence state |
|---|---|---|
| Stability and shelf life, new `stability` Skill | Done, scoped fixed-batch linear Q1E | Partial: agreement/facts done; four calibration misses; printed shelf life not reproduced; live scenarios pending |
| RP across runs, new `potency-assay` Skill | Done, replicated determinations with random run | Partial: agreement/calibration/facts done; published example unavailable; live scenarios pending |
| Comparability / biosimilarity | Not started | No gate claimed |
| Specifications, tolerance intervals, Cpk/Ppk | Not started | No gate claimed |
| Oldest six modules' interpretation facts; Arrhenius | Not started | No gate claimed |

## Methods and boundaries

- Stability: each batch has a linear regression; ANCOVA tests slopes first at
  alpha .25 and intercepts only if slopes can pool. Selected pooled/common-slope/
  separate lines use a common residual mean square. Shelf life is the earliest
  mean confidence-bound crossing: one-sided 95% lower/upper, or two-sided 95%.
  The minimum batch result governs; source-linked direction/specification/storage
  are mandatory. Batch-specific curvature at .05 withholds a linear result.
- The Q1E statistical-analysis decision-tree branches require explicit supporting
  data, continued pattern and commitment-study declarations. Room/no accelerated
  change: min(2X,X+12); refrigerated/no change: min(1.5X,X+6); room/change at 3–6
  months and no intermediate change: min(1.5X,X+6). Other cases are capped at X,
  the shortest observed batch duration. No future-batch tolerance guarantee.
- Potency: natural-log RP = fixed mean + random run + residual, or fixed intercept
  and slope on log nominal RP. Requires ≥3 runs and ≥2 independent determinations
  per run/level; one RP per run cannot identify both components. No shared-reference
  estimates, technical-well pseudoreplication or heterogeneous known-SE meta-analysis.
  REML and MLS/MOVER reuse the unchanged variance primitives. GLS t intervals use
  runs−1 df. Intermediate geometric CV = 100 sqrt(exp(total log variance)−1).
- `dose_runs` imports verified `dose_fit.py` / `equivalence.py` comparisons and
  copies source artifacts. Source-linked suitability is per determination; any
  failure withholds combination. Prespecified single Grubbs screening flags and
  retains every run. No automatic exclusions. Validation criteria are sourced;
  bias/precision point rules and statistical interval supplements are labelled.
  Range needs ≥3 contiguous passing tested levels plus declared global linearity
  equivalence. No interpolation across failing levels or extrapolation.
- Every new analysis writes facts from saved results before the manifest, including
  must-mention evidence and failing items. Rendering does not refit.

All new runtime numerics are **Python**. R is optional and used only as a new-method
validation oracle; no new `r_bridge` method. Existing optional MMRM KR is unchanged.
No software or regulatory equivalence is claimed.

## Agreement and published examples

[Benchmark record](cmc_benchmarks_0.10.0.json): R 4.6.0, lme4 2.0.6, emmeans 2.0.4.
308 regenerated R fields pass: stability `lm`/`anova`, predictions/SE and `emmeans`;
potency REML coefficients/covariance/components via lme4 and independent R MLS/MOVER
quadratic forms. Maximum relative difference 6.895e-7, absolute 3.039e-8; stability
alone 3.492e-13 relative. Fixed tolerances: regression 1e-8 relative, REML/MLS 2e-5;
absolute 1e-8 below reference magnitude 1e-4. No old-method R benchmarks were rerun.

Koleva et al. (2016), *Science, Engineering & Education* 1(1):106–112: 15 data
values transcribed from Table 1; publisher PDF and official Q1E PDF hashes pinned
in [source manifest](../fixtures/stability/source_manifest.json). Ten printed
regression/ANCOVA/SE fields pass half-last-decimal tolerances. Slope F=4.0862069,
p=.05461982; separate batches govern. Q1E mean-bound shelf life is 5.587563 months;
the paper prints 3.16423 under a different interval construction. **Printed shelf-life
reproduction is not met**, and the discrepancy is retained. No matching published
replicated random-run potency worked example was obtained; **that gate is not met**.
Synthetic fixtures are never labelled published reproductions.

## Calibration

[All rows](cmc_calibration_0.10.0.json). Script, bounds, seed 20261012 and scenarios
were committed in `6d3db15` before running; 1000 datasets/row, all evaluable, no errors.
Truth is calculated directly from generating parameters. Coverage floor 93.6216%;
size bounds 22.2614–27.7386% at alpha .25; TOST ceiling 6.3784% at alpha .05.
Shelf rows evaluate the selected regression before curvature-based withholding;
non-overestimation uses the uncapped first crossing, so a duration cap cannot make
coverage pass trivially. The intercept-size row is unconditional on the slope screen.

| Scenario | Metric | Rate % | MCSE, pp | Result |
|---|---|---:|---:|---|
| pooled | Mean lower bound coverage, month 12 | 92.2 | 0.85 | **Miss** |
| pooled | Shelf life ≤ true | 96.4 | 0.59 | Pass |
| pooled | Slope test size | 24.9 | 1.37 | Pass |
| pooled | Intercept test size, unconditional | 23.7 | 1.34 | Pass |
| unpoolable | Mean lower bound coverage, month 12 | 94.5 | 0.72 | Pass |
| unpoolable | Shelf life ≤ true | 94.4 | 0.73 | Pass |
| common_slope | Mean lower bound coverage, month 12 | 93.8 | 0.76 | Pass |
| common_slope | Shelf life ≤ true | 93.5 | 0.78 | **Miss** |
| near_poolability | Mean lower bound coverage, month 12 | 92.1 | 0.85 | **Miss** |
| near_poolability | Shelf life ≤ true | 90.7 | 0.92 | **Miss** |
| balanced_ratio_0 | Combined log RP coverage | 98.4 | 0.40 | Pass |
| balanced_ratio_0 | Intermediate precision interval coverage | 95.2 | 0.68 | Pass |
| balanced_ratio_05 | Combined log RP coverage | 96.3 | 0.60 | Pass |
| balanced_ratio_05 | Intermediate precision interval coverage | 95.7 | 0.64 | Pass |
| balanced_ratio_2 | Combined log RP coverage | 95.4 | 0.66 | Pass |
| balanced_ratio_2 | Intermediate precision interval coverage | 95.5 | 0.66 | Pass |
| unbalanced | Combined log RP coverage | 95.6 | 0.65 | Pass |
| unbalanced | Intermediate precision interval coverage | 95.0 | 0.69 | Pass |
| small_three_runs | Combined log RP coverage | 99.4 | 0.24 | Pass |
| small_three_runs | Intermediate precision interval coverage | 95.8 | 0.63 | Pass |
| intercept_upper_margin | TOST type I at margin | 5.7 | 0.73 | Pass |
| slope_upper_margin | TOST type I at margin | 2.0 | 0.44 | Pass |

18/22 pass. Four stability misses remain scientific limitations of the declared
poolability-selection procedure; no alternate procedure was substituted to erase
them. Each is disclosed in `stability/SKILL.md`, saved results and facts. Potency
mean intervals are conservative at a zero component and with three runs.

## Live-agent scenarios and defects

[Scenario record](agent-scenarios/0.10.0/RESULTS.md): all five **pending reviewer run**.
CMC1/2/4/5 reached Claude CLI but returned `Not logged in` (exit 1); no workaround.
CMC3 preparation exposed a NumPy-boolean serialization defect, fixed with a workflow
regression test. A separate missing-RP path could retain NaN instead of JSON null;
fixed and tested. Report download attributes were aligned with the install-smoke
contract. No unresolved coding defects were observed in completed checks.
Authentication failures and backend tests do not count as live-agent passes.
**Update (reviewer, 2026-09-28):** the reviewer then ran all five scenarios in an authenticated
environment (attempt `claude-run`); all five passed. See the scenario record.

## Preservation, tests, build and install

- [Legacy check](legacy_fixture_checks_0.10.0.json): 76 existing configurations,
  680 byte-identical artifacts against dev/0.9.3, zero differences. No old calibration
  reruns; shared REML, MLS/MOVER, dose-fit and equivalence numerics are unchanged.
- Full pytest **499 passed**, 43 new CMC tests; eight existing pandas warnings.
  [Release checks](release_checks_0.10.0.json): wheel built from `84b89f3`.
- Both new and both changed Skills pass skill-creator quick validation.
- Wheel built with `uv build` from a clean `git archive`; existing wheels retained.
  All 13 new fixtures match source exactly; 79 facts files match byte for byte.
  Isolated install smoke passed 87 configurations (optional R unavailable there). [Clean environment](clean_environment_0.10.0.json), [install smoke](install_smoke_0.10.0.json).
- Scope is this macOS host; no cross-platform claim. Figures visually checked;
  offline HTML structure and links checked, browser visual QA not completed.
