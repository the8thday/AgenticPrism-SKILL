# 0.13.5 — multivariable regression and GLM

Development baseline: `cc76b2f` (0.13.4 plus optional consultation).
Implementation and schema scope: plan `plans/RELEASE_0.13.5_PLAN.md` (prespecified in `fb28807`; development repository only),
[Skill](../skills/multivariable-regression/SKILL.md),
[contract](../skills/multivariable-regression/references/input-and-model.md).

Independent-unit Gaussian OLS/HC3, binary or grouped-binomial logit, Poisson log
with exposure/offset, and explicit Pearson-scale quasi-Poisson. Declared
continuous/factor predictors, references and pairwise interactions; conditional
Wald tests and supported mean predictions. No automatic variable/family selection.

## Five evidence gates

| Gate | Observed scope |
|---|---|
| Established implementation | Base-R `lm/glm` and explicit HC3 sandwich: 11/11 fixtures passed, including independent design coding, coefficient/SE/CI, covariance, term tests, fitted means and supported predictions. Includes audit-only fits; see [record](regression_benchmarks_0.13.5_initial.json) |
| Published worked example | Star98 official statsmodels example: 21 printed six-decimal coefficients and deviance matched within printed rounding tolerances. Fixed-binomial inference remains withheld for extra-binomial dispersion. Gaussian/HC3 and Poisson/quasi worked-example gates UNMET. [Record](regression_published_0.13.5.json) |
| Calibration | Prespecified in `fb28807`: seed 131350001, 1000 draws for each of 12 scenarios, all draws retained, no execution errors; 15/15 bounded criteria passed. Stress limitations below remain part of evidence. [Record](regression_calibration_0.13.5.json) |
| Live-agent misuse | 6/6 fresh-context forward tests passed, implementer-scored: Gaussian/Logistic/Poisson execution, separation, overdispersion and unknown percentage denominators. Actual responses and inspected artifacts retained in [record](agent-scenarios/0.13.5/RESULTS.md); no independent-certification claim |
| Interpretation facts | Results-only extraction before manifest, source hashes, parent/prediction withholding and required disclosures; tested without refitting |

## Calibration limits retained

Correct-design reportable coverage: Gaussian 0.93662 (994 reported), Gaussian
HC3 0.959, heteroscedastic HC3 0.958, binary 0.955, grouped binomial 0.94461
(993 reported), Poisson 0.94439 (989 reported), proportional-variance quasi-Poisson
0.948. These are conditional observations, not uniform post-selection guarantees.

- Heteroscedastic classical stress withheld 997/1000, but the remaining three
  fits covered only 1/3 (MCSE 0.272). The variance gate is not a validity certificate.
- Negative-binomial stress: fixed Poisson withheld every fit; quasi-Poisson
  covered 0.936 (MCSE 0.007740), with no registered stress bound. It is not a
  general substitute for a negative-binomial model.
- Rare binary stress reported 937/1000, with conditional coverage 0.98399; retain
  all 63 withheld fits and unconditional results. No sparse-Wald guarantee.
- Perfect separation withheld 1000/1000. This calibration does not cover every
  quasi-separation, leverage, collinearity or misspecified-mean configuration.

Influence and sparse-information cutoffs are review heuristics, not automatic
exclusions. Dispersion/Breusch-Pagan p < .01 gates are approximate, especially
for sparse counts. HC3/quasi options require explicit choice and cannot fix
dependence, post-treatment conditioning, or causal confounding.

## Engineering checks

Initial targeted tests: 25 passed, 7 failed because the test harness passed argv
to a zero-argument CLI main, and demanded byte equality after a CSV round trip
of log exposure. Corrected the harness invocation and used 1e-12 numerical
tolerance for exposure/log-offset equivalence; 32/32 then passed. No regression
estimator change was needed for those failures.

A subsequent schema review added rejection of `Intercept` predictor names and
ambiguous encoded factor/interaction names. The final targeted suite passed
33/33; all 13 valid fixture compute results were exactly unchanged from the
calibration/agent-frozen implementation ([record](regression_schema_guard_0.13.5.json)).
Original calibration provenance and fresh-agent input hashes are retained;
new evidence status text is generated after actual behavioral scoring.

All 180 legacy configurations matched the baseline across 1,343 scientific
artifacts ([record](legacy_fixture_checks_0.13.5.json)); unchanged methods were
not recalibrated. The three changed/new Skills pass structural validation
([record](skill_checks_0.13.5.json)).

The full suite passed 831 tests with 8 warnings in 665.16 seconds. It began
before the final name-collision guard; the final 33-test targeted suite includes
that extra test. See [test scope](pytest_0.13.5.json). A Gaussian HC3 coefficient
figure was also inspected for readable labels/intervals.

Clean git-archive build from implementation commit `fadb6a4` passed: all 104
Python source files match the wheel ([build record](archive_build_0.13.5.json)).
Initial offline cache failures and the local-cache-only repair are preserved in
[build prerequisites](build_prerequisites_0.13.5.json); no source change or network
download was needed.

An isolated wheel-only environment outside the repository, with locked
dependencies and no PYTHONPATH, passed 23 configurations with 215 scientific
artifacts exactly equal to source. Six estimated/withheld states passed report
re-rendering and hash verification without scientific changes
([clean install](clean_install_0.13.5.json)). These checks concern this macOS host;
they do not establish cross-platform equivalence. Final verification records
and documentation follow the build commit; runtime source remains unchanged.
Runtime evidence is generated by `scripts/update_evidence_0135.py`; it copies
the original rows and notes, without refitting or replacing inconvenient results.

No Prism numerical comparison, regulatory-system validation, predictive external
validation or cross-platform claim. PK/PD remains deferred. This release stays local.
