# Release evidence — local development 0.9.2

Branch `dev/0.9.2` starts from `dev/0.9.1`, reviewer commit `237de01`.
Scope stops at 0.9.2. No push, publication, main merge or 0.10 work.

## Scope and gates

| Item | Status / boundary |
|---|---|
| 1. Precision intervals | Done within the stated scope. Orthogonal mean-square MLS; unbalanced correlated MOVER; positive sums; separate component upper bounds. Satterthwaite remains available and the compatibility default. |
| 2. ADA lower cut-point bounds | Done within the stated scope. Parametric mixed-panel moment approximation and independent-pair order statistics; all three tiers. |
| 3. Full method validation | Not started. No accuracy/total error, profiles, linearity/hook, matrix parallelism, selectivity or stability acceptance workflow. |
| 4. ADA sensitivity / drug tolerance | Not started. |
| 5. Dynamic cut point / older facts | Not started. New interval results have facts. |

## Precision: method and evidence

`variance_intervals.py` separates REML point estimates from mean-square interval
centers. Balanced orthogonal strata use integer df and independent mean squares.
Graybill–Wang MLS coverage is approximate even when mean squares are exact.
Unbalanced/fixed-adjusted designs use ordered orthogonal quadratic forms,
E(MS)=M v, REML plug-in moments 2 tr(A_i V A_j V), moment-matched df and
correlated MOVER. This is our application of a published construction, not an
existing VCA unbalanced-MLS algorithm. Signed MS coefficients use signed MOVER.
The declared term order is preserved; no response-driven ordering. Positive
component sums and all one-sided upper bounds are explicit artifacts. CV uses
its declared fixed denominator. No USP <1033> compliance claim is made.

Sources: [Graybill–Wang 1980](https://doi.org/10.1080/01621459.1980.10477565),
[Burdick–Graybill 1984](https://www.stat.cmu.edu/technometrics/80-89/VOL-26-02/v2602131.pdf),
[Zou et al. 2009](https://publish.uwo.ca/~gzou2/10MOVERs.pdf).
The 1984 one-way unweighted/shortest-interval example checks the MLS primitive;
it does not establish the validity of our general unbalanced extension.

[R estimates](variance_benchmarks_0.9.2.json): VCA 1.5.2 / lme4 2.0.6, eight
public/synthetic fixtures. Maximum absolute variance difference 9.2832e-7
(VCA), 6.1770e-7 (lme4); existing successful fits are numerically unchanged.
[New interval benchmarks](intervals_benchmarks_0.9.2.json) regenerate independent
R matrix/quantile calculations and tolerance 3.0.0 outputs. Mean-square total
interval maximum relative error is below 1e-13 on these eight fixtures;
raw near-zero correlations can have relative error 40.70 (absolute below 1e-14).
Normal-bound/tolerance-factor maximum relative error is below 5e-13; order
statistics and ranks agree exactly, including public rADA data and fixtures.
This is matched arithmetic evidence, not blanket software equivalence.

[Worked examples and hashes](../fixtures/intervals_092/source_manifest.json):
Burdick–Graybill bottle 95% shortest MLS interval [0.0019308170, 0.0215855852],
versus printed [0.00193, 0.02159], max absolute difference 4.41481e-6.
Young 2016 milk lower tolerance limit 0.9610332832 versus 0.9610333,
absolute difference 1.67503e-8. This checks the normal-percentile primitive,
not an ADA clinical example. The tolerance 3.0.0 help's logistic example
(seed 100, n=200) is regenerated and reproduced in Python exactly; its executed
`side=1` call is used despite a two-sided label in the help prose.

## Precision calibration: every design and component upper bound

[Detailed record](variance_calibration_0.9.2.json): 1000 datasets per row,
seed 20260929; fitted-model parametric bootstrap seed 923, 1000 per fixture.
95% coverage bound 93.6216%; failures count against all-generated coverage.
Values below are coverage % ± MCSE in percentage points. Total is two-sided;
component columns are **one-sided upper** coverage. Dash means absent component.

| Design | Total | Lot upper | Run upper | Analyst upper | Residual upper |
|---|---:|---:|---:|---:|---:|
| small_nested | 95.1 ± 0.6826 | 95.7 ± 0.6415 | 94.7 ± 0.7085 | — | 96.2 ± 0.6046 |
| small_unbalanced_nested | 95.3 ± 0.6693 | 95.1 ± 0.6826 | 95.7 ± 0.6415 | — | 96.3 ± 0.5969 |
| unbalanced_nested_crossed_fixed | 95.9 ± 0.6270 | 92.8 ± 0.8174 **MISS** | 95.4 ± 0.6624 | 95.2 ± 0.6760 | 96.0 ± 0.6197 |
| zero_lot_boundary | 95.1 ± 0.6826 | 100.0 ± 0.0000 | 94.7 ± 0.7085 | — | 96.2 ± 0.6046 |
| six_lots | 95.5 ± 0.6556 | 94.9 ± 0.6957 | 95.1 ± 0.6826 | — | 94.4 ± 0.7271 |
| bootstrap_unbalanced | 95.3 ± 0.6693 | 95.5 ± 0.6556 | 94.5 ± 0.7209 | — | 95.1 ± 0.6826 |
| bootstrap_mixed_fixed | 95.7 ± 0.6415 | 95.6 ± 0.6486 | 94.9 ± 0.6957 | 95.2 ± 0.6760 | 95.0 ± 0.6892 |

Only nested-crossed lot upper coverage misses; disclose 92.8% ± 0.8174 points.
It is an approximation limitation, not corrected by hiding the result. Total
intervals may be very wide: balanced three-lot mean variance width 169.8584
(MCSE 5.3758) for true total 6. Zero-lot coverage 100% has plug-in MCSE zero;
it does not show that a zero estimate proves absence of variation.

Reviewer seed 11 check reproduced **84.6% ± 1.1414 points** classical
Satterthwaite and **95.6% ± 0.6486** MLS. Historical
[0.9.1 calibration](variance_calibration_0.9.1.json) is unchanged. New evidence
supersedes it for MLS/MOVER only; Satterthwaite misses still apply to that option.
We retain its default to avoid silently changing previous configurations;
Skill guidance recommends explicitly choosing MLS with justified applicability.

## ADA: confidence direction and calibration

[Hoffman–Berger 2011](https://doi.org/10.1016/j.jim.2011.08.019) distinguishes
average and confidence-level cut points; [Shen et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25356783/)
explicitly uses a lower percentile limit to assure FPR **at least** target.
This corrects the opposite-goal framing in 0.9.1 without editing its history.
`ada_bounds.py` uses crossed ANOVA total variance, grand-mean variance, effective
n/df and a noncentral-t lower percentile limit (mixed-panel approximation).
Nonparametric mode uses the largest admissible binomial order statistic on
prespecified pairs with distinct subjects AND runs. No pooled-cell independence.
Both target a marginal new subject/random run, not each realized run. This
separate estimand retains run diagnostics but does not apply point-mode dynamic
selection; floating still requires NC tracking. No IQR handling or exclusions.

[ADA calibration](ada_bounds_calibration_0.9.2.json): seed 20260930, 1000 panels
per row, 80 subjects × 6 runs, subject variance .09, residual .01, run variance
0/.04. Exact conditional normal FPR; no noisy future Bernoulli draw. Every
panel is reportable. Target confidence .90; acceptance bound 88.1026%.
Mean FPR and confidence-attainment MCSE are shown in percentage points.

| Tier | Method | Run variance | Confidence attainment % ± MCSE | Mean FPR % ± MCSE | Result |
|---|---|---:|---:|---:|---|
| screening | parametric | 0.0 | 90.8 ± 0.9140 | 7.62868 ± 0.06735 | pass |
| screening | parametric | 0.04 | 91.9 ± 0.8628 | 10.34778 ± 0.14465 | pass |
| screening | nonparametric | 0.0 | 96.6 ± 0.5731 | 28.56191 ± 0.50049 | pass |
| screening | nonparametric | 0.04 | 97.3 ± 0.5126 | 28.51938 ± 0.50801 | pass |
| confirmatory | parametric | 0.0 | 91.2 ± 0.8959 | 2.03840 ± 0.02913 | pass |
| confirmatory | parametric | 0.04 | 92.5 ± 0.8329 | 3.14539 ± 0.06495 | pass |
| confirmatory | nonparametric | 0.0 | 93.3 ± 0.7906 | 13.84852 ± 0.39129 | pass |
| confirmatory | nonparametric | 0.04 | 93.8 ± 0.7626 | 13.97255 ± 0.39868 | pass |
| titer | parametric | 0.0 | 91.5 ± 0.8819 | 0.33055 ± 0.00766 | pass |
| titer | parametric | 0.04 | 93.1 ± 0.8015 | 0.61627 ± 0.01961 | pass |
| titer | nonparametric | 0.0 | 99.5 ± 0.2230 | 13.84852 ± 0.39129 | pass |
| titer | nonparametric | 0.04 | 99.2 ± 0.2817 | 13.97255 ± 0.39868 | pass |

Targets are screening 5%, confirmatory 1%, titer 0.1%, not silently inferred
from tier. Six-pair nonparametric titer mean FPR about 14% is conservative but
potentially operationally unsuitable. Passing confidence attainment is not
closeness to the target. Point-mode numerical outputs and its titer bias remain.
No new runtime R: Python implements both interval families. R packages are
validation-only in ignored `.r-lib/`; only existing MMRM KR uses `r_bridge`.

## Defects, live agents and release checks

- Six-lot simulation exposed 16 failed REML score checks. Centered retry only
  on failed fits fixed them, preserving successful legacy numerics. A bootstrap
  trial covariance also failed Cholesky; rejecting invalid line-search trials
  fixed it. Regression tests and unchanged-seed reruns have zero fit failures.
  [Earlier attempts](attempt_logs_0.9.2/) retain the failures. Bounds never changed.
- Initial live-agent attempt exposed a harness timing defect: the release
  version changed while agents used a copied collection. Version is now frozen;
  fresh attempts rerun the scenarios. Guidance also clarifies approximate MLS
  coverage, matching calibration truth, conditional FPR wording, admissible ranks
  and shared-run cell correlation. Earlier failed attempts remain inspectable.
  [Live-agent record](agent-scenarios/0.9.2/RESULTS.md): V3/V4/A6 and final A7 pass;
  no pending executions. Developing-assistant scoring; failed attempts are retained.
- [Release checks](release_checks_0.9.2.json): 422 tests pass (8 warnings), three
  Skills validate, 65 isolated install cases pass, 55 facts match source exactly.
  R-visible/hidden checks pass; 51 wheel source files and 34 preserved build files match.
- [Legacy check](legacy_fixture_checks_0.9.2.json): 50 configurations, 492
  byte-identical artifacts. All scientific artifacts for 34 unchanged-module
  configs match; old ADA/precision numerics exactly equal. New inert defaults,
  dependent hashes and scoped point-mode limitation wording are the only differences.
- No edits under runs/, build/ or to historical validation. Build uses a temporary
  source copy. Same-host macOS checks do not establish cross-platform support.
