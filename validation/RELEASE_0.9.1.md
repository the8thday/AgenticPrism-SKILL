# Release evidence — local development 0.9.1

Branch `dev/0.9.1` starts from `dev/0.9.0` at `0e87d50`.
Scope stops at 0.9.1; no push, publication, main merge or 0.10 work.

## Scope and five-gate status

| Requested item | Status and boundary |
|---|---|
| 1. ADA recalibration | Done: exact conditional FPR, reviewer reproduction, corrected Skill/contract/facts; existing ADA numerics unchanged. |
| 2. General variance components | Partial overall (live-agent gate pending); numerical scope implemented: Python REML independent random intercepts, explicit nested/crossed terms, unbalanced data, optional fixed categorical/numeric adjustments, SD/CV and Satterthwaite intervals. Calibration misses below; live-agent gate pending approval. |
| 3. Cut-point confidence bounds | Not started. A lower threshold raises FPR; the requested guarantee that FPR stays below target requires an upper tolerance-type bound. Neither direction is implemented. |
| 4. Bioanalytical method validation | Not started. Precision components alone do not establish accuracy, total error, quantitation limits or SOP acceptance. |
| 5. ADA sensitivity / drug tolerance | Not started. |
| 6. Dynamic cut point / older facts | Not started. New precision results and corrected ADA evidence have facts. |

## ADA: why the historical record is superseded

[0.9.0 calibration](ada_calibration_0.9.0.json) is preserved unchanged.
One future negative per panel measured FPR with avoidable Bernoulli noise;
the floating/IQR screening failure labels were not reliable method diagnoses.
The new script evaluates `norm.sf((log(threshold)-(4+shift))/sqrt(.1))`,
or the corresponding normal inhibition tail, on each reportable panel.
Withheld panels are excluded explicitly, never counted as successes.
Paired latent draws allow fixed/floating comparison with exactly known NC drift.
Criterion frozen before execution: mean conditional FPR <= 1.2 target + 2 MCSE.
The 20% relative tolerance permits modest point-percentile bias; it is not a
confidence guarantee for individual panels. Share above twice target is descriptive.

Reviewer seed 7, 400 panels: fixed/floating **5.4645477%**, MCSE **0.0947606
percentage points**; IQR **5.4743027%**, MCSE **0.0947718 points**.
All have 381 evaluable panels. Finding A is reproduced to its printed precision.
The mean remains above 5%; recalibration does not remove percentile bias.

[Detailed ADA record](ada_calibration_0.9.1.json): seed 20260927, 1000 panels
per row, 80 subjects × 6 runs × 2 symmetric technical wells. Normal log effects,
subject variance .09 and cell variance .01; noiseless NC. Values below are %;
MCSE is in percentage points. Every row has 953 evaluable / 47 withheld panels.

| Scenario | Mean FPR ± MCSE | Target | Bound | Share FPR>2 target ± MCSE | Result |
|---|---:|---:|---:|---:|---|
| fixed_screening_log | 5.40183 ± 0.05828 | 5 | 6.11656 | 1.3641 ± 0.3757 | pass |
| floating_screening_log | 5.40183 ± 0.05828 | 5 | 6.11656 | 1.3641 ± 0.3757 | pass |
| fixed_screening_nonparametric | 5.58355 ± 0.06926 | 5 | 6.13851 | 3.4627 ± 0.5923 | pass |
| floating_screening_nonparametric | 5.58355 ± 0.06926 | 5 | 6.13851 | 3.4627 ± 0.5923 | pass |
| fixed_confirmatory | 1.19834 ± 0.02067 | 1 | 1.24134 | 11.2277 ± 1.0227 | pass |
| fixed_titer_log | 0.14886 ± 0.00413 | 0.1 | 0.12825 | 23.9244 ± 1.3820 | **MISS** |
| fixed_screening_boxcox_lambda_zero | 5.40183 ± 0.05828 | 5 | 6.11656 | 1.3641 ± 0.3757 | pass |
| fixed_screening_biological_iqr_exclude | 5.40290 ± 0.05827 | 5 | 6.11655 | 1.3641 ± 0.3757 | pass |

Separate future-shift demonstration (+0.6 log units): FPR **60.4496% ± 0.2044 points**; all reportable panels exceed twice target (share MCSE 0). No method pass/miss label.

## General variance components: method and R agreement

`variance_components.fit_components` extends the shared primitive; the old
`crossed_anova` function remains byte-identical. `precision.py` validates DEFAULTS,
literal-true gates and design rationale. Nested terms use composite factor IDs;
fixed effects adjust means. No random slopes or correlated residuals; dense
solver limited to 1000 observations. No silent missing-row deletion.
CV uses a declared positive reference treated as fixed, with rationale.

Python was retained: expected REML information `trace(P K_j P K_k)/2`, its
inverse on active components (VCA Giesbrecht–Burns convention), and df
`2*v²/Var(v)` are directly implementable. Sums include component covariances.
Existing `mixed_inference._hessian` supplies final optimizer score refinement.
Boundary component intervals are unavailable, not [0,0]; total intervals
condition on the active set. Null upper endpoints within an interval are unbounded.
Intervals are marginal approximations, not simultaneous confidence statements.

[Regenerated benchmark](variance_benchmarks_0.9.1.json): eight fixtures, VCA 1.5.2
and lme4 2.0.6, R 4.6.0. All meet frozen tolerance 5e-4 × max(1,abs(reference)).
Component variance maximum absolute difference 9.2832e-7; vs tight lme4 6.1770e-7.
Maximum relative differences for references >=1e-6: variance 1.0128e-4,
df 2.0256e-4, lower CI 2.3054e-7, upper CI 2.4836e-7,
component covariance 4.5461e-7, fixed coefficients 4.4973e-9.
Raw maximum relative errors are also retained: variance 1 at near-zero boundaries,
and covariance 6.5726 for numerical near-zero entries. No blanket relative-
precision or Prism equivalence claim is made. VCA and lme4 share a fitting engine;
Python is the independent numerical implementation.

Published worked example: [VCA vignette](https://cran.r-project.org/web/packages/VCA/vignettes/VCA_package_vignette.html),
`dataEP05A2_3`: total/day/day:run/error variances 35.546313/12.231099/7.031193/
16.284021 reproduced with maximum absolute difference **4.84969e-7**, within
5.1e-7 printed precision. The published example is balanced ANOVA; interior REML
has the same point estimates. Its printed intervals are not asserted equivalent.
EP05-A2 datasets 1–3 and Penicillin are included. Source/archive/data hashes,
public attribution and exact export/generator scripts are in
[the source manifest](../fixtures/variance_reml/source_manifest.json).
R VCA is validation-only in ignored `.r-lib/`; no runtime R dependency was added.
Existing optional MMRM KR remains the only R runtime method, for its adjusted
covariance/df machinery. The bridge and installer are unchanged.

## Variance interval calibration: every row

[Detailed record](variance_calibration_0.9.1.json), seed 20260929, 1000 datasets
per row. Three lots × four runs × three replicates, or a 2–4 replicate unbalanced
variant; crossed scenario adds four analyst levels and a fixed covariate.
Bound: .95 − 2 sqrt(.95×.05/evaluable). Coverage and MCSE below are %/percentage
points; unavailable intervals are excluded. All-generated proportions and
estimate means with MCSE are retained in JSON. Zero fit failures after correction.

| Design | Component | Evaluable | Coverage ± MCSE | Lower bound | Result |
|---|---|---:|---:|---:|---|
| small_nested | lot | 926 | 96.3283 ± 0.6180 | 93.5676 | pass |
| small_nested | lot:run | 984 | 95.3252 ± 0.6730 | 93.6104 | pass |
| small_nested | residual | 1000 | 95.9000 ± 0.6270 | 93.6216 | pass |
| small_nested | total | 1000 | 85.1000 ± 1.1261 | 93.6216 | **MISS** |
| small_unbalanced_nested | lot | 928 | 96.6595 ± 0.5899 | 93.5691 | pass |
| small_unbalanced_nested | lot:run | 972 | 94.9588 ± 0.7018 | 93.6019 | pass |
| small_unbalanced_nested | residual | 1000 | 94.6000 ± 0.7147 | 93.6216 | pass |
| small_unbalanced_nested | total | 1000 | 84.9000 ± 1.1322 | 93.6216 | **MISS** |
| unbalanced_nested_crossed_fixed | lot | 827 | 93.5913 ± 0.8516 | 93.4843 | pass |
| unbalanced_nested_crossed_fixed | lot:run | 959 | 95.0991 ± 0.6971 | 93.5924 | pass |
| unbalanced_nested_crossed_fixed | analyst | 892 | 94.9552 ± 0.7328 | 93.5405 | pass |
| unbalanced_nested_crossed_fixed | residual | 1000 | 94.7000 ± 0.7085 | 93.6216 | pass |
| unbalanced_nested_crossed_fixed | total | 1000 | 83.4000 ± 1.1766 | 93.6216 | **MISS** |
| zero_lot_boundary | lot | 387 | 0.0000 ± 0.0000 | 92.7842 | **MISS** |
| zero_lot_boundary | lot:run | 989 | 97.6744 ± 0.4792 | 93.6140 | pass |
| zero_lot_boundary | residual | 1000 | 96.0000 ± 0.6197 | 93.6216 | pass |
| zero_lot_boundary | total | 1000 | 96.6000 ± 0.5731 | 93.6216 | pass |

Total coverage misses and true-zero lot failure remain unresolved statistical
limitations of this VCA-matched approximation, disclosed in Skill and facts.
We did not tune bounds, scenarios or seeds or claim nominal small-sample coverage.
Boundary censoring of individual intervals must accompany conditional coverage.

## Defects, guidance and checks

- Initial printed-total comparison missed by 5.78154e-7: relative-objective
  optimizer stopping was too early. Score refinement fixed it without widening
  tolerance. Three initial simulation fits failed the constrained score check;
  a constrained optimizer fallback fixed all three on the unchanged scenarios.
  Initial outputs and code are retained in [attempt logs](attempt_logs_0.9.1/).
- A tiny-df chi-square underflow produced an unrepresentable lower endpoint,
  incorrectly treated as zero in calibration. Such intervals now have explicit
  unavailable status; the same-seed rerun gives zero-lot coverage 0/387, not 1/388.
  New regression test covers this public EP05-A2 boundary case. Zero observed
  binomial MCSE is not a confidence guarantee. Pre-fix records are retained.
- Validation-only VCA source compilation initially failed on a macOS deployment
  target; a valid local target allowed installation. No global configuration changed.
- New Skill stops on undeclared nesting, CV denominator or acceptance limits,
  separates zero estimates from absence of variation and discloses coverage misses.
- Live Claude scenarios V1/V2/A5 are prepared. Automatic approval review rejected
  the external Claude call over repository-data export; explicit approval requested.
  See [scenario status](agent-scenarios/0.9.1/RESULTS.md). Not claimed executed.
- Full pytest: **385 passed**, four existing pandas warnings, 185.55 seconds;
  36 new-module tests. All three changed Skills pass quick_validate.
- Legacy check against dev/0.9.0: **42 configurations, 484 byte-identical
  artifacts**. Eight ADA configurations change only calibration evidence,
  associated limitations and dependent hashes; all other JSON content is equal.
- **341 historical validation files and 687 tracked run files** unchanged.
  Both pre-existing untracked run directories remain untouched. Build's 34 files
  are unchanged; `uv build` ran in a temporary source copy.
- Wheel/source: **48 packaged source/resource files** identical. Isolated clean
  environment passed **39 exact facts configurations**, including all eight new
  fixtures. R-hidden/R-present precision and ADA results are byte-identical;
  old KR refuses without R and runs with R. No fresh runtime R install.
- Initial smoke check stalled rebuilding fonts with an unwritable home cache;
  it was interrupted. A shared temporary writable Matplotlib cache was added
  to the check script. Corrected-wheel smoke: **49 configurations passed**,
  including analyze/render/verify, immutable facts/results and plate import.
- [Release checks](release_checks_0.9.1.json) contain commands, log hashes, wheel
  hash and scope. This is same-host macOS validation, not cross-platform evidence.
