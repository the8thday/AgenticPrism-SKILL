# Release evidence — local development 0.8.0

0.8.0 adds three specialists and extends one. The collection now has nine
Skills and one package. No equilibrium, kinetic, dose-response, ELISA or
independent-group numerical model was changed. Runtime dependency added:
statsmodels 0.15.0 (pinned). R is used only as a validation oracle and is not
needed to run any Skill. Work remains local.

The release was built in two stages. The first commit (61e29a5, Codex) added
one-factor repeated measures with Wald and parametric-bootstrap inference; its
record is kept unchanged in the [appendix](#appendix--first-080-commit-record-historical).
A review followed, and the later commits (872fa08 to the release commit) added
Satterthwaite inference, arm × condition designs, time-to-event and tumor
growth. Statements in the appendix that no Satterthwaite df, treatment-by-time
effects or random slopes are provided are superseded by the sections below.
Its simulation (100 experiments, coarse smoke bounds) is superseded by the
1,000-experiment calibration below.

## Implemented scope

| Skill | Added in 0.8.0 |
|---|---|
| repeated-measures | One within-unit factor (RM ANOVA with GG; random-intercept REML), and between-unit arm × within-unit condition: split-plot ANOVA (type III via the multivariate model on orthonormal within contrasts, GG epsilon from the pooled within-arm covariance) and a two-way random-intercept model with type III F tests. Inference: `satterthwaite` (t/F with Satterthwaite df and the Fai–Cornelius multi-df combination, as lmerTest), `parametric_bootstrap` (one-factor only) or `wald_asymptotic`. Families: arm vs control at each condition, or condition vs first within each arm, with Holm p and Bonferroni simultaneous limits. Satterthwaite is withheld when the random-intercept variance is on its zero boundary. |
| time-to-event | Kaplan–Meier with Greenwood SE and log or log-log limits; medians with limits (R midpoint rule); landmark survival; log-rank on k−1 df (generalized inverse with rank df when some arm shares no risk set with an event, as R `survdiff`) or Monte Carlo permutation; pairwise vs-control Holm families; Cox (Efron ties) hazard ratios with Bonferroni limits and optional prespecified baseline covariates; `cox.zph` (km transform) PH score tests. Event coding and censoring reasons are declared before testing. |
| tumor-growth | log(V + offset) random-slope mixed model (animal random intercept and slope, unstructured 2×2 covariance, REML by Woodbury sufficient statistics with Newton refinement); per-arm growth rates and doubling times; equal-rates F; rate differences vs control (Holm p, Bonferroni limits, Satterthwaite df); model geometric-mean T/C at the analysis day. Secondary observed TGI% and T/C% at the analysis day with Fieller limits (Welch df, Bonferroni family), flagged as survivor-biased whenever an animal was removed earlier. |

Shared implementation: `mixed_inference.py` holds REML deviance functions for
random-intercept and random-slope models and Satterthwaite df on the
unconstrained variance parameters. Derivatives use Richardson extrapolation
(numDeriv defaults, as lmerTest). A fixed-step version was off by 9e-4 in
relative df from lmerTest on the tumor fixture, where a Cholesky off-diagonal is
tiny; Richardson brought it to 1.7e-7.

## Numerical checks against R

R 4.6.0 with nlme 3.1.169, lme4 2.0.6, lmerTest 3.2.1, afex 1.5.1, emmeans
2.0.4, pbkrtest 0.5.5 (Kenward–Roger df recorded for reference only) and
survival 3.8.6. Public datasets are exported by R at validation time. Values are
maximum relative differences unless marked absolute.

| Check | Data | Observed agreement |
|---|---|---|
| RM ANOVA F, GG epsilon | Orthodont (complete) | 7.8e-16, 2.2e-16 |
| Random-intercept REML vs nlme: β, SE, τ, σ, REML log-lik | Orthodont complete / incomplete; synthetic fixtures | β ≤ 8.8e-11; SE ≤ 1.4e-8; τ ≤ 7.1e-8; σ ≤ 1.6e-7; log-lik ≤ 5.4e-13 (absolute) |
| One-factor Satterthwaite contrasts vs lmerTest: df, t | Orthodont complete / incomplete; synthetic complete / incomplete | df ≤ 2.3e-8; t ≤ 1.2e-8 |
| Split-plot ANOVA vs afex (GG) | BodyWeight, Orthodont | ≤ 1.4e-13 |
| Two-way mixed type III F and ddf vs lmerTest | BodyWeight, ChickWeight, Orthodont, fixture | ≤ 9.1e-7 (ChickWeight); others ≤ 4.6e-7 |
| Two-way contrasts vs emmeans/lmerTest: estimate (absolute), df, t | same | estimate ≤ 2.8e-8; df ≤ 1.4e-7; t ≤ 3.8e-8 |
| KM survival and limits (log, log-log), medians with limits (absolute) | lung, veteran | ≤ 7.8e-16; medians identical |
| Log-rank χ² vs `survdiff` | lung, veteran | ≤ 8.7e-15 |
| Reduced-rank log-rank (one arm censored before any event) | 3-arm unit-test case | χ² 5.051660516605 on 1 df, p 0.024602349954, equal to `survdiff` within 1e-10 |
| Cox β, SE, LR/Wald/score, log-lik vs `coxph` (Efron) | lung, veteran | ≤ 1.6e-14 |
| `cox.zph` (km) per-term and global χ² | lung, veteran | ≤ 3.7e-15 |
| Random-slope growth rates, SE, df; equal-rates F and ddf; rate-difference df vs lmerTest | BodyWeight, ChickWeight, tumor fixture | rates ≤ 4.8e-10; SE ≤ 5.8e-8; df ≤ 1.9e-7; F ≤ 1.1e-7 |

Records: [repeated_benchmarks.json](repeated_benchmarks.json),
[survival_benchmarks.json](survival_benchmarks.json),
[tumor_growth_benchmarks.json](tumor_growth_benchmarks.json). Reproduce with
`scripts/validate_repeated_benchmarks.py`, `validate_survival_benchmarks.py`
and `validate_tumor_growth_benchmarks.py` on a machine with R and these
packages. Agreement with R covers the stated models and data; it is not a claim
of equivalence with Prism, SAS or any other lmerTest/survival option.

## Calibration simulations

All bounds were fixed in the scripts before running: rejection ≤ 0.05 + 2
binomial SE, coverage ≥ 0.95 − 2 SE, over evaluable simulations. A miss is
reported as a miss; nothing was retuned afterwards. 1,000 simulations per row.

### Repeated measures ([repeated_calibration.json](repeated_calibration.json), seed 20260929)

One-factor: Gaussian random intercept (subject SD 2 or 0.3, residual SD 1), four
conditions, 15% MCAR omissions; family = conditions vs the first. Bootstrap used
199 refits per experiment (the Skill default is 999).

| Scenario | Method | Reportable | Omnibus / effect rejection | Family rejection | Simultaneous coverage | Bound | Outcome |
|---|---|---|---|---|---|---|---|
| 8 units, subject SD 2 | Satterthwaite | 999 | 4.7% | 4.0% | 96.0% | ≤ 6.38% / ≥ 93.62% | passed |
| | bootstrap | 1000 | 3.9% | 3.7% | 96.0% | same | passed |
| | Wald (comparator) | 1000 | 9.1% | 6.1% | 93.9% | none | liberal |
| 12 units, subject SD 2 | Satterthwaite | 1000 | 5.0% | 5.1% | 94.9% | ≤ 6.38% / ≥ 93.62% | passed |
| | bootstrap | 1000 | 4.6% | 6.0% | **93.4%** | same | **missed coverage** |
| | Wald (comparator) | 1000 | 8.4% | 6.7% | 93.3% | none | liberal |
| 12 units, subject SD 0.3 (weak clustering) | Satterthwaite | 673 (327 withheld at the variance boundary) | 5.5% | 4.6% | 95.4% | ≤ 6.68% / ≥ 93.32% | passed |
| | bootstrap | 1000 | 4.1% | 4.3% | 95.4% | ≤ 6.38% / ≥ 93.62% | passed |
| | Wald (comparator) | 1000 | 6.4% | 5.6% | 94.4% | none | — |
| Split-plot, 3 arms × 5 times × 8, residual SD 0.5→3.5 (nonspherical), complete | GG | 1000 | arm 5.2%, time 4.7%, arm × time 4.7% | 3.0% | 97.0% | ≤ 6.38% / ≥ 93.62% | passed |
| Two-way random intercept, 3 × 5 × 8, 15% MCAR | Satterthwaite type III | 1000 | arm 4.5%, time 5.2%, arm × time 5.0% | 2.8% | 97.2% | same | passed |

One miss: the bootstrap's simultaneous coverage at 12 units was 93.4%, below the
93.6% bound, although its rejection rates were within bounds. With 199 refits the
max-|t| quantile is estimated from few draws; whether 999 refits (the Skill
default) fixes this was not simulated. The Skill therefore recommends
Satterthwaite and states the miss. In the weak-clustering scenario Satterthwaite
was withheld in 33% of datasets because the random-intercept variance was
estimated at zero; the bootstrap still gave results there and held its bounds.
The two-way families are conservative (2.8–3.0%), as expected for Bonferroni
over ten correlated contrasts.

This run used the code after the REML polish and before the switch to
Richardson derivatives. For random-intercept models the two derivative methods
agree to about 1e-7 in relative df, far below the Monte Carlo error here.

### Time to event ([survival_calibration.json](survival_calibration.json), seed 20261001)

Exponential event times with independent uniform censoring.

| Check | Setting | Result | Bound | Outcome |
|---|---|---|---|---|
| Log-rank type I error, permutation (999 draws) | 2 arms × 6 / 8 / 10 | 4.9% / 4.3% / 4.5% | ≤ 6.38% | passed |
| Log-rank type I error, asymptotic (comparator) | same | 6.6% / 5.6% / 5.1% | none | liberal at n = 6 |
| KM log-log limits cover S = 0.5 at the true median | n = 10 / 20 | 96.4% / 96.3% | ≥ 93.62% | passed |
| KM log limits (comparator) | n = 10 / 20 | 91.9% / 93.3% | none | undercover |
| Cox Wald HR limits, true HR 0.5 | 2 arms × 20 (× 10 reported) | 95.1% (96.9%) | ≥ 93.62% | passed |
| `cox.zph` global type I error under PH | 2 arms × 20 | 6.2% | ≤ 6.38% | passed |

These results set two Skill defaults: permutation log-rank for small arms and
log-log limits.

### Tumor growth ([tumor_growth_calibration.json](tumor_growth_calibration.json), seed 20261003)

log V = log 100 + u0 + (rate + u1)·day + e with u0 ~ N(0, 0.15²), u1 ~ N(0, 0.015²),
e ~ N(0, 0.12²); 3 arms × 10 animals; days 0–28; an animal is removed after its
first volume above a threshold (MAR). Analysis day 21.

| Scenario | Rates per day | Animals removed before day 21 (mean) | Result | Bound | Outcome |
|---|---|---|---|---|---|
| A null | 0.08 / 0.08 / 0.08, threshold 2000 mm³ | 0.0 | equal-rates type I 5.8%; rate-difference family coverage 94.2% | ≤ 6.38%; ≥ 93.62% | passed |
| B early dropout | 0.12 / 0.08 / 0.04, threshold 1000 mm³ | 2.1 | rate-difference family coverage 95.1%; observed-TGI Fieller family coverage 92.2% (comparator) | ≥ 93.62% for the model | passed |
| C no dropout | same rates, no threshold | 0 | rate-difference coverage 94.5%; observed-TGI coverage 95.4% | ≥ 93.62% each | passed |

In scenario A, tumors reached about 540 mm³ by day 21, so no animal was removed;
it checks the null without dropout. Dropout under the null was not simulated.
The TGI truth (61.5% and 88.1%) is E[ΔV] by Monte Carlo on the generating model.
Scenario B shows why the Skill leads with the growth model: under MAR removal
the model kept nominal coverage while survivor-based TGI did not.

## Live-agent sessions

[agent-scenarios/0.8.0/RESULTS.md](agent-scenarios/0.8.0/RESULTS.md):
three new scenarios (R1 arm × day body weight with missing cells, T1 humane
euthanasias coded as censored, G1 day-28 TGI after vehicle dropout) and the nine
0.7.0 misuse scenarios rerun against the 0.8.0 router. All twelve passed, one
session each, scored by the developing assistant. In R1 and T1 the agent stopped
before modelling to ask for the missingness reason or event definition. In G1 it
led with the growth model and reported the survivor bias of observed TGI.

## Tests, package and installation

- Full suite: 226 passed, 0 failed, 0 skipped
  ([release_checks_0.8.0.json](release_checks_0.8.0.json); the JUnit XML is kept
  locally). Tests for the new code: `test_repeated.py` (25, including
  Satterthwaite), `test_repeated_two_way.py` (12), `test_time_to_event.py` (15,
  including the R-matched reduced-rank log-rank) and `test_tumor_growth.py` (11).
- All nine Skill folders pass skill-creator `quick_validate.py`.
- Wheel sha256 1d83f731…: installed into a separate Python 3.13 environment with
  the locked dependencies. For all seven 0.8.0 fixtures (RM ANOVA, one-factor
  mixed with bootstrap, split-plot, two-way mixed, survival permutation, adjusted
  Cox, tumor growth), results from the wheel equal results from the source
  exactly. Earlier modules reproduce as before
  ([clean_environment_0.8.0.json](clean_environment_0.8.0.json)).
- Install smoke from an unrelated working directory: analyze, verify, render
  and verify again for 17 fixtures, with render leaving `results.json` unchanged
  ([install_smoke_0.8.0.json](install_smoke_0.8.0.json)). Repeated-measures
  installed-wheel check including both two-way designs
  ([installed_repeated_0.8.0.json](installed_repeated_0.8.0.json)).
- macOS arm64 only; Linux and Windows were not tested.

## Remaining boundaries

- Only synthetic fixtures and public R datasets were used; no team-measured
  study has been validated.
- Repeated measures: compound symmetry only (a random intercept does not remedy
  nonsphericity; use split-plot GG for complete data); no MMRM/unstructured or
  AR(1) covariance, Kenward–Roger df, nested or crossed effects, covariates or
  non-Gaussian outcomes. The two-way parametric bootstrap is not implemented.
- Tumor growth: exponential growth on the log scale only (no Gompertz/logistic),
  no cage effects or heteroscedastic measurement error, and MAR removal only.
  Informative (MNAR) removal is not modelled; the Skill routes time to endpoint
  to time-to-event instead.
- Time to event: no competing risks, recurrent events, time-varying covariates,
  frailty/clustered subjects or interval censoring. Calibration covered
  exponential hazards with independent censoring.
- Calibration used one generating model per module and small animal-study sizes;
  it does not establish nominal behavior under other designs.
- Agreement with R is not Prism equivalence; no Prism project benchmark was run.

---

## Appendix — first 0.8.0 commit record (historical)

Kept as written for commit 61e29a5. Superseded where the sections above say so.
Its test count (185), installed-wheel check (wheel sha256 9adcf04e…, moved to
`dist/superseded/`) and Skill count (seven) describe that commit only.


Added a separate repeated-measures specialist and shared numerical/reporting
implementation. This release supports one within-unit factor and a Gaussian
random intercept, with explicit model, missingness and contrast declarations.
No existing kinetic, equilibrium, dose-response, ELISA or independent-group
numerical model was changed. Runtime dependency: statsmodels 0.15.0, pinned with
its dependencies. The collection contains seven Skills. Work remains local.

### Implemented scope

- Complete one-factor repeated-measures ANOVA, always using Greenhouse–Geisser
  corrected df/p; raw F and uncorrected p retained. Paired comparison families
  have Holm p and Bonferroni simultaneous t intervals.
- Condition cell means plus a unit random intercept, fitted by REML using
  statsmodels MixedLM. The fixed-effect covariance is explicit plug-in GLS.
  Missing observations are retained as an audited pattern, not imputed or
  silently removed as complete cases. Available-case use requires a stated
  conditional MAR rationale.
- Explicit `wald_asymptotic` (z/chi-square, no small-sample df correction), or
  centered parametric bootstrap Wald/max-|t| inference with refitted REML
  variances on every simulated subject vector. This is not an ML
  likelihood-ratio bootstrap, KR or Satterthwaite method. Bootstrap inference is
  withheld if fewer than 199 or 95% of requested refits succeed.
- Saved omnibus, observed summaries, estimated condition means, contrasts,
  missing cells, marginal/conditional residuals, optimizer attempts and all
  bootstrap pivots. Original inputs, exclusions, configuration, source and
  dependency hashes are preserved. Render does not refit or resample.

### Numerical checks

| Check | Observed result | Scope |
|---|---|---|
| RM F vs statsmodels AnovaRM | Relative agreement within 1e-12 | Complete synthetic 16-unit, four-condition design |
| GG epsilon | Matches a separate orthonormal Helmert/eigenvalue computation within 1e-12 | Same design; F-tail calculation checked independently |
| Paired contrasts | Match scipy paired t, hand Holm step-down and Bonferroni t endpoints | Prespecified control family |
| Balanced REML | Matches closed-form ANOVA variance components and contrast SE | Positive random-intercept variance; numerical tolerance 2e-4 |
| Incomplete-data GLS | Matches a separately constructed dense V matrix solution within 1e-10 | Explicit observed design, no imputation |
| Independent R nlme | Maximum absolute condition-mean difference 1.45e-8; maximum relative plug-in SE difference 2.11e-6; conditional fitted values differ by at most 1.43e-6 | Complete and incomplete synthetic fixtures; R 4.6.0, nlme 3.1.169; variance components also within 3e-4 relative tolerance. **No R p-value/df or bootstrap equivalence claim** |

Reproduce the R benchmark with
`.venv/bin/python scripts/validate_repeated_benchmarks.py` on a machine with
R/nlme installed. R is used only for validation, not for running the Skill.
Full values and source hashes: [repeated_benchmarks.json](repeated_benchmarks.json).

### Prespecified simulation checks

Seed 20260928; rejection bounds were set before running to
0.05 + 3 sqrt(0.05 × 0.95 / N), with coverage lower bound equal to one minus
that bound. Bounds were not tuned to observed results.

| Scenario | Results | Interpretation |
|---|---|---|
| 2,000 null Gaussian RM experiments, n=12, k=4, unequal condition residual SDs 0.5/1/2/3 plus subject SD 2 | Uncorrected omnibus rejection 0.0735; GG rejection 0.0460; paired-family rejection 0.0395; simultaneous coverage 0.9605 | Corrected checks passed the rejection bound 0.06462 and coverage bound 0.93538. Uncorrected omnibus is a misspecified comparator |
| 100 null random-intercept experiments, n=12, k=4, subject SD 2, residual SD 1, 15% MCAR omissions, 199 bootstrap refits per experiment | Bootstrap omnibus rejection 0.06; family rejection 0.07; simultaneous coverage 0.93; 100/100 reportable. Wald comparator: omnibus 0.08, family 0.07, coverage 0.93 | Bootstrap checks passed **coarse smoke bounds** 0.11538 / 0.88462. This is preliminary evidence, not precise confirmation of 5% error or 95% coverage. Wald has no acceptance bound |

Each mixed simulation retained only designs meeting the declared minimum
observed counts; selection used the missingness pattern, not outcomes. No
boundary original fits occurred in this simulation. Boundary handling is covered
by a separate deterministic regression test, not by a coverage study.

Reproduce with `.venv/bin/python scripts/validate_repeated_simulation.py`.
[Full simulation record](repeated_simulation.json) includes per-experiment
results. Bootstrap repeats calibrate conditional on the fitted covariance and
observed missingness pattern. More repeats do not fix covariance misspecification
or informative dropout.

### Workflow, package and report checks

- [Tests](release_checks_0.8.0.json) record the observed full-suite and focused
  rerun results. New tests cover formula agreement, pairing, missingness,
  duplicate technical replicates, disconnected designs, scale/row invariance,
  variance boundary, deterministic bootstrap, failed-refit withholding, literal
  IDs (`NA`, `001`), empty table headers, failed model serialization, provenance
  and render-only hash invariance.
- The new specialist, router and group-comparison Skills pass the skill-creator
  `quick_validate.py` checks. No independent live-agent behavioral evaluation
  was run for this release.
- [Installed wheel check](installed_repeated_0.8.0.json): a separate virtual
  environment on the same macOS host, no source import path, locked dependencies;
  RM ANOVA, mixed Wald and mixed bootstrap results/pivots compared with source;
  analyze, render and verify checked. This is not a Windows/Linux validation.
- Example outputs are generated locally using the commands below. Example PNGs for
  individual trajectories and mixed-model residual/QQ panels were visually
  inspected. Offline embedded SVG and render invariance are checked by tests;
  browser interaction QA was not performed.

From the repository root, generate the example reports with:

```sh
.venv/bin/agentic-prism analyze --config fixtures/repeated_synthetic/config_rm.json --output runs/repeated-anova-0.8.0-final
.venv/bin/agentic-prism analyze --config fixtures/repeated_synthetic/config_mixed.json --output runs/repeated-mixed-0.8.0-final
```

Each command creates `report.html` and the accompanying analysis artifacts in
its output directory. Choose a new output directory if it already exists.
Generated reports, figures, per-draw bootstrap outputs and the full pytest XML
log are kept locally rather than tracked in Git. The fixture inputs/configs,
tests, validation scripts and compact validation records remain versioned.

### Remaining boundaries

Only synthetic data were used. No team-measured antibody experiment has been
validated. The model does not support treatment-by-time effects, random slopes,
multiple/nested/crossed grouping factors, covariates, heterogeneous residual
variance, serial correlation, unstructured MMRM covariance or non-Gaussian
outcomes. A random intercept does not generally remedy nonsphericity. No missing
not at random (MNAR) sensitivity analysis, KR/Satterthwaite df, variance-component
confidence intervals or random-effect prediction intervals are provided.

Fewer than 30 units is flagged as a caution, not a proven validity cutoff.
Convergence and positive covariance do not establish scientific suitability.
At a near-zero random-intercept variance boundary, estimates are retained and
flagged; inference remains approximate. Neither the numerical benchmark nor
these simulations establish equivalence to Prism or to every R mixed-model
inference method.
