# Release evidence — local development 0.9.0

Branch `dev/0.9.0` starts from `dev/0.8.2` at
`de8e579ac78093e3f8786347d684654c3d97b6f4`. This is a deliberately scoped ADA
release, not completion of the full roadmap 0.9 section. No 0.10 work was
started. No push, publication, merge into main or public-tree build was performed.
The roadmap's Python-default and five-gate rules apply to the delivered scope.

## Requested items and delivered boundaries

| Requested item | Status | Delivered / deliberately deferred |
|---|---|---|
| ADA screening, confirmatory and titer cut points | **Partial overall; scoped core implemented** | Point thresholds, NC ratio/difference normalization, log or prespecified Box-Cox, normality/skew diagnostics, prespecified analytical flags and single-pass biological IQR handling, fixed/floating gates and subject/run variance components. Complete balanced repeated negative panels only. |
| Dynamic cut points | **Not started as a deployable method** | Heterogeneity produces a dynamic recommendation and withheld common threshold. Per-run audit percentiles are not a future-run algorithm. |
| Cut-point lower confidence bound; general analyst × day × plate REML | **Not started** | No cut-point interval is reported. Shared primitive separates two crossed factors and residual only. No fitted Box-Cox lambda, robust-MAD method, automatic analytical deletion, incomplete panels or confirmatory floating/log-ratio model. |
| Shared variance-component primitive | **Done for balanced additive two-factor crossed ANOVA** | Pure Python, raw and nonnegative component estimates, negative-component flags, conservative Gaussian component intervals. Reusable API in `src/agentic_prism/variance_components.py`; this is not the general nested/crossed REML engine requested for all future designs. |
| ADA sensitivity and drug-tolerance matrix | **Not started** | No positive-control dilution interpolation or drug-tolerance inference. |
| Bioanalytical method validation | **Not started** | Accuracy/precision, total error, precision-profile LLOQ/ULOQ, dilution/hook, parallelism, selectivity/specificity and stability remain outside scope. Existing ELISA plate QC was not relabelled as complete validation. |
| Interpretation facts | **Done for ADA; older-module retrofit not started** | New ADA facts include saved thresholds/reportability, diagnostics, source pointers, component intervals and all calibration rows/misses. Six older modules were not retrofitted. |

The main implementation uses a complete panel of at least 50 independent
subjects measured in every one of at least 3 runs, with equal numbers of at
least 2 technical replicates per cell. These minimums are software applicability
gates, not proven adequacy criteria. A single drug-naive representative negative
population and one reagent lot are required. Disease-matched drug-naive samples
can be appropriate; disease status itself is not an exclusion criterion.

Technical wells are averaged first. Confirmatory responses are
`100*(1-mean(inhibited)/mean(uninhibited))`, retaining negative inhibition.
Other tiers use raw, signal/NC or signal-NC values and the declared transformation.
Parametric thresholds use the sample mean plus `qnorm(1-FPR)` times sample SD;
nonparametric thresholds use R type-7 interpolation, followed by the inverse
transformation. The FPR and rationale must be supplied; neither 99% nor 99.9%
is silently declared a universal confirmatory requirement.

Raw and normalized run diagnostics are saved. Run means use subject-blocked
ANOVA; median Levene and pooled Shapiro p values are explicitly **descriptive**
with repeated subjects. Non-significance does not establish equality. Fixed /
floating selection is a scoped decision rule; remaining normalized mean or
variance differences withhold a common threshold. Floating additionally requires
raw drift and NC tracking. No hypothesis test proves future assay suitability.
Analytical IQR flags withhold for investigation. Biological exclusion is once
at the whole-subject level under the declared policy; no iterative trimming.

The variance model is `Y_ij = mu + subject_i + run_j + residual_ij` on the
analysis scale after cell averaging. The residual includes interaction and
cell-mean measurement noise; it is not pure technical repeatability. Negative
unconstrained ANOVA estimates remain visible; clipped values are not constrained
REML. Three simultaneous chi-square mean-square intervals (Bonferroni) propagate
to conservative component bounds. Normal independent additive effects and common
residual variance are required. The interval method is supplemental to the
published point-estimate examples, with its own calibration below.

## Why Python; R scope

All new runtime numerics are Python/NumPy/SciPy. Closed-form balanced ANOVA,
percentiles and diagnostics can be checked directly against R, so no additional
runtime R method is justified. R is only a validation oracle for this release,
using the existing ignored `.r-lib/`. The existing optional MMRM Kenward–Roger
bridge and package pins are unchanged. Python ADA also works with R hidden.
Observed oracle versions: R 4.6.0, lme4 2.0.6, car 3.1.5, jsonlite 2.0.0. R warns
that this lme4 build was compiled under R 4.6.1; the recorded comparisons ran
successfully under the observed R 4.6.0 host.

## Five gates: evidence and limits

### 1. Regenerated R comparisons

`scripts/validate_ada_benchmarks.py` runs `scripts/ada_oracle.R` afresh on eight
ADA fixture/configuration matrices and the public Penicillin matrix. Base R
`lm`/ANOVA, `qnorm`, `quantile(type=7)`, `shapiro.test`, `car::leveneTest`,
chi-square interval arithmetic and `lme4::lmer` provide the oracles. All nine
rows meet the predeclared comparison tolerances. Positive balanced ANOVA
components alone are compared to constrained REML; negative components are
retained without an equality claim.

| Quantity | Maximum relative difference | Maximum absolute difference |
|---|---:|---:|
| ANOVA mean squares | 4.6649889e-14 | 2.13162821e-13 |
| Raw variance components | 2.61033218e-13 | 8.43769499e-15 |
| Conservative component interval endpoints | 5.07200382e-14 | 7.81597009e-14 |
| Run F and p | 1.75404228e-13 | 1.13686838e-12 |
| Median Levene F and p | 8.12022787e-14 | 3.90798505e-14 |
| Shapiro W and p | 2.13686764e-12 | 1.9095836e-13 |
| Parametric percentile | 2.844994e-16 | 8.8817842e-16 |
| Nonparametric percentile | 1.36278595e-16 | 4.4408921e-16 |
| Positive components vs lme4 REML | 8.07224498e-07 | 2.69029576e-07 |

[Full record](ada_benchmarks_0.9.0.json), including actual R outputs and all
per-field differences. No Prism, SAS, full rADA workflow or general REML
software-equivalence claim is made. Bounds were not relaxed after comparison.

### 2. Published worked examples and source hashes

- [rADA 1.1.9 vignette](https://cran.r-project.org/web/packages/rADA/vignettes/rada_vignette.html):
  seven published nonparametric cut points (six day/operator combinations plus
  pooled) reproduced after the vignette's CV>=20% cell filter. Maximum absolute
  difference from printed values is **4.103217655e-6** (printed-rounding tolerance
  5.1e-6). This is a scalar percentile reproduction on public **simulated** data,
  not adoption of rADA's replicate independence or its cut-point intervals.
  The product instead aggregates wells and requires a complete panel. The
  original unfiltered public matrix is also used in the end-to-end fixture/R
  comparison, and correctly yields withheld common deployment. Original data,
  canonical CSV, acquisition URLs, transformations and archive/vignette hashes
  are in [the source manifest](../fixtures/ada_public/source_manifest.json).
- [Douglas Bates, 2011, section 4, printed page 11](https://pages.cs.wisc.edu/~st850-1/text/1SimpleH.pdf):
  Penicillin crossed plate/sample random effects. Python estimates are
  **0.716908212560386 / 3.7309178743961393 / 0.30241545893719807** for plate,
  sample and residual, matching printed 0.7169 / 3.7309 / 0.3024 at rounding
  precision. Current direct lme4 output is independently regenerated. Dataset,
  attribution and source hashes are in
  [the Penicillin manifest](../fixtures/variance_penicillin/source_manifest.json).
  This checks the shared primitive, not an ADA application.
- A separately consulted Johns Hopkins lecture prints sample variance 3.7311.
  Its difference is **0.00018212560386077215**, above the original 5.1e-5
  printed-precision bound. That **reference discrepancy remains a miss** in
  the record; it was not hidden or used to widen a tolerance. The primary
  author's printed example and current tight lme4 agree with the closed-form
  result. No runtime or simulation parameter was changed because of this.

The broader Shankar/Devanarayan procedures motivated the design, but their full
worked workflows were not reproduced. No team-measured ADA dataset was tested.

### 3. Calibration — every prespecified row

Seed **20261002**, **1,000 datasets per row**. Scenarios, target rates and
bounds were frozen in `scripts/validate_ada_calibration.py` before execution.
Cut-point rows use 80 subjects × 6 runs × 2 technical wells; log-scale subject
variance .09 and residual variance .01. Equal +/-1% technical wells have the
intended arithmetic cell mean. Confirmatory inhibition is `10 + 5*latent`.
Known NC drift is removed exactly in floating rows; noisy estimated NC is not
calibrated. One independent new negative subject per dataset gives a Bernoulli
false-positive outcome. Different rows use different fixed draws.

FPR bound = target + 2 sqrt(target*(1-target)/evaluable). This equals the roadmap
0.05 rule for screening and is stricter for the declared confirmatory/titer
rates. The machine record also retains the generic 0.05 bound. Withheld panels
are excluded from the evaluable denominator and are not successes. No execution
exceptions remained in the final simulation.

| Scenario | Evaluable / withheld | Positives | Observed FPR | Target | Upper bound | Result |
|---|---:|---:|---:|---:|---:|---|
| fixed_screening_log | 957 / 43 | 47 | 4.9112% | 5% | 6.4090% | pass |
| floating_screening_log | 951 / 49 | 65 | 6.8349% | 5% | 6.4135% | **MISS** |
| fixed_screening_nonparametric | 951 / 49 | 58 | 6.0988% | 5% | 6.4135% | pass |
| floating_screening_nonparametric | 951 / 49 | 44 | 4.6267% | 5% | 6.4135% | pass |
| fixed_confirmatory | 949 / 51 | 12 | 1.2645% | 1% | 1.6460% | pass |
| fixed_titer_log | 949 / 51 | 1 | 0.1054% | 0.1% | 0.3052% | pass |
| fixed_screening_boxcox_lambda_zero | 946 / 54 | 45 | 4.7569% | 5% | 6.4172% | pass |
| fixed_screening_biological_iqr_exclude | 950 / 50 | 68 | 7.1579% | 5% | 6.4142% | **MISS** |
| fixed_unanticipated_future_shift | 947 / 53 | 583 | 61.5628% | 5% | 6.4165% | **MISS** |

**Floating parametric screening, the biological-IQR-exclusion row, and the
unanticipated future-shift row miss their bounds.** These remain in the Skill,
results and interpretation facts. The IQR row cannot be compared causally with
the no-exclusion row because their draws differ. The extra future shift is +0.6
on the log scale; it demonstrates a fixed threshold's failure under new drift,
not a mechanism for selecting a more favorable threshold after the fact.
Passing rows do not establish performance under general Box-Cox lambda (only
lambda=0 was calibrated), non-Gaussian effects, noisy NC, heterogeneous runs,
other panel sizes or cut-point uncertainty. No cut-point interval exists here.

Variance-component intervals: coverage lower bound **0.9362159512**, including
simultaneous coverage, for each row's 1,000 evaluable datasets. Gaussian crossed
random effects are generated directly; negative raw component datasets remain.

| Scenario; levels A×B; true variances A/B/error | Component coverage A / B / error | Simultaneous coverage | Negative raw component datasets | Result |
|---|---|---:|---:|---|
| positive_components; 60×6; 0.09/0.04/0.01 | 99.4% / 98.2% / 98.3% | 95.9% | 0 | pass |
| zero_run_component; 60×6; 0.09/0.0/0.01 | 98.3% / 99.9% / 97.9% | 96.1% | 588 | pass |
| small_crossed_design; 8×4; 0.5/0.25/1.0 | 99.1% / 98.8% / 99.0% | 96.9% | 252 | pass |

[Machine record](ada_calibration_0.9.0.json) includes every cut-point draw and
all denominators, bounds and source hashes. An initial run stopped while
serializing a NumPy boolean in the progress JSON. Only that bool-to-built-in
conversion was corrected; seeds, scenarios and bounds were identical in the
complete rerun. The initial script/log are preserved in
[attempt logs](attempt_logs_0.9.0/). No old method's calibration or R benchmark
was rerun.

### 4. Skill interpretation and live misuse scenarios

The new Skill and router specify applicability, estimand, transformation,
normalization, exclusions, diagnostics, stopping conditions and calibration
misses. Claude Code **2.1.283**, model **claude-opus-5-5**, ran twelve sessions:
four roadmap misuse scenarios under three preserved attempt labels.
All executions completed; the final core stopping guards passed in all four.
A1/A3 complete scientific explanations retain overgeneralizations after two
focused guidance corrections. In particular A3 still says floating adjusts
location but not width, despite the contract distinguishing ratio-based
multiplicative drift and difference-based additive drift. This residual is
**not fixed at the agent-narrative level**. Runtime withholding is correct;
scientific review remains necessary. Repeating a prompt is not proof of general
narrative reliability. See [the English scenario assessment](agent-scenarios/0.9.0/RESULTS.md).

Other observed tier, outcome-driven exclusion and nested-lot interpretation
problems were addressed in guidance. Final A4 correctly distinguishes a possible
prespecified nested model from this implementation's unsupported scope. No
independent human re-score or real-team end-to-end agent study was performed.

### 5. Interpretation facts and immutable artifacts

ADA saves raw input, resolved config, analysis cells, per-run audit percentiles,
preprocessing/exclusion logs, results, diagnostics and source-linked facts before
the manifest. Facts copy the saved numbers, reportability and calibration record;
there is no refitting in extraction. Render verifies hashes and writes only
presentation artifacts. Tests cover withheld/flagged states, source-value
agreement, config gates, duplicate/missing cells, technical aggregation,
transform domains, negative inhibition, extreme empirical tails, whole-subject
exclusion, variance boundaries, row/scale invariance and render immutability.

## Preservation, tests and packaging

- Full `.venv/bin/python -m pytest -q`: **349 passed**, zero failures/skips,
  **181.29 seconds**. This adds 43 ADA/shared-component tests to the 306 baseline.
  Four pandas future keyword-only warnings are recorded; they do not change the
  pinned-environment results.
- Both changed Skills (`ada-cut-point`, router) passed `quick_validate.py`.
- `validate_legacy_090.py`: **34 configurations, 428 scientific artifacts**
  byte-identical to dev/0.8.2. No exceptions or explanatory changes. Thirty-five
  pre-existing Python files are unchanged; only version and dispatch changed in
  the three remaining old Python files. All old R resources are unchanged.
- **294 historical validation files** and **687 tracked run files** unchanged;
  older validation-index entry text retained. The two pre-existing untracked
  repeated-measures run directories were not modified or staged.
- `uv build` produced the local 0.9.0 wheel and sdist. **42 Python files and
  3 packaged R resources** match source bytes. Wheel SHA-256:
  `9ba0bed5f7a14a910bf204d6ff3b6c9d3c0e16912081d79cadd4ccd10d73bd90`.
- `validate_clean_environment.py` passed, including **31 byte-identical
  wheel/source facts configurations** (8 new ADA configurations), exact new
  results, historical checks at their original tolerances and plate import.
  See [clean-environment record](clean_environment_0.9.0.json).

- `validate_install_smoke.py` passed from the isolated wheel-installed Python:
  **41 configurations**, analyze/render/verify, immutable results/facts and plate
  import, from an unrelated temporary cwd with no source import path. Existing
  optional-R KR configurations were included.
- With R hidden from PATH, ADA results are byte-identical to the R-present run,
  doctor reports R unavailable, and existing KR refuses with exit 2 and a clear
  message. With R present, KR and ADA run. No new R runtime method or installer
  behavior was added; the R-present check reuses the pinned library. See
  [optional-R checks](optional_r_checks_0.9.0.json).

Commands and log hashes are in [release checks](release_checks_0.9.0.json).
Validation is on the same macOS
host in isolated Python environments; no Linux/Windows execution, fresh R
compiler installation, matched Prism project, full regulatory-validation claim
or biological generalization is made.

## Reproduction

```sh
.venv/bin/python scripts/validate_ada_benchmarks.py
.venv/bin/python scripts/validate_ada_calibration.py
.venv/bin/python scripts/validate_legacy_090.py
.venv/bin/python -m pytest -q
uv build
.venv/bin/python scripts/validate_clean_environment.py
# Use the wheel-installed Python reported by clean_environment_0.9.0.json:
/path/to/isolated/venv/bin/python scripts/validate_install_smoke.py
.venv/bin/python scripts/validate_optional_r_090.py
```

The public canonical inputs are committed. For regenerating them, download the
pinned rADA 1.1.9 archive, verify its source-manifest SHA-256, load its
`data/lognormAssay.RData` in R and `write.csv(..., row.names=FALSE)` to the
committed raw CSV path. The fixture generator records its archive/vignette
hashes from `/tmp/rADA.tar.gz` and `/tmp/rada_vignette.html`; these acquisition
files are not a runtime dependency. No generated reports are committed.
