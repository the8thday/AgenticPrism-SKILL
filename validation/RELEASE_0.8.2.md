# Release evidence — local development 0.8.2

Branch `dev/0.8.2` starts from `dev/0.8.1` (`536e99f`). This release implements
only the requested 0.8 follow-ups. No 0.9 work, push or publication was performed.
The roadmap's Python-default rule and five gates were used. Numerical agreement,
simulation calibration and live-agent interpretation are reported separately.

## Scope and implementation

| Requested item | Status | Delivered scope |
|---|---|---|
| Runtime discovery / F1 | Done | Resolve the Skill's real path, locate the collection, check its `.venv/bin/agentic-prism`, run doctor, then consider PATH. Failure of `which` alone cannot justify installation. New attempts preserve the 0.8.1 failure. |
| Nonparametric tests | Done, within the declared contract | Mann–Whitney and paired signed-rank, HL shift/pseudomedian and interval, exact small untied inference; normal approximation for ties; KW with Dunn Holm/Bonferroni; complete-block Friedman. |
| MMRM / KR | Done, with calibration misses | Between-arm × categorical-visit marginal REML, common US or ordered-visit AR(1), numeric baseline covariates, Python Satterthwaite and optional R KR. US small-sample simulation rows miss bounds; see below. |
| One-factor bootstrap B=999 | Done | Recalibration of the existing implementation; no numerical change to the legacy bootstrap. |
| Interpretation facts | Done for both new analysis types; partial overall | Both new methods save source-linked results, reportability, intervals, diagnostics and limits. Other older modules beyond the three already covered in 0.8.1 were not retrofitted. |

The two existing specialist Skills were extended with separate config contracts:
`analysis_type=nonparametric` under group-comparison and `analysis_type=mmrm`
under repeated-measures. All new options live in `DEFAULTS`; resolvers reject
unknown settings and require literal-true applicability declarations plus
rationales. Unit identity, pairing, visit order, missingness and contrast family
are explicit. Technical replicates are not silently promoted to independent units.

### Numerical choices and limits

- **Rank tests:** exact distribution counts use integer dynamic programming;
  exact rank inversion gives the small-sample interval. Automatic exact inference
  requires fewer than 50 observations per sample, no ties and (signed-rank) no
  zero differences. Forced exact inference with ties refuses. R 4.6 can perform
  additional conditional exact tied calculations; those are not implemented.
  Tied data use R's `exact=FALSE` convention. The product exposes two-sided tests.
- HL is the median of pairwise differences or Walsh averages. For normal
  approximation, the R-compatible root estimate is also saved, separately from
  sample HL. Its inverted interval concerns the same population location
  parameter. Neither is a generic difference of medians. Unbounded and reduced
  confidence intervals are explicit; a nominal normal-approximation confidence
  level is not empirical coverage. Root tolerances limit precision near zero.
- **Dunn:** pooled ranks with tie correction, two-sided normal p values, then
  canonical Holm or Bonferroni. `dunn.test` 1.4.1's native Holm step values may
  be nonmonotone and require its stopping rule. The oracle therefore compares
  its Z/raw p with `stats::p.adjust` adjusted p, retaining native step values in
  the record. This is not a claim that every native adjusted-p field is identical.
  Friedman provides an omnibus test only, without post-hoc comparisons.
- **MMRM:** cell means plus declared numeric baseline columns; no additional
  random intercept. US uses a Cholesky covariance; AR(1) uses variance and a
  transformed correlation. Python REML and the existing derivative-based
  Satterthwaite machinery are used by default. Visit lag is the declared visit
  order, not elapsed time. Covariance is shared across arms. Available cases
  require MAR conditional on the model. Input-count gates allow fitting; they
  do not establish adequate information or power.
- The MMRM primary test is arm × visit interaction. Visit-specific arm-minus-
  control contrasts use Holm p values and Bonferroni simultaneous intervals.
  Numerical failures withhold inference. This does not add KR to the old
  random-intercept contract, a one-arm MMRM, two-way bootstrap, heterogeneous
  arm covariance, nested effects or MNAR correction.

## Why KR uses R

Only MMRM Kenward–Roger uses R at runtime. Its adjusted coefficient covariance,
second-order derivatives and F scaling are delegated to the pinned mature
`mmrm` implementation rather than introducing a separate, unvalidated Python
KR implementation. Python remains the default for the new rank tests and
MMRM REML/Satterthwaite, and all pre-existing methods.

One shared `r_bridge.py` sends JSON files to `Rscript --vanilla` and receives
JSON. Fixed R scripts ship inside the package; user strings never construct R
source or formulas. `install.py --with-r` opts into pinned packages in `.r-lib/`;
R itself must already be installed. `doctor` reports availability and versions.
The run manifest includes R, pinned/loaded package versions and script hashes,
including when a fitted R model is withheld. Missing R refuses only KR clearly.

Observed versions include R 4.6.0, mmrm 0.3.18, jsonlite 2.0.0, dunn.test 1.4.1,
nlme 3.1.169 and pbkrtest 0.5.5. The complete installed pin set is recorded in
[the optional-install check](optional_r_install_0.8.2.json); validation-only
packages remain in the ignored `.r-lib/`.

## Five gates and numerical evidence

### 1–2. R agreement and published worked examples

`scripts/validate_nonparametric_benchmarks.py` regenerates R `wilcox.test`,
`kruskal.test`, `friedman.test` and Dunn outputs on public examples and tied
synthetic fixtures: **14 configuration/variant rows passed their predeclared
tolerances**. Exact inference is checked independently by enumeration in tests.

| Rank-test comparison | Observed maximum difference |
|---|---|
| Published untied MW: exact statistic, p, HL and interval | zero |
| Published untied signed-rank: exact outputs | zero |
| KW/Dunn compared statistics and p values | absolute 1.1102231e-16; relative 2.8291666e-16 |
| Friedman, public and tied fixture | absolute 1.0658141e-14; relative 5.4022351e-15 |
| Approximate HL/root estimates and interval endpoints | absolute 2.4132113e-5 |
| Near-zero tied MW lower endpoint | absolute 1.0662337e-5; **relative 0.7049680** |

The last row is not high relative precision: both implementations use a 1e-4
root tolerance; the predeclared absolute endpoint tolerance is 3e-4. Relative
errors are retained, not hidden behind a blanket agreement claim. Full fields,
tolerances and output hashes are in
[nonparametric_benchmarks_0.8.2.json](nonparametric_benchmarks_0.8.2.json).


A later boundary check found cancellation in `1 - CDF` for a completely
separated 49-versus-49 sample: before correction Python returned
1.3322676295501878e-15 instead of the combinatorial probability
7.850029192963324e-29. Both exact null distributions are symmetric, so the
implementation now evaluates the opposite tail directly. Three regression
cases (n=12, 20, 49 per group) compare with `2 / choose(2n,n)` and regenerated
R exact probabilities, using a relative tolerance of 1e-12. Their observed
maximum relative error is 9.853679946e-15. The rank benchmark and all 20
nonparametric calibration rows were rerun after this numerical fix; every
calibration count/rate was unchanged. The fix affects only the new exact
rank-test implementation.

The public examples reproduce Hollander and Wolfe (1973) datasets as worked in
the official R help: placental permeability MW, Hamilton paired depression
scores, mucociliary KW and rounding-time Friedman. The original MW one-sided
`greater` p is 0.1272061272061272; the Python exact-tail identity gives
0.1272061272061272 (absolute difference zero). Two-sided intervals
and Dunn are explicitly supplemental calculations. The Friedman help says
18 players, but its matrix contains 22 rows; all 22 published rows were retained.
Source URLs, help snapshots, generator and input hashes are in the
[rank-test source manifest](../fixtures/nonparametric/source_manifest.json).

`scripts/validate_mmrm_benchmarks.py` regenerates direct R `mmrm` and independent
`nlme::gls` outputs for **eight configurations** (public/synthetic × US/AR1 ×
Satterthwaite/KR). It compares coefficients, covariance matrices, interaction
and contrast inference. Scaled tolerance is 2e-3 × max(1, |reference|) for Python
and independent fits; the KR/direct-package tolerance is 1e-10 on that scale.

| MMRM comparison | Observed maximum relative difference |
|---|---:|
| Python Satterthwaite vs direct mmrm, public AR1 | 8.5164103e-9 |
| Python Satterthwaite vs direct mmrm, public US | 1.8109159e-6 |
| Python Satterthwaite vs direct mmrm, synthetic AR1 | 1.9716840e-9 |
| Python Satterthwaite vs direct mmrm, synthetic US | 8.6674808e-7 |
| KR bridge vs direct mmrm, all four datasets/covariances | 5.9969275e-15 |
| Independent nlme GLS fit checks, all eight configurations | 5.0421591e-5 |
| Separate complete-CS check: mmrm KR-Linear vs pbkrtest F / denominator df | 5.7880792e-7 / 6.0447366e-7 |
| Same CS p value | 2.0369174e-5 (absolute 6.0830640e-20) |

The pbkrtest check uses complete Orthodont data and the matched
**Kenward–Roger–Linear** covariance adjustment in mmrm. It is an independent
special-case check, not evidence of pbkrtest equivalence for the product's
US/AR1 full-KR adjustment. Independent nlme checks validate fits, not KR df.
See [mmrm_benchmarks_0.8.2.json](mmrm_benchmarks_0.8.2.json).

The worked example is the [mmrm 0.3.18 introduction](https://openpharma.github.io/mmrm/v0.3.18/articles/introduction.html):
`FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID)`.
Equivalent cell-means coding plus baseline indicators reproduces its model.
The public FEV dataset is simulated clinical example data, not antibody-study
evidence. Of 800 source rows, 537 are complete on the published model columns;
the 263 omissions are explicitly recorded. AR1 and KR are supplemental variants.
Source, data, R generator, original default fit and tighter fit hashes are in the
[MMRM source manifest](../fixtures/mmrm/source_manifest.json).

The original default R optimizer gave a public US covariance mismatch above
the preset tolerance. Matching BFGS convergence at `reltol=1e-12` resolved it;
no benchmark tolerance was widened. Both the original default worked output
and the tighter reference are preserved. Their maximum coefficient absolute
change is 0.0001839390 (maximum relative 0.00390645); covariance absolute change
is 0.0087751093 (relative 0.00511780). Thus the tighter computation is not claimed
to reproduce every printed default-optimizer digit. The initial CS comparison
also mixed full and linear KR adjustments (F relative difference 0.00352749);
matching the adjustment, without changing the bound, fixed that comparison.

No matched Prism or SAS benchmark was run, and no equivalence with them is claimed.

### 4–5. Guidance, live misuse checks and facts

The router and group/repeated specialists explain runtime discovery, estimands,
ties and zeroes, interval precision, missing blocks, visit spacing, missingness,
withholding and calibration limits. New scenarios N1, N2 and M1 are added to
`scripts/run_agent_scenarios.py`; F1 is rerun under new labels. Full transcripts,
input/source hashes, model/version and results are preserved in
[the scenario record](agent-scenarios/0.8.2/RESULTS.md). Reviews here were performed
by the developing agent; there was no independent human re-score.

Both new workflows write `interpretation_facts.json` before the manifest. Facts
copy saved numbers with source pointers; extraction performs no new inference.
Tests cover exact-distribution enumeration, pairing/applicability gates,
unbounded or reduced-confidence intervals, fit failure, missing R, source-value
agreement and scientific-file immutability after rendering. Fitted covariance,
R environment and raw audit values remain available without overriding
reportability. Other older modules' facts remain outside this release's scope.

## Calibration: every prespecified row

All rows generated 1,000 datasets. Bounds were fixed before simulation: rejection ≤ 0.05 + 2√(0.05×0.95/n), coverage ≥ 0.95 − 2√(0.05×0.95/n), using the declared evaluable denominator. At n=1,000 these are 6.378405% and 93.621595%. No seed, scenario or bound was changed after results were seen. A passing conditional interval row does not mean a requested 95% interval was available in every dataset.

### Rank tests — seed 20260928

Independent standard-normal nulls or symmetric paired normal differences; ties were generated by rounding to integers. KW and Friedman use three groups/conditions. All tests were evaluable in all 1,000 datasets in each row. Family rejection refers to all-pairs Dunn; HL coverage concerns a zero shift/pseudomedian.

| Method | n per group / blocks | Ties | Adjustment | Test rejection | Family rejection | HL interval coverage | 95% intervals / reduced level | Result |
|---|---:|---|---|---:|---:|---:|---|---|
| mann_whitney | 8 | no | none | 4.5000% | — | 95.5000% | 1000 / 0 | pass |
| mann_whitney | 8 | yes | none | 3.9000% | — | 96.8000% | 1000 / 0 | pass |
| mann_whitney | 12 | no | none | 4.9000% | — | 95.1000% | 1000 / 0 | pass |
| mann_whitney | 12 | yes | none | 3.4000% | — | 97.1000% | 1000 / 0 | pass |
| wilcoxon_signed_rank | 8 | no | none | 4.1000% | — | 95.9000% | 1000 / 0 | pass |
| wilcoxon_signed_rank | 8 | yes | none | 2.4000% | — | 98.3766% | 616 / 384 | pass |
| wilcoxon_signed_rank | 12 | no | none | 3.0000% | — | 97.0000% | 1000 / 0 | pass |
| wilcoxon_signed_rank | 12 | yes | none | 3.9000% | — | 97.6264% | 969 / 31 | pass |
| kruskal_wallis | 8 | no | holm | 4.0000% | 2.8000% | — | — | pass |
| kruskal_wallis | 8 | no | bonferroni | 4.7000% | 4.8000% | — | — | pass |
| kruskal_wallis | 8 | yes | holm | 4.0000% | 3.6000% | — | — | pass |
| kruskal_wallis | 8 | yes | bonferroni | 5.5000% | 4.4000% | — | — | pass |
| kruskal_wallis | 12 | no | holm | 4.2000% | 3.9000% | — | — | pass |
| kruskal_wallis | 12 | no | bonferroni | 5.2000% | 5.0000% | — | — | pass |
| kruskal_wallis | 12 | yes | holm | 3.7000% | 2.5000% | — | — | pass |
| kruskal_wallis | 12 | yes | bonferroni | 5.3000% | 4.9000% | — | — | pass |
| friedman | 8 | no | none | 4.5000% | — | — | — | pass |
| friedman | 8 | yes | none | 5.5000% | — | — | — | pass |
| friedman | 12 | no | none | 6.2000% | — | — | — | pass |
| friedman | 12 | yes | none | 4.7000% | — | — | — | pass |

The tied signed-rank rows had 384/1,000 (n=8) and 31/1,000 (n=12) intervals whose saved confidence level was below 95%. They were excluded from the 95% coverage denominator and were not counted as successes. No interval was unavailable in these rows. Coverage bounds for the remaining 616 and 969 intervals were 93.2437% and 93.5997%. These normal/rounded-normal nulls do not validate skewed paired differences or unequal-shape location-shift interpretations.

[Machine record](nonparametric_calibration_0.8.2.json); reproduce with `scripts/validate_nonparametric_calibration.py`.

### MMRM — seed 20260930

Two arms × 12 units × three visits, zero cell means, matching Gaussian covariance and 15% MCAR omissions. US truth is [[1,.4,.2],[.4,1.5,.5],[.2,.5,2]]; AR(1) has variance 1 and ρ=.5. Missingness alone is redrawn until each arm/visit has six observations. The primary test is arm × visit interaction; the family is arm minus control at each visit. The four rows use different fixed draws, so differences between methods are not paired simulation comparisons.

| Covariance | Inference | Reportable / withheld | Interaction rejection | Family rejection | Simultaneous coverage | Result |
|---|---|---|---:|---:|---:|---|
| unstructured | satterthwaite | 1000 / 0 | 7.0000% | 5.5000% | 94.5000% | **MISS** |
| unstructured | kenward_roger | 1000 / 0 | 7.1000% | 6.7000% | 93.3000% | **MISS** |
| ar1 | satterthwaite | 1000 / 0 | 5.7000% | 4.3000% | 95.7000% | pass |
| ar1 | kenward_roger | 1000 / 0 | 5.3000% | 4.8000% | 95.2000% | pass |

**US/Satterthwaite misses the interaction rejection bound. US/KR misses both rejection bounds and the coverage bound.** These misses are retained. KR is computed by the pinned R package itself, and the fixture comparisons establish numerical agreement only. No change to the covariance model, cutoffs, seed or exclusion rule was made to make these rows pass. Skills and interpretation facts disclose the US limits. No claim is made for MNAR, covariance misspecification, unequal arm covariances or covariate calibration.

[Machine record](mmrm_calibration_0.8.2.json); reproduce with `scripts/validate_mmrm_calibration.py --work NEW_DIRECTORY`.

### One-factor bootstrap B=999 — seed 20260929

The same fixed data draws, per-dataset seeds and bounds as the original one-factor script were used. Only bootstrap replications changed from 199 to 999. Four conditions, Gaussian subject intercept, residual SD=1, 15% MCAR; contrasts are conditions versus the first. Each row generated 1,000 datasets.

| Units | Subject SD | Reportable / withheld | Omnibus rejection | Family rejection | Simultaneous coverage | Result |
|---:|---:|---|---:|---:|---:|---|
| 8 | 2 | 1000 / 0 | 4.4% | 4.4% | 95.5% | pass |
| 12 | 2 | 1000 / 0 | 4.8% | 5.9% | 94.0% | pass |
| 12 | 0.3 | 1000 / 0 | 4.5% | 4.4% | 95.6% | pass |

The formerly missed 12-unit/SD=2 row has B=999 coverage 94.0%, compared with the historical B=199 value 93.4%; its fixed lower bound is 93.621595%. This does not establish coverage under a different covariance, non-Gaussian data, MNAR dropout, or the current user's experiment.

The run was resumed from preserved checkpoints while worker counts changed (6 → 12 → 6 → 10) to manage host load. Saved `(scenario,index)` records were retained, and the final checkpoint has 3,000 unique rows. No seed, draw, bound or scenario was selected or changed in response to the inference results. The JSON seconds field measures the last continuation only. Checkpoint hashes and prefix-preservation checks are in the release checks; prior scheduling logs are preserved under `attempt_logs_0.8.2/`.

[Machine record](bootstrap_999_calibration_0.8.2.json); reproduce with `scripts/validate_bootstrap_999.py --work NEW_DIRECTORY --workers 6`. The `--resume DRAW_FILE` option continues an interrupted checkpoint without changing the generated data or random seeds.

## Preservation and defects found during validation

`validate_legacy_082.py` compares the 0.8.1 source with the new source in one
Python environment, without rerunning old R benchmarks or old calibration
scripts. **22 configurations reproduce 295 scientific artifacts byte for byte.**
One additional artifact changes only its bootstrap-evidence limitation text;
all estimates, tests, pointers and source hashes in that facts file are unchanged.
Thirteen legacy numerical source files are byte-identical. The record is
[legacy_fixture_checks_0.8.2.json](legacy_fixture_checks_0.8.2.json).

- The runtime-location failure was fixed by an ordered procedure and verified
  in fresh F1 attempts. Initial explanation defects about fixed-effect SEs,
  within-unit clustering, HL confidence and incomplete blocks were addressed
  in guidance and rerun; complete live explanations still have the scoped
  residual problems documented in the scenario record.
- The tiny exact-MW upper-tail cancellation was fixed numerically and checked
  against R and an independent combinatorial probability; all affected new-method
  checks were rerun. No legacy model was changed.
- The public US oracle's loose optimizer and the mismatched full/linear KR
  comparison were corrected as described above. Original failed logs are kept in
  [attempt_logs_0.8.2](attempt_logs_0.8.2/).
- The first MMRM calibration attempt stopped after its initial Python segment
  when the R runner was being edited and R read an incomplete file. This was a
  validation orchestration error. The fixed runner was held unchanged for the
  complete rerun, using the same seed, scenarios and bounds. The initial error
  log is preserved; it is not counted as an inferential rejection or withheld
  model in the final 4,000-dataset run.
- Claude CLI returned `subtype=success` and exit zero for two API 429 quota
  errors. The recorder now retains the explicit error and terminal fields.
  Those attempts are not passes. The subsequent valid attempts are separate.
- US calibration misses were not hidden, converted into a new sample-count
  cutoff, or repaired by changing the simulation. The method agrees with its
  stated numerical oracle but does not meet those small-sample error bounds.

The historical validation index text is retained. **212 historical validation
files and 687 tracked run files remain unchanged** against `dev/0.8.1`. The two
pre-existing untracked repeated-measures run directories were neither edited nor
staged. No historical benchmark record was copied to a new version just to
repeat an unchanged check. Only new-method and B=999 evidence is added.

## Tests, packaging and installation

- The router, group-comparison and repeated-measures Skills passed the final
  `quick_validate.py` checks after the bootstrap and disclosure guidance updates.
- Final full pytest: **306 passed**, zero failures/errors/skips, console elapsed
  time **395.92 seconds**. This includes 34 nonparametric and 22 MMRM tests;
  the 0.8.1 baseline contained 250 tests.
- `uv build` succeeded. **38 Python files and three packaged R resources**
  match the source byte for byte. The final numerical wheel SHA-256 is
  `dfaff7f762a9bd5faa585a11b03cbc189b88574a569493f83edb972076179686`.
  Earlier same-version wheels and sdists are preserved in `dist/superseded/`.
- The final wheel passed `validate_clean_environment.py`: **23 configurations**
  have byte-identical interpretation facts between wheel and source, including
  eight rank-test and four synthetic MMRM configurations. Results also match the
  source exactly; historical-run comparisons retain their original tolerances.
  See [clean_environment_0.8.2.json](clean_environment_0.8.2.json).
- `validate_install_smoke.py` passed from the isolated wheel-installed interpreter:
  **33 configurations**, with analysis, render, verify, results/facts immutability
  and the plate import check. The interpreter ran from an unrelated temporary
  working directory with no source `PYTHONPATH`. KR was available and both KR
  fixture covariances were included.
- Optional-R installation checks passed in two fresh Python virtual environments.
  Without R on PATH, doctor reports optional R unavailable, KR exits 2 with a
  clear message, and Python MMRM runs. With `--with-r`, doctor reports the pinned
  versions and KR/Python MMRM run. The enabled case reuses this host's already
  installed pinned `.r-lib/` through a symlink: it is not a clean R download or
  compiler test. See [optional_r_install_0.8.2.json](optional_r_install_0.8.2.json).

No independent team data or matched Prism/SAS analysis, Linux/Windows run,
independent human scenario re-score or fresh R compiler installation was
performed. Two-way bootstrap, single-arm MMRM, KR for the legacy random-intercept
contract, other older modules' facts, tumor/survival follow-ups and 0.9 remain
outside this release. The numerical core, new fixtures and checks are locally
committed in logical steps; no push, publication or public-tree build was run.

The final live-agent core guards passed in F1/N1/N2/M1. F1's complete narrative
remains partial because it still overgeneralizes covariance reliability at the
variance boundary, despite the explicit narrower Skill rule. This remaining
interpretation limitation is recorded, not treated as a passing full narrative.
The final M1 response stopped before fitting and did not propose US, so that
attempt does not exercise every prospective US disclosure. See the full
[scenario assessment](agent-scenarios/0.8.2/RESULTS.md).

Commands, logs, wheel/source hashes, optional-R scope, preserved checkpoints and
file-preservation checks are recorded in
[release_checks_0.8.2.json](release_checks_0.8.2.json).
