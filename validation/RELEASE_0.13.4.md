# Release 0.13.4 — reliable 4PL reporting and everyday endpoints

Local development on `dev/0.13.4`, baseline `dev/0.13.3` (`1762b91180b38d46ee08cc6845ebfe895e7208d2`). No publication is part of this release.

## Scope and decisions

- Default 4PL: approximate replicate pure-error lack-of-fit F (unweighted only) and ordered dose-mean residual-sign runs diagnostic, each at 0.005. Failure withholds endpoints, summaries and relative potency. The fitted residual signs do not justify an exact-test interpretation. With no replicate pure error, or with relative weighting, that F check is ineligible. Passing either check cannot establish mechanism or rule out a hook.
- Relative ECx/ICx: opt-in declared percentages including 10/20/50/80/90; profile the requested endpoint while reoptimizing Hill and plateaus. IC90 is 90% inhibition through the fitted span, including for a decreasing response. Pointwise intervals; symmetric 4PL only. No supported extrapolated point estimate, no automatic alternative-model selection.
- Paired AUC: one trapezoid per unit and condition on the same interval, matched IDs, explicit incomplete-pair policy, paired t contrasts, Holm p values and Bonferroni intervals over the full declared family (including withheld contrasts). All retained/excluded/incomplete unit-condition curves are audited. Complete-pair inference is not a correction for informative dropout.
- `location_test`: sourced-null ordinary one-sample t, one-sample log-scale t, paired log-ratio t and independent log-scale Welch. Ratio estimates concern geometric means. Minimum three independent units/pairs; positive inputs for ratios, no fabricated controls or pseudocounts.
- Existing independent AUC and ordinary two-group methods are unchanged. The shared 4PL fitter is explicitly isolated inside ELISA and standard-calibration initialization to preserve their separate QC and numerical methods.

## Five evidence gates

| Capability | Established implementation | Published worked example | Registered calibration | Live-agent misuse | Interpretation facts |
|---|---|---|---|---|---|
| 4PL adequacy | Base-R replicate F agrees, 3/3 fixtures; combinatorial runs independently enumerated in tests | UNMET | Five bounded gate rows passed; unreplicated hook limitation retained | See behavioral record | Saved diagnostics and observed gate rows |
| ECx/ICx | Base-R reparameterized nls/confint: 10/10 reportable endpoint cases pass; 5 audit-only hook endpoint comparisons miss | UNMET | 20/20 bounded conditional-coverage rows pass; narrow-range stress remains poor | See behavioral record | Pointwise/relative definitions, range/parent withholding, full observations |
| Paired AUC | Independent R trapezoids and paired t: 3/3 contrasts pass | UNMET | 4/4 bounded family coverage/FWER rows pass | See behavioral record | Pair losses, declared family, missingness limits and evidence |
| One-sample / ratio | Base-R t.test: 4/4 methods pass | UNMET | 7/8 bounded rows pass; one-sample log-ratio coverage misses | See behavioral record | Sourced null, scale, independent units, exact miss |

No matched Prism project benchmark or general Prism equivalence is claimed. Runtime evidence is copied by `scripts/update_evidence_0134.py` from saved validation records, with no fitting in facts or rendering.

## Numerical reference checks

Only base R is used (`scripts/everyday_0134_oracle.R`). R constructs its own starting values and profiles an explicit log endpoint parameter. The initial oracle with `nls` tolerance 1e-8 did not converge on some curves (8/15 checks pass); its JSON and source are retained as `everyday_benchmarks_0.13.4.json` and `scripts/everyday_0134_oracle_initial.R`. The recheck uses standard 1e-5 convergence tolerance and fixes a read.csv type warning; it does not change the Python fitter or comparison tolerances.

Recheck `everyday_benchmarks_0.13.4_recheck.json`: **20/25 checks pass**. All ten valid endpoint cases, three lack-of-fit checks, three AUC contrasts and four location tests pass. The five audit-only endpoint estimates from the invalid hook fixture differ by 2e-5 to 1.2e-4 relative, exceeding the fixed 1e-5 point tolerance despite nearly identical SSE. All are withheld by Python; some R intervals fail to profile. These are recorded numerical misses, not discarded observations or a claim of full agreement.

The original comparison criteria were retained: location/AUC rtol 1e-9 and atol 1e-12, nonlinear estimate rtol 1e-5, SSE rtol 1e-7, interval-end differences at most 0.003 of the R width.

## Registered calibration

Script `scripts/validate_everyday_0134_calibration.py` and plan committed in **b93b8f8** before execution. Seed **131340001**, 1000 simulations per scenario, 13 scenarios, zero execution errors. JSON retains all 13000 draws, all descriptive metrics, bounds and Monte Carlo SE. **36/37 bounded metrics pass**. No calibration rerun replaced a miss.

The one-sample log-ratio coverage is 0.936, below 0.936216 (MCSE 0.007740); its type-I rate 0.064 remains within the separately prespecified upper bound. Normal/lognormal simulations do not validate arbitrary distributions.

Narrow-dose stress (descriptive, not a bounded pass): C80 is reportable in 11/1000 draws; only 1/11 intervals covers the true C80 (0.091, MCSE 0.087). C90 is withheld in all 1000. This is disclosed in every endpoint run: selecting an apparently in-range endpoint does not guarantee post-selection coverage. Extend the dose range for high-effect endpoints. In the unreplicated hook scenario the new shape checks detect none, but other existing gates withhold every fit. Do not claim the new runs diagnostic detects unreplicated hooks reliably.

| Scenario | Metric | n | Rate | MCSE | Bound | Passed |
|---|---|---:|---:|---:|---:|---|
| dose_increasing | shape_withheld | 1000 | 0.004000 | 0.001996 | 0.019439 | True |
| dose_increasing | C10_reportable_coverage | 996 | 0.950803 | 0.006853 | 0.936188 | True |
| dose_increasing | C20_reportable_coverage | 996 | 0.949799 | 0.006919 | 0.936188 | True |
| dose_increasing | C50_reportable_coverage | 996 | 0.944779 | 0.007237 | 0.936188 | True |
| dose_increasing | C80_reportable_coverage | 996 | 0.950803 | 0.006853 | 0.936188 | True |
| dose_increasing | C90_reportable_coverage | 996 | 0.945783 | 0.007175 | 0.936188 | True |
| dose_decreasing | shape_withheld | 1000 | 0.003000 | 0.001729 | 0.019439 | True |
| dose_decreasing | C10_reportable_coverage | 997 | 0.939819 | 0.007532 | 0.936195 | True |
| dose_decreasing | C20_reportable_coverage | 997 | 0.940822 | 0.007473 | 0.936195 | True |
| dose_decreasing | C50_reportable_coverage | 997 | 0.949850 | 0.006912 | 0.936195 | True |
| dose_decreasing | C80_reportable_coverage | 997 | 0.955868 | 0.006505 | 0.936195 | True |
| dose_decreasing | C90_reportable_coverage | 997 | 0.955868 | 0.006505 | 0.936195 | True |
| dose_relative | shape_withheld | 1000 | 0.000000 | 0.000000 | 0.019439 | True |
| dose_relative | C10_reportable_coverage | 1000 | 0.940000 | 0.007510 | 0.936216 | True |
| dose_relative | C20_reportable_coverage | 1000 | 0.942000 | 0.007392 | 0.936216 | True |
| dose_relative | C50_reportable_coverage | 1000 | 0.942000 | 0.007392 | 0.936216 | True |
| dose_relative | C80_reportable_coverage | 1000 | 0.944000 | 0.007271 | 0.936216 | True |
| dose_relative | C90_reportable_coverage | 1000 | 0.941000 | 0.007451 | 0.936216 | True |
| dose_unreplicated | shape_withheld | 1000 | 0.000000 | 0.000000 | 0.019439 | True |
| dose_unreplicated | C10_reportable_coverage | 1000 | 0.952000 | 0.006760 | 0.936216 | True |
| dose_unreplicated | C20_reportable_coverage | 1000 | 0.949000 | 0.006957 | 0.936216 | True |
| dose_unreplicated | C50_reportable_coverage | 1000 | 0.945000 | 0.007209 | 0.936216 | True |
| dose_unreplicated | C80_reportable_coverage | 1000 | 0.943000 | 0.007332 | 0.936216 | True |
| dose_unreplicated | C90_reportable_coverage | 1000 | 0.948000 | 0.007021 | 0.936216 | True |
| hook_replicated | shape_withheld | 1000 | 1.000000 | 0.000000 | 0.871540 | True |
| one_sample_t | coverage | 1000 | 0.956000 | 0.006486 | 0.936216 | True |
| one_sample_t | type_I_or_FWER | 1000 | 0.044000 | 0.006486 | 0.070676 | True |
| one_sample_ratio_t | coverage | 1000 | 0.936000 | 0.007740 | 0.936216 | False |
| one_sample_ratio_t | type_I_or_FWER | 1000 | 0.064000 | 0.007740 | 0.070676 | True |
| paired_ratio_t | coverage | 1000 | 0.956000 | 0.006486 | 0.936216 | True |
| paired_ratio_t | type_I_or_FWER | 1000 | 0.044000 | 0.006486 | 0.070676 | True |
| independent_ratio_welch | coverage | 1000 | 0.950000 | 0.006892 | 0.936216 | True |
| independent_ratio_welch | type_I_or_FWER | 1000 | 0.050000 | 0.006892 | 0.070676 | True |
| paired_auc | coverage | 1000 | 0.952000 | 0.006760 | 0.936216 | True |
| paired_auc | type_I_or_FWER | 1000 | 0.048000 | 0.006760 | 0.070676 | True |
| paired_auc_mcar | coverage | 1000 | 0.939000 | 0.007568 | 0.936216 | True |
| paired_auc_mcar | type_I_or_FWER | 1000 | 0.061000 | 0.007568 | 0.070676 | True |

## Compatibility and extraction

`legacy_fixture_checks_0.13.4.json`: **172 configurations**, **1270 byte-identical scientific artifacts** against the release baseline. Unchanged methods are not recalibrated. Default 4PL is intentionally changed in diagnostics/status, not optimizer estimates, predictions, original C50 intervals or resolved default configuration. Its saved grids, predictions and numeric fields compare exactly.

The public `dose_public` Compound-A is now withheld for approximate lack of fit (F=7.11998, df=4/16, p=0.00171618). Compound-B remains reportable. The independent four-parameter numerical comparison is preserved in tests; an obsolete unconditional reportability assertion was updated to assert the new state. The old default-4PL absence-of-evidence assertion was likewise updated because default 4PL now has its own evidence record. The initial full test run (793 passed, 2 such assertion failures) is retained here as observed history.

Unit tests exercise parent-fit propagation into RP, withheld estimate removal in facts, endpoint-specific profiling, direction and units, row-order-independent pairing, zero/duplicate/incomplete inputs, entirely excluded AUC curves, multiplicity, and immutable scientific hashes after rendering. Both paired AUC and location facts receive only their applicable 0.13.4 evidence; earlier independent-AUC evidence is not repurposed.

## Behavioral review

Criteria fixed before execution: `agent-scenarios/0.13.4/CRITERIA.md`. Six initial Claude CLI attempts could not run: `Not logged in · Please run /login`, zero input/output tokens. They are execution failures, not scored behavioral failures or passes; metadata is retained under `initial/`.

A separate independent forward pass uses six fresh evaluator contexts and only the realistic request, raw data, Skill and runtime. Evaluators do not receive expected answers, fixes or prior conclusions. Scoring is by the implementer, not independent validation of the whole release. See `agent-scenarios/0.13.4/results.json` and `RESULTS.md` for actual outcomes and artifact hashes. Earlier completed runs may say evidence pending because the extraction record was finalized later; their historical artifacts are not rewritten.

## Release checks

Final tests and clean archive / isolated wheel-install checks passed as recorded below. The forward pass found a CLI completion bug for location_test after valid artifacts had been written; the dispatch was fixed and all four CLI methods received subprocess regression tests (20/20 targeted tests pass). Build/install records are generated only from committed source. No `runs/`, `build/` or pre-0.13.4 validation records are changed.

- Full suite after expectation updates: **795 passed**, 8 pre-existing pandas deprecation warnings, 342.30 s. After the CLI dispatch correction: **20 targeted tests passed**, including the four added subprocess cases, 21.87 s. The original evaluator scenario was re-executed into a new directory: CLI exit 0, stdout primary equals saved results, artifact hashes verified (`agent-scenarios/0.13.4/cli_recheck.json`).
- Six independent forward behavioral scenarios passed (implementer-scored); original CLI/login failures remain recorded.
- Clean archive build: **passed**, commit `4d61c31348ab761025417e7f214b647a4ba2fe2a`; all 102 Python source files in the wheel match the archive exactly. Wheel and sdist hashes are in `archive_build_0.13.4.json`. Early offline attempts lacked setuptools in the dedicated temporary cache (one retry occurred before cache preparation completed; another required restoring uv archive links). Existing local cached distributions were copied with their archive layout, then the unchanged build command succeeded. No new statistical run or source adjustment was used to remedy this environment issue.
- Isolated wheel install: **passed**, 15 configurations and 163 scientific artifacts exactly equal to source (`clean_install_0.13.4.json`). A separate wheel-only environment uses the lockfile, no editable install, no PYTHONPATH; import location is asserted inside that environment. ICx, hook, paired AUC and one-sample log-ratio reports render and verify successfully, with unchanged science hashes after another style render. Same macOS host only.
- All four changed Skill folders pass `skill-creator` quick validation. Existing `runs/`, `build/` and historical validation remain untouched.
