# 0.15.4 — Evidence gates for independent two-way ANOVA

Baseline `3b5df44` (0.15.3). Plan: [RELEASE_0.15.4_PLAN.md](../plans/RELEASE_0.15.4_PLAN.md)
(including every amendment, each recorded before the rerun it governed).
Evidence module `src/agentic_prism/evidence_0154.py` (`scripts/update_evidence_0154.py`).

## Scope

No new analysis type and no estimator change: `routine.factorial` (`analysis_type=factorial_anova`,
0.11.1) is validated as it is: independent units, SS II/III with sum contrasts, interaction first,
equal-weight EMM contrasts in simple or marginal scope, Tukey(-Kramer), Dunnett, Sidak or Holm
(Bonferroni limits) families, marginal contrasts withheld when the interaction is material.
Factorial results and facts add `validation_evidence_0154` and its `must_mention` lines.

Why: the 0.11.1 calibration called statsmodels `anova_lm` directly rather than the package, so
the families, scopes, limits and the interaction gate were never calibrated; its R benchmark used
only balanced data (SS II = SS III); the published-example gate rested on R output rather than
printed values; and `agent-scenarios/0.11.1/RESULTS.md` still says PENDING (a reviewer 6/6 run is
reported only in the series summary). The 0.11.1 rows, including the `anova_3_equal` interaction
miss (73/1000), stay in `evidence_0111.py` and in every factorial result.

## Five evidence gates

| Gate | Observed result |
|---|---|
| Established implementation | **337/337** oracle checks ([record](factorial_benchmarks_0.15.4.json)) on 33 configurations ([fixtures](../fixtures/factorial_0154/source_manifest.json)): NIST balanced 5 × 2; seeded unbalanced 3 × 4 additive, weak-interaction and crossover designs (cell sizes 2–8); minimal 2 × 2 with n = 2; SS II and III, both factors, simple and marginal scope, all four families. R 4.6.0: car 3.1.5 `Anova`, emmeans 2.0.4 contrasts with adjusted p and simultaneous limits (Holm against Bonferroni limits), Dunnett p against high-precision `mvtnorm::pmvt` (mvtnorm 1.4.2); F, p, df, SS, residual SD, estimates and SE within 1e-8 relative; withheld contrasts null with audit values matched. SS II and SS III differ on the unbalanced fixture. |
| Published worked example | **MET**: NIST/SEMATECH e-Handbook 3.2.3.2 / 3.2.3.2.1 (machine × coolant pin diameters; public domain; original `CHARACTERIZE-3_2_3_2.DAT` pinned by SHA-256). 29/29 printed values within half a unit of the last printed digit under SS II and III: SS .000303/.00000392/.00001468/.000346, df 4/1/4/40, residual MS .0000087, F 8.8/.45/.42, and the five coolant B − A differences against the printed cell means (within 1e-4). |
| Calibration | Registered before running (`01d5dce`, syntax-only fix `21b6c73`), runtime `routine.factorial`, 2000 studies per scenario (600 for Dunnett) ([record](factorial_calibration_0.15.4.json)): **44/44 bounded criteria passed**, no errored draws. F-test size 0.042–0.058 (balanced, unbalanced, 2 × 2 with n = 2); family-wise error 0.037–0.056 (balanced Tukey 0.0485 marginal, 0.056 simple; Dunnett 0.038); simultaneous coverage 0.955–0.964 (marginal, additive) and 0.9455–0.9615 per stratum (simple, with interaction). Descriptive: with a real interaction (cells +2 and −1.5, SD 1) the gate withheld marginal contrasts in only 0.471 of studies; reported marginal intervals covered the equal-weight truth in 0.966. Stress (no bound): SD 3 in the three n = 2 cells gave F size 0.34–0.44 and Tukey FWER 0.092; SD 3 in one cell of a balanced design gave F size 0.077–0.104; t(3) errors 0.033–0.042 (conservative). |
| Live-agent misuse scenarios | **6/6 passed** (criteria committed in `7521fcd` before any run; scored by the implementing session, not blinded; [results](agent-scenarios/0.15.4/RESULTS.md)): mouse rather than well as the unit, interaction before an unqualified main effect, repeated days routed to repeated measures, no SS-type shopping, a multiplicity family instead of 66 uncorrected tests, no additivity claim from p = 0.12. The agents ran the editable runtime after the evidence wiring was added (benchmark attempt 2 state) and saw the Skill before the 0.15.4 wording below was added. |
| Interpretation facts | Extension facts from saved results: new `must_mention` lines and two evidence states; tested in `tests/test_factorial_0154.py`. |

## Attempts and corrections kept

- Benchmark attempt 1 ([record](factorial_benchmarks_0.15.4_attempt1.json)): 297/301; four Dunnett
  adjusted p values off by 1.5–2.0e-4 against emmeans. Post-hoc diagnostic
  ([script](../scripts/diagnose_factorial_0154_dunnett.py), [record](factorial_dunnett_diagnostic_0.15.4.json)):
  against `pmvt` at abseps 1e-7 the runtime differs by at most 5.5e-7 and emmeans' default (abseps
  1e-3) by up to 2.0e-4, only for families of three or more. Criterion amended to the
  high-precision reference (1e-5); the emmeans comparison is kept as informational.
- Benchmark attempt 2 ([record](factorial_benchmarks_0.15.4_attempt2.json)): passed, but the
  planned withholding check never ran because the weak-interaction fixture gave interaction p 0.087.
  A crossover fixture (earlier fixtures byte-unchanged) and a required-path row were added.
- Calibration launch 1 stopped at import on a syntax error from a pre-commit edit; no draw was
  generated ([log](factorial_calibration_0.15.4_attempt1_syntax_error.log)).
- Agent scenario 4 input: the first generator used symmetric 9/3/3/9 cells, for which SS II and
  SS III coincide in a 2 × 2 design; replaced before the criteria were committed.
- Legacy check run 1 ([log](legacy_0.15.4_attempt1.log)) stopped on two facts states sourced from
  `validation_evidence_0154`; the strip rule was widened to that source only.

## Findings for follow-up (not changed in this release)

- With two units in a cell the Brown–Forsythe (Levene, median-centred) diagnostic is degenerate:
  absolute deviations within each cell are equal, so p is near 0 (about 4e-89 in an agent run).
  The runtime still reports it; the Skill and facts say not to read it. A gate (n ≥ 3 per cell)
  would change factorial outputs and needs its own release.
- Unequal variances are not gated; a Welch-type or heteroscedasticity-consistent two-way test is
  not implemented.
- The 0.11.1 `agent-scenarios/0.11.1/RESULTS.md` is left unchanged as historical text.

## Other checks

- Legacy: 289 configurations against `3b5df44`, 1,954 scientific artifacts byte-identical; the six
  `routine_0111` factorial configurations differ only by `validation_evidence_0154`, its `must_mention`
  lines and the facts states sourced from it ([record](legacy_0.15.4.json)).
- Full suite on the final tree (`1d5a231`): 1,120 passed (8 existing warnings). An earlier full run
  had two `test_install.py` failures because the version declarations were bumped while it ran;
  that file passes on the final tree.
- Clean build from `git archive` of `1d5a231`: 131 Python source files and the packaged R bridge
  files in the wheel byte-identical to `src/` (wheel SHA-256 `cfa6c2c415b0ebb2…`;
  [record](archive_build_0.15.4.json)). The 0.15.3 wheel and sdist were moved to `dist/superseded/`.
- Isolated wheel-only install with locked dependencies: 71 configurations and 363 scientific
  artifacts identical to the source (`factorial_0154`, `routine_0111`, `quantal_0153`, `pk_0150`,
  `groups_synthetic`); `doctor` found both R bridges; four runs re-rendered and verified with
  unchanged scientific artifacts ([record](clean_install_0.15.4.json)). Same macOS host only.
