# 0.14.3 — follow-ups to 0.14.2

Baseline `c5ad2a3` (0.14.2, local). Plan: [plans/RELEASE_0.14.3_PLAN.md](../plans/RELEASE_0.14.3_PLAN.md).
Both additions are opt-in; existing outputs are checked against 0.14.2 at the end.

## Item 1 — RO exposure-occupancy EC50: `exposure_response.interval: cluster_robust_t`

[Skill](../skills/receptor-occupancy/SKILL.md), [contract](../skills/receptor-occupancy/references/input-and-model.md).
Bell–McCaffrey (CR2) cluster-robust variance of log EC50 (and log h) at the least-squares
estimate, subject as cluster, Satterthwaite degrees of freedom (Pustejovsky–Tipton working
independence). The default stays the 0.14.2 percentile bootstrap so older configurations
reproduce; the Skill recommends the new interval. Prespecified in `0934d73`.

| Gate | Observed scope |
|---|---|
| Established implementation | R `clubSandwich` 0.7.0 (`vcovCR` type CR2, `coef_test` Satterthwaite) on the linearized design after a base-R fit: 5/5 fixtures within 1e-6 relative (largest 2.7e-8), including h estimated and a withheld case. [Record](receptor_occupancy_benchmarks_0.14.3.json) |
| Published worked example | UNMET. |
| Calibration | Disclosed before registration: a 300-draw exploratory pilot on two 0.14.2 scenarios chose CR2 with Satterthwaite df over CR1/CR2 with t(G−1) and z. Registered: seed 141501 (0.14.2 draws reused for matched scenarios), 1000 draws per scenario: **3/3 bounded criteria passed** — EC50 coverage among reportable fits 0.953 (free/baseline, 8 subjects), 0.963 (bound/saturation), 0.961 (dose 1, 50 % bracketed). On the same draws the 0.14.2 bootstrap covered 0.878 and 0.888. Descriptive: 6 subjects 0.948 (median df 4.7); h estimated: EC50 0.941, h 0.945; 12 subjects 0.958. Coverage counting withheld fits as not covered (the 0.14.2 registered metric) 0.836 and 0.832, because only 0.877 / 0.864 of dose-5 datasets were reportable (occupancy often stayed above 50 %); at dose 1 all were reportable. [Record](receptor_occupancy_calibration_0.14.3.json) |
| Live-agent misuse | 2/2 passed (robust interval reported over the requested bootstrap; animal kept as the unit instead of 40 independent samples). [Record](agent-scenarios/0.14.3/RESULTS.md) |
| Interpretation facts | Extraction-only; the run states which interval and its df |

## Item 2 — enzyme inhibition, `tight_binding` (Morrison)

[Skill](../skills/enzyme-inhibition/SKILL.md), [contract](../skills/enzyme-inhibition/references/input-and-model.md).
Inhibitor titrations at one substrate concentration per experiment; Morrison quadratic in a
cancellation-free form; declared total enzyme fixed or active enzyme estimated; profile-F
intervals; Ki from the declared mechanism with declared Km (mixed: Ki_app only); range,
depletion, titration-regime and active-enzyme flags. Prespecified in `2113192`.

| Gate | Observed scope |
|---|---|
| Established implementation | Base R `nls` (port) on log parameters with profile-F limits found exactly (uniroot on constrained nls): 13/13 reportable experiment fits within 1e-4 relative (largest 5.0e-6); two experiments withheld by this implementation were not compared (R's nls failed on both). **First attempt kept:** comparing limits with `stats::confint`, which interpolates the profile, failed (limits differing by up to 4.7e-4 or NA in R). [Record](enzyme_morrison_benchmarks_0.14.3.json), [attempt 1](enzyme_morrison_benchmarks_attempt1_0.14.3.json) |
| Published worked example | UNMET: no printed Morrison dataset with fitted values located (GraphPad's guide gives the equation only). |
| Calibration | Seed 143201, 1000 draws per scenario (three experiments each), all kept: **7/8 bounded criteria passed.** Ki coverage among reportable fits 0.951 (Ki_app = 0.1 Et), 0.943 (= Et), 0.953 (= 10 Et) with enzyme fixed; 0.959 with active enzyme estimated (0.6 of declared); summary Ki coverage 0.950, 0.949, 0.954. **Miss retained:** active-enzyme coverage 0.822 (bound 0.936), open limits counted as not covered. Descriptive: 40 % inactive enzyme fixed as active gave Ki coverage 0.192 (log10 error −0.066); a classical fit without depletion at Ki_app = Et gave log10 Ki error +0.146; titration regime (Ki_app = Et/100, enzyme estimated) reportable in 0.273 with coverage 0.912 and log10 error +0.297 among those. [Record](enzyme_morrison_calibration_0.14.3.json) |
| Calibration, post hoc | Same draws regenerated from the seeds (reportability identical): the lower active-enzyme limit was open in 0.147 of reportable fits and two-sided intervals covered 0.964 — the miss comes from open limits, not interval width. [Record](enzyme_morrison_active_enzyme_posthoc_0.14.3.json) |
| Live-agent misuse | 4/4 passed (Cheng–Prusoff refused for Morrison; active enzyme estimated instead of an untitrated label; Ki_app only when the mechanism is unknown; titration-regime Ki reported as an upper bound). One weaker point: an agent paraphrased the active-enzyme miss as narrow intervals. [Record](agent-scenarios/0.14.3/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results |

## Checks across 0.14.3

- Rapid-equilibrium enzyme outputs: three sentences that said Morrison analysis was not
  implemented now point to the `tight_binding` option (declared in the plan before the check).
- Unchanged methods: 257 configurations against 0.14.2 (`c5ad2a3`); 1,790 scientific artifacts
  byte-identical; the five 0.14.2 enzyme fixtures differ only in that wording (and the facts'
  record of the results hash), every number identical ([record](legacy_fixture_checks_0.14.3.json)).
  The first run of the check stopped on that results-hash field; the script now allows exactly
  that path.
- Full suite: `.venv/bin/python -m pytest -q` 971 passed (2026-10-17).
- Live-agent misuse scenarios 6/6 (criteria `8c8a1f3`, not blinded).
- R validation package added to the git-ignored `.r-lib`: clubSandwich 0.7.0 (CRAN 2026-05-04,
  maintainer J. Pustejovsky).
- Clean build from `git archive` of `a831764`: 118 Python source files in the wheel are
  byte-identical to `src/` (wheel SHA-256 `86fa57d244849dfc…`; [record](archive_build_0.14.3.json)).
  The offline build first failed because the local uv cache had been cleared; the cache was
  refilled online and the build rerun from the same archive.
- Isolated wheel-only install with locked dependencies: 34 configurations and 179 scientific
  artifacts identical to the source (all 0.14.1–0.14.3 fixture folders plus `groups_synthetic`);
  six runs re-rendered and verified with unchanged scientific artifacts
  ([record](clean_install_0.14.3.json)). This also covers 0.14.1 and 0.14.2, which had no
  separate build. Same macOS host only; not a cross-platform claim.
