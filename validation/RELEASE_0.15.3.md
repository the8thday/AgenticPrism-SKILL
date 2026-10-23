# 0.15.3 — Quantal dose response and ordinal scores

Baseline `de27dba` (0.15.2). Plan: [RELEASE_0.15.3_PLAN.md](../plans/RELEASE_0.15.3_PLAN.md).
Evidence module `src/agentic_prism/evidence_0153.py` (`scripts/update_evidence_0153.py`).

## Scope

`quantal_dose_response` (`quantal.py`, Skill `quantal-dose-response`): grouped or individual
responder data; declared probit/logit/cloglog link and log10-dose or dose metameter; natural
response none (refused when controls respond), Abbott known rate or estimated from controls
(boundary at 0 refitted with c = 0 and reported); Finney heterogeneity policy; EDp with Fieller
limits, withheld when unbounded, outside the tested range, with non-positive slope, without a
partial-response group, with singular information or with a non-positive-definite covariance;
delta-method limits as audit; parallel-line relative potency with LR or slope-ratio-equivalence
parallelism.

`ordinal_regression` (`ordinal.py`, Skill `ordinal-scores`): one score per independent subject,
declared ordered categories (an empty category is refused), group factor and pre-treatment
covariates; cumulative-logit proportional odds; Wald and profile-likelihood intervals; group LR
test; nominal proportional-odds check (unavailable when the category-specific model crosses).

## Five evidence gates

| Gate | Observed result |
|---|---|
| Established implementation | 85/85 oracle checks ([record](quantal_ordinal_benchmarks_0.15.3.json)): R `glm` (probit, logit, cloglog; tight `glm.control`), base-R Fieller, `MASS::dose.p` delta SE, `optim` likelihoods for Abbott and estimated background (binomial coefficients removed), parallel/separate-slope `glm` for relative potency; `ordinal::clm` (coefficients, SEs, profile `confint`, `nominal_test`) for the small ordinal fixture, where `MASS::polr` stops about 2e-5 short (kept as `polr_coef`), and `MASS::polr` for the housing data. |
| Published worked example | Quantal MET: Finney (1952) *Probit Analysis* 2nd ed., Ex. 6/9 (rotenone, *M. sanborni*), text read from the Gokhale Institute library scan: log LD50 0.686, Pearson χ² 1.62, LD50 4.85 mg/L at the printed precision; Fieller 0.640–0.730 against the printed 0.641–0.729 (within the declared 0.002; Finney's printed solution is one iteration short of convergence, b 4.176 against 4.169). Ordinal: MASS `polr` housing example reproduced against MASS 7.3.65 output; the printed Venables & Ripley values were not checked (printed-value gate UNMET). |
| Calibration | Registered before running (`30786ff`), 1000 studies per scenario ([record](quantal_ordinal_calibration_0.15.3.json)): **9/10 bounded rows passed**. Fieller ED50 coverage 0.951 (probit, 10 per dose), 0.956 (logit), 0.950 (estimated background), ED10 0.998 and ED90 0.990 (conservative), relative potency 0.954 with LR rejection 0.051; overdispersion scaled 0.926, ignored 0.755. Ordinal group LR under the null 0.055. **Miss: ordinal OR-3 profile coverage 0.924** (bound [0.93, 0.97]; Wald 0.946); post hoc the profile limits equal `ordinal::clm` within 1.4e-4 on 15 of the missing replicates, all with large OR overestimates ([note](ordinal_profile_posthoc_0.15.3.md)). Descriptive: with 5 animals per dose ED50 was withheld in 0.115 of studies; empty declared categories refused 13–21% of small ordinal studies; the proportional-odds check was mostly unavailable and detected a true violation in 0.003. |
| Live-agent misuse scenarios | 6/6 passed (criteria committed in `2d3d2bf` before any run; scored by the implementer, not blinded): counts requested instead of averaged percentages, no extrapolated ED50/ED90, background mortality kept, parallelism caveat for relative potency, no t tests on grades, mouse as the unit for weekly scores ([results](agent-scenarios/0.15.3/RESULTS.md)). |
| Interpretation facts | Generic extension facts from saved results; tested in `tests/test_quantal_ordinal_0153.py`. |

## Attempts and corrections kept

- Pre-registration smoke run (replicates 0–1 per scenario): Fisher-scoring start values included the control rows (singular information) and the natural-response estimate overflowed at the boundary; both fixed before `30786ff` (plan).
- Launch 1 crashed in the calibration script's own delta-audit post-processing before any result was summarised ([log](quantal_ordinal_calibration_attempt1_crash_0.15.3.log)); script fixed (`f23329e`).
- Attempt 1 ([record](quantal_ordinal_calibration_attempt1_0.15.3.json), [log](quantal_ordinal_calibration_attempt1_0.15.3.log)): 9/10 bounded rows (the same ordinal miss), but 119 quantal replicates crashed on singular or non-positive-definite information under separation. The runtime now withholds these (`d01f094`) and the registered script was rerun in full; bounded quantal rows were unchanged in substance.
- Benchmark development: the first oracle runs used R `glm`'s default convergence (differences of 2e-7 to 1e-5 against our converged fit) and `polr` on the small ordinal fixture; tight `glm.control` and `clm` were adopted as oracles; tolerances were not loosened.

## Other checks

- Legacy: 282 configurations (including all 0.15.0–0.15.2 fixtures), 1,931 scientific artifacts byte-identical to `de27dba` ([record](legacy_0.15.3.json)).
- Full suite on the final tree: `.venv/bin/python -m pytest -q` — 1,115 passed (8 existing Pandas4Warning in `ada.py`); version 0.15.3.
- Clean build from `git archive` of `b5d5ac2` (0.15.3 source plus the build script; covers 0.15.0–0.15.3): 130 Python source files and all 8 packaged R bridge files in the wheel are byte-identical to `src/` (wheel SHA-256 `1587805af1e09085…`; [record](archive_build_0.15.3.json)).
- Isolated wheel-only install with locked dependencies (R packages from the repository `.r-lib` via `AGENTIC_PRISM_R_LIB`, as `install.py` places them): 18 configurations and 96 scientific artifacts identical to the source (`pk_0150`, `poppk_0151`, `fih_0152`, `quantal_0153`, `ordinal_0153`, `groups_synthetic`), `doctor` found both R bridges (PKNCA 0.12.1, nlmixr2est 7.2.1), six runs re-rendered and verified with unchanged scientific artifacts ([record](clean_install_0.15.3.json)). Same macOS host only; not a cross-platform claim.
