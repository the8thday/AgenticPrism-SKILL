# 0.15.1 — Nonlinear antibody PK: `population_pk` with R nlmixr2

Branch `dev/0.15.0` (continued) from `2e1ec73` (0.15.0). Plan:
[RELEASE_0.15.1_PLAN.md](../plans/RELEASE_0.15.1_PLAN.md). Evidence module
`src/agentic_prism/evidence_0151.py` (`scripts/update_evidence_0151.py`).

## Scope

New analysis type `population_pk` (Skill `population-pk`, upgraded from the 2026-10-13 guidance
prototype): one population fitted across dose levels with two-compartment linear, linear plus
Michaelis–Menten, or quasi-steady-state TMDD elimination, after IV bolus, IV infusion or SC
dosing (F fixed with a source). R nlmixr2est 7.2.1 / rxode2 5.1.8 (FOCEi) estimate every value
through `r_scripts_pop/pop_runner.R` (fixed templates; JSON parameter names whitelisted and
turned into symbols only through `bquote`) with 67 exact pins (`install.py --with-r-pop`,
installer `install_pop.R`, which compares versions with `package_version`). The MMRM and NCA
bridge folders are unchanged. Gates withhold parameters for failed covariance, |correlation| >
0.95, saturation outside the observed concentration range, natural-scale RSE > 100% (jointly for
nonlinear parameters), and individual estimates for shrinkage > 30%.

Not covered: covariates, interoccasion variability, PD, full kon/koff TMDD, M3 BLQ, VPC, bootstrap.

## Five evidence gates

| Gate | Observed result |
|---|---|
| Established implementation | nlmixr2 (FOCEi). Independent SciPy re-solution of 326 IPRED values (linear infusion, MM IV, QSS SC fixtures): maximum relative difference 1.7e-3, which **misses the prespecified 1e-3 tolerance** at low concentrations ([record](poppk_benchmarks_0.15.1.json)). Post-hoc diagnostic ([script](../scripts/poppk_0151_ipred_tolerance_diagnostic.R), [output](poppk_ipred_tolerance_diagnostic_0.15.1.json)): rxode2 solving the same equations at tolerance 1e-12 reproduces the SciPy values, while nlmixr2's reported IPRED differs from rxode2 by 0.16-0.17% at every rxode2 tolerance tried; the model equations and mapping agree, and the cause inside nlmixr2's reported IPRED is not established. |
| Published worked example | UNMET for printed estimates. `nlmixr2data::Bolus_2CPTMM` (Schoemaker et al. 2019, *CPT:PSP* 8:923, CC BY-NC; 120 subjects, multiple doses): the published simulation truth lies inside the 95% CI for Vc (70.98; 67.2–75.0), Q, Vmax (976; 901–1058) and Km (230; 205–257), but not Vp (truth 40, estimate 49.8, CI 45.2–54.8); the model's extra linear CL (truth 0) is correctly withheld as not identified. The paper prints estimates only in figures. No published QSS-TMDD example with data was located. |
| Calibration | Registered before running (`807d15c`), rerun in full after attempt 1 (see below), 100 simulated studies per scenario ([record](poppk_calibration_0.15.1.json)): **9/9 bounded rows passed**, no fit errors. MM IV (1/3/10 mg/kg × 4): coverage CL 0.92, V1 0.94, Vmax 0.89, Km 0.89 among reportable fits; every parameter reportable in 0.79 of studies. Linear-truth data fitted with MM: Vmax and Km both passed the gates in 0.02 (bound ≤ 0.10). Linear infusion: coverage CL 0.96, V1 0.94, Q 0.90, V2 0.90. QSS SC (descriptive): KSS reportable in 0.17 of fits, every parameter in 0.13; coverage CL 0.91 (74 fits), R0 0.76 (17 fits). PK data alone rarely identify QSS target parameters; Vmax/Km coverage is below nominal 0.95 with 12 animals. |
| Live-agent misuse scenarios | 6/6 passed (criteria committed in `bd5cbb0` before any run; scored by the implementer, not blinded): one population fit instead of per-dose fits, declared model with descriptive comparisons, TMDD not claimed from one dose, total/free asked, no absolute CL or F from SC-only data, no ranking under 61% shrinkage ([results](agent-scenarios/0.15.1/RESULTS.md)). |
| Interpretation facts | Generic extension facts from saved results: parameter states and reasons, IIV and shrinkage, comparisons, evidence; tested in `tests/test_population_pk_0151.py`. |

## Attempts and corrections kept

- Pre-registration: a smoke run of replicate 0 per scenario showed a linear-truth fit whose Vmax CI spanned 1.9e-15 to 2.9e7 yet passed the earlier gates (RSE was only a flag); RSE withholding (> 100%) and joint nonlinear withholding were added before the registration commit `807d15c` (recorded in the plan). Those smoke draws are part of the registered set.
- During development: nlmixr2 reports the sandwich covariance as `"r,s (full)"` and internally fixed thetas carry no SE, which first marked every fit as failed covariance; `f$shrink["var", ]` is the EBE variance, not shrinkage (shrinkage is `Shrink(SD)%` in `parFixedDf`); nlmixr2's `%RSE` divides by the log-scale estimate, so a natural-scale RSE is used. All fixed before registration.
- Calibration attempt 1 ([kept](poppk_calibration_attempt1_0.15.1.json), [log](poppk_calibration_attempt1_0.15.1.log)): 9/9 bounded rows passed, but 32 replicates (17 linear-null, 15 QSS) crashed on an exp overflow when nlmixr2 returned a huge log-scale SE for an unidentified parameter. Fixed (`06ce843`: the parameter is withheld as not identified) and the unchanged registered script was rerun in full.
- Pins: R normalises `15.6.0-1` to `15.6.0.1`, so the string comparison used by `r_scripts/install.R` would fail for such pins; the population installer and runner compare with `package_version`, and the MMRM installer was not changed.

## Other checks

- Legacy: 277 configurations (including the four 0.15.0 `pk_nca` fixtures through the unchanged PKNCA bridge), 1,908 scientific artifacts byte-identical to `2e1ec73` ([record](legacy_0.15.1.json)); no unchanged benchmark or calibration rerun.
- Full suite on the final tree: `.venv/bin/python -m pytest -q` — 1,092 passed (8 existing Pandas4Warning in `ada.py`), with PKNCA and nlmixr2 pins installed; version 0.15.1.
