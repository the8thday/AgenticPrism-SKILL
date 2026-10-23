# 0.15.0 — PK execution: `pk_nca` with R PKNCA

Branch `dev/0.15.0` from `b3a49dc` (dev/0.14.4 with Prism extraction). Plan:
[RELEASE_0.15.0_PLAN.md](../plans/RELEASE_0.15.0_PLAN.md); roadmap:
[ROADMAP_0.15-0.16.md](../plans/ROADMAP_0.15-0.16.md). Evidence module
`src/agentic_prism/evidence_0150.py` (`scripts/update_evidence_0150.py`).

## Scope

New analysis type `pk_nca` (Skill `pk-nca`, upgraded from the 2026-10-13 guidance prototype):
single-dose NCA for IV bolus, IV infusion and extravascular dosing. R PKNCA 0.12.1 computes every
NCA value through a fixed packaged script (`src/agentic_prism/r_scripts_pk/pk_runner.R`) with
exact pins (`r_scripts_pk/packages.csv`, 25 packages, `install.py --with-r-pk`). The PK bridge is
separate from the MMRM bridge, whose script folder, pins and recorded hashes are unchanged.
Python checks the contract, applies declared reportability, summarises reportable subject
values and fits the Smith et al. (2000) power model. IV bolus uses two PKNCA passes because
PKNCA's `cl.obs` uses `aucinf.obs` without the back-extrapolated C0; pass 2 inserts C0 and its
AUCinf must equal pass 1 `aucivinf.obs` when both choose the same lambda-z.

Not covered: multiple dose / steady state, partial AUCs, urine, sparse or destructive sampling,
bioequivalence, unit conversion. Population PK remains guidance (0.15.1 adds runtime models).

## Five evidence gates

| Gate | Observed result |
|---|---|
| Established implementation | PKNCA 0.12.1 (CRAN 2025-08-19). Independent formula recheck in `scripts/validate_pk_0150_benchmarks.py`: 1,034 reported values on Theoph (12 subjects, regenerated from R `datasets`) and four synthetic designs agree within 7.0e-15 relative; dose-proportionality slope and interval equal R `lm`/`confint` within 2.3e-15 ([record](pk_benchmarks_0.15.0.json)). |
| Published worked example | MET. Kim et al. 2018, *Transl Clin Pharmacol* 26(1):10-15, Table 3 (Phoenix WinNonlin, Theoph subject 8, 320 mg, linear trapezoidal): Cmax 7.56, Tmax 2.02, AUClast 88.55995, t½ 8.510037883, AUCinf 103.906687, CL/F 3.079686301, Vz/F 37.81050811, all reproduced at the printed precision. With lin-up/log-down the same subject gives AUClast 86.81, so the AUC method is a required declaration. |
| Calibration | Registered before running (commit `6969b79`), 2000 draws per scenario: 6/6 bounded rows passed; 90% CI coverage for b 0.894-0.913 ([record](pk_calibration_0.15.0.json)). Descriptive power row: with exactly proportional exposure, 6 per dose, CV 20% over a 10-fold range, 0.201 of studies concluded proportional within 0.80-1.25 (3 per dose, CV 30%: 0.004). Deterministic NCA itself has no calibration. |
| Live-agent misuse scenarios | 6/6 passed in `rerun1` (criteria committed in `67f3136` before any run; scored by the implementer, not blinded): half-life shopping refused, BLQ rule requested, CL/F and Vz/F labelled, destructive sampling declined as individual NCA, AUC0-8 not relabelled AUCinf, dose proportionality reported inconclusive ([results](agent-scenarios/0.15.0/RESULTS.md)). The first attempt ended in HTTP 429 errors and is kept. |
| Interpretation facts | Generic extension facts from saved results, with per-subject parameter states, withholding reasons, the dose-proportionality conclusion and this evidence copied into `must_mention`; tested in `tests/test_pk_nca_0150.py`. |

## Attempts and corrections kept

- Pre-freeze check: PKNCA reproduced the printed WinNonlin values with the linear method before the contract was frozen.
- The first scratch run of the independent recheck disagreed on IV bolus lambda-z (up to 0.4% in t½) because the oracle used points after time 0 rather than after the analysed profile's Tmax (PKNCA's declared rule; with noisy early samples Tmax can be the second sample). The oracle was corrected; the runtime was not changed.
- After the prespecification commit, testing showed PKNCA fits manually declared lambda-z points regardless of `min.hl.points` (2 declared points gave a half-life). The runtime now withholds declared fits with fewer than `min_points`. The registered calibration tests only the power model, which is unchanged.
- First live-agent attempt (`validation/agent-scenarios/0.15.0/`): all six sessions ended with API rate-limit errors (HTTP 429) before analysis; kept unchanged and not scored.
- Population-PK template: `nlmixr2::nlmixr2` is not exported in nlmixr2 7.0.1, so the 2026-10-13 template failed; it now calls `nlmixr2est::nlmixr2` and was executed (SAEM, seed 20261013: ka 1.58, CL/F 2.76, V/F 31.5).

## Other checks

- Legacy: 273 configurations, 1,888 scientific artifacts byte-identical to `b3a49dc` ([record](legacy_0.15.0.json)); no unchanged benchmark or calibration rerun.
- Live-agent harness now symlinks the repository `.r-lib` beside the copied runtime, as `install.py` places it.
- Full suite on the final tree: `.venv/bin/python -m pytest -q` — 1,075 passed (8 existing Pandas4Warning in `ada.py`), with PKNCA pins installed; version 0.15.0.
