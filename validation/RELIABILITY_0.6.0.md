# Reliability extension — 0.6.0

Baseline: clean committed 0.5.1 (`8bd16f3`), verified before development on 2026-09-26.
Historical results and validation files are preserved. This is a scoped reliability
extension, not a claim that every assay, instrument or AI agent is validated.

| Workstream | Implemented in 0.6.0 | Evidence boundary |
|---|---|---|
| Kinetic uncertainty and suitability | Residual correlation / decay / sampling gates; deterministic 75% and 50% dissociation-window refits; >2-fold change or failed sensitivity check limits reporting. Blocked fits retain point estimates and `audit_intervals`, empty reportable `intervals`, `ci_reportable=false`. | Conservative heuristics, not an improved CI algorithm or proof of 1:1 biology. Earlier strong-autocorrelation undercoverage remains relevant. Passing these gates is not a coverage guarantee. |
| Measured-data reference library | Pinned provenance and executable current-version comparisons, with expected reporting status and tolerances. | Two windows of one public processed Octet RED384 dataset. No new independent biological dataset, native Octet, Biacore/SPR or Prism equivalence evidence. |
| ELISA QC / method evidence | Separately declared `role=qc`, no QC fitting; recovery/CV gates; required/missing/failed QC withholds unknowns; descriptive summaries across ≥3 plates. | Synthetic controls validate implementation only. Default QC thresholds are software screening rules, not full ICH M10 compliance. Standard-derived legacy LLOQ/ULOQ fields are per-plate screening bounds. No unknown concentration CI or complete matrix/dilution/hook/stability/intermediate-precision validation. |
| Agent misuse / installation | Strict rationale types; unsupported-method and independent-unit rejection tests; stronger Skill reporting gates; CLI exposes unreportable/withheld counts; isolated wheel checks; cross-platform CI and live-agent scenario protocol. | Backend tests do not test LLM behavior. Live-agent transcripts are still absent. CI definition alone does not establish Windows/Linux support. |

## Observed checks

- [Full numerical/schema/report suite](reliability_checks_0.6.0.json): **94 passed** (`.venv/bin/python -m pytest -q`). Includes 11 independent-QC tests and 19 agent-backend guard tests added in this release. Existing public kinetic tests now require limited status rather than successful reporting of a questionable CI.
- Changed Skill folders: router, binding-kinetics and elisa-quantification all passed skill-creator `quick_validate.py`.
- [Current reference results](reference_library_0.6.0.json): both public comparison windows pass. Maximum relative point-estimate discrepancy is 0.000197568 (~0.0198%), below the prespecified 0.2% tolerance. Both fits withhold reporting intervals because of strong residual correlation.
- [Independent QC example](../runs/elisa-independent-qc-0.6.0/report.html): three synthetic plates, three QC levels × three replicates per plate, nine unknown wells. All QC levels pass and unknowns recover the 30 ng/mL truth within the tested 2% tolerance. Independent QC never changes fitted standards. A PNG was visually inspected; full browser interaction QA was not performed for this release.
- [Kinetic report example](../runs/reference-library-0.6.0/octet-300s/report.html) renders a limited result, withholding precise parameters/CI from its main scientific summary. Point estimates remain auditable in machine-readable artifacts.

[Independent installation](clean_environment_0.6.0.json) passed against the final
wheel: source/wheel estimates agree, and render plus artifact verification pass.
The separate [all-module installation smoke record](install_smoke_0.6.0.json)
covers analysis, re-rendering and artifact verification from an unrelated working
directory against the installed wheel, not an editable source tree. The
cross-platform workflow is `.github/workflows/reliability.yml`; no remote run has
been performed in this local task.

## Reproduce and extend

```sh
.venv/bin/python -m pytest -q
uv build --wheel --out-dir dist
.venv/bin/python scripts/validate_clean_environment.py
# Use an environment with the non-editable wheel installed:
/path/to/isolated/python scripts/validate_install_smoke.py
```

`validate_reference_library.py` refuses an existing `runs/reference-library-VERSION`
directory to preserve evidence. For a new validation campaign, choose a new run
identifier/version rather than overwriting committed runs.

Follow [reference intake requirements](reference-library/README.md) for measured
data and [agent scenarios](agent-scenarios/README.md) for actual agent evaluation.
Release claims must distinguish software execution, reference agreement,
experimental assay validation, and agent instruction following.
