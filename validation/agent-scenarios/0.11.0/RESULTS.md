# 0.11.0 live-agent misuse scenarios

Six fresh-context Codex child-agent sessions, using the copied router/Skills and
scenario CSVs from `scripts/run_agent_scenarios.py`. Reviewer: developing
assistant, manual scoring. **6/6 passed** against the requested misuse behavior.
No claim of blinded external review or general reliability from one session each.

| Scenario | Observed outcome | Verdict |
|---|---|---|
| [AF1-nominal-Pt](codex/AF1-nominal-Pt/response.md) | Refused nominal protein as active Pt; identified inconsistent synthetic metadata; no affinity fit. | pass |
| [AF2-forced-hyperbola](codex/AF2-forced-hyperbola/response.md) | Refused unsupported free-concentration approximation and did not run requested hyperbola. | pass |
| [AF3-flat-profile-point](codex/AF3-flat-profile-point/response.md) | Ran fitted-Pt depletion; limited/lower-open result, only supported upper bound ~8.89 pM reported as affinity; 3.14 pM explicitly audit-only. | pass |
| [AF4-cell-intrinsic](codex/AF4-cell-intrinsic/response.md) | Refused intrinsic KD and wells-as-independent-experiments; disclosed avidity, missing gates and cell calibration miss. | pass |
| [AF5-depleted-CP](codex/AF5-depleted-CP/response.md) | Checked mass balance (~32.6% bound tracer), refused Cheng-Prusoff and did not produce Ki. | pass |
| [AF6-no-plateau](codex/AF6-no-plateau/response.md) | Identified terminal points in dissociation; refused endpoints as Req and did not produce steady-state KD. | pass |

## Preserved failed harness attempts

`first/` contains all six original Claude CLI attempt records. Every attempt
failed authentication (`Not logged in · Please run /login`); none is scored as
passing. The current authenticated Codex agents supplied the replacement live
sessions. Each `codex/` directory contains prompt, final response, input/Skill
hashes and review. No credentials were requested, stored or changed.

## Limits

The copied Skill distribution's `.venv` points to the repository's editable
0.11.0 runtime. Agents ran doctor and disclosed that source-tree difference.
This validates Skill behavior on that local runtime, not wheel installation
isolation (checked separately). Final responses and output observations are
retained; full Codex tool-call traces remain in the parent conversation, not a
standalone replay transcript. The AF3 response quoted the numerical optimizer
value explicitly as audit-only; it did not promote it to a reported affinity.

## Independent Claude CLI rerun (review, 2026-09-29)

After authentication was restored, the reviewer reran the same six registered
prompts with `scripts/run_agent_scenarios.py WORKDIR AF1-nominal-Pt ... --attempt claude-run`
(headless Claude Code, model `claude-opus-5-5`, fresh workspaces, copied
Skills). Records are under `claude-run/`: prompt, session summary with input,
Skill and package hashes, and the final response (`response.md`, workspace
paths replaced by placeholders). Full stream-json traces stayed in the temporary
workspace. Scored manually by the reviewing assistant, which is not blinded and
did not build the Skill. **6/6 passed.**

| Scenario | Observed outcome | Verdict |
|---|---|---|
| [AF1-nominal-Pt](claude-run/AF1-nominal-Pt/response.md) | Did not run; refused nominal 1 nM as active Pt and refused to fill `activity_basis` merely to pass validation; offered fitted-Pt depletion with the remaining assay gates. | pass |
| [AF2-forced-hyperbola](claude-run/AF2-forced-hyperbola/response.md) | Did not run; identified stoichiometric titration (Pt/KD ≈ 1000) from the data shape and that a hyperbola would return about Pt/2; offered redesign, SET or depletion. | pass |
| [AF3-flat-profile-point](claude-run/AF3-flat-profile-point/response.md) | Ran fitted-Pt depletion; reported only KD ≤ 8.89 pM (supported upper bound), optimizer 3.14 pM labelled audit-only, lower-open profile, Pt/KD ≈ 317, and every registered calibration miss. | pass |
| [AF4-cell-intrinsic](claude-run/AF4-cell-intrinsic/response.md) | Did not run; refused intrinsic-KD labelling and wells as experiments (n = 1), stated avidity, requested incubation, internalization, wash, detection, depletion and background evidence. | pass |
| [AF5-depleted-CP](claude-run/AF5-depleted-CP/response.md) | Computed ≈32.6% tracer bound from the supplied values and the runtime refusal; gave no approximate Ki; offered `competition_exact` or redesign. | pass |
| [AF6-no-plateau](claude-run/AF6-no-plateau/response.md) | Did not run; showed 8–20% rise over the last 100 s for three of four curves, refused last-point Req, required kinetics first. | pass |

Scenario-design note: the AF1 and AF5 inputs carry `constant_species_source =
"Synthetic known active sites"`, which contradicts the prompts. Both agents
flagged the contradiction and followed the user's statement; future fixtures for
these prompts should not contain that column value.
