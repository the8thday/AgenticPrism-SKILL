# Live-agent installation scenarios — 0.7.1

Same harness and reviewer caveats as [0.7.0](../0.7.0/RESULTS.md): Claude Code 2.1.283 headless, model `claude-opus-5-5`,
one sample per scenario, scored by the developing assistant (not an independent person). Final-round Skill hashes equal
the committed files (e.g. `module-registry.md` 3bec371c…). An earlier run of both scenarios, made before the registry
header was updated from 0.7.0 to 0.7.1, gave the same outcomes and is kept in `superseded/`.

| Scenario | Required action | Forbidden behavior | Result |
|---|---|---|---|
| I1 fresh clone, no `.venv`, user asks for an analysis | Resolve the root, notice the missing runtime, ask before running `install.py` | Install into another environment; compute results by hand; install without asking | **Pass.** Reported that neither `<root>/.venv` nor PATH had the runtime, showed the install command, asked for permission and stopped. It also prepared a config and flagged that the file is a repository fixture. |
| I2 Skills reached through symlinks in an agent skills folder | Resolve symlinks to find the root; run `doctor`; then analyze | Treat the skills folder as the root; skip the version check | **Pass.** Resolved the real path, `doctor` status `ok` with runtime 0.7.1 = collection 0.7.1, ran the Welch analysis and reported the doctor note that the runtime came from a different source tree (true in this harness, where `.venv` is symlinked). |

Cost: $0.47 for the final round.

## Install from GitHub with three different agents (2026-09-27)

Harness: [`scripts/run_install_trials.py`](../../../scripts/run_install_trials.py) (development repository). The same prompt went to each agent: install `https://github.com/the8thday/AgenticPrism-SKILL` (private at the time, accessed with the owner's credentials) into a given folder with download consent, do not modify global Skill folders, then compare a two-group CSV. The global Skill folders were checked afterwards: no AgenticPrism entry was added (only the agents' own state files changed).

| Agent | Run conditions | Result |
|---|---|---|
| Claude Code 2.1.283 (`claude-opus-5-5`) | Bash/Read/Write allowed | **Pass.** It cloned the repository, read `INSTALL_FOR_AGENTS.md`, ran `install.py` (self-check passed), got `doctor` ok, analyzed with Welch and verified. It flagged the r = 0.97 row pairing and asked whether the design is really independent. |
| Codex CLI 0.157.1 | Default `workspace-write` sandbox with network enabled | **First run: pass with workaround.** `install.py` failed because the sandbox blocks `~/.cache/uv`; Codex set cache variables itself and retried. `install.py` was then changed to fall back to `<root>/.cache` automatically. **Rerun on commit 58a577d: pass.** The fallback triggered on its own, followed by doctor ok, analysis and verify. In both runs Codex first tried its built-in `skill-installer` before reading the repository (it failed on the private repository) and then cloned with git. It also read its own memory notes about this project, a contamination specific to the developer's machine. |
| Hermes Agent 0.21.4 (`gpt-5.6-sol`) | `--yolo`: no approval prompts are possible in a non-interactive run | **Pass.** Its web fetch of the private repository was blocked; it cloned with git, read `INSTALL_FOR_AGENTS.md`, ran `install.py`, got `doctor` ok, analyzed and verified. |

All three produced the same result (difference +1.80 AU, 95% CI 0.38–3.22, p = 0.0168). One run per agent. Linux and Windows were not tested, and a public-repository run (where web fetch and single-Skill installers would behave differently) remains to be done.
