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
