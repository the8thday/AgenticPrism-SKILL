# Release evidence — development 0.7.1 (installable distribution)

Changes since 0.7.0 are packaging only; no numerical code changed. The analysis
modules are unchanged; see [RELEASE_0.7.0.md](RELEASE_0.7.0.md).

| Change | Evidence | Observed result / boundary |
|---|---|---|
| `install.py` (standard library only) creates `<clone>/.venv` with uv, or venv + pip, installs the package editable with dependencies constrained by `requirements-lock.txt`, then self-checks by running doctor, analyze, render and verify | Fresh copies of the repository without `.venv`, on macOS | uv path: Python 3.13.9 fetched by uv, self-check passed, `doctor` status `ok` with no pin mismatches (≈50 s). pip path (`--no-uv`) on the system Python 3.14.7: self-check passed (≈80 s). Five reference fixtures (multi-group, 5PL ELISA, single-cycle kinetics, equivalence potency, public KD) gave **identical** `results.json` values on 3.14 and 3.13 (maximum relative difference 0). Linux and Windows not executed. |
| `agentic-prism doctor [--collection ROOT]` reports version, interpreter, dependency-pin mismatches and whether runtime and collection versions match; exits 1 on error | `tests/test_install.py` | Matches on this repository; errors on a folder without `skills/`; the CLI exit code is 1 |
| `--link-skills DIR` symlinks Skill folders for agents that discover Skills from a directory, and never overwrites existing entries | `tests/test_install.py` | Symlinks resolve to the collection; a pre-existing folder is left untouched; a rerun reports links already in place |
| Shared runtime instructions (`skills/agentic-prism/references/runtime.md`) used by all six Skills: resolve symlinks, use `<root>/.venv` (with the Windows path), ask before running `install.py`, run `doctor` | [agent-scenarios/0.7.1](agent-scenarios/0.7.1/RESULTS.md) | Live agent: fresh clone → asked permission to install (I1); symlinked Skills → resolved the root and ran doctor (I2). A personal absolute path was removed from the equilibrium Skill. |
| `requires-python` corrected from `>=3.11` to `>=3.12` | Package metadata of the pinned wheels | The pinned numpy 2.5.3, scipy 1.18.1 and contourpy 1.4.0 require Python ≥ 3.12 |
| Install from GitHub with Claude Code, Codex and Hermes | [agent-scenarios/0.7.1](agent-scenarios/0.7.1/RESULTS.md#install-from-github-with-three-different-agents-2026-09-27) | All three cloned the repository, installed, ran doctor, analyzed and verified with identical results. Codex's sandbox blocked the home uv cache, so `install.py` now falls back to `<root>/.cache` and `<root>/.python` (checked again with a read-only `HOME` for both uv and pip, and in a Codex rerun). |
| CI workflow now also runs the installer path (install, doctor, pytest) on three OSes | `.github/workflows/reliability.yml` | Definition only; there is no git remote, so it has not run |

Suite: [162 passed](release_checks_0.7.1.json) (158 + 4 installer tests) at the first 0.7.1 build; the cache-fallback test added afterwards brings it to 163.
[Isolated wheel](clean_environment_0.7.1.json) and [installed-runtime smoke](install_smoke_0.7.1.json)
checks passed against the same 0.7.1 wheel on macOS.
