# Locating and checking the runtime

All Skills call the same versioned runtime. Follow this order even when only
interpreting saved results, and even if `which agentic-prism` returns nothing.

1. Resolve the **real absolute path of the Skill file being used**, including
   symlinks (`realpath /path/to/SKILL.md`, or Python `Path(...).resolve()`). From
   `<root>/skills/<specialist>/SKILL.md`, the collection root is three parent
   directories above the file. Confirm `<root>/pyproject.toml` and `<root>/skills`
   exist. Do not derive the root from the current working directory.
2. Check `<root>/.venv/bin/agentic-prism` on macOS/Linux, or
   `<root>\.venv\Scripts\agentic-prism.exe` on Windows. If present, use this
   absolute executable and run `doctor --collection <root>` before any analysis
   or saved-run verification. Record the resolved executable and doctor result.
3. Only if the collection executable is absent, check PATH (`command -v
   agentic-prism`, or `where agentic-prism`). If found, run that absolute
   executable with `doctor --collection <root>` against the same collection.
4. Doctor must have a non-error status and matching `package_version` and
   `collection_version`. A warning means dependency versions differ from the
   lock file; report it. An existing executable that fails or has a version
   mismatch is an unusable/mismatched runtime, not evidence of no installation.
5. Conclude **not installed** only after resolving the real Skill path, checking
   the collection executable, and checking the PATH fallback. State which paths
   were checked. For an absent or unusable runtime, explain the specific finding
   and ask before `python3 <root>/install.py` (Windows: `py <root>\install.py`),
   unless installation has already been authorized. A copied standalone Skill
   with no collection root is an incomplete collection; request its location.

A failed PATH lookup alone never justifies an installation request. Commands in
Skills use `agentic-prism` as shorthand for the resolved absolute executable.
Use that executable for `verify` on saved runs before reading interpretation facts.

In user-facing rerun commands, expand that shorthand to the actual resolved
absolute executable. Do not give a bare `agentic-prism` command after discovering
that only the collection executable is available. Saved `config.resolved.json`
uses the run's copied input; retain that relative layout and use a new output.
