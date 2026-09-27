# Locating and checking the runtime

All Skills in this collection call one versioned command-line runtime,
`agentic-prism`, that lives next to them. Every analysis number comes from it.

1. **Find the collection root.** Take the real path of the Skill file you are
   reading, resolving symlinks (for example `realpath SKILL.md`). The root is the
   folder that contains `skills/`: two levels above the Skill's own folder. A Skill
   folder copied on its own has no root and cannot work.
2. **Use the root's runtime.** macOS/Linux: `<root>/.venv/bin/agentic-prism`;
   Windows: `<root>\.venv\Scripts\agentic-prism.exe`. Do not rely on the current
   working directory.
3. **If that file does not exist,** the collection has not been installed. Tell the
   user, and with their permission run `python3 <root>/install.py` (Windows:
   `py <root>\install.py`). It creates `<root>/.venv`, downloads the pinned
   dependencies (network access needed) and runs a self-check. Do not install
   packages into another environment instead.
4. **Check it once per session:** `<runtime> doctor --collection <root>`. `status`
   must not be `error`, and `package_version` must equal `collection_version`. On a
   version mismatch (for example after `git pull`), ask the user to rerun
   `install.py`. `warning` means dependency versions differ from the release lock
   file; say so when reporting results. An `agentic-prism` found on PATH is
   acceptable only if this check passes against the same root.

Commands in the Skills write `agentic-prism` for the resolved runtime.
