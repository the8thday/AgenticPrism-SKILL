# Installing AgenticPrism — instructions for AI agents

You are an AI agent (Claude Code, Codex, Hermes, Cursor or similar) and a user has
asked you to install or use this Skill collection. Follow these steps in order.
Human-oriented instructions are in [README.md](README.md).

## What you are installing

Six Agent Skills (`skills/`) plus one versioned Python runtime (`src/`,
`install.py`). They are **one unit**: every Skill finds its runtime at
`<root>/.venv`, where `<root>` is the cloned repository. Therefore:

- Do **not** install individual Skill folders with a skill registry or hub
  command (for example `hermes skills install <url>`), and do not copy Skill
  folders elsewhere. A Skill separated from `<root>` cannot find its runtime.
- Do **not** `pip install` into the user's global or project Python. The
  runtime lives only in `<root>/.venv`.

## 1. Ask the user once, before changing anything

Tell the user what will happen and get agreement on:

- **Location** for the clone. Propose `~/AgenticPrism` unless they prefer another.
- **Downloads**: the git clone (about 10 MB) and Python packages (about
  220–300 MB in `<root>/.venv`). Internet access is required.
- **Skill discovery** (optional): whether to link the Skills into your agent's
  Skills folder (step 4), so future sessions find them automatically.

## 2. Clone

```sh
git clone https://github.com/the8thday/AgenticPrism-SKILL.git ~/AgenticPrism
```

If `<root>` already exists and is this repository, update it with
`git -C <root> pull` instead.

## 3. Install the runtime

```sh
cd <root>
python3 install.py        # Windows: py install.py
```

- It needs [uv](https://docs.astral.sh/uv/) (preferred; it downloads the
  validated Python 3.13 itself) **or** Python 3.12+. If neither is available,
  tell the user and stop; do not try other package managers.
- It needs network access. In a sandboxed agent (for example Codex's default
  sandbox), request permission for network access and for writing to `<root>`.
  Do not bypass the sandbox on your own.
- Success ends with `self-check passed: analyze, render and verify work`. On
  failure, show the user the error text verbatim.

## 4. Make the Skills discoverable (optional, with the user's consent)

```sh
python3 install.py --link-skills <skills folder>
```

| Agent | User Skills folder (confirm against your agent's documentation) |
|---|---|
| Claude Code | `~/.claude/skills` |
| Codex | `~/.codex/skills` |
| Hermes Agent | `~/.hermes/skills` |
| Other agents | Their documented user Skills folder. If you do not know it, skip linking and use the Skills by path (step 6). |

The command creates symlinks only and **never overwrites** an existing entry.
If it reports `skipped … something else already exists there`, the user already
has a Skill with that name: tell them, and do not delete or rename their Skill.
If symlinks cannot be created (some Windows setups), skip this step and use
the Skills by path. A new agent session may be needed before linked Skills are
discovered.

## 5. Verify

```sh
<root>/.venv/bin/agentic-prism doctor --collection <root>
# Windows: <root>\.venv\Scripts\agentic-prism.exe doctor --collection <root>
```

`status` must be `ok` or `warning`, and `package_version` must equal
`collection_version`. Report any `problems` to the user.

## 6. Use

Tell the user the collection is ready, and how to start:

> Use AgenticPrism (`<root>/skills/agentic-prism/SKILL.md`) to analyze my data.

Then read `<root>/skills/agentic-prism/SKILL.md` and follow it. The Skills
contain their own scientific checks; do not invent missing experimental facts
such as assay applicability, independent units or equivalence margins. Ask the
user for them.

## Update and uninstall

- **Update**: `git -C <root> pull`, then `python3 <root>/install.py`.
- **Uninstall**: remove only the symlinks in the Skills folder that point into
  `<root>/skills`, then delete `<root>` (which contains `.venv`).

## Troubleshooting

| Symptom | Action |
|---|---|
| `Python 3.12+ is required` | Ask the user to install uv or a newer Python; rerun step 3. |
| Download or network errors | Network is blocked. Ask the user to allow network access for this step. |
| `doctor` reports a version mismatch | The code was updated without reinstalling. Rerun `python3 install.py`. |
| A Skill says the runtime is missing | Step 3 did not complete in this `<root>`. Rerun it. |
