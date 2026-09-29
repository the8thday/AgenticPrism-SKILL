# 0.12.0 live-agent gate — PENDING

No Claude CLI or substitute agents were run. Reviewer command:

```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0120 --release 0.12.0 --attempt reviewer-1
```

Five prompts and their fixtures are fixed in `fixtures/agent_scenarios_0120/scenarios.json`. The reviewer must retain raw traces and score scientific refusal/rationale.

Run on `dev/0.12.0`. `--release` selects the prompt set, not a historical
runtime checkout. Reviewing a later branch does not retrospectively validate
this release's frozen code.
