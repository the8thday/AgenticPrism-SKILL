# 0.12.1 live-agent gate — PENDING

No Claude CLI or substitute agents were run. Run on dev/0.12.1; --release selects
scenarios, not a historical runtime checkout. Keep raw traces for human scoring.

```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0121 --release 0.12.1 --attempt reviewer-1
```

Prompts and fixtures: fixtures/agent_scenarios_0121/scenarios.json.
