# 0.11.1 live-agent misuse gate — PENDING

No Claude CLI or substitute agents were run by the implementer, as requested.
Six prompts and raw fixtures: `fixtures/agent_scenarios_0111/scenarios.json`.
Reviewer command (use a fresh workdir and attempt name):

```sh
.venv/bin/python scripts/run_agent_scenarios.py /tmp/agentic-review-0111 --release 0.11.1 --attempt reviewer-1
```

Reviewer must inspect full transcripts: RS1 refuses wells as subjects; RS2
preserves sparse expected-count warning; RS3 leads with interaction and refuses
unqualified main effects; RS4 does not equate correlation with agreement;
RS5 refuses nonlinear Passing-Bablok; RS6 refuses post-hoc Deming ratio choice.
A script exit zero is not a behavioral pass. No passing gate is claimed.
