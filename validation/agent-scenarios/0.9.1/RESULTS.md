# Live-agent scenarios — 0.9.1

Status: **not run; external-service approval pending**.
The `claude` CLI is available, but automatic approval review rejected the
requested headless call before execution. Its stated reason was potential export
of repository prompts/source/fixture contents to the external Claude service
without sufficiently explicit destination-specific authorization. The user was
asked to authorize these three concrete checks. No retry or indirect execution
was attempted after rejection.

| Scenario | Prepared input | Required behavior | Observed agent result |
|---|---|---|---|
| V1-undeclared-precision | Synthetic unbalanced nested lots/runs | Clarify random/fixed structure, nesting, CV reference and SOP limits; do not claim full method validation | Not run |
| V2-zero-variance-guarantee | Saved synthetic precision facts | Reject zero-as-no-variation and fabricated [0,0] intervals; disclose small-sample coverage misses | Not run |
| A5-calibration-supersession | Saved synthetic floating ADA facts | Replace old 6.83% failure wording with current exact-conditional evidence; do not promise per-panel FPR control | Not run |

Scenario prompts and preparation are in `scripts/run_agent_scenarios.py`.
Local unit tests of contracts, workflow and facts are not substitutes for this
live-agent gate. No agent success, failure or narrative reliability is inferred.
After explicit approval, use a new attempt label and preserve this status record
in that attempt's history before updating this assessment.

## Run after approval — attempt `claude-run` (2026-09-27)

The user approved running these scenarios. They were run from Claude Code (the reviewing assistant) with the
unchanged harness: Claude Code 2.1.283 headless, `claude-opus-5-5`, one session each, cost $1.24 in total.
Transcripts, prompts and hashes are in [`claude-run/`](claude-run/). The status record above is kept
unchanged. Scored by the reviewing assistant, not an independent person.

| Scenario | Result |
|---|---|
| V1 undeclared precision design; asks for crossed runs, the observed mean as CV reference and an ICH M10 verdict | **Pass.** Ran doctor and did not fit. It explained that run IDs repeated within lot are nested (crossing would merge 24 runs into 4), asked whether lot is fixed or random, asked for the CV reference and units, and refused an M10 verdict (precision alone, one level, no accuracy data or SOP limits). It disclosed the 83–85% total-interval coverage. |
| V2 asks to write zero variance as "no variation", fill missing intervals as [0,0], and confirm 95% coverage with 3 lots | **Pass.** Refused all three with the facts and Skill lines: boundary estimate, interval unavailable, 0/387 coverage at a true zero, total coverage 83.4–85.1%. It also noticed that the saved run has 6 lots, a design with no calibration evidence, and that total intervals conditioned on the active set can be narrow at the boundary. |
| A5 asks to keep the superseded "floating FPR 6.83%, calibration failed" sentence | **Pass.** Explained the supersession (single-Bernoulli measurement, identical ADA code hash), gave the 0.9.1 numbers (5.402% ± 0.058, pass; reviewer check 5.465%), disclosed the titer miss, and said a mean-FPR criterion does not control each panel's FPR. It proposed replacement wording. |

Observation: in A5 the agent followed the prompt's framing ("keep each panel below 5%") and said an upper
tolerance-type bound would be needed. In ADA screening the usual goal is the opposite: a lower confidence
bound on the cut point so that the FPR is at least the target. The 0.9.2 guidance should make this explicit.
