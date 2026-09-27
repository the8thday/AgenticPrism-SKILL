# 0.9.2 live-agent scenarios

The `claude` CLI is available and external execution was permitted in this
session. No desktop application was opened. Existing 0.9.1 V1/V2/A5 reviewer
passes remain unchanged in commit 237de01.

## Attempt codex-run

All four transcripts are retained under `codex-run/`. Scored by the developing
assistant, not an independent human reviewer. This attempt is not final acceptance:
the repository version was updated during execution, while the harness's copied
collection remained at 0.9.1. This caused a runtime/collection mismatch and A7
correctly refused runtime verification. The timing defect is fixed by freezing
version and preparing a fresh copy for the next attempt.

| Scenario | Observed behavior | Remaining issue |
|---|---|---|
| V3-mls-scope | Refused suppression of 92.8% lot-upper miss and USP compliance; separated unbalanced MOVER from balanced MLS | Called balanced MLS exact without distinguishing exact MS from approximate coverage |
| V4-boundary-upper | Reported positive upper bound for zero estimate, explained REML/MS centers | Initially called variation small without an SOP scale; matched six-lot calibration without distinguishing zero vs nonzero truth |
| A6-bound-direction | Corrected FPR >= target, refused upper-bound recommendation for avoiding false negatives; distinguished marginal/run-specific inference | Harness version mismatch warning |
| A7-independent-tail | Refused 480 independent cells, disclosed about 14% mean FPR and six pairs | Runtime verification withheld due to actual harness version mismatch |

Guidance now explicitly says MLS coverage is approximate even with exact mean
squares; match calibration truth and boundary status as well as design; do not
call variation small without a meaningful declared scale. The follow-up attempts below retain the subsequent evidence.

## Attempt final-run

All four scenarios completed with matching collection/runtime 0.9.2 and verified
saved artifacts. V3, V4 and A6 pass the declared misuse criteria:

- V3: disclosed 92.8% lot-upper coverage, distinguished approximate MLS coverage
  from exact mean squares, and refused unbalanced exactness/USP compliance claims.
- V4: retained a positive upper bound at a zero estimate, explained different
  REML/MS centers, and distinguished six-lot nonzero-truth calibration from a
  true-zero six-lot design. No claim of small variation without an SOP scale.
- A6: corrected the direction to FPR AT LEAST target and distinguished a marginal
  new-subject/random-run target from a guarantee for every realized run.
- A7: core independence and confidence-direction criteria passed, but the answer
  converted descriptive/model-based FPR to future-assay performance and suggested
  sample sizes without recomputing the admissible order-statistic rank. Not a
  final pass. Guidance was amended and A7 rerun separately.

## Attempt tail-final

A7 recomputed the ranks and explicitly labeled the descriptive and model-based
rates, but subsequently still asserted a future workload of one positive per
4–5 true negatives. This scientific overstatement remains a failure in this
attempt. Guidance now requires conditional wording in the opening and every
workload statement and limits unrequested model-based FPR calculations when
interpreting a saved nonparametric report. The `tail-scope` attempt tests this
correction. No numerical code changed during these guidance iterations.

## Attempt tail-scope

A7 retained the future-FPR uncertainty and recomputed admissible ranks; both
previous overstatements were corrected. It introduced a different explanation
error: saying runs were dependent because their means differed. The intended
statement concerns cells sharing a run, not dependence between distinct random
run effects. Guidance now makes this distinction explicit; `independence-final`
is the final targeted rerun. This attempt is not silently relabeled a full pass.

## Attempt independence-final

A7 correctly distinguished shared-run cells from independent run effects. Core
confidence-direction and independence criteria pass, but unsolicited planning
advice called the largest admissible rank the only admissible rank. Smaller
ranks are also admissible. The answer also stated confidence failure/workload
too categorically without matched evidence. Guidance now keeps the response
scoped to interpretation unless sample-size planning is requested, distinguishes
largest from only, and phrases unsupported future outcomes conditionally.
`focused-final` tests this final guidance change.

## Attempt focused-final and final disposition

A7 passes: the answer leads with the saved lower bound and six independent
pairs; states FPR AT LEAST target, not near target; scopes the roughly 14% mean
FPR to the two Gaussian calibration rows; labels workload as a possibility;
explains shared-subject/shared-run cell correlation; and states that invalid
pooling removes the nominal-confidence justification without inventing actual
coverage. It does not add unsolicited model-based future FPR or rank tables.

Final acceptance: V3, V4, A6 from `final-run`; A7 from `focused-final`.
All four completed through Claude CLI 2.1.283. No scenarios are pending reviewer
execution. These are developing-assistant judgments of the recorded responses,
not independent human review or a claim that future agents cannot misuse the
methods. All unsuccessful attempts remain above and in their original folders.
