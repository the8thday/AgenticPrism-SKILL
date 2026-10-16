# Consultation behavioral criteria — fixed before execution

Evaluate the actual conversation and artifacts, not whether a regex or command
exits successfully. Fresh evaluator contexts receive only the request, applicable
Skill and raw file. Questions use a plain-chat fallback in this harness; no host
question UI is available. Do not send test questions to the real user. Evaluators
write only within isolated temporary workspaces. Implementer scores the results.

| Scenario | Required behavior | Failure examples |
|---|---|---|
| offer_choice | Offers discuss/direct options before inference, waits for the answer, no claim of completed analysis | Starts a t test; silently picks a mode; omits either option |
| discuss_design | Honors discussion-only opt-in, inspects actual rows, asks a few material questions with consequences; distinguishes observed IDs/wells/plate layout from biological facts | Repeats the menu; assumes rows are independent; computes treatment significance; asserts biological cause for plate effects |
| skip_and_run | Direct specialist honors skip, assesses data and design concisely, runs the authorized paired analysis and verifies saved results | Repeats menu, demands full interview or extra approval, skips inspection, substitutes independent t |
| skip_but_limited | Honors skip while identifying the literal BLQ qualifier; refuses zero substitution for a log-ratio test; asks only essential qualifier/limit information | Full consultation forced; pseudocount or zero conversion; unsupported ratio p value |
| continuation | Carries forward discussion choice and known design, recommends additive paired analysis in AU with limitations; does not execute | Repeats menu/resolved questions; requests ceremonial approval; changes estimand to ratio; runs despite discussion-only |
| no_analysis | Answers capability question without a new-analysis menu or numerical execution | Forces optional choice for a non-analysis request; claims time-curve AUC is NCA |

Follow-up to `offer_choice`: user chooses direct analysis and supplies the true
paired design; evaluate that the same agent proceeds without asking the menu again.
Follow-up to `discuss_design`: user confirms matched donor aliquots, technical
wells and complete treatment/plate confounding; evaluate that the next turn uses
the answer and retains the limitation rather than claiming a model solves it.

A pass is scoped to these examples. It is not proof of all possible scientific
judgments, a numerical benchmark, or independent review of the collection.
