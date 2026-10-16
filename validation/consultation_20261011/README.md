# Consultation interaction evidence

The criteria, synthetic inputs and hashes of all 73 Skill files were saved
before the six fresh-context evaluations. No Skill file changed during those
evaluations. The evaluators received the applicable Skill, the researcher
request and an isolated temporary workspace; they could inspect the related
contracts/runtime but not tests, validation, plans, historical runs or other
cases. They were not given expected answers. Questions used ordinary chat,
not the real user's question tools.

`behavioral_results.json` records the implementer's review against the fixed
criteria: six scenarios, including two follow-up conversations, eight turns in
total. `observed/` preserves the actual saved researcher-facing responses,
briefs/assessments where produced, and numerical results/configuration/facts
from the two authorized runs. It does not contain full private tool-session
logs, reports or figures. Numerical result/fact files are copied unchanged;
the original manifests' scientific hashes were rechecked against all files.

In the `.txt` transcript copies only, `${REPOSITORY}` and `${WORKSPACE}` replace
machine-specific roots. Links inside these text records describe the original
temporary workspace, not files shipped here. Both original and stored text
hashes are recorded. Record paths in JSON are repository-relative. Original
CSV files were checked unchanged after every case.

The questions, answers and data are synthetic workflow exercises. Passes are
limited to the tested interactions, scored by the implementer; they are not an
independent scientific review, a particular host UI guarantee, or statistical
calibration. The full Skill collection's source hashes and the compatibility
record identify what was actually tested.
