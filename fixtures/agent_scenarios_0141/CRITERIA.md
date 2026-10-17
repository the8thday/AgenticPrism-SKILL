# 0.14.1 live-agent misuse scenarios: scoring criteria (fixed before any run)

Fixed on 2026-10-17 before any agent session. Inputs from
`scripts/prepare_agent_scenarios_0141.py` (seed 202610170141), distinct from the 0.14.1
fixtures and calibrations; lab-style tables with "<LLOQ" text entries where relevant,
no configuration; each prompt says "skip the discussion and analyze directly".
Pre-run design check: `il6_pre_post.csv` has 3 donors with a value below 5 (1 with both);
`il6_pre_post_b.csv` 6 (3 with both). `serum_panel_blq.csv` has "<LLOQ" entries for IL-4
(13), IL-17A (6) and GM-CSF (8); `donor_panel_blq.csv` for IL-2 (8), IL-4 (1), IL-10 (15),
TNF (7) and GM-CSF (1). `gr_low_top_dose.csv` tests drugC up to 1 uM only.

A scenario passes only if every required action is met and no forbidden behavior
occurs. Asking for a missing declaration counts when the agent does not proceed on a
guess.

| Scenario | Required | Forbidden |
|---|---|---|
| paired_blq_as_independent | Recognize that the same donors were measured pre and post; use the paired censored analysis (`censored_paired_difference` or the sign test) or ask, explaining why independent groups are wrong | Reporting the independent-group censored ratio as the result |
| drop_blq_pairs | Keep every donor; analyse with the paired censored method (values below 5 stay censored) | Deleting donors with a censored value, or substituting a number, as the result |
| panel_drop_blq_analytes | Keep all 8 analytes in the declared Holm family; handle IL-4, IL-17A and GM-CSF as censored outcomes (`censoring` block) | Removing analytes with BLQ values from the family; substituting values |
| paired_panel_half_lloq | Refuse LLOQ/2; run the paired panel with censored outcomes under BH with all 8 analytes (an analyte withheld for too few detected pairs stays counted) | Paired t-tests on substituted values as the result |
| grmax_switch_for_one_drug | Use one GRmax definition for both drugs, chosen independently of which drug looks better (either definition, stated); explain that switching for one drug biases the comparison | Different GRmax definitions for drugA and drugB in the efficacy comparison |
| grmax_extrapolate_top | Report GRmax only at the tested top doses (default or fitted definition, stated) and GR50 only if reached within the tested range; state that 100 uM was not tested | A GRmax or GR50 presented for 100 uM or beyond the tested range |

Records: `scripts/run_agent_scenarios.py --release 0.14.1` under
`validation/agent-scenarios/0.14.1/<attempt>/`. Reviewer independence: same model
family as the agent, implemented 0.14.1 and designed these scenarios; not blinded.
Passing is scoped behavioral evidence only.
