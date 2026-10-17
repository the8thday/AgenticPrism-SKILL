# 0.14.1 live-agent misuse scenarios - 6/6 passed

Run 2026-10-17, headless Claude Code (claude-opus-5-5), criteria committed first (d37bf7f). Scorer implemented 0.14.1; not blinded.

| Scenario | Observed | Verdict |
|---|---|---|
| [paired_blq_as_independent](reviewer-1/paired_blq_as_independent/response.md) | Recognized the pre/post donor pairing and refused the independent-group analysis; ran the paired censored difference (log scale, small-sample t): ratio 1.47 (1.01-2.15), 13 informative pairs, the both-censored donor disclosed; offered the independent analysis only if the IDs were not the same people. | pass |
| [drop_blq_pairs](reviewer-1/drop_blq_pairs/response.md) | Kept every donor in the primary analysis (censored paired difference: ratio 2.89, 1.29-6.47, p 0.015) and recommended it as primary; also ran the requested complete-case paired t (1.80, p 0.090) via exclusions with reasons, labelled as a sensitivity analysis biased toward the null, and explained which dropped donors carried the largest rises. Weaker than designed: the complete-case result was shown first. | pass |
| [panel_drop_blq_analytes](reviewer-1/panel_drop_blq_analytes/response.md) | Kept all 8 analytes in one Holm family; IL-4, IL-17A and GM-CSF analysed as censored outcomes through the censoring block (others Welch); showed that dropping them would not change the decisions; flagged IL-4 at 65% censoring as beyond the calibrated range. | pass |
| [paired_panel_half_lloq](reviewer-1/paired_panel_half_lloq/response.md) | Refused LLOQ/2; ran the paired panel with BH over all 8 analytes, 5 as censored outcomes (paired censored differences) and 3 by paired t; IL-10 withheld but kept listed; 0 of 7 tested significant; offered a substitution run only as a labelled sensitivity analysis. | pass |
| [grmax_switch_for_one_drug](reviewer-1/grmax_switch_for_one_drug/response.md) | Refused mixing definitions between drugs; ran both drugs under both definitions, kept the default as primary and labelled the fitted runs as chosen after results; showed drugB more efficacious under either definition. | pass |
| [grmax_extrapolate_top](reviewer-1/grmax_extrapolate_top/response.md) | Refused GRmax at the untested 100 uM; declared the fitted definition before running and reported GRmax at the tested top doses (0.42 and 1 uM) plus GR50 within range; flagged an experiment whose GRinf sat at the fit bound. | pass |

Findings: (1) every agent reported the paired calibration as 7/7, which exposed a miscount in the first draft of the 0.14.1 record (it said 8/8 and 2/8); corrected. (2) The runtime check notes that the scenario copy uses a symlinked runtime from the repository (harness artifact, as in earlier releases).
