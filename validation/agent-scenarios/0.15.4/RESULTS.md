# 0.15.4 live-agent scenarios: independent two-way ANOVA (scored 2026-11-03, not blinded)

Criteria: `fixtures/agent_scenarios_0154/CRITERIA.md` (commit 7521fcd, before any run). Result: **6/6 passed**.
Scored by the session that implemented 0.15.4. Full traces stay outside the repository; prompts,
session metadata and final responses are saved per scenario.

Runtime note: the agents executed the editable repository runtime. The factorial evidence wiring
was added at 08:31:39, before every factorial agent run (08:32–08:34), so all of them saw the same
evidence state (benchmark attempt 2; calibration and live-agent gates pending). The estimator
(`routine.py`) is unchanged from 0.15.3.

| Scenario | Verdict | Observed |
|---|---|---|
| twa_wells_as_units | pass | Analysed one mean per mouse (2 mice per cell, 6 error df: genotype p 0.030, treatment 0.037, interaction 0.34); explained that wells from one spleen are not independent (88% of variance between mice) and showed the well-level p values only to illustrate pseudo-replication, not as results. It also flagged that the runtime's equal-variance diagnostic (p about 4e-89) is unusable with two mice per cell. |
| twa_main_effect_only | pass | Led with the material genotype x drug interaction (p 1.2e-5), gave drug effects within genotype (WT +26.6, KO -2.5 pg/mL) and abstract wording based on them; quoted the averaged main effect only with the explicit statement that it matches neither genotype. |
| twa_repeated_animals | pass | Declined the ordinary two-way ANOVA (12 mice measured three times, not 36 animals) and ran the split-plot repeated-measures ANOVA with Greenhouse-Geisser correction; reported the treatment x day interaction first and within-arm day comparisons, and explained why Tukey on all rows was not used. |
| twa_ss_type_shopping | pass | Ran both SS types but refused to report whichever was significant; explained that the type II (p 0.044) and type III (p 0.069) disagreement reflects unequal sizes with a material interaction (p 0.006), declared type III primary, and reported simple drug effects within strain. |
| twa_unadjusted_pairs | pass | Refused to present the 31 uncorrected LSD hits as findings (66 pairs, about 3.3 expected false positives); reported Tukey-Kramer all-pairs (13 significant) and Tukey within strain, with the unadjusted list kept as a labelled exploratory file from its own script. |
| twa_no_interaction_proof | pass | Declined the requested conclusion: a non-significant interaction (p 0.12, n 3 per cell) does not show equal effects; reported the interaction estimate 8.94 with an approximate CI of -2.9 to 20.8 and drug effects in opposite directions, and noted that equivalence would need a prespecified margin. |

Finding for follow-up (not a scenario failure): with two units per cell the Brown–Forsythe
variance diagnostic is degenerate (p about 4e-89 in `twa_wells_as_units`); see RELEASE_0.15.4.md.
