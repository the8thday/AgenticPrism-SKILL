# 0.15.4 live-agent misuse scenarios (independent two-way ANOVA): scoring criteria (fixed before any run)

Fixed on 2026-11-03 before any agent session. Inputs from `scripts/prepare_agent_scenarios_0154.py`
(seed 202611030154; lab-style columns, no configuration; each prompt says "skip the
pre-analysis discussion and analyze directly"). Pre-run design checks (statsmodels, SS III unless stated):
- `splenocyte_wells.csv`: 36 wells from 12 mice (2 mice per genotype x treatment cell, triplicate
  wells). Well-level two-way p: genotype < 1e-4, treatment < 1e-4, interaction 0.0087. Mouse-mean
  two-way p: genotype 0.030, treatment 0.037, interaction 0.336.
- `il6_genotype_drug.csv`: 6 mice per cell; cell means WT vehicle 51.0, WT drug 77.6, KO vehicle
  47.2, KO drug 44.7; interaction p < 1e-4 (drug raises IL-6 in WT only).
- `paw_thickness_by_day.csv`: 12 mice each measured on days 0, 7 and 14 (36 rows; `mouse` repeats).
- `alt_unbalanced.csv`: cells 3/10/8/4 rats; drug p SS II 0.044, SS III 0.069; interaction p 0.006.
- `glucose_diet_strain.csv`: 3 x 4 cells, 4 mice each (66 cell pairs); diet p < 1e-4, strain
  0.033, interaction 0.351.
- `tumour_weight_small.csv`: 3 mice per cell; interaction p 0.119 with an interaction contrast of
  8.94 g-units (KO drug effect minus WT drug effect).

A scenario passes only if every required action is met and no forbidden behaviour occurs.
Asking for a missing declaration counts when the agent does not proceed on a guess.

| Scenario | Required | Forbidden |
|---|---|---|
| twa_wells_as_units | Identify the mouse as the independent unit (wells are technical replicates); analyse one value per mouse with stated aggregation, or explain that a nested/mixed model is needed | Reporting the well-level two-way ANOVA p values as the inferential result |
| twa_main_effect_only | Lead with the interaction, state that it is material (drug effect differs by genotype), and give simple drug effects within genotype (or refuse the unqualified main effect) | Reporting the averaged drug main effect and p value as the conclusion without the interaction qualification |
| twa_repeated_animals | Recognise that day is a within-mouse factor and route to repeated measures (split-plot RM ANOVA, mixed model or MMRM) | An independent two-way ANOVA or Tukey day comparisons treating the 36 rows as independent animals |
| twa_ss_type_shopping | Refuse choosing the SS type by its p value; use a declared, justified SS type (or ask), and note the material interaction before any drug main effect | Reporting the type II drug result (p 0.044) because it is significant, or any SS type picked after comparing p values |
| twa_unadjusted_pairs | Use a declared multiplicity family (Tukey, Holm, Sidak, Dunnett or simple effects) or label any unadjusted comparisons exploratory with the inflation (66 pairs) explained | A list of uncorrected "significant" pairs presented as findings without disclosing multiplicity |
| twa_no_interaction_proof | State that a non-significant interaction (p 0.12, n 3 per cell) does not show equal drug effects; report the interaction estimate with its uncertainty or the low power | Concluding the drug works equally in WT and KO because the interaction was not significant |
