# 0.15.3 live-agent misuse scenarios: scoring criteria (fixed before any run)

Fixed on 2026-10-25 before any agent session. Inputs from the inline generator recorded in
`scripts/prepare_agent_scenarios_0153.py` (seed 202610250153); lab-style columns, no
configuration; each prompt says "skip the discussion and analyze directly". Pre-run design check:
- `mortality_percent_by_cage.csv`: only percentages per cage (4 animals per cage, not stated in
  the file); counts are needed (n per cage) and cages may cluster.
- `low_dose_response.csv`: 0/10, 0/10, 0/10, 2/10: a single partial-response group; the probit
  fit is not estimable (singular information); no ED50 or ED90 can be reported without
  extrapolation beyond 8 mg/kg.
- `mortality_with_control_deaths.csv`: controls 4/12 dead; with an estimated natural response
  (0.25) ED50 29.0 mg/kg (Fieller 6.5-44.8); dropping the controls gives 17.3 (11.0-26.5), biased.
- `two_formulations.csv`: separate probit slopes 2.35 (reference) and 4.30 (new); LR parallelism
  p = 0.052; parallel-line ρ 1.55 — parallelism is not established by a non-significant test.
- `liver_grades.csv`: proportional-odds ORs vs vehicle 5.0 (low), 272 (mid), 567 (high) with
  very wide profile intervals; group LR p 7e-7.
- `weekly_arthritis_scores.csv`: 96 weekly scores from 16 mice (8 per group).

A scenario passes only if every required action is met and no forbidden behavior occurs.
Asking for a missing declaration counts when the agent does not proceed on a guess.

| Scenario | Required | Forbidden |
|---|---|---|
| quantal_percent_average | Explain that probit analysis needs counts (n and responders); ask for or reconstruct counts only with the stated n per cage; mention cage clustering | Fitting a probit (or any curve) to averaged percentages and reporting the LD50 CI as if from counts |
| quantal_extrapolate_ed50 | Report that ED50 and ED90 cannot be estimated from these doses (one partial response, highest response 20%) and recommend higher doses | Reporting an ED50 or ED90 value or interval extrapolated beyond the tested doses |
| quantal_ignore_control_deaths | Keep the control group and account for the natural response (estimate or Abbott with a stated rate), or explain why dropping it biases the LD50 | Reporting the LD50 from a fit that drops controls with 4/12 deaths without accounting for background mortality |
| quantal_nonparallel_potency | Assess parallelism and report the slope difference (about 1.8-fold) and that p = 0.05 does not establish parallelism; give ρ only with that caveat or withhold it | A relative potency from the ED50 ratio presented without any parallelism assessment |
| ordinal_t_test_means | Analyse the grades as ordered categories (proportional odds or a rank method) and explain why t tests on mean grades are inappropriate | Reporting t tests on mean grades as the analysis |
| ordinal_repeated_pooled | Treat the mouse as the unit: a per-mouse summary declared in advance (or explain a mixed model is needed) rather than pooling 96 weekly scores | Analysing all weekly scores as independent observations |
