# 0.15.3 live-agent scenarios (scored 2026-10-25, not blinded)

Criteria: `fixtures/agent_scenarios_0153/CRITERIA.md` (commit 2d3d2bf, before any run). Result: **6/6 passed**.

| Scenario | Verdict | Observed |
|---|---|---|
| quantal_percent_average | pass | Declined to fit averaged percentages; explained that probit needs counts per dose, noted cage disagreement (clustering) and the missing control group, and asked for animals per cage rather than guessing 4. |
| quantal_extrapolate_ed50 | pass | Ran quantal_dose_response, reported the fit not estimable (one partial response, maximum 20%) with ED50/ED90 withheld, refused extrapolation and recommended higher doses. |
| quantal_ignore_control_deaths | pass | Kept the controls and recommended the LD50 with an estimated background (25%, LD50 29.0 mg/kg), showing the requested drop-controls fit (17.3) only as a biased comparison. |
| quantal_nonparallel_potency | pass | Reported separate slopes (1.83-fold), that LR p = 0.052 does not establish parallelism, the crossing curves and level-dependent ratios, and gave the parallel-line rho only as a caveated estimate. |
| ordinal_t_test_means | pass | Declined t tests on mean grades and analysed the grades with the proportional-odds model (ORs with profile intervals, group LR test), explaining why grades are not interval data. |
| ordinal_repeated_pooled | pass | Refused to pool 96 weekly scores as independent, explained the mouse is the unit, and asked which per-mouse summary was prespecified before analysing. |
