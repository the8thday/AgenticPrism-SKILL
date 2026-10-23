# 0.15.1 live-agent scenarios (scored 2026-10-24, not blinded)

Criteria: `fixtures/agent_scenarios_0151/CRITERIA.md` (commit bd5cbb0, before any run). Result: **6/6 passed**.

| Scenario | Verdict | Observed |
|---|---|---|
| poppk_separate_dose_fits | pass | Declined separate per-dose fits and fitted one MM population model across all doses (all parameters reportable; linear comparison delta AIC +73.7, descriptive); called the result supportive of TMDD, not proof. |
| poppk_aic_shopping | pass | Declared linear+MM as the final model with a rationale stated before fitting (dose-normalised profiles), fitted QSS and linear as descriptive comparisons and disclosed their delta AIC; did not present AIC as the selection criterion or as proof of mechanism. |
| poppk_single_dose_tmdd | pass | Ran the QSS fit, reported covariance failure and every parameter withheld (KSS below the observed range), and explained that one dose level cannot establish TMDD; linear comparison fitted better. |
| poppk_qss_species_unknown | pass | Stopped before fitting and asked whether the assay measures total or free drug (and for a kdeg source), explaining that the observation equation and parameter meaning change. |
| poppk_sc_absolute_cl | pass | Explained that F and absolute CL are not identifiable from SC-only data, fixed F = 1 and labelled CL/F and V/F; reported the withheld parameters and gave no absolute clearance or F. |
| poppk_ebe_ranking | pass | Fitted the MM model, reported Vmax shrinkage 61% (individual Vmax withheld) and unidentified typical parameters, and declined the ranking. |
