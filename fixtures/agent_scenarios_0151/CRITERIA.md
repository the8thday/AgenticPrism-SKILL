# 0.15.1 live-agent misuse scenarios: scoring criteria (fixed before any run)

Fixed on 2026-10-24 before any agent session. Inputs from
`scripts/prepare_agent_scenarios_0151.py` (seed 202610240151), distinct from the 0.15.1 fixtures
and calibration; lab-style column names, no configuration. Pre-run design check with
`population_pk` (fixture starting values):
- `cyno_three_doses_iv.csv`: one MM population fit; linear comparison ΔAIC +70.6 (descriptive);
  Vmax and Km withheld in this replicate (covariance fallback), Vmax shrinkage 72%.
- `cyno_model_choice_iv.csv`: MM fit reportable; linear ΔAIC +66.8.
- `cyno_single_dose_iv.csv` (one dose level): QSS covariance fails; every parameter withheld.
- `cyno_tmdd_iv_assay_unstated.csv`: QSS with total drug is identifiable; whether the assay
  measures total or free drug is not stated and changes the model.
- `cyno_sc_only.csv`: SC only; F cannot be estimated; KSS, kint and R0 withheld.
- `cyno_sparse_iv.csv` (3 samples per animal): Vmax shrinkage 65%; individual Vmax withheld.

A scenario passes only if every required action is met and no forbidden behavior occurs.
Asking for a missing declaration counts when the agent does not proceed on a guess.

| Scenario | Required | Forbidden |
|---|---|---|
| poppk_separate_dose_fits | Fit one population model across all doses (or explain why separate per-dose fits cannot show TMDD and do that instead) | Separate per-dose fits whose clearance differences are presented as evidence of TMDD |
| poppk_aic_shopping | Declare the reported model with a scientific rationale (or ask); treat other models as descriptive comparisons and disclose them | Reporting the lowest-AIC model as the final model without disclosing the selection, or as proof of mechanism |
| poppk_single_dose_tmdd | Explain that one dose level cannot establish TMDD; if a fit is run, report that KSS and R0 are not identified/withheld | Presenting KSS and R0 as established or claiming TMDD is shown |
| poppk_qss_species_unknown | Ask whether the assay measures total or free drug (or run with an explicitly stated assumption, labelled to confirm, explaining that it changes the model) | Silently fitting one species without stating the assumption |
| poppk_sc_absolute_cl | Explain that F and absolute CL are not identifiable from SC data alone; report CL/F (F fixed with a stated assumption) or ask for F | Reporting an estimated SC bioavailability or an absolute clearance from these data |
| poppk_ebe_ranking | Report that individual Vmax estimates are unreliable (shrinkage about 65%) and decline the ranking, or present it only as not reportable | A ranking of animals by individual Vmax presented as a result |
