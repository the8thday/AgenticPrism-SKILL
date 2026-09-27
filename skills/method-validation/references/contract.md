# Method-validation input and numerical contract (0.9.3)

`analysis_type: "method_validation"`, `schema_version: 1`, one `experiment` per run.

## Config

```json
{
  "analysis_type": "method_validation",
  "input": "ap.csv",
  "source": "Validation study VS-012, runs 1-6",
  "experiment": "accuracy_precision",
  "assay": {"platform": "Sandwich ELISA", "matrix": "human serum", "concentration_unit": "ng/mL",
            "lloq": 1.0, "uloq": 100.0, "concentrations_back_calculated": true, "independent_runs": true,
            "rationale": "Each run has its own calibration curve, day and analyst session; QCs independent of standards."},
  "criteria": {"source": "SOP BA-004 v3 citing ICH M10 (2022) 4.2.4.2",
               "accuracy_percent": 20, "accuracy_percent_edge": 25,
               "precision_cv_percent": 20, "precision_cv_percent_edge": 25,
               "total_error_percent": 30, "total_error_percent_edge": 40},
  "statistics": {"confidence_level": 0.9, "cv_denominator": "observed_mean",
                 "tolerance_beta": 0.8, "profile_limit_percent": 30}
}
```

Criteria keys that do not apply to the experiment must be omitted (the tool
refuses them). `_edge` limits apply to LLOQ/ULOQ roles.

| Experiment | Required criteria | Required assay fields | Required columns |
|---|---|---|---|
| accuracy_precision | accuracy, precision, total error (and `_edge`) | — | level, role (lloq/qc/uloq), run_id, nominal |
| dilution_linearity | accuracy_percent, precision_cv_percent | lloq, uloq | series_id, dilution_factor, nominal (neat) |
| parallelism | parallelism_cv_percent | — | sample_id, dilution_factor |
| selectivity | accuracy_percent (high), accuracy_percent_edge (LLOQ), required_pass_fraction, blank_pass_fraction | lloq | source_id, role (blank/lloq/high), nominal |
| specificity | accuracy_percent_edge, required_pass_fraction, blank_pass_fraction | lloq | source_id, role (blank/lloq/uloq), nominal, interferent |
| stability | accuracy_percent | — | condition, level, nominal |

Every file also has `observation_id`, `value`, `status`
(`quantified`/`below_lloq`/`above_uloq`), `exclude` (`true`/`false`) and
`exclusion_reason`. `value` is empty unless quantified: never substitute the
LLOQ or ULOQ number. `nominal` may be empty only for blank rows.

## Numerics

- **Accuracy/precision** per level: one-way random-run ANOVA (classical
  unbalanced n0 = (N − Σnᵢ²/N)/(p−1), as `VCA::anovaVCA`); negative between-run
  estimates are set to zero for the CVs and noted. Bias % = 100(mean/nominal −
  1) over all results. Within-run CV = √MS_W / D; between-run CV = √(σ²_B + MS_W)
  / D with the declared denominator D; total error = |bias %| + between-run CV.
  Per-run within-run accuracy and CV are also reported.
- Supplements: t interval on the mean (df = runs − 1, variance Σnᵢ²/N² σ²_B +
  MS_W/N); MLS (Graybill–Wang) interval for MS_B/n0 + (1 − 1/n0)MS_W; Mee (1984)
  beta-expectation interval mean ± t(ν, (1+β)/2) · √(σ²_B+σ²_W) · √(1 + 1/(p n0 B²)),
  B² = (R+1)/(n0R+1), ν = (R+1)²/((R+1/n0)²/(p−1) + (1−1/n0)/(p n0)).
  The accuracy-profile range interpolates, on log concentration, where the
  interval crosses ±`profile_limit_percent`; it never extrapolates.
- **Dilution linearity**: corrected = value × dilution factor; per factor with
  all rows quantified, mean accuracy vs neat nominal and CV across series. Hook
  effect: any row whose expected in-well concentration exceeds the ULOQ but is
  not reported `above_uloq`. Headline `passes_without_hook`.
- **Parallelism**: per sample with ≥3 quantified dilutions, CV of corrected
  concentrations and the log-log slope; common slope with sample intercepts and
  its interval; per-sample slopes beyond `slope_margin` are flagged.
- **Selectivity/specificity**: replicates averaged per source; blanks pass when
  below LLOQ; spiked rows pass when |accuracy| ≤ limit; pass fraction per role
  with a Clopper–Pearson interval.
- **Stability**: mean accuracy per condition and level; t interval as a
  supplement.

Minimum-design checks (6 runs, 3 replicates, 5 levels, 3 dilutions, 3 series,
10 matrices, 3 aliquots) produce diagnostics, not refusals.
