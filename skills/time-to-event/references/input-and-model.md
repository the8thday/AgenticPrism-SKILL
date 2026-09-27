# Time-to-event contract (0.8.0)

## CSV

Required: `subject_id,arm,time,event`. One row per subject. `time` is finite and
≥ 0 in `study.time_unit` from `study.time_origin`; `event` is 1 (event observed)
or 0 (censored at `time`). Optional numeric covariate columns named in
`comparison.covariates` must be present for every included subject (no silent
case deletion). Optional `exclude=true|false` with `exclusion_reason`. Arms in
the file must equal the declared arms; every arm needs ≥ 2 included subjects.

## JSON

```json
{
  "analysis_type": "time_to_event",
  "input": "survival.csv",
  "source": "Study identifier and data preparation",
  "study": {"endpoint": "Humane endpoint: tumor volume >= 2000 mm3 (event)",
            "time_origin": "randomization", "time_unit": "day",
            "censoring_rationale": "Alive without endpoint at day-60 study end",
            "rationale": "Ten mice randomized per arm; comparisons vs vehicle declared in the protocol"},
  "comparison": {"arms": ["vehicle", "mAb-low", "mAb-high"], "control_arm": "vehicle",
                 "post_hoc": "arm_vs_control", "logrank_inference": "permutation",
                 "permutation_draws": 9999, "seed": 20260930, "covariates": [],
                 "landmarks": [28, 42], "conf_type": "log-log", "confidence_level": 0.95}
}
```

`time_unit`: hour, day, week, month or year. `post_hoc`: `arm_vs_control`,
`all_pairs` or `none`. `permutation_draws`: 999–100000. `conf_type`: `log` or
`log-log`. `report.show_confidence_bands` toggles KM bands (display only).

## Methods

- **Kaplan–Meier:** S(t) = Π (1 − d_j/n_j). Standard error of the cumulative
  hazard (Greenwood): sqrt(Σ d_j / (n_j (n_j − d_j))). `log` limits:
  exp(log S ± z·se), upper capped at 1. `log-log` limits:
  exp(−exp(log(−log S) ∓ z·se/log S)). Where S = 1 the limits are 1; where S = 0
  they are undefined. Identical to R `survfit` (survival 3.8.6).
- **Median and landmark survival:** the median is the first time S ≤ 0.5 (if S
  equals 0.5 exactly over an interval, the midpoint to the next drop, as in R
  `quantile.survfit`); its limits apply the same rule to the limit curves. "Not
  reached" when the curve (or a limit) never reaches 0.5. Landmark survival is the
  step-function value; beyond the last follow-up it is withheld unless the curve
  has already reached 0.
- **Log-rank:** O − E per arm with the hypergeometric covariance; chi-square on
  k − 1 df (R `survdiff`, rho = 0). Permutation option: arm labels permuted
  (Monte Carlo, fixed seed), p = (1 + #{χ²* ≥ χ²}) / (B + 1). Pairwise tests use
  the two arms' subjects only; Holm across the declared family.
- **Cox model:** Efron ties, Newton–Raphson from β = 0 with step halving;
  covariates centred. Wald, likelihood-ratio and score tests; HR = exp(β) with
  Wald limits (pairwise HRs: Bonferroni over the family). Reference arm: the
  control arm, else the first declared arm. Withheld when an arm has no events or
  the fit does not converge (monotone likelihood). Identical to R `coxph`.
- **Proportional-hazards test:** score test for X·g(t) added to the model at
  (β̂, 0) with g = 1 − KM(t−) centred on event times; per term (arm indicators
  together) and global. Identical to R `cox.zph(transform = "km")`.

Assumptions: independent subjects, non-informative censoring, common time zero,
and for Cox proportional hazards. Not supported: competing risks, recurrent
events, time-varying covariates, stratified or frailty models, interval
censoring, left truncation.
