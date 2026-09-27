# Tumor-growth contract (0.8.0)

## CSV

Required: `animal_id,arm,day,volume`. One row per animal and day; `volume` ≥ 0 in
`study.volume_unit` (mm3 or cm3); `day` numeric in `study.time_unit` (day or
week). Each animal belongs to one arm and has a `baseline_day` measurement and at
least two measurement days; each arm has ≥ 3 animals. Optional
`exclude=true|false` with `exclusion_reason`.

## JSON

```json
{
  "analysis_type": "tumor_growth",
  "input": "tumors.csv",
  "source": "Study identifier and caliper data preparation",
  "study": {"volume_unit": "mm3", "time_unit": "day", "time_origin": "first dose (day 0)",
            "removal_rule": "Euthanasia after the first volume above 2000 mm3",
            "dropout_rationale": "Removal depends only on the last observed volume (MAR)",
            "rationale": "Randomized by baseline volume; comparisons vs vehicle and the day-21 readout prespecified"},
  "analysis": {"arms": ["vehicle", "mAb-low", "mAb-high"], "control_arm": "vehicle",
               "baseline_day": 0, "analysis_day": 21, "log_offset": 1.0, "confidence_level": 0.95}
}
```

## Growth model

y = log(V + offset). Fixed effects: one intercept and one slope per arm (cell
means). Random effects: animal intercept and slope with an unstructured 2 × 2
covariance G; independent residuals with variance σ². REML by direct optimization
of the deviance over the Cholesky factor of G (log diagonal) and log σ², from
several starts, then polished. Satterthwaite t and F use a Richardson-
extrapolated Hessian of the REML deviance (as lmerTest/numDeriv) and gradients of
c'C(u)c; they are withheld when the slope variance or the random-effect
correlation is at a boundary. Outputs: per-arm slope (per time unit) with
pointwise limits and doubling time ln2/slope (limits from the slope limits when
the lower limit is positive); equal-slopes F (k − 1 df); slope differences vs
control (Holm p, Bonferroni limits); model T/C = exp[(a_T − a_C) + (b_T − b_C) ×
analysis_day] with Bonferroni limits.

## Observed TGI% and T/C%

Animals with baseline and analysis-day measurements only. ΔV = V(analysis) −
V(baseline); TGI% = 100 (1 − mean ΔT / mean ΔC); T/C% = 100 mean V_T / mean V_C.
Limits by Fieller's theorem for a ratio of independent means with
Welch–Satterthwaite df evaluated at the estimate, at the Bonferroni-adjusted level
for the arm-vs-control family. When the control mean is not distinguishable from
zero the Fieller set is unbounded and limits are reported as missing. TGI is not
defined when the control did not grow on average.

Assumptions: independent animals; exponential growth (linear log volume) over the
analysed window; Gaussian random effects and residuals on the log scale; MAR
dropout for the model. Not supported: Gompertz/logistic growth curves, regrowth
after regression modelled as change points, cage effects, AR(1) residuals.
