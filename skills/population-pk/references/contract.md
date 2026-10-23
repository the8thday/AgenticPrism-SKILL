# population_pk input and model contract (0.15.1)

## Input CSV

Long table with dose rows and observation rows.

| Column | Meaning |
|---|---|
| `subject` | Animal or subject ID |
| `time` | Time in `design.time_unit` (≥ 0), one clock per subject |
| `dose` | Dose amount on dose rows (empty on observations); multiple doses allowed |
| `concentration` | Observation (empty on dose rows; empty when `blq` is true) |
| `infusion_duration` | IV infusion only: duration on every dose row |
| `group`, `blq`, `exclude`, `exclusion_reason` | Optional; every excluded row needs a reason |

Observations at or before a subject's first dose (predose samples) must be
excluded with a reason. Doses are into the central compartment (IV) or a depot (SC).

## Configuration

```json
{"schema_version": 1, "analysis_type": "population_pk", "input": "pk.csv", "source": "Study X report",
 "design": {"analyte": "mAb-X total", "matrix": "serum", "species": "cynomolgus monkey", "subject_unit": "animal",
            "route": "iv_bolus", "dose_unit": "mg/kg", "concentration_unit": "ug/mL", "time_unit": "h"},
 "model": {"structural": "linear_mm_2cmt", "rationale": "Expected target-mediated elimination (declared before fitting)",
           "initial_estimates": {"cl": 3e-4, "v1": 0.05, "q": 5e-4, "v2": 0.04, "vmax": 1e-3, "km": 2},
           "initial_estimates_source": "Typical cynomolgus IgG values; NCA",
           "random_effects": ["cl", "v1", "vmax"], "residual": "proportional", "residual_initial": {"proportional": 0.2}},
 "estimation": {"method": "focei", "seed": 15101, "max_outer_iterations": 5000},
 "blq": {"rule": "drop", "rationale": "SAP: BLQ excluded (M1)"},
 "comparison_models": [{"structural": "linear_2cmt", "additional_initial_estimates": {}}]}
```

- Parameters: `linear_2cmt` cl, v1, q, v2; `linear_mm_2cmt` adds vmax, km (parallel
  linear and Michaelis-Menten elimination); `qss_tmdd_2cmt` adds kss, kint, kdeg, r0
  (Gibiansky quasi-steady-state, target turnover from kdeg and baseline r0), with
  `observed_species` total or free. SC adds ka and needs `bioavailability`
  `{value, source}` (fixed).
- `fixed` (with `fixed_source`) fixes parameters; `initial_estimates` must cover
  every other parameter. `random_effect_initial_variance` defaults to 0.1.
- Residual `proportional` or `combined` (`residual_initial` has `additive` too).
- `comparison_models` are optional and descriptive; each lists starting values for
  its extra parameters.

Units are composed, never converted: volumes `dose_unit/concentration_unit`,
clearances per time unit, Vmax `dose_unit/time`, Km, KSS and R0 in the
concentration unit.

## Output (`results.json` → `primary`)

- `parameters[p]`: estimate, log estimate and SE, natural-scale RSE, 95% CI,
  unit, `status` (`reportable`, `withheld`, `fixed`), `reasons`, `flags`.
- `iiv[p]`: omega variance, CV%, shrinkage (nlmixr2 Shrink(SD)%),
  `individual_estimates` (`withheld` above 30% shrinkage).
- `residual`, `covariance` (method, `sandwich`/`fallback`/`failed`),
  `high_correlation_pairs`, `nonlinearity_check`, `fit` (objective, AIC, BIC, n),
  `comparisons` (ΔAIC, descriptive), `observed_concentration_range`, BLQ dropped,
  excluded rows, nlmixr2 versions and warnings, the exact model code.
- `fit_table` (DV, PRED, IPRED, IWRES, individual parameters) drives the
  goodness-of-fit and individual-fit figures.

Gates: failed covariance withholds every interval; |correlation| > 0.95 withholds
both parameters; Km or KSS outside the observed concentration range withholds the
nonlinear parameters; natural-scale RSE > 100% (or an overflowing log-scale SE)
withholds a parameter as not identified, and any withheld nonlinear parameter
withholds all of them; 50–100% RSE is a flag.
