# pk_nca input and output contract (0.15.0)

## Input CSV

One long table, one row per subject and sample time, single dose at time 0.

| Column | Meaning |
|---|---|
| `subject` | Subject or animal ID (the unit of analysis) |
| `time` | Time after the dose (≥ 0) in `design.time_unit` |
| `concentration` | Measured value; empty when `blq` is true |
| `dose` | Dose in `design.dose_unit`, constant within subject |
| `group` | Optional dose group label for summaries (otherwise grouped by dose) |
| `blq` | Optional true/false; BLQ rows must leave `concentration` empty |
| `exclude`, `exclusion_reason` | Optional declared exclusions; every excluded row needs a reason |
| `infusion_duration` | Required for `intravascular_infusion`, constant within subject |
| `lambda_z_point` | Required for `lambda_z.selection: declared_points` |

Other columns are refused. Duplicate subject-time rows, negative values and doses
that vary within a subject are refused.

Time zero: for infusion and extravascular dosing declare `time_zero.rule`
`observed` (a predose time-0 row exists; BLQ is fine) or `insert_zero` (a 0 at
time 0 is added and logged). For IV bolus a time-0 row, if present, must be the
predose sample (BLQ or 0); C0 is always back-extrapolated.

## Configuration

```json
{"schema_version": 1, "analysis_type": "pk_nca", "input": "pk.csv", "source": "Study X bioanalytical report v2",
 "design": {"analyte": "mAb-X total", "matrix": "serum", "species": "cynomolgus monkey", "subject_unit": "animal",
            "route": "intravascular_bolus", "dose_unit": "mg/kg", "concentration_unit": "ug/mL",
            "time_unit": "h", "time_basis": "nominal"},
 "blq": {"first": "keep_as_zero", "middle": "drop", "last": "drop", "rationale": "SAP section 9.2"},
 "nca": {"auc_method": "lin_up_log_down", "auc_method_rationale": "SAP",
         "lambda_z": {"selection": "best_adjusted_r2", "min_points": 3, "adj_r2_factor": 0.0001}},
 "reportability": {"min_adj_r2": 0.9, "min_span_ratio": 2, "max_aucinf_extrap_percent": 20, "source": "SAP section 9.4"},
 "dose_proportionality": {"parameters": ["cmax", "aucinf"], "acceptance": [0.8, 1.25], "confidence_level": 0.9,
                          "source": "Smith et al. 2000", "rationale": "Exposure scaling across 1-10 mg/kg"}}
```

`blq` is required when any post-dose row is BLQ. `reportability` may instead be
`{"descriptive_only": true, "rationale": "..."}`: terminal parameters are then
`not_assessed` and left out of summaries (and `aucinf` dose proportionality is
refused). `dose_proportionality` is optional.

## Calculation (PKNCA)

All values come from PKNCA's `pk.nca` with the declared options
(`auc.method`, `conc.blq`, `min.hl.points`, `adj.r.squared.factor`,
`allow.tmax.in.half.life = FALSE`). IV bolus runs two PKNCA passes: pass 1 on the
observed samples gives C0 (PKNCA order c0, logslope, c1, cmin, set0), Cmax and Tmax;
pass 2 inserts C0 at time 0 for AUC, AUMC, lambda-z, CL, Vz, MRT and Vss. When both
passes choose the same lambda-z, pass-2 AUCinf must equal pass-1 `aucivinf.obs`.

## Output (`results.json` → `primary`)

- `subjects[]`: `values` (c0, cmax, tmax, tlast, clast, auclast, aumclast,
  lambda_z and its points/times/R²/adjusted R²/span ratio, half_life, aucinf,
  aucinf_extrap_percent, aumcinf, cl, vz, mrt, vss), `status` per parameter
  (`reportable`, `withheld`, `not_estimated`, `not_assessed`, `not_applicable`),
  `terminal_reasons`, `aucinf_reasons`, PKNCA exclusion text and `diagnostics`.
- `parameter_labels` (CL/F, Vz/F for extravascular; MRT_iv for intravascular) and
  composed `units` (no conversion).
- `summaries[]` per group: n, n_not_reportable, geometric mean, geometric CV%,
  mean, SD (Tmax/Tlast: median, min, max), reportable values only.
- `dose_proportionality.parameters[k]`: slope b, SE, df, CI, dose ratio r,
  acceptance interval for b, conclusion (`proportional_within_acceptance`,
  `disproportionate`, `inconclusive`) and r^(b−1) with its CI; or `withheld`.
- `pknca` (version, options, warnings), `preprocessing_log`, `excluded_samples`,
  `bolus_two_pass_checks`; `mapped_input` is the exact table passed to PKNCA.

Withholding: lambda-z not estimated, a declared criterion failed, or fewer declared
points than `min_points` withholds t½ and every AUCinf-dependent parameter;
%AUCinf extrapolated above its limit withholds the AUCinf-dependent ones. Cmax,
Tmax and AUClast stay reportable.
