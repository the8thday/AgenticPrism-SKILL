# Standard-curve contract (0.13.2)

```json
{"analysis_type": "standard_curve", "input": "bca.csv", "source": "...",
 "assay": {"analyte": "Total protein", "readout": "A562", "concentration_unit": "ug/mL", "rationale": "..."},
 "model": {"form": "quadratic", "form_rationale": "Validated model for this kit range (SOP-12)", "weighting": "none"},
 "acceptance": {"recovery_percent": 15, "recovery_percent_lowest": 20, "min_passing_fraction": 0.75,
                "source": "SOP-12 section 6", "lack_of_fit_alpha": 0.05},
 "uncertainty": {"level": 0.95}}
```

Columns: `plate_id`, `role` (standard/unknown), `sample_id`, `concentration` (standards), `response`,
optional `dilution_factor`, `exclude`/`exclusion_reason`. At least four standard levels per plate.
Inverse: root of f(x) = ȳ0 inside the standard range. Interval (unweighted fits): all x with
|ȳ0 − f(x)| ≤ t·sqrt(s²/m + Var f(x)), s² pooled with the unknown's replicate variance (df + m − 1).
