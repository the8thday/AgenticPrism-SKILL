# quantal_dose_response input contract (0.15.3)

Grouped CSV: `dose`, `n`, `responders` (optional `preparation`, `exclude`,
`exclusion_reason`). Individual CSV (`data_format: individual`): `dose`, `response`
(0/1), optional unique `subject`; rows are aggregated by preparation and dose. Dose 0
rows are controls (used only by `natural_response: estimate`, refused for `none`
when they contain responders).

```json
{"schema_version": 1, "analysis_type": "quantal_dose_response", "input": "mortality.csv", "source": "Study X",
 "design": {"endpoint": "mortality at 14 days", "response_definition": "dead or moribund", "subject_unit": "animal",
            "dose_unit": "mg/kg", "independence_rationale": "singly housed, randomised"},
 "data_format": "grouped",
 "model": {"link": "probit", "metameter": "log10_dose", "rationale": "Protocol"},
 "natural_response": {"method": "none", "rationale": "Vehicle controls 0/10"},
 "heterogeneity": {"policy": "finney_if_significant", "alpha": 0.05, "rationale": "Finney (1952)"},
 "effective_doses": {"percent": [10, 50, 90], "confidence_level": 0.95},
 "relative_potency": {"reference": "reference", "parallelism": {"method": "lr_test", "alpha": 0.05}, "rationale": "Same moiety"}}
```

Model: P = c + (1 − c) F(a + b x); ML by Fisher scoring with expected information.
EDp: x_p = (F⁻¹(p) − a)/b with Fieller limits (g = t² V(b)/b²; normal quantile, or t
on k − parameters df when the covariance is scaled by χ²/df). Withheld when g ≥ 1,
outside the tested range, slope ≤ 0, no partial-response group, singular
information or a non-positive-definite covariance. `parallelism` may instead be
`{"method": "slope_ratio_equivalence", "margin": [0.8, 1.25], "confidence_level": 0.9, "source": "..."}`.
Relative potency ρ = D_reference/D_test at equal effect (> 1: test more potent).
`relative_potency` is optional.
