# ordinal_regression input contract (0.15.3)

CSV: `subject` (unique), `score` (one of the declared categories), the group column,
declared numeric covariates, optional `exclude` + `exclusion_reason`.

```json
{"schema_version": 1, "analysis_type": "ordinal_regression", "input": "histopathology.csv", "source": "Study X",
 "design": {"outcome": "hypertrophy grade", "categories": ["0", "1", "2", "3", "4"], "subject_unit": "rat",
            "independence_rationale": "one terminal score per rat"},
 "group": {"column": "group", "reference": "vehicle"},
 "covariates": [{"name": "baseline_weight_g", "rationale": "Protocol pre-treatment covariate"}],
 "proportional_odds_check": {"method": "nominal_lr", "alpha": 0.05},
 "confidence_level": 0.95}
```

Model: logit P(Y ≤ j) = ζ_j − η (MASS::polr parameterisation; covariates centred);
ML with analytic gradient, numerical Hessian for Wald SEs; profile-likelihood
intervals; group LR test; nominal LR test with category-specific group effects
(`unavailable` when that model gives crossing cumulative probabilities). |log OR|
> 10 is withheld as separation.
