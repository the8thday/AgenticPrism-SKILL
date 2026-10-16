# Censored-comparison input contract (0.14.0)

```json
{"analysis_type": "censored_comparison", "method": "censored_lognormal", "input": "il6.csv", "source": "...",
 "design": {"unit": "donor", "rationale": "One serum sample per donor", "outcome": "IL-6", "outcome_unit": "pg/mL",
            "limit_source": "Sample LLOQ x dilution from validation report VR-12"},
 "comparison": {"groups": ["vehicle", "mAb"], "reference": "vehicle", "confidence_level": 0.95}}
```

Columns: `unit_id` (unique), `group`, `value`, `censored` (`none`, `left`,
`right`); optional `exclude` with `exclusion_reason`. For a censored row `value`
is the limit itself. Text values ("<5", "BLQ", "ND") are refused.

| Method | Model | Output |
|---|---|---|
| `censored_lognormal` | ln(value) ~ Normal(mean per group, common SD); likelihood uses the density for detected values, Φ for left-censored and 1 − Φ for right-censored (R `survreg`, `Surv(type = "interval2")`, gaussian) | Geometric means and ratios versus the reference with Wald intervals from the observed information; likelihood-ratio test of equal means; SD on the log scale |
| `peto_peto` | Peto-Peto k-sample test on flipped values M − value (R `survdiff(rho = 1)`); left-censoring only | Chi-square, df, p value |

Withheld: a group with fewer than 3 detected values (parametric); a fit whose
gradient does not vanish or whose information matrix is singular. Flags: censored
fraction above 0.5, fewer than 5 units in a group.
