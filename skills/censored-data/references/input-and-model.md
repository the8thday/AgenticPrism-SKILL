# Censored-comparison input contract (0.14.0; paired designs 0.14.1)

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

## Paired designs (0.14.1)

```json
{"analysis_type": "censored_comparison", "method": "censored_paired_difference", "input": "il6_pre_post.csv", "source": "...",
 "design": {"unit": "donor", "rationale": "Pre- and post-dose serum from each donor", "outcome": "IL-6", "outcome_unit": "pg/mL",
            "limit_source": "Sample LLOQ x dilution from validation report VR-12"},
 "comparison": {"groups": ["pre", "post"], "reference": "pre"},
 "paired": {"scale": "log", "scale_rationale": "Multiplicative effects, right-skewed values (declared before results)",
            "pairing_rationale": "Same donor before and after dosing", "inference": "small_sample_t"}}
```

Same columns; `unit_id` appears once in each of exactly two groups. Incomplete
pairs are refused; excluding a value requires excluding both rows of that unit.
Left censoring only. Difference = other group − reference.

| Method | Model | Output |
|---|---|---|
| `censored_paired_difference` | Each pair gives an interval for the difference; a censored value spans [0, limit] (`scale: original`, as NADA2 `cen_paired`) or (−∞, ln limit] (`scale: log`). Normal ML with exact, one-sided and interval-censored differences (R `survreg`, `interval2`, gaussian). | Mean difference (or ratio of paired geometric means) with interval and p value. `inference: small_sample_t` (default): ML SE × √(k/(k−1)) with t(k−1) over k informative pairs, equal to the paired t test without censoring. `wald_z`: NADA2's z. |
| `censored_sign_test` | NADA2 `cen_signtest`: ties are pairs whose direction cannot be determined | n+, n−, ties, exact binomial p (primary), Fong tie-corrected p |

Withheld: fewer than 3 pairs with both values detected; non-converged likelihood
or singular information; a sign test with no untied pair. Flags: censored
fraction above 0.5 in a group; fewer than 10 pairs; on the log scale, pairs with
both values censored (uninformative, counted).
