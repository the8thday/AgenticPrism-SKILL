# Panel-comparison input contract (0.14.0)

```json
{"analysis_type": "panel_comparison", "input": "panel.csv", "source": "...",
 "design": {"unit": "mouse", "rationale": "One serum sample per mouse", "pairing": "independent",
            "scale": "log", "scale_rationale": "Cytokine concentrations are approximately log-normal"},
 "comparison": {"group_a": "vehicle", "group_b": "mAb", "confidence_level": 0.95, "alpha": 0.05},
 "families": {"cytokines": ["IL2", "IL6", "IFNG"], "chemokines": ["CXCL9", "CXCL10"]},
 "multiplicity": {"procedure": "holm", "rationale": "Claims will be made about individual cytokines (plan section 4)"}}
```

Columns: `unit_id`, `group`, `outcome`, `value`, `unit` (one unit per outcome);
optional `exclude` with `exclusion_reason`. Every measured outcome must be in
exactly one family and every declared outcome must be measured.

| Step | Method |
|---|---|
| Per outcome | Welch t (independent) or paired t on complete pairs (incomplete pairs withhold the outcome); log scale gives ratios of geometric means |
| Holm | Holm adjusted p within the family; Bonferroni simultaneous intervals at 1 − (1 − level)/m |
| Benjamini-Hochberg | Step-up adjusted p within the family (R `p.adjust("BH")`); intervals unadjusted |

Withheld outcomes (fewer than 3 units or pairs, zero variance, incomplete pairs)
are listed and excluded from the family size.

## Censored outcomes (0.14.1)

Add `"censoring": {"method": "censored_lognormal", "limit_source": "Plate LLOQ per analyte, validation report VR-7"}`
with `design.scale: log`, and a `censored` column (`none`, `left`, `right`; the limit as
`value`). A `censored` column without the block is refused.

| Design | Outcome with any censored value |
|---|---|
| Independent | Tobit: ln value ~ Normal(group mean, common SD), left/right censoring (R `survreg`, `interval2`); ratio with SE × √(N/(N−2)) and N − 2 df (Student's pooled t without censoring) |
| Paired | Interval-censored paired ln differences (as `censored-data` paired); SE × √(k/(k−1)), k − 1 df over k informative pairs (paired t without censoring); left censoring only; pairs with both values censored are uninformative |

Withheld: fewer than 3 detected values in a group (independent), fewer than 3 pairs
with both values detected (paired), non-converged likelihood. Outcomes without
censored values use Welch or paired t as above. The family's Holm/BH adjustment and
simultaneous intervals use each outcome's own df.
