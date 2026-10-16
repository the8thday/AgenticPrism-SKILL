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
