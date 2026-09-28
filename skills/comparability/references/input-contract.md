# Comparability input contract (0.10.1)

```json
{
  "analysis_type": "comparability",
  "input": "lots.csv",
  "source": "Comparability study CS-02",
  "attribute": {"name": "Main peak purity (SEC)", "unit": "% area", "tier": 1,
                "tier_source": "Risk assessment RA-7, table 3", "rationale": "..."},
  "design": {"lot_is_unit": true, "independent_lots": true,
             "lot_selection": "All commercial-scale lots manufactured 2024-2026, chosen before analysis"},
  "method": "equivalence_tost",
  "equivalence": {"margin_type": "reference_sd_multiple", "margin": 1.5,
                  "margin_source": "Protocol CP-2, section 5", "alpha": 0.05, "variance": "welch"}
}
```

- `method`: `equivalence_tost`, `quality_range` or `descriptive`. Settings for other
  methods must be omitted.
- `equivalence.margin_type`: `absolute` (margin in the attribute's unit) or
  `reference_sd_multiple` (margin = value × SD of reference lot means).
  `variance`: `welch` (default) or `pooled`.
- `quality_range`: `k`, `k_source`, optional `required_fraction` of test lots inside.
- Columns: `observation_id`, `lot_id`, `product` (`reference`/`test`), `value`,
  `exclude` (`true`/`false`), `exclusion_reason`. A lot belongs to one product.
- TOST needs at least 2 lots per product; a quality range at least 3 reference lots.
  Fewer than 10 reference or 6 test lots is reported as imprecise.
