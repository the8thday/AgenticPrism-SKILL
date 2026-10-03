# Sample-size input and method contract (0.13.1)

`analysis_type: "sample_size"`, `schema_version: 1`, no `input` file.

```json
{
  "analysis_type": "sample_size",
  "source": "Planning memo PM-07",
  "purpose": "Two-arm xenograft study, tumor volume day 21",
  "design": "two_sample_t",
  "solve_for": "n",
  "alpha": 0.05,
  "sidedness": "two_sided",
  "target_power": 0.8,
  "assumptions": {
    "difference": 150, "difference_source": "Smallest difference that changes the go/no-go decision",
    "sd": 180, "sd_source": "Pooled vehicle-arm SD, studies X-12 to X-15",
    "allocation_ratio": 1
  },
  "sensitivity": {"assumption": "sd", "values": [140, 180, 220, 260]}
}
```

`solve_for: "power"` takes an integer `n` (per group; pairs for `paired_t`; total
subjects for `logrank`) instead of `target_power`.

| Design | Assumptions (each with `_source`) | Optional | Method |
|---|---|---|---|
| `two_sample_t` | difference, sd | allocation_ratio (n2/n1) | exact noncentral t, pooled SD |
| `paired_t` | difference, sd_difference | — | exact noncentral t |
| `one_way_anova` | group_means (≥3), sd | — | exact noncentral F, λ = k n f² |
| `two_proportions` | p1, p2; `method` arcsine or normal_approximation | allocation_ratio (arcsine only) | Cohen h (pwr.2p2n.test) or pooled normal (power.prop.test) |
| `logrank` | hazard_ratio, event_probability | allocation_ratio | Schoenfeld events |
| `tost_two_means` | difference (true), sd, margin [low, high] | allocation_ratio | exact TOST, alpha per one-sided test |

`n2 = ceil(allocation_ratio × n1)`. The solved n is the smallest integer whose power
reaches the target. Results add Cohen's effect size, a 40-point `power_curve` and
`sensitivity` rows; an unattainable target below n = 100000 is a failing item.
