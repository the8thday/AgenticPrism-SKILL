# Curve AUC contract (0.13.2)

```json
{"analysis_type": "curve_auc", "input": "killing.csv", "source": "...",
 "design": {"unit": "donor", "unit_rationale": "...", "outcome": "Killing", "outcome_unit": "%", "time_unit": "hour"},
 "auc": {"interval": [0, 72], "baseline": "first_value", "baseline_value": null,
         "incomplete_policy": "withhold_unit", "policy_rationale": "Wells lost to contamination, unrelated to killing"},
 "comparison": {"groups": ["isotype", "TCE"], "post_hoc": "vs_control", "control_group": "isotype", "confidence_level": 0.95}}
```

Columns: `unit_id`, `group`, `time`, `value`, optional `exclude`/`exclusion_reason`. Each unit must have
observations exactly at both interval ends (no extrapolation); interior gaps are joined linearly.
`common_interval` shortens every unit to the earliest last time. Output: per-unit AUC and status, group
mean AUC with t intervals, Welch contrasts (Holm p, Bonferroni CI).
