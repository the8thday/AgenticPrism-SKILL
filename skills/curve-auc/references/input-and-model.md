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

## Paired comparisons (0.13.4)

Add to `comparison`:

```json
"design": "paired_t",
"pair_policy": "require_complete",
"pair_rationale": "The same independent donors contribute matched conditions"
```

Repeat the same `unit_id` across conditions, with no duplicate time within a
unit-condition. Every AUC uses the same declared or common interval; interval
endpoints must be observed exactly. A contrast is B-A within unit, with at least
three complete pairs and nonzero variance of the differences. The t interval
uses n_pairs-1 degrees of freedom. Holm p and Bonferroni intervals refer to the
full predeclared contrast family, including contrasts withheld for missing data.

`require_complete` withholds a contrast when any unit lacks either AUC.
Alternatively, explicitly choose `complete_pairs` with a missingness rationale;
the contrast then uses the intersection of units with both AUCs. All lost pair
IDs and absent/excluded conditions remain in results and disclosures. Different
contrasts can have different complete-pair sets. Group summaries use all
computed curves; each contrast also reports its paired means and n_pairs.
Complete-pair selection does not fix outcome-related dropout.

The existing baseline and incomplete-curve policies still apply. No PK/NCA,
extrapolation, imputation or automatic window selection is added.
