# Nested and ANCOVA contracts (0.13.2)

## nested_comparison

```json
{"analysis_type": "nested_comparison", "input": "wells.csv", "source": "...",
 "design": {"unit": "mouse", "unit_rationale": "Each mouse randomized and dosed independently",
            "replicate": "well", "replicate_rationale": "Organoid wells derived from that mouse's tumor",
            "outcome": "Viability", "outcome_unit": "%"},
 "comparison": {"groups": ["vehicle", "mAb"], "post_hoc": "vs_control", "control_group": "vehicle", "confidence_level": 0.95}}
```

Columns: `group`, `unit_id`, `value`, optional `well_id`, `exclude`, `exclusion_reason`. Each unit in one group,
at least two units per group, at least one unit with replicates. Output: per-group units and observations,
REML variance components, ICC, design effect 1 + (m̄ − 1)·ICC and effective n, Satterthwaite omnibus F and
contrasts (Holm p, Bonferroni CI), and the unit-means analysis (Welch t for two groups, pooled-variance
contrasts after one-way ANOVA otherwise).

## ancova

```json
{"analysis_type": "ancova", "input": "study.csv", "source": "...",
 "design": {"unit": "mouse", "unit_rationale": "...", "outcome": "Day-21 volume", "outcome_unit": "mm3",
            "covariates": ["baseline_volume_mm3"], "covariate_timing": "before_treatment",
            "covariate_rationale": "Measured at randomization, prespecified in the plan"},
 "comparison": {"groups": ["vehicle", "mAb"], "post_hoc": "vs_control", "control_group": "vehicle",
                "confidence_level": 0.95, "slope_homogeneity_alpha": 0.05}}
```

Columns: `unit_id`, `group`, `outcome`, each covariate. One row per unit. Model: cell-means OLS with
centered covariates and a common slope; homogeneity by the extra-sum-of-squares F for all group ×
covariate terms; adjusted means at the overall covariate means (emmeans); adjusted contrasts (Holm,
Bonferroni) withheld when homogeneity is rejected or an adjustment point lies outside a group's range.
