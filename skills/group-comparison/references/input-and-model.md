# Group comparison input and design contract

CSV columns: `observation_id,independent_unit_id,group,outcome,value,unit`. Optional `exclude=true|false` and `exclusion_reason` preserve exclusions. Input contains one outcome in one physical unit and exactly the declared groups. Repeated wells, cells or technical reads must be reduced upstream to one prespecified unit-level value with provenance; this module never silently averages them. Nonfinite values and duplicate observation IDs fail.

## Two groups (`analysis_type=group_comparison`)

Config example: [synthetic paired comparison](../../../fixtures/groups_synthetic/paired_config.json). Set `input`, and `comparison.design` to `paired_t` or `independent_welch`. Also declare `group_a`, `group_b`, `outcome`, `unit`, an experiment-specific `rationale`, and optionally `confidence_level` (default 0.95). A second [synthetic independent comparison](../../../fixtures/groups_synthetic/independent_config.json) demonstrates the other design. Each included independent unit occurs once in each paired condition, or once in total in an independent design.

The estimand is `mean(group B)−mean(group A)`. The paired design analyzes within-unit differences and uses `n_pairs−1` degrees of freedom. The independent design uses Welch's standard error and Satterthwaite degrees of freedom. Both use a two-sided Student-t interval for the mean difference.

## Three or more independent groups (`analysis_type=multi_group_comparison`)

Config examples: [Dunnett vs isotype](../../../fixtures/multigroup_synthetic/config_dunnett.json) and [Games-Howell all pairs](../../../fixtures/multigroup_synthetic/config_games_howell.json).

```json
{"analysis_type": "multi_group_comparison", "input": "data.csv",
 "comparison": {"design": "one_way_anova", "post_hoc": "dunnett_vs_control", "control": "isotype",
                "groups": ["isotype", "mAb-A", "mAb-B", "mAb-C"], "outcome": "tumor_volume_change",
                "unit": "fold", "rationale": "One value per independently randomized animal",
                "confidence_level": 0.95, "sd_ratio_warning": 3}}
```

`groups` (≥3, declared order) must equal the groups in the file; every unit ID appears once. `control` is required for `dunnett_vs_control` and `holm_welch_vs_control` and forbidden otherwise.

| design | post_hoc | method |
|---|---|---|
| `one_way_anova` | `dunnett_vs_control` | Classic F; Dunnett many-to-one with pooled variance (SciPy multivariate-t, fixed QMC seed recorded in the manifest) |
| `one_way_anova` | `tukey_all_pairs` | Classic F; Tukey-Kramer studentized range, pooled variance |
| `welch_anova` | `games_howell_all_pairs` | Welch F; Games-Howell with per-pair Welch df |
| `welch_anova` | `holm_welch_vs_control` | Welch F; Welch t vs control, Holm step-down p, Bonferroni simultaneous intervals |
| either | `none` | Omnibus only |

Contrast direction is later-declared group minus earlier group (for vs-control families: group minus control). Intervals are simultaneous at `confidence_level` for the family; p values are family-adjusted. Contrasts are not gated on the omnibus test. `eta_squared` is given for the classic design. Outputs: `omnibus.csv`, `group_summary.csv` (n, mean, SD, SEM), `contrasts.csv`.

Assumptions: independent units, approximately normal outcome within groups; the classic design also assumes equal variances. Diagnostics flag an SD ratio above `sd_ratio_warning` under the classic design and Welch designs with a group below six units. Repeated measures and scoped random-intercept models use the separate [repeated-measures contract](../../repeated-measures/references/input-and-model.md). Covariates, general clustering and nonparametric tests are outside this contract. A small p value alone is not an effect-size or equivalence claim; plotting style does not enter the calculations.
