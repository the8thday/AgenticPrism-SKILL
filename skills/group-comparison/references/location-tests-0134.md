# Single-sample and ratio tests (0.13.4)

Use `analysis_type=location_test`. This extends the group-comparison specialist;
existing `group_comparison` Welch/paired differences remain unchanged.

| method | Estimand | Analysis |
|---|---|---|
| one_sample_t | Arithmetic mean versus declared null | One-sample Student t |
| one_sample_ratio_t | Geometric mean versus positive null | One-sample t on natural logs; exponentiate interval |
| paired_ratio_t | Geometric mean of within-unit B/A | One-sample t on log(B)-log(A); exponentiate |
| independent_ratio_welch | Geometric mean B / geometric mean A | Welch t on log observations; exponentiate |

```json
{
  "analysis_type": "location_test", "input": "responses.csv", "source": "Study ID",
  "method": "paired_ratio_t",
  "design": {"unit": "donor", "rationale": "Matched aliquots from independent donors",
             "outcome": "IL2", "outcome_unit": "pg/mL"},
  "comparison": {"group_a": "control", "group_b": "treated", "null_value": 1,
                 "null_source": "Prespecified no-change ratio", "confidence_level": 0.95}
}
```

CSV: `observation_id,independent_unit_id,group,outcome,value,unit`; optional
`exclude,exclusion_reason`. One included row per independent unit and group.
For a one-sample test set group_b=null and provide only group_a. Declare the null
and its source before analysis (for example 100% recovery or fold change 1).
No imaginary repeated control group is created to implement a one-sample test.

All ratio inputs and the ratio null must be strictly positive: do not add
pseudocounts, drop zeroes, or equate below-quantification values to zero. Pairing
is by unit ID, not row order; complete pairs are required after explicit
exclusions. Independent groups must have disjoint IDs. Technical wells do not
count as donors. Minimum three units per group or three pairs; zero or
numerically constant analysis-scale variance withholds inference.

The normality assumption concerns the analysis scale, especially within-unit
log ratios in a paired test. Report original units for a single geometric mean
and dimensionless B/A for a ratio comparison, together with the interval and
two-sided p value. A geometric-mean ratio is not a ratio of arithmetic means.
No data-driven choice between difference/ratio tests, automatic exclusions,
equivalence or assay-acceptance conclusion from nonsignificance.

Run analyze and verify as usual, then read saved facts. Rendering uses saved
results and does not repeat tests. Examples: `fixtures/everyday_0134/`.
Evidence: [0.13.4](../../../validation/RELEASE_0.13.4.md).
