# Outlier-screen input contract (0.13.6)

```json
{"analysis_type": "outlier_screen", "method": "generalized_esd", "input": "replicates.csv", "source": "...",
 "design": {"unit": "independent determination", "rationale": "Six determinations per lot on separate days",
            "outcome": "Potency", "outcome_unit": "% of reference",
            "purpose": "SOP screen before reporting; flags are reviewed for assignable causes",
            "prespecified_source": "SOP QC-123 section 7.4"},
 "test": {"alpha": 0.05, "max_outliers": 2, "scale": "original"}}
```

Columns: `unit_id`, `group`, `value`; optional pre-existing `exclude` (true/false)
with `exclusion_reason`. Every group is screened separately.

| Method | Settings | Statistic |
|---|---|---|
| `grubbs` | `alpha`, `sides`: `two_sided`, `max` or `min` | G = max deviation / SD; critical (n−1)/√n · √(t²/(n−2+t²)), t at α/(2n) two-sided or α/n one-sided, n−2 df; Bonferroni-type p; n ≥ 3 |
| `generalized_esd` | `alpha`, `max_outliers` r | Rosner R_i and λ_i for i = 1..r; outliers = largest i with R_i > λ_i; n − r ≥ 3 |
| `rout` | `q` (false discovery level, e.g. 0.01) | Lorentzian robust location; RSDR = 68.27th percentile of \|residuals\| × √(N/(N−1)); t = \|res\|/RSDR, N−1 df; α_i = Q(N−i+1)/N over ranks int(0.7N)..N; first significant rank and all larger residuals flagged; n ≥ 4 |

`scale: log` (positive values only) requires `scale_rationale`. Output: test
tables, flagged units, and `summary_all` / `summary_without_flagged` per group.
