# Potency-ranking input contract (0.14.0)

```json
{"analysis_type": "potency_ranking", "input": "ec50_table.csv", "source": "...",
 "design": {"unit": "independent experiment", "rationale": "Benchmark on every plate; separate thaws",
            "quantity": "EC50", "quantity_unit": "nM", "estimate_source": "Reportable runs in runs/dose_*"},
 "comparison": {"reference": "BENCH", "confidence_level": 0.95,
                "fold_threshold": 3, "fold_threshold_source": "Triage rule TR-7", "minimum_experiments": 3}}
```

Columns: `candidate`, `experiment_id`, `ec50` (blank when not estimated), `unit`,
`status` (`estimated`, `above_range`, `below_range`, `failed`); optional
`source_run`, `note`. One row per candidate and experiment.

| Step | Method |
|---|---|
| Per candidate | d_e = log10 EC50(candidate) − log10 EC50(reference) in each experiment with an estimated reference; mean, SE, Student t with n − 1 df |
| Relative potency | RP = 10^−mean(d); intervals back-transformed |
| Family | Holm adjusted p values over the m reportable candidates; Bonferroni simultaneous intervals at level 1 − (1 − level)/m; unadjusted intervals for context |
| Threshold tier | `at_least_Fx_more_potent` if the simultaneous lower RP limit ≥ F; `at_least_Fx_less_potent` if the upper limit ≤ 1/F; otherwise `not_resolved_against_threshold` |

Withheld: any non-estimable run with an estimated reference; fewer than
`minimum_experiments` shared experiments; zero between-experiment variance.
