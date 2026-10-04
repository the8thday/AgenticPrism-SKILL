# Nonlinear model contract (0.13.3)

```json
{"analysis_type": "nonlinear_fit", "input": "serum_stability.csv", "source": "...",
 "model": {"name": "one_phase_decay", "rationale": "First-order loss of intact antibody in prior serum studies", "weighting": "none"},
 "design": {"x_label": "Time in serum", "x_unit": "day", "y_label": "Intact antibody", "y_unit": "% of day 0", "rationale": "..."},
 "replicates": {"independent_unit": "experiment_id", "independence_source": "Separate serum lots and days"},
 "uncertainty": {"level": 0.95}}
```

Columns: `sample_id`, `experiment_id`, `curve_id`, `x`, `y`, optional `exclude`/`exclusion_reason`.

| Model | Equation | Parameters | Derived | Gate |
|---|---|---|---|---|
| one_phase_decay | plateau + (y0 − plateau)·e^(−kx) | y0, plateau, k | half_life = ln2/k | plateau withheld if max x < 3 half-lives |
| one_phase_association | y0 + (plateau − y0)(1 − e^(−kx)) | y0, plateau, k | half_time = ln2/k | same |
| two_phase_decay | plateau + span_fast·e^(−k_fast x) + span_slow·e^(−k_slow x), k_fast > k_slow | 5 | half-lives of both phases | all withheld unless two phases beat one (F test, p < 0.01) |
| exponential_growth | y0·e^(kx) | y0, k | doubling_time = ln2/k | — |
| michaelis_menten | vmax·x/(km + x) | vmax, km | — | both withheld if max x < km |
| logistic_growth | asym / (1 + e^((xmid − x)/scal)) | asym, xmid, scal | rate = 1/scal | asym withheld if max x < xmid + 2·scal |

Rates are optimized on the log scale. Each parameter has a profile-F interval: all values with
S ≤ S_min·(1 + F(1, df)/df), other parameters re-optimized (the same set as R `confint.nls`). Derived
times use the mapped rate interval. A Wald–Wolfowitz runs test on residual signs is reported (descriptive).
