# GR input contract (0.14.0)

```json
{"analysis_type": "growth_rate_inhibition", "input": "counts.csv", "source": "...",
 "design": {"cell_line_rationale": "Why these lines", "readout": "Nuclei count (imaging)", "concentration_unit": "uM",
            "time_unit": "h", "independent_unit": "experiment (separate seeding day)", "independence_source": "Lab notebook"}}
```

Columns (authors' "input 1" layout plus `experiment`): `cell_line`, `agent`, `time`,
`experiment`, `concentration` (> 0), `cell_count`, `cell_count__ctrl` (untreated
at the endpoint), `cell_count__time0`; optional `exclude` with `exclusion_reason`.
One row per curve and concentration.

| Quantity | Definition |
|---|---|
| GR value | 2^(log2(x/x0) / log2(x_ctrl/x0)) − 1, ratios floored at 1e-6 |
| Curve | GR(c) = GRinf + (1 − GRinf)/(1 + (c/GEC50)^h); bounded least squares (GRinf −1..1, GEC50 within 100× of the tested range, h 0.1..5) |
| Flat test | F test of sigmoid versus the mean, alpha 0.05, df (2, n − 3 + 1) as in the reference code; flat curves have no GEC50/GR50 |
| GR50 | Fitted concentration where GR = 0.5; withheld if not reached, flat or outside the tested range |
| GRmax | Default (authors): minimum observed GR at the two highest concentrations. Opt-in `fit.grmax_definition: fitted_minimum_top_two` (0.14.1): minimum of the selected fit (sigmoid or flat level) at the two highest concentrations; observed value kept as `GRmax_observed_audit` |
| GR_AOC | Trapezoid area of 1 − GR over log10 concentration ÷ log10 range |
| Summaries | t intervals for GRmax and GR_AOC; geometric t interval for GR50 (withheld if any GR50 not estimated) |

`gr_values.csv` stores every GR value. Withheld curves: untreated control grew
less than half a doubling (< √2-fold), fewer than 5 concentrations.
