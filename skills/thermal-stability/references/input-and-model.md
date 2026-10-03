# Thermal-unfolding input and model contract (0.13.1)

```json
{"analysis_type": "thermal_unfolding", "input": "melt.csv", "source": "Prometheus run 2026-10-01",
 "assay": {"technique": "nanodsf_ratio", "readout": "350/330 nm ratio", "scan_rate_c_per_min": 1.0,
           "buffer": "20 mM histidine pH 6.0", "reversibility": "irreversible_or_unknown", "rationale": "..."},
 "model": {"transitions": 3, "transitions_source": "IgG1: CH2, Fab, CH3", "fit_range_c": null, "dh_bounds_kj_mol": [50, 5000]},
 "derivative": {"window_c": 3.0, "polyorder": 2},
 "replicates": {"independent_unit": "replicate_id", "independence_source": "Separate dilutions on separate days"},
 "comparisons": [{"id": "v1", "reference_sample": "ref", "test_sample": "variant"}],
 "uncertainty": {"level": 0.95}}
```

Columns: `sample_id`, `replicate_id`, `curve_id` (one per sample replicate),
`temperature_c`, `signal`, optional `exclude`/`exclusion_reason`.

Model: f_i(T) = 1/(1 + exp(−(dH_i/R)(1/Tm_i − 1/T))), T in K;
y = a + b(T − T0) + Σ A_i f_i + c(T − T0) f_k, T0 the mean fitted temperature. Tm_1, log
gaps between successive Tm and log10 dH_i are optimized (variable projection for the
k + 3 linear terms); each Tm has a profile-F interval with all other parameters
re-optimized. Gates: at least 5 points beyond Tm_1 − width and Tm_k + width (10–90 %
width 2 ln 9 R Tm²/dH), adjacent Tm further apart than their mean width, amplitude
above 3 × RMSE, closed interval, no bound hit. Derivative: Savitzky–Golay first
derivative on a uniform grid, extrema above 20 % of the largest, top k.
