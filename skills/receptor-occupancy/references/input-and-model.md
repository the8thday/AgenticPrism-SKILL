# Receptor-occupancy input contract (0.14.2; EC50 interval option 0.14.3)

```json
{"analysis_type": "receptor_occupancy", "input": "ro.csv", "source": "...",
 "design": {"receptor": "PD-1", "drug": "anti-PD-1 IgG4", "cell_population": "CD3+ T cells", "subject_unit": "patient", "matrix": "whole blood within 4 h"},
 "method": {"name": "free_over_total", "rationale": "Assay plan section 3", "reagent_evidence": "Free clone competes (report RO-07); total clone does not (RO-08)"},
 "summary": {"confidence_level": 0.95},
 "exposure_response": {"concentration_unit": "ug/mL", "rationale": "One-site binding expected", "hill": 1, "interval": "cluster_robust_t"}}
```

Columns: `subject`, `group`, `timepoint`; `baseline` (true for the one predose sample per
subject; free formats); format columns below; `drug_concentration` when
`exposure_response` is declared; optional `exclude` with `exclusion_reason`. One sample per
subject and timepoint; each subject in one group.

| Format | Columns | RO (%) |
|---|---|---|
| `free_baseline` | `free_mfi`, `free_background_mfi` | 100 (1 − F/F_predose), F = free − background |
| `free_over_total` | plus `total_mfi`, `total_background_mfi` | 100 (1 − (F/T)/(F/T)_predose) |
| `bound_saturation` | `bound_mfi`, `bound_saturated_mfi`, `bound_background_mfi` | 100 B/B_sat |

Not computable (listed): denominator not above background. No clipping. Summaries per
group and timepoint: mean with t interval over subjects (withheld with fewer than 3
subjects or zero variance), median, counts below 0 and above 100.

Exposure-occupancy: RO = 100 C^h/(C^h + EC50^h) (Emax 100, baseline 0; h fixed or
`"estimate"`), unweighted least squares over post-dose samples; percentile interval from a
subject-cluster bootstrap. Withheld: fewer than 6 subjects or 8 samples, EC50 outside the
observed concentrations, 50% occupancy not bracketed by the data, more than 5% failed
bootstrap fits.

`interval` (0.14.3): `cluster_robust_t` (recommended) gives exp(log EC50 ± t(df) SE) with the
Bell–McCaffrey CR2 cluster-robust SE (subject as cluster) and Satterthwaite df, which depend on
the number of subjects; it is withheld additionally when one subject determines the fit (singular
I − H_gg). `bootstrap_replicates` and `seed` are refused with it. `cluster_bootstrap` (the default,
kept so 0.14.2 configurations reproduce) is the 0.14.2 percentile bootstrap, which undercovered with
6–8 subjects. With `hill: "estimate"` the robust interval for h is reported too.
