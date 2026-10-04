# qPCR relative-quantification contract (0.13.2)

```json
{"analysis_type": "qpcr_relative", "input": "cq.csv", "source": "...",
 "assay": {"chemistry": "SYBR Green with melt curve", "instrument": "...", "cq_method": "Auto threshold", "rationale": "..."},
 "genes": {"targets": ["IL6"], "references": ["GAPDH", "ACTB"], "reference_rationale": "geNorm M < 0.5 in the pilot"},
 "efficiency": {"mode": "declared", "values": {"IL6": 1.93, "GAPDH": 1.95, "ACTB": 2.0}, "source": "Standard curves, run R-07"},
 "design": {"conditions": ["vehicle", "mAb"], "control_condition": "vehicle", "pairing": "independent", "pairing_rationale": "..."},
 "qc": {"cq_max": 35, "technical_sd_flag": 0.3, "reference_shift_flag_cycles": 0.5, "source": "Lab SOP-4"}}
```

Columns: `sample_id` (biological replicate), `condition`, `gene`, `cq` (number, empty or "Undetermined"),
optional `role` (sample/ntc), `pair_id` (paired designs), `exclude`/`exclusion_reason`.
log2 NRQ = −Cq_target·log2 E_target + mean over references of Cq_ref·log2 E_ref (E = 2 gives −ΔCq).
Fold change = 2^(mean log2 NRQ in a condition − mean in the control), with Welch or paired t intervals.
