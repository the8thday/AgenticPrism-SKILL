# fih_dose input contract (0.15.2)

Design-only: `"input": null` (omit it). Every number is declared with a source.

```json
{"schema_version": 1, "analysis_type": "fih_dose", "source": "Programme X nonclinical summary v3",
 "molecule": {"name": "mAb-X", "modality": "protein_gt_100kda", "molecular_weight_kda": 150},
 "human": {"body_weight_kg": 70, "source": "Typical adult weight"},
 "noael": {"studies": [{"species": "monkey", "noael_mg_per_kg": 10, "source": "4-week GLP study"}],
           "scaling": "mg_per_kg", "species_selection": "monkey", "selection_rationale": "Only cross-reactive species",
           "safety_factor": 10, "safety_factor_rationale": "Guidance default"},
 "mabel": {"targets": [{"basis": "receptor_occupancy", "name": "10% RO", "kd_nM": 1.0, "target_occupancy_percent": 10, "source": "SPR"},
                       {"basis": "effect_concentration", "name": "EC20 PBMC", "value": 0.05, "unit": "nM", "source": "PBMC assay"}],
           "central_volume_l_per_kg": 0.04, "central_volume_source": "Human plasma volume"},
 "allometric_pk": {"species": "cynomolgus monkey", "animal_body_weight_kg": 3.5,
                   "parameters_per_kg": {"cl": 2.1e-4, "v1": 0.04, "q": 6e-4, "v2": 0.03},
                   "units": {"cl": "L/h/kg", "v": "L/kg", "vmax": "mg/h/kg", "km": "ug/mL"},
                   "exponents": {"cl": 0.85, "q": 0.85, "v1": 1.0, "v2": 1.0}, "exponents_source": "Deng et al. 2011",
                   "candidate_doses_mg_per_kg": [0.001, 0.01, 0.1], "route": "iv_infusion", "infusion_duration_h": 1,
                   "horizon_h": 672, "kd_nM": 1.0, "source": "Cyno population PK"}}
```

## Calculations

- `noael.scaling`: `bsa_table` HED = NOAEL / divisor (FDA 2005 Table 1: mouse 12.3, hamster
  7.4, rat 6.2, ferret 5.3, guinea_pig 4.6, rabbit 3.1, dog 1.8, monkey 3.1, marmoset 6.2,
  squirrel_monkey 5.3, baboon 1.8, micro_pig 1.4, mini_pig 1.1; divisor = 37/Km);
  `body_weight_exponent` HED = NOAEL × (W_animal/W_human)^0.33 (`animal_body_weight_kg` per
  study); `mg_per_kg` HED = NOAEL (required for `protein_gt_100kda`). `species_selection`
  `most_sensitive` or a named species with `selection_rationale`; MRSD = HED / safety factor.
  A selected species that is not the most sensitive, and a safety factor below 10, are listed
  as failing items for review.
- MABEL: receptor occupancy RO gives total drug = KD·RO/(1−RO) + RO·R_T (1:1; R_T optional);
  effect concentration in nM or ug/mL; ug/mL = nM × MW(kDa) × 1e-3; dose (mg/kg) =
  concentration (ug/mL = mg/L) × central volume (L/kg). Sensitivity: KD ÷3 and ×3.
- Allometry: human per-kg parameter = animal per-kg × W_a × (W_h/W_a)^b / W_h (Km unchanged);
  2-compartment simulation (optional Michaelis-Menten) on a 2001-point grid: Cmax, concentration
  at the horizon, AUC to the horizon, AUCinf and terminal half-life (linear model), hours above
  each MABEL concentration (grid approximation) and RO at Cmax and the horizon when `kd_nM` is
  given. Sensitivity: human CL for CL exponents 0.75, 0.85, 1.0.
