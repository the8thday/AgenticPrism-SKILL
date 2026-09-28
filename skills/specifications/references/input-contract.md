# Specification input contract (0.10.1)

```json
{
  "analysis_type": "specification",
  "input": "release.csv",
  "source": "Release data, lots 2023-2026",
  "attribute": {"name": "Protein concentration", "unit": "mg/mL", "rationale": "..."},
  "design": {"independent_units": true, "unit_description": "one release result per lot"},
  "tolerance": {"method": "normal", "sides": "two_sided", "content": 0.99, "confidence": 0.95, "normal_factor": "exact"},
  "capability": {"lsl": 47.5, "usl": 52.5, "limits_source": "Specification SP-4 v2", "confidence": 0.95, "subgroup_column": null}
}
```

- `tolerance.method`: `normal` or `nonparametric` (or omit the section).
  `sides`: `two_sided`, `lower`, `upper`. `normal_factor` (two-sided normal only):
  `exact` (default; Odeh integral) or `howe` (NIST/SEMATECH form). One-sided normal
  intervals use the exact noncentral-t factor.
- `capability`: `lsl` and/or `usl` with `limits_source`; optional `subgroup_column`
  for within-subgroup Cp/Cpk (Pp/Ppk always use the overall SD).
- Columns: `observation_id`, `value`, `exclude`, `exclusion_reason`, and the
  subgroup column when declared.
