# Input and configuration contract

Input is a UTF-8 CSV long table. Required columns:
`sample_id, experiment_id, curve_id, concentration, concentration_unit, response, response_unit`.
Curve IDs are globally unique within a run and map to one sample, experiment and
response unit. Preserve author IDs separately if needed. Concentrations are
nonnegative and finite; responses are finite and may be nonpositive only in
linear mode. Supported concentration units: M, mM, uM/µM/μM, nM, pM. Mass
concentrations are not automatically converted. Response values are not silently
normalized or background subtracted.

Optional columns: `observation_id` (unique), `technical_replicate_id`,
`response_sd`, `exclude` (`true`/`false`), `exclusion_reason`, and source metadata.
Missing observation IDs are stable source-row IDs. Blank experiment_id is allowed
only with `replicates.independent_unit: none`; it does not imply independence.
Exclusions require reasons and remain in saved data and plots. Duplicate
concentrations are observations, not automatically duplicate records.

`column_map` maps original names to canonical names. Units always come from the
CSV; there is deliberately no competing configuration unit override. Unknown
response units are preserved as labels, never invented. Extra columns are retained.

This version fits every supplied included row. If rows are pre-aggregated means,
state this in `source`; provide the standard error of that mean as `response_sd`
when weighting means. The column name denotes the SD of the supplied observation,
not necessarily the raw individual measurements. It does not infer n or recover
replicates from SD/SEM. `inverse_sd` uses residual/SD (squared weight 1/SD²), only
in linear mode; positive finite SD is required. These are relative noise scales;
the overall residual variance is estimated. Cross-experiment SD is not a
within-curve weight. Unweighted is the default.

Configuration is JSON; unknown keys are errors. Example:

```json
{
  "input": "titration.csv",
  "source": "Describe the actual experiment and data layer",
  "assay": {
    "equilibrium_supported": true,
    "signal_proportional": true,
    "single_site_supported": true,
    "concentration_basis": "free",
    "rationale": "Replace with actual evidence; these flags are not a template approval"
  },
  "fit": {"baseline_mode": "fitted", "residual_scale": "linear", "weighting": "unweighted"},
  "uncertainty": {"parameter_ci": "profile_f", "level": 0.95},
  "replicates": {"independent_unit": "none"},
  "report": {"plot_style": "prism_like"}
}
```

All defaults and supported options are defined in `src/agentic_prism/schema.py`
and written into `config.resolved.json`. Input paths resolve relative to the
configuration file. The resolved run config points to its copied `input.csv`.
Only Chinese HTML, 85/180 mm figure widths, 300/600 dpi PNG and the three export
formats SVG/PDF/PNG are supported in this release.
