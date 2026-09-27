# 4PL input and model contract

The CSV requires `sample_id,experiment_id,curve_id,concentration,concentration_unit,response,response_unit`. Optional `observation_id`, `replicate_id`, `exclude`, and `exclusion_reason` preserve provenance. Every curve belongs to one sample, experiment, concentration unit, and response unit; repeated rows remain distinct observations. Excluded rows require a reason. At least six distinct positive concentrations and ten included observations are required.

Set `dose_scale=linear` for physical concentrations, including genuine zero-dose controls. Set `dose_scale=log10` only when the input concentration column is already the base-10 logarithm of its declared unit; a numeric zero then means concentration **one unit**, not a zero-dose control. The resolved run keeps original X and canonical X. Supported molar units are M, mM, uM/µM, nM, pM; mass/volume units are mg/mL, ug/mL/µg/mL, ng/mL. Units are normalized only within the molar or mass family. Reports and figures show each curve in its input unit. `source_unit` means a public example omitted a physical unit; do not treat its half-response value as a measured biological potency.

## Configuration

The config requires `assay.endpoint` (`EC50` or `IC50`), `assay.direction` (`increasing` or `decreasing`), `tested_agent`, `response_definition`, `relative_half_response_supported=true`, and an experiment-specific rationale. Endpoint and observed readout direction are separate decisions. If the data trend significantly the other way (Spearman on positive doses, p<0.05), the curve gets `observed_direction_opposite_to_declared` and is withheld.

```json
{
  "analysis_type": "dose_response_4pl",
  "input": "potency.csv",
  "assay": {"endpoint": "EC50", "direction": "increasing", "tested_agent": "...",
            "response_definition": "OD450 binding signal", "relative_half_response_supported": true,
            "rationale": "Replace with the actual assay justification"},
  "fit": {"weighting": "relative", "fixed_bottom": null, "fixed_top": null},
  "replicates": {"independent_unit": "experiment_id", "conditions_comparable": true},
  "comparisons": [{"id": "Test-vs-Ref-E1", "reference_curve": "Ref-E1", "test_curve": "Test-E1"}],
  "comparison_settings": {"parallelism_alpha": 0.05}
}
```

All defaults and options are in `src/agentic_prism/dose_schema.py`; unknown keys are errors. See [the potency fixture](../../../fixtures/dose_synthetic/config_potency.json).

## Model and fitting

`Y = Bottom + (Top−Bottom)/(1 + 10^(−s × h × (log10 C−log10 C50)))`, where `s=+1` for increasing and `s=−1` for decreasing responses, `h>0`, Top ≥ Bottom (Top is the higher plateau), and `C50` is the midpoint between the plateaus. A true linear-input zero is evaluated at its model limit and drawn in a separate linear panel. The output names the declared `EC50` or `IC50`. It is not an absolute Y=50 interpolation (unless plateaus are fixed at 0/100), an equilibrium KD, or a competition Ki.

Only log10 C50 and log10 h are optimized numerically (multistart). For each proposal, the free plateaus are solved as nuisance parameters with Top−Bottom ≥ 0.

- `fixed_bottom` / `fixed_top` (response units) remove those parameters. Residual df is `n − (2 + number of free plateaus)`. `fixed_top` must exceed `fixed_bottom`. Fix plateaus only from independent control evidence; fixing to hide a missing plateau makes the estimate depend on that assumption.
- `weighting=unweighted` minimizes Σ(Y−Ŷ)². `weighting=relative` minimizes Σ((Y−Ŷ)/Ŷ)² with Ŷ the curve value, Prism's relative (1/Y²) weighting. It requires all included responses and fixed plateaus to be positive. In simulation with 8% CV noise, relative weighting gave 95.5% coverage and unweighted 93.5% (validation record).

The profile-F interval for C50 reoptimizes Hill and plateaus at each candidate C50, with threshold `SS_min·(1 + F(1, df)/df)` on the (weighted) objective. It is conditional on the model, constraints and error assumptions, and is not a confidence band for the curve.

## Replicate summary

`sample_summary.csv` gives, per sample, the geometric mean C50 after averaging log10 C50 of technical curves within each experiment and weighting experiments equally, with a Student-t interval on the log scale. This requires `independent_unit=experiment_id`, `conditions_comparable=true`, a single dose-unit family, and every curve of the sample reportable. Otherwise the status says why it was withheld.

## Curve comparison and relative potency

Each comparison fits three models to the reference and test curves:

1. separate fits (the individual 4PL fits);
2. shared C50 with separate plateaus and Hill: an extra-sum-of-squares F test answers "does C50 differ";
3. the parallel model, with shared plateaus and Hill but separate C50: RP = C50_ref/C50_test with a profile-F interval on log10 RP.

The parallelism F test compares (3) with (1) using `1 + free plateaus` numerator df. RP is reportable only if parallelism is not rejected at `parallelism_alpha`, both single-curve fits are reportable, the parallel model is identifiable and away from bounds, and the interval closes. Curves from different `experiment_id` values are flagged, because the RP then includes between-experiment variation. `potency_summary.csv` gives the geometric-mean RP across independent experiments (one same-experiment comparison per experiment, all reportable), with a log-scale t interval.

The F test detects non-parallelism in proportion to precision. USP <1032>/<1034> recommend equivalence testing against pre-set margins for release assays; use `parallelism_method=equivalence` for that (0.7.0, below).

## Equivalence parallelism and RP acceptance (0.7.0)

```json
"comparison_settings": {
  "parallelism_method": "equivalence",
  "equivalence": {"confidence_level": 0.90, "hill_ratio_limits": [0.8, 1.25],
                  "bottom_difference_limits": [-0.05, 0.05], "top_difference_limits": [-0.3, 0.3],
                  "rationale": "Source of the margins, e.g. historical reference-vs-reference runs"},
  "rp_acceptance_limits": [0.8, 1.25],
  "rp_acceptance_rationale": "Specification or protocol that defines these limits"
}
```

The unrestricted model is the two separate fits. For each curve, `Cov = s²(JᵀJ)⁻¹` of the weighted residual Jacobian over its free parameters (plateaus, log10 C50, log10 h) uses the pooled `s² = (SSE_ref + SSE_test)/(n − 2k)`. The Hill ratio test/reference has a Wald interval on log10 scale; Bottom and Top differences (test − reference, response units) are tested only when their limits are given (a plateau fixed in both fits is reported as not tested). Parallel = every tested interval inside its margin, a two one-sided tests decision at (1 − level)/2 per side. When `equivalence` is selected the F test is still reported but is informational; `parallelism_equivalence_not_demonstrated` withholds RP. Outputs: `parallelism_equivalence.csv`, `parallelism_equivalent` and `rp_acceptance_within` columns in `comparisons.csv`, and `rp_acceptance_within` in `potency_summary.csv`. Hill-ratio limits must bracket 1, plateau limits must bracket 0, and RP acceptance limits must bracket 1; each set needs a rationale. The [equivalence example](../../../fixtures/dose_synthetic/config_potency_equivalence.json) uses illustrative margins only. Simulated operating characteristics (pass rate at the margin) are in the validation record. The software cannot judge whether margins are appropriate.

## Examples

Use [the synthetic EC50](../../../fixtures/dose_synthetic/config_ec50.json), [IC50](../../../fixtures/dose_synthetic/config_ic50.json) and [potency](../../../fixtures/dose_synthetic/config_potency.json) configs for schema examples only. The [public example](../../../fixtures/dose_public/config.json) retains its upstream unit as `source_unit`; its original application uses a fixed-response-50 IC50, so this repository's relative midpoint is a different estimand. Source hashes and license are in [the manifest](../../../fixtures/dose_public/source_manifest.json). For model selection and plateau cautions see the [GraphPad equation guide](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_choosing_a_dr_equation.htm), the [normalization guide](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_pros_and_cons_of_normalizing.htm) and [relative weighting](https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_weighting_tab.htm).
