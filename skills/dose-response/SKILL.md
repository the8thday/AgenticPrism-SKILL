---
name: dose-response
description: Fit empirical 4PL concentration-response curves for relative EC50 or IC50, with optional fixed plateaus and relative (1/Y²) weighting, summaries across independent experiments, and reference-vs-test comparison (relative potency, parallelism by F test or by predeclared equivalence margins, optional RP acceptance limits, shared-C50 F test), plus diagnostics and offline reports. Use for ELISA/FACS binding EC50, blocking or neutralization IC50, reporter, ADCC/CDC and other dose-response data; not equilibrium KD, Ki, or ELISA unknown interpolation.
---

# Dose response (4PL)

Use this specialist for a measured response across a concentration series when a four-parameter, symmetric sigmoid is scientifically appropriate. Read the [input and model contract](references/input-and-model.md) before configuring new data. This module reports **relative EC50 or IC50**, the concentration halfway between fitted (or fixed) Bottom and Top. It is the absolute Y=50 concentration only when the plateaus are fixed at 0 and 100 on a normalized scale.

1. Establish the tested agent, assay endpoint, what the response measures, dose units, whether X contains linear or log10 concentrations, the observed curve direction, and whether a plateau midpoint is the intended quantity. `endpoint=IC50` does not itself imply that the measured response falls: an inhibition-percent readout may rise with dose. Do not copy fixture applicability flags into real assays without checking them.
2. Decide the fit options from the assay, not from which result looks better:
   - `fit.fixed_bottom` / `fit.fixed_top` only when the plateau is known independently, typically 0 and 100 after normalization to real controls. Top is always the higher plateau.
   - `fit.weighting: "relative"` when noise grows in proportion to signal (common for OD, luminescence, fluorescence). It needs positive responses. Check the residual plot either way.
   - `replicates.independent_unit: "experiment_id"` with `conditions_comparable: true` only for truly independent, comparable experiments (separate plates or days with fresh dilutions). Never infer independence from dates or well replicates.
   - `comparisons` for reference-vs-test relative potency, one per experiment, pairing curves from the same plate.
3. Preserve individual observations and exclusions in the [canonical CSV](references/input-and-model.md). Plate-reader grids with a plate map can be converted by `agentic-prism import-plate` with `"target": "dose_response_4pl"` (layers or constants for `sample_id`, `experiment_id`, `curve_id`, `concentration`, units; the `curve_id` layer selects used wells); see the ELISA contract's plate-import section for the grid format. No normalization, averaging, outlier deletion, or reference subtraction happens automatically. Do not convert mass/volume concentration to molarity without a justified molecular weight and explicit upstream conversion.
4. Locate and check the runtime per the [runtime instructions](../agentic-prism/references/runtime.md) (collection root two folders above this Skill after resolving symlinks; `<root>/.venv`; `doctor` once; ask before running `install.py`). Run `agentic-prism analyze --config <config.json> --output <new-directory>`, then `agentic-prism verify --run <directory>`. Inspect `results.json`, `diagnostics.json`, `sample_summary.csv`, `comparisons.csv`, `potency_summary.csv`, residuals, and the offline report. `render` switches visual style without refitting.
5. Report the declared endpoint in the input unit, fitted or fixed Bottom/Top, signed Hill slope, weighting, and the profile-F interval when closed. Withhold a precise estimate when the report does: half-effect outside the tested doses, a numerical bound hit, an open interval, or `observed_direction_opposite_to_declared`, which usually means the configured direction is wrong. A Hill slope is a descriptive shape parameter, not proof of cooperativity.
6. For comparisons, report RP = C50_reference / C50_test (>1 means the test article is more potent), its profile-F interval, and the parallelism decision. RP is withheld when parallelism is not accepted, either curve is not reportable, or the interval is open. The geometric-mean RP across independent experiments is the quantity to quote; a single-plate RP interval does not include between-plate variability.
7. Parallelism method. `comparison_settings.parallelism_method=f_test` (default) rejects at `parallelism_alpha`; it is a significance test. For quality decisions (lot release, bridging, biosimilarity-style comparisons) USP <1032>/<1034> prefer `equivalence`: intervals (default 90%) for the Hill-slope ratio and, when declared, Bottom/Top differences must lie inside **predeclared margins**. Margins must come from the laboratory's historical reference-vs-reference data or an approved protocol, recorded in `equivalence.rationale`. Never derive, widen or choose margins from the data being analyzed, and never borrow the fixture's illustrative margins; if the user has none, stop and ask. `rp_acceptance_limits` (bracketing 1, with a rationale) adds a pass/fail check that the whole RP interval lies inside the limits, per comparison and for the cross-experiment geometric mean. An open or missing interval is never a pass. Report the decision together with the interval and margins, not as a bare pass/fail.

Current scope: per-curve symmetric 4PL; unweighted or relative weighting; optional fixed plateaus; per-sample summaries on log C50 across declared independent experiments; pairwise parallel-line relative potency with F-test or equivalence parallelism and optional RP acceptance limits. No 5PL for EC50/IC50 (5PL is only in ELISA calibration), biphasic model, weighting by supplied SDs, automatic plate normalization, mixed-model cross-plate potency, or Prism numerical equivalence is validated. See [validation](../../validation/README.md).

The [AgenticPrism router](../agentic-prism/SKILL.md) selects this specialist by experimental purpose. Use it directly when the task is clearly EC50/IC50 dose response or relative potency.

## Across-run potency (0.10.0)

For a random-run combination of independently fitted RP determinations or
validation across nominal RP levels, use [potency-assay](../potency-assay/SKILL.md).
It imports verified dose-response comparisons using equivalence parallelism,
retains per-run suitability, and withholds combination when a supplied run fails.
At least two independent determinations per run are needed to separate run and
residual variance; shared reference fits cannot count as independent.


0.13.0 adds saved-result `interpretation_facts.json`. Read all states and required
disclosures before interpretation; incomplete summaries, independence, QC and
interval limits remain visible. Facts and rendering never refit the analysis.
