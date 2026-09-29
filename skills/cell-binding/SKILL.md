---
name: cell-binding
description: Analyze cell-surface binding titrations for apparent KD with measured nonspecific controls, explicit background and receptor depletion, profile intervals and independent-experiment summaries. Use for justified equilibrium flow or cell MFI binding; never label cell EC50 as intrinsic KD.
---

# Cell binding

Use the packaged Python runtime and read [the input contract](references/input-contract.md).
Establish the estimand before fitting: this specialist reports **apparent KD**.
Bivalent IgG binding on cells is avidity-influenced. A functional cell EC50 is not
an intrinsic KD and cannot be converted by changing a label.

1. Inspect cell line, receptor-density stratum, valency, incubation time,
   temperature and internalization, detection and wash protocol. Require an
   experiment-specific rationale for equilibrium and negligible wash dissociation.
   Secondary antibody detection needs evidence that detection does not distort
   equilibrium occupancy. Do not invent applicability declarations.
2. Establish the measured nonspecific series (isotype, antigen-negative cells or
   excess competitor), its comparability to the total series, and explicit
   unstained/secondary-only background handling. The nonspecific slope is fitted
   jointly, never silently subtracted from another experiment. Declare
   `nonspecific.baseline`: `shared` only for controls on the same cells;
   antigen-negative cells require `separate` (own control intercept).
3. Declare negligible depletion with evidence or quadratic depletion with sourced
   cells/mL and quantified receptors/cell. Do not infer calibrated receptor counts
   from RNA, ordinary MFI or an unsourced nominal number.
4. Locate `<collection>/.venv/bin/agentic-prism` using the
   [runtime instructions](../agentic-prism/references/runtime.md). Run `doctor`.
   Use `analyze --config CONFIG --output NEW_DIRECTORY`, then `verify --run DIR`.
5. Lead with `interpretation_facts.json`: reportable apparent KD or withholding
   reasons, profile status, every calibration miss and unfulfilled evidence gate.
   Audit optima are not reportable points. Read all `must_mention` and failing items.
6. Wells and technical curves are not independent experiments. Use declared
   experiment IDs; summarize equally on log KD only when conditions are comparable.
   Any incomplete estimate withholds the entire sample summary. Do not pool cell
   lines or receptor-density strata.
7. Return report.html, fit_results.csv, sample_summary.csv, facts, resolved config
   and rerun.txt. Explain assay-specific limitations; a fit does not validate biology.

`render --run DIR --style standard` or `prism_like` uses saved artifacts without
refitting. New evidence is tracked in [0.11.0](../../validation/RELEASE_0.11.0.md).
Published worked-example gate is UNMET; synthetic fixtures do not replace it.

Registered 0.11.0 calibration: unweighted hyperbolic apparent-KD coverage with
fitted nonspecific slope was **93.4%**, below the **93.62%** bound (1,000
simulations; MCSE 0.785 percentage points). Preserve this miss. Relative weighting
and quadratic cell depletion have numerical comparisons but no separate coverage
calibration; do not extend the unweighted coverage evidence to them.
