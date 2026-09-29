---
name: equilibrium-binding
description: Analyze single-site equilibrium binding titrations for KD with explicit assay checks, profile intervals, replicate handling, diagnostics and offline Prism-like reports. Use for supported concentration-response binding CSVs or re-rendering saved results; not ELISA calibration, generic EC50/IC50 or SPR/BLI rate fitting.
---

# Equilibrium binding

Use the versioned `agentic-prism` runtime shipped with this skill collection.
Do not generate ad hoc numerical scripts for production analysis.

## New analysis

1. Inspect the actual files. Read [input contract](references/input-schema.md)
   and [scientific rules](references/equilibrium-binding.md) before configuring a
   new assay. Establish what is titrated, the response, equilibrium support,
   single-site applicability, concentration basis and experimental independence.
   A good curve fit alone does not establish KD or intrinsic affinity. If these
   facts are missing, request them; do not set applicability booleans to true to
   make the command run. ELISA response, bivalent binding and cell binding require
   their own assay justification, not an automatic apparent-KD label.
2. Create a JSON config beside the input or use a provided compatible config.
   Keep choices explicit and include the evidence rationale. See the distribution
   `fixtures/synthetic/config.json` for schema, not biological defaults. Unknown
   methods/options must fail; do not silently replace them.
3. Locate and check the runtime as described in the
   [runtime instructions](../agentic-prism/references/runtime.md): resolve this file's real path; the collection root is two folders above
   this Skill; use `<root>/.venv` and run `doctor` once. If the runtime is not
   installed, ask the user before running `install.py`.
4. Run `agentic-prism analyze --config <config.json> --output <new-run-dir>`.
   Output must be a new directory. Exit 0 means the workflow ran without hard
   curve failures, not that every KD is reliable. Exit 3 means partial/failed
   curves with saved artifacts; inspect `results.json` and `diagnostics.json`.
   Exit 2 means a configuration/input/runtime error; inspect `failure.json` when
   present. Do not erase a failed run to conceal errors.
5. Inspect both diagnostics and the rendered report. Follow
   [report rules](references/report-and-plotting.md). In particular, do not promote
   audit-only out-of-range estimates to precise reported KD or statistical bounds.
6. Return `report.html`, `fit_results.csv`, `sample_summary.csv`, the resolved
   configuration and rerun command. Explain material assumptions, failed curves,
   baseline-conditioned intervals and unconfirmed experiment relationships.

## Re-render only

For style changes on an existing run, call
`agentic-prism render --run <run-dir> --style prism_like` (or `standard`). This
verifies immutable scientific-artifact hashes and never calls the fitter. Read
[report rules](references/report-and-plotting.md); preserve the original analysis
config and record new rendering provenance separately.

## Unsupported requests

Do not fit kon/koff from endpoint data, report functional EC50 as KD, or add
unvalidated multi-site models. Use [the parent registry](../agentic-prism/references/module-registry.md)
for scope. Model or statistical changes belong in a separately validated release.

## 0.11.0 opt-in models

Read [affinity-depth contract](references/affinity-depth.md) for
`one_site_depletion`, `solution_equilibrium_titration` or `competition_exact`.
Use [cell-binding](../cell-binding/SKILL.md) for cell-surface apparent KD.
All equilibrium runs now save `interpretation_facts.json`; lead with reportable
facts, required disclosures and withheld items. Old hyperbola science is unchanged.
