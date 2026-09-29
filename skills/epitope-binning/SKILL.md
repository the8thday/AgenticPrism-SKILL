---
name: epitope-binning
description: Analyze directed epitope competition matrices from declared tandem, premix or sandwich assays with matched reference and self-block controls, prespecified thresholds, asymmetric-pair diagnostics, hierarchical bootstrap stability and reciprocal-block graph communities.
---

# Epitope binning

Read [input contract](references/input-and-model.md) and the shared
[runtime](../agentic-prism/references/runtime.md). Establish the assay format,
what the detector measures, valency, background and whether every response and
control uses one scale. Require declared independent experiments and threshold
sources before analysis. Do not tune thresholds after seeing a matrix.

Use `analysis_type=epitope_binning` and canonical ordered-pair CSV. Vendor files
must first use a verified importer or an explicit documented conversion. Never
infer an unsupported export layout. Keep missing/ambiguous cells and both pair
directions; asymmetric blocking must not become a non-competing classification.

Self-block controls must demonstrate actual suppression relative to reference;
a self-normalized value of one alone is insufficient. Failed or missing
controls withhold affected pairs. No imputation precedes clustering.

Hierarchical clusters describe outgoing blocking profiles. Bootstrap support
resamples panel features and is conditional on this panel; it is not biological
confidence or pvclust multiscale AU. Graph communities require reciprocal
blocking edges; a connected component can contain non-competing members. Never
merge bins by eye or conclude structural epitope identity from competition.

Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`, then read
`interpretation_facts.json` and all failing items before interpreting results.
The directed heatmap, community table and branch stability must be considered
together. See [release evidence](../../validation/RELEASE_0.12.1.md).
