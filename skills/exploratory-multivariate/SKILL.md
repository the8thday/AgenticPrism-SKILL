---
name: exploratory-multivariate
description: Exploratory principal component analysis (PCA) and clustering of independent units measured on several variables - cytokine or biomarker panels, assay profiles, sample characterization - with declared variables, transform and scaling, parallel-analysis component retention, k-means or hierarchical clustering with gap-statistic k selection, bootstrap cluster stability and silhouettes, and batch/group annotations shown descriptively. Use to describe structure and generate hypotheses; not for testing group differences (group-comparison, panel-comparison), supervised prediction (multivariable-regression) or single-cell pipelines.
---

# Exploratory PCA and clustering

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use `analysis_type: pca`
or `analysis_type: clustering`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Exploration, not inference.** PCA and clustering describe these data. They
   give no p-values and cannot confirm groups, mechanisms or subtypes. Say so when
   the user wants a confirmatory answer, and point to the analysis that tests it.
2. **Units and variables.** One row per independent unit (donor, animal, lot);
   average technical replicates first by a declared rule. Declare the variables
   before looking at results. Columns such as batch, run, site or treatment arm are
   `annotations`: they are never used to build components or clusters and are shown
   against them afterwards. Do not drop or add variables to make a plot look better.
3. **Transform and scale.** Declare log10 for skewed concentrations. Standardize
   when variables have different units or ranges; unscaled (covariance) PCA lets the
   largest-variance variable dominate. Values below a quantification limit are not
   numbers: do not substitute zero or LLOQ/2; say the data need another approach or
   exclude those variables with a stated reason.
4. **Missing values.** Rows with a missing variable are excluded and listed. Do not
   impute with zeros or means to keep units.
5. **Clustering choices before results.** Declare the method (k-means for compact
   groups; hierarchical with Ward, average or complete linkage), the distance and
   whether k is fixed (with a reason) or chosen by the gap statistic. Do not try
   several methods or k values and keep the most convincing.

## Run and report

PCA: report the variance explained, the components retained by parallel analysis
(or the declared number), loadings of retained components, and flags (close
eigenvalues, few units per variable, unscaled dominance). Component signs are
arbitrary. If a component lines up with batch or run in `annotation_summary`,
report that as a batch effect, not biology.

Clustering: report k and how it was chosen. If the gap statistic gives k = 1, say
the data support no cluster structure and do not present groups. Report cluster
sizes, profiles, Jaccard stability (below 0.6 is unstable) and silhouette width.
Never test the clustering variables between clusters (t test, ANOVA): the clusters
were built to differ on them. A cluster that matches an annotation (batch, site)
reflects that factor. Clusters are hypotheses until confirmed in independent data.

Respect `must_mention` and failing items in `interpretation_facts.json`.

## Evidence (0.14.4)

See [the 0.14.4 record](../../validation/RELEASE_0.14.4.md): ISLR's printed
USArrests PCA reproduced, base-R / `cluster` oracle, registered calibration, and the
unmet clustering published-example gate.
