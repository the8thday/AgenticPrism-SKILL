# Exploratory PCA / clustering input contract (0.14.4)

```json
{"analysis_type": "pca", "input": "cytokines.csv", "source": "...",
 "design": {"unit": "donor", "unit_rationale": "One serum sample per donor", "purpose_rationale": "Explore marker structure before a confirmatory study"},
 "data": {"id_column": "donor", "variables": ["IL6", "IL8", "TNF", "MCP1"], "annotations": ["batch", "arm"],
          "log10": ["IL6", "IL8", "TNF", "MCP1"], "log_rationale": "Right-skewed concentrations"},
 "scaling": {"method": "standardize", "rationale": "Different units and ranges"},
 "components": {"retain": "parallel_analysis", "replicates": 1000, "seed": 14401, "quantile": 0.95}}
```

```json
{"analysis_type": "clustering", "input": "profiles.csv", "source": "...",
 "design": {"unit": "lot", "unit_rationale": "...", "purpose_rationale": "..."},
 "data": {"id_column": "lot", "variables": ["m1", "m2", "m3"], "annotations": ["site"]},
 "scaling": {"method": "standardize", "rationale": "..."},
 "clustering": {"method": "kmeans", "rationale": "Compact groups expected", "starts": 100, "seed": 14402},
 "k": {"select": "gap", "max": 8, "reference_replicates": 100},
 "stability": {"bootstrap_replicates": 100}}
```

Every column must be the id, a variable, an annotation, or `exclude`/`exclusion_reason`.
Variables must be numeric (text such as "<LLOQ" is refused); log10 needs positive values.
Rows with a missing variable are excluded and listed.

| Item | Rule |
|---|---|
| Processing | Optional log10; centre; `standardize` divides by the SD (n − 1), `none` keeps units |
| PCA | SVD; eigenvalue = s²/(n − 1); unit loadings (largest absolute loading positive); scores = U S |
| Retention | Declared integer, or parallel analysis: column permutations, keep component j while its eigenvalue exceeds the declared quantile of permuted eigenvalue j |
| PCA flags | Neighbouring eigenvalues within North's rule (λ√(2/n)); fewer than 5 units per variable; one variable over half the variance unscaled; nothing retained |
| Clustering | `kmeans` (k-means++ starts, Lloyd, Hartigan transfers; best of `starts`) or `hierarchical` (`linkage` ward/average/complete; `distance` euclidean or correlation = 1 − r; Ward needs Euclidean); `space` variables or `pca_scores` with `pca_components` |
| k | `fixed` (≥ 2) or `select: gap`: W = within-cluster SS; references uniform in the PC-aligned box; smallest k with gap(k) ≥ gap(k+1) − s(k+1); k = 1 reports no clusters |
| Descriptive | Hennig bootstrap Jaccard per cluster (< 0.6 flagged), silhouette widths, cluster means/SD, annotation cross-tabs, PC–annotation eta² or r |

Artifacts: `pca_scores.csv`, `pca_loadings.csv` or `cluster_assignments.csv`, with
`results.json` and `interpretation_facts.json`.
