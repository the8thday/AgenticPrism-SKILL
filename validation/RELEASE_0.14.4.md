# 0.14.4 — exploratory PCA and clustering

Baseline `b54523e` (0.14.3, published as public `0f0bffb`). Plan:
[plans/RELEASE_0.14.4_PLAN.md](../plans/RELEASE_0.14.4_PLAN.md). New analysis types only; earlier
outputs are checked byte for byte at the end.

[Skill](../skills/exploratory-multivariate/SKILL.md), [contract](../skills/exploratory-multivariate/references/input-and-model.md).
`pca`: declared variables, optional log10, standardize or none, SVD with a fixed sign convention,
parallel-analysis or declared retention, North's-rule and dominance flags, annotations described
against scores. `clustering`: k-means (k-means++ starts, Lloyd, Hartigan transfers) or
hierarchical (Ward / average / complete; Euclidean or correlation distance), on variables or PC
scores; declared k or gap statistic (k = 1 allowed); Hennig bootstrap Jaccard, silhouettes,
descriptive profiles; the output states that testing clustering variables between clusters is
invalid. Prespecified in `f1fa634`; method amended before the calibration rerun in `8b09e70`.

| Gate | Observed scope |
|---|---|
| Published worked example | PCA: ISLR (James, Witten, Hastie, Tibshirani) PCA lab, `prcomp(USArrests, scale = TRUE)`: **all printed loadings (16), centres (4) and scales (4) match** at printed precision (columns sign-aligned; signs are arbitrary). USArrests regenerated from R's `datasets` at validation time. Clustering: UNMET (no printed result with its data located). [Record](multivariate_benchmarks_0.14.4.json) |
| Established implementation | R `prcomp`, `kmeans` (100 starts), `hclust`, `cluster::silhouette` (cluster 2.1.8.2) and gap log W recomputed in R on this implementation's reference data: 6/6 fixtures agree (eigenvalues, loadings and scores to machine precision; identical partitions and within-cluster SS; merge heights; silhouette; gap values within 5.4e-5 and the same k). **First attempt kept:** gap values differed by up to 0.02 because Lloyd-only k-means stopped short of R's Hartigan–Wong optimum on uniform reference data; k-means now refines each start with Hartigan single-point transfers and the default starts rose to 100 (0.7 % of reference fits still short of R's 1000-start optimum, by at most 0.05 %). [Attempt 1](multivariate_benchmarks_attempt1_0.14.4.json) |
| Calibration | Seed 144401, 1000 draws per scenario, all kept. **First run stopped** after 5 of 12 scenarios when the k-means defect was found ([log](multivariate_calibration_attempt1_stopped_0.14.4.log)); the unchanged script was rerun in full after the fix and the four PCA rows reproduced exactly. **5/5 bounded criteria passed:** parallel analysis retained a component in 0.064 of independent-noise datasets (bound 0.071) and exactly two in 1.000 of two-factor datasets; gap chose k = 1 in 1.000 of structureless datasets and the true k = 3 in 0.938 (k-means) and 0.985 (Ward). Descriptive: n = 30 two-factor 0.898; an unscaled variable × 100 took 1.000 of PC1's loading; **the gap statistic is conservative** — with two informative variables out of five it never chose k = 3 (k = 1 in 1.000), and for two elongated groups it chose k = 2 in 0.014 (k-means) and 0.050 (average linkage); forcing k = 3 on structureless data flagged 0.444 of clusters unstable (mean Jaccard 0.61) while true clusters had Jaccard ≥ 0.75 in 1.000. [Record](multivariate_calibration_0.14.4.json) |
| Live-agent misuse | 6/6 passed (unscaled dominance shown as a units artifact; plate, not drug, on PC1; no zero imputation; no between-cluster t tests; no subtypes from structureless data; no silhouette shopping). Finding: on data with three real groups the gap statistic chose k = 1, consistent with its conservatism. [Record](agent-scenarios/0.14.4/RESULTS.md) |
| Interpretation facts | Extraction-only; facts equal saved results |

## Checks across 0.14.4

- Pre-freeze note (in the plan): on the first fixture design (two informative variables of five)
  both this code and `cluster::clusGap` chose k = 1 with the same log W; the fixture was changed and
  the design kept as a descriptive calibration row.
- Unchanged methods: 267 configurations, 1,850 byte-identical scientific artifacts against
  0.14.3 (`b54523e`) ([record](legacy_fixture_checks_0.14.4.json)).
- Full suite: `.venv/bin/python -m pytest -q` 979 passed (2026-10-17).
- Live-agent misuse scenarios 6/6 (criteria `833a821`, not blinded).
- Follow-ups: model-based clustering (Gaussian mixtures with BIC) as an alternative to the
  conservative gap rule; missing-data PCA.
- Clean build from `git archive` of `57f5350`: 120 Python source files in the wheel are
  byte-identical to `src/` (wheel SHA-256 `38959caeb7c30e64…`; [record](archive_build_0.14.4.json)).
- Isolated wheel-only install with locked dependencies: 18 configurations and 106 scientific
  artifacts identical to the source (`multivariate_0144`, the 0.14.3 fixture folders and
  `groups_synthetic`); six runs re-rendered and verified with unchanged scientific artifacts
  ([record](clean_install_0.14.4.json)). Same macOS host only; not a cross-platform claim.
