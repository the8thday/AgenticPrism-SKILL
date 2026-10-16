# Rater-agreement input contract (0.13.6)

Long table: `unit_id`, `rater`, `rating`; optional `exclude` (true/false) with
`exclusion_reason` (whole units only). No other columns.

## Cohen's kappa

```json
{"analysis_type": "rater_agreement", "method": "cohen_kappa", "input": "scores.csv", "source": "...",
 "design": {"unit": "slide", "rationale": "One section per tumor", "raters": "Two pathologists",
            "rater_rationale": "Blinded independent reads", "outcome": "HER2 IHC score", "outcome_unit": "category"},
 "kappa": {"categories": ["0", "1+", "2+", "3+"], "ordered": true,
           "weights": "linear", "weights_source": "Study protocol 4.3"}}
```

`weights`: `none` (default), `linear`, `quadratic` (ordered categories only) or a
square agreement-weight matrix in category order (diagonal 1, entries 0–1)
with `weights_source`. Disagreement weights v convert as 1 − v / max(v).

| Quantity | Method |
|---|---|
| Kappa | (p_o − p_e) / (1 − p_e) with agreement weights |
| Interval | Fleiss, Cohen & Everitt (1969) non-null large-sample variance, Wald, truncated to [−1, 1] and flagged |
| Test vs 0 | z with the FCE null variance |
| Context | Observed/expected agreement, maximum kappa given marginals (unweighted); prevalence and bias indices and PABAK (2 × 2) |

Expected agreement 1 (a single category used throughout) withholds kappa.

## ICC

```json
{"analysis_type": "rater_agreement", "method": "icc", "input": "titers.csv", "source": "...",
 "design": {"unit": "serum sample", "rationale": "...", "raters": "Three operators",
            "rater_rationale": "Operators sampled from the lab's trained pool", "outcome": "titer", "outcome_unit": "log2"},
 "icc": {"model": "two_way_random", "type": "agreement", "unit": "single",
         "rationale": "Future samples will be read by one operator from the pool; absolute titers matter"}}
```

| Form | Estimate | Interval |
|---|---|---|
| one_way_random (ICC(1)) | (MSR − MSW) / (MSR + (k−1) MSW) | Exact F |
| consistency (ICC(C,1) = Shrout-Fleiss ICC(3,1)) | (MSR − MSE) / (MSR + (k−1) MSE) | Exact F |
| agreement (ICC(A,1) = Shrout-Fleiss ICC(2,1)) | (MSR − MSE) / (MSR + (k−1) MSE + k (MSC − MSE)/n) | McGraw-Wong approximate F |
| average of k | Spearman-Brown of the single form | Spearman-Brown of the single limits |

Two-way models need the same raters on every unit; `two_way_random` and
`two_way_mixed` share formulas and differ in the population of raters the
conclusion covers. Zero residual or between-unit mean square withholds.
