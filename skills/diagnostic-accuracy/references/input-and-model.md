# ROC input contract (0.13.6)

```json
{"analysis_type": "roc_analysis", "input": "markers.csv", "source": "Study and file provenance",
 "design": {"unit": "patient", "rationale": "One serum sample per patient at diagnosis",
            "reference_standard": "Biopsy histology", "reference_source": "Protocol section 6.2"},
 "markers": [{"column": "marker_a", "label": "Marker A", "unit": "ng/mL",
              "direction": "higher_indicates_condition", "direction_source": "Marker biology, prior study X"}],
 "comparison": null,
 "thresholds": [{"marker": "marker_a", "cutoff": 2.0, "source": "Assay manufacturer cut-off"}],
 "exploratory_youden": false,
 "confidence_level": 0.95}
```

Columns: `unit_id` (unique), `condition` (1 present / 0 absent by the reference
standard), one or two declared marker columns, optional `exclude` (true/false)
with `exclusion_reason`. No other columns.

| Quantity | Method |
|---|---|
| AUC | Mann-Whitney probability P(marker_condition > marker_no_condition) + ½ P(tie), in the declared direction |
| AUC interval | DeLong (1988) placement-value variance; Wald interval on logit(AUC), back-transformed |
| Test vs 0.5 | Mann-Whitney U, two-sided, normal approximation with tie and continuity correction (R `wilcox.test(exact = FALSE)`) |
| Paired comparison | DeLong covariance of two markers' placement values; difference marker 2 − marker 1, Wald interval and z test |
| Cut-off | Called positive if marker ≥ cut-off (higher direction) or ≤ cut-off (lower); Clopper-Pearson exact intervals |
| Youden (optional) | max(sensitivity + specificity − 1) over observed cut-offs; exploratory, not reportable |

`roc_curve.csv` lists every empirical operating point. Withholding: AUC 0 or 1
(complete separation) or zero DeLong variance withholds the AUC interval; a
non-positive paired variance withholds the comparison.
