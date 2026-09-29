# Epitope binning contract prepared for 0.12.1

Analysis type epitope_binning; a separate specialist Skill. Explicit format:
tandem, premix or sandwich. Canonical response data, not a guessed instrument
parser. Each row declares independent experiment ID, first antibody, second
antibody, response, matched no-block reference, matched self-block control and
background. Config declares readout/valency, matched control scale, independence,
background handling and sources, assay rationale and thresholds fixed before
viewing results.

Normalized blocking is (reference-response)/(reference-self), preserving the
ordered pair. No silent clipping or averaging reciprocal observations. Thresholds
have an explicit ambiguous zone. Failed self-controls, insufficient control span,
missing pairs or missing experiments remain visible. Experiments get equal
weights; wells do not become experiments.

Hierarchical clustering uses directed response profiles with average linkage and
Euclidean distance. Complete profiles only; if the matrix is incomplete, retain
the matrix and withhold clustering rather than impute missing competition.
Bootstrap branch stability resamples profile features, as pvclust does. This is
conditional panel stability and must not be presented as independent biological
replication or approximately unbiased (AU) support.

Graph bins use a predeclared reciprocal-block requirement. Both directions must
be classified blocked to create an undirected graph edge. Asymmetric pairs stay
flagged and have no inferred reciprocal edge. Connected communities can contain
non-cliques; membership does not assert identical epitopes or structural contact.
Thresholds and graph rules are never selected after viewing the matrix.

Self-control validation also checks the residual self-block response relative to
the matched no-block reference after declared background handling. A normalized
diagonal of1 alone is tautological when the same self response is its anchor;
it does not establish an effective self-block. Residual self signal must be at
most1-self_block_min (default20%) of the positive reference. This gate is separate
from control-span and observed diagonal checks.

## CSV columns

| Column | Meaning |
|---|---|
| `experiment_id` | Declared independent experiment; exactly one row per ordered pair within it |
| `first`, `second` | Antibody identifiers in experimental order |
| `response` | Format-matched measured endpoint |
| `reference_response` | Matched no-block reference on the same scale |
| `self_response` | Matched self-block endpoint on the same scale |
| `background` | Explicit measured background, or zero when already corrected |

Config requires `source`, `format`, `readout`, `valency`, `assay_rationale`,
`controls_share_response_scale=true`, `independent_experiments` and
`independence_source`; `background_handling` is `subtract_declared` or
`already_corrected`, with `background_source`. Declare
`thresholds_predeclared=true`, `threshold_source`, `control_span_source` and
`clustering.cut_source`. Defaults are blocked >=0.8, nonblocked <=0.2,
self-control >=0.8, minimum control span1, asymmetry difference0.3;
clustering average/Euclidean cut0.5 and1000 feature resamples, seed121260930.
These defaults are recorded choices to justify before analysis, not inferred
assay-specific thresholds. Experiment intervals use equal weights and a95%
Student t interval; fewer than two experiments yield no experiment interval.

A fully specified synthetic example is
[config_complete.json](../../../fixtures/epitope_0121/config_complete.json).
Its source and independence declarations apply to the simulation only.
