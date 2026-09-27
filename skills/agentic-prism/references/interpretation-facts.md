# Interpretation facts, schema 1

Implemented in 0.8.1 for repeated-measures, time-to-event and tumor-growth runs.
Read this artifact after `verify`, before drafting the scientific interpretation.
Earlier runs and the other specialists do not contain it.

## Contract

`interpretation_facts.json` is a deterministic view of saved `results.json` and
`config.resolved.json`. It does not fit a model, calculate a new interval, change
a test, or decide whether the declared experimental assumptions are true.
The run manifest hashes it. `source_sha256` identifies its result, configuration,
diagnostic and normalized-input files; `file.json#/path/to/value` references use
JSON Pointer escaping. Rendering preserves it byte for byte.

- `analysis_type`, `schema_version`: dispatch and artifact schema, independently
  of the package version.
- `model`, `estimand`, `declared_design`: the fitted model, the question it can
  address, and the user's declarations (not independently verified evidence).
- `primary_result_ids`: narrative order, chosen by design, never by p value.
- `results`: individually identified tests or estimates. Each has an estimand,
  context, source references, `reportability`, and reasons. `reportability` is
  `reportable`, `limited` (only with the stated caveats), or `withheld`.
  A withheld result has null report-facing estimate, test and interval values;
  its `source` points to the audit result. Do not recover an audit number to
  bypass a withheld state. A missing point estimate is distinct from a missing
  interval. Null limits are never zero, infinity or a search bound.
- `interval`: confidence level, method, multiplicity, bounds and availability.
  `not_applicable` means this is a test, with no parameter interval.
  `unavailable` / `partial` mean the saved result does not provide both limits;
  they do not alone establish a mathematically unbounded confidence set.
- `required_disclosures`: every saved diagnostic, reportability flag and status,
  with its source. A flag can describe only a submodel: a successful KM curve
  does not make a withheld Cox model reportable. Median-not-reached and other
  limitations expressed only through null values also appear on result rows.
- `limitations`: model assumptions and interpretation boundaries to retain in
  the narrative. They are not evidence that those assumptions hold.

## Narrative rules

Lead with the estimand and the primary estimates and their intervals, alongside
the prespecified tests. Include reasons for limited or withheld results and all
diagnostics material to the claimed conclusion. Do not equate non-significance
with equivalence. Never choose a model or contrast after seeing its p value.
The facts cannot validate endpoint definitions, independence, MAR, censoring or
the rationale in a config: ask the user when those are missing or contradictory.

Repeated measures lead with the corrected omnibus / interaction and declared
mean-difference family; fixed-effect means retained after failed inference are
audit information. Survival leads with KM summaries and log-rank; medians not
reached stay unavailable, Cox is separate, and a rank-zero log-rank supplies no
test of arm equality. Tumor growth leads with the growth model, then observed
readout-day ratios explicitly restricted to measured animals. A positive offset
makes model T/C a ratio on the V + offset scale and doubling time a time for
V + offset to double. Neither is an unqualified raw-volume quantity.

## Extension rule

Add a module-specific extractor, source/state tests, and Skill guidance before
claiming another module is covered. Do not infer new scientific conclusions
from an unfamiliar field name. Config defaults and numerical results remain
owned by the statistical modules. New inference methods require their own
numerical agreement, worked-example and calibration evidence.
