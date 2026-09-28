---
name: stability
description: Estimate shelf life for a quantitative quality attribute from long-term batch stability data using ICH Q1E linear regression, ordered poolability tests and mean confidence bounds. Requires declared specifications, direction, storage and applicability. Does not cover bioanalytical sample stability or accelerated Arrhenius models.
---

# Stability and shelf life

Follow [runtime discovery](../agentic-prism/references/runtime.md), then read the
[input contract](references/input-contract.md). This Skill belongs to AgenticPrism 0.10.0.
Use the shared CLI with `analysis_type=stability` and a new output directory.

## Gate before fitting

Establish the attribute, original unit, months since baseline, at least three
batches, and the long-term storage condition. Require the user to supply lower,
upper or both specification limits with a source, the direction with its source,
and the storage condition with its source. Do not infer specifications from the
observed range, infer storage conditions, or borrow another product's limits.
Confirm linearity, independent Gaussian errors and a common residual variance
with a scientific rationale. One attribute and storage condition per run; no
silent transformation, exclusions, smoothing, or merging storage conditions.
Repeated aliquots with correlated errors need another model.

## Analysis and interpretation

- ANCOVA tests slopes at alpha 0.25 first. Only if slopes can be pooled does the
  intercept test follow at 0.25. Never skip these tests to get a pooled result.
- Report the chosen pooled/common-slope/separate-batch model and both test
  states, including an intercept test that was not performed.
- Shelf life uses the first crossing of the mean bound: lower one-sided 95% for
  decreasing attributes, upper one-sided 95% for increasing attributes, or
  two-sided 95% for unknown/bidirectional change. These are not prediction or
  tolerance limits for individual units or future production batches.
- Separate batches use the pooled residual mean square under declared common
  variance; the minimum batch estimate governs. The overall product proposal
  must also respect every other critical attribute.
- Report supported months and governing batch. A crossing beyond the observed
  range is an audit quantity, not permission to extrapolate. Without complete
  user declarations, report the capped estimate and say it is not extrapolated.
- Q1E decision-tree declarations require their source and supporting evidence.
  This implementation covers statistical-analysis branches only: room storage
  without accelerated change by six months allows at most min(2X, X+12);
  refrigerated allows min(1.5X, X+6). Room storage with change at 3–6 months
  and no significant intermediate change allows min(1.5X, X+6). Early
  accelerated change, refrigerated accelerated change, frozen storage, or
  incomplete evidence does not allow extrapolation. X is the shortest observed
  batch duration. Supporting data, continued change pattern and a commitment
  study must be declared. The agent cannot declare them on the user's behalf.
- Inspect curvature and observed specification failures. Significant
  batch-specific curvature at 0.05 withholds the linear shelf-life result.
  A successful command does not establish an acceptable shelf life.

Read verified `interpretation_facts.json` before writing a narrative. Lead with
its reportable/withheld state, supported shelf life and whether it is extrapolated.
Include every `must_mention` item and relevant `failing_items`; cite saved artifacts.
Rendering changes style only.

## Evidence limits that must be disclosed

In 1,000-dataset Gaussian simulations, poolability-selected mean-bound coverage
was 92.2% for identical batches and 92.1% for near-poolable slopes; both missed
93.62%. Shelf-life non-overestimation was 93.5% for common slopes with different
intercepts and 90.7% for near-poolable slopes, also misses. Do not describe this
procedure as providing uniform 95% coverage after model selection. MCSE and all
rows are in [release evidence](../../validation/RELEASE_0.10.0.md).
The published Koleva example reproduces regression, ANCOVA and standard errors;
its printed shelf life uses a different interval and was not reproduced.
No Arrhenius, nonlinear degradation, random-batch prediction, or regulatory
approval is implemented.
