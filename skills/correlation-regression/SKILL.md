---
name: correlation-regression
description: Analyze independent paired measurements using Pearson or Spearman correlation, linear regression, Deming, Passing-Bablok, or Bland-Altman agreement. Distinguishes association, constant bias, proportional bias and agreement; requires declared measurement-error and acceptance assumptions.
---

# Correlation, regression and method comparison

Read [input contract](references/input-and-model.md), then use the collection
runtime described in [runtime](../agentic-prism/references/runtime.md).

Establish one pair per independent subject, the measurement scale, readout,
range and intended estimand. Wells are not subjects. Do not turn correlation or
R squared into agreement. Do not remove outliers automatically.

Deming needs Var(error Y)/Var(error X), declared before analysis with a source;
never tune it to agreement. Passing-Bablok requires a plausible positive linear
relation and nonnegative measurements; lead with a failed CUSUM linearity
check and withhold its coefficients. For agreement, declare acceptance limits
with a source before seeing results. Bland-Altman normal quantile intervals
assume independent normal differences with constant spread; a trend against
pair mean is a diagnostic, not a correction.

Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`.
Read `interpretation_facts.json`: reportability takes precedence over optimizer
or workflow completion. Report constant and proportional bias separately;
carry forward every calibration miss and unmet evidence gate. Rendering only
reads saved results. Consult [0.11.1 evidence](../../validation/RELEASE_0.11.1.md).
