---
name: nonlinear-models
description: Fit a declared common nonlinear model per curve - one-phase decay or association (half-life, half-time), two-phase decay, exponential growth (doubling time), Michaelis-Menten (Vmax, Km) or logistic growth - with profile-F intervals for every parameter, runs test, and gates for plateau support, two- vs one-phase structure, saturation and asymptote support, plus summaries across independent experiments. Use for serum or thermal stability decay, internalization or uptake kinetics, enzyme or cleavage kinetics and growth curves. Not for sigmoid dose-response (dose-response), binding KD or kon/koff (binding Skills), or PK (not covered).
---

# Common nonlinear models

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[contract](references/input-and-model.md). Implemented in 0.13.3.

## Before fitting

1. **Declare one model** with a rationale from mechanism or prior data (`model.rationale`). Never fit
   several models and keep the best-looking one; if the user wants a model comparison, it must be a
   predeclared question (two- vs one-phase is tested automatically and only withholds).
2. **Route elsewhere when another Skill owns the question:** sigmoid concentration-response to
   dose-response, equilibrium titrations to equilibrium-binding, sensorgrams to binding-kinetics,
   melting curves to thermal-stability.
3. **Weighting**: none, or 1/y² for roughly constant relative error (declared, fixed from observed y).
4. **Design support**: decay and association curves should run for at least three half-lives or the
   plateau is withheld; Michaelis-Menten substrate concentrations should go well beyond Km or Vmax and
   Km are withheld (and, rarely, still underestimated); logistic asymptotes need data past the inflection.
5. **Replicates**: `replicates.independent_unit: experiment_id` only for independent experiments;
   technical repeats of one sample are averaged within the experiment before the between-experiment summary.

## Run and report

Report each curve's parameters with profile-F intervals, derived half-life / doubling time with the
interval mapped from the rate, the gates that withheld anything, and the runs-test flag. Summaries
across independent experiments use geometric means for rates and derived times. Withheld parameters are
reported as withheld, not as numbers.

## Evidence (0.13.3)

Checked against base R only (stats::nls with R's own starting values - selfStart models, a log-linear
start or an R-side grid - and confint profile intervals), following the 0.13.3 R-package audit: 18/18
datasets agree (SSE about 1e-14 relative, estimates 1.3e-7, interval ends within 1.5e-4 of the width).
Calibration (1000 per row): 11/11 bounded coverage rows passed (0.942-0.967), including 1/y² weighting
under multiplicative error. Stress rows: plateaus of curves observed for under one half-life were always
withheld; with substrate below Km, 99.8% of fits were withheld but 2 of 1000 slipped through with Km
underestimated. Published worked example and live-agent gates are pending. See
[the 0.13.3 record](../../validation/RELEASE_0.13.3.md).
