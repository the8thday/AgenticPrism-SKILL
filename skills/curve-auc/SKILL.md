---
name: curve-auc
description: Compute area under time curves per biological unit over a declared interval and baseline, with independent Welch or paired within-unit comparisons, Holm p values and Bonferroni intervals. Use for killing, confluence, cytokine, body-weight or tumor-volume time courses. Requires explicit dropout/pair policies; not PK/NCA or a correction for outcome-related dropout.
---

# Curve AUC

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[contract](references/input-and-model.md). Implemented in 0.13.2.

## Before computing

1. **Unit.** One curve per independent unit (animal, donor, independent culture). Technical wells of one
   culture are averaged upstream or analysed as a nested design, not entered as units.
2. **Interval and baseline** are declared before looking at the curves: `auc.interval`, and `baseline`
   none, first_value (area relative to the starting value) or declared_constant. Net area: segments below
   the baseline count negative. Never pick the interval that maximizes a difference.
3. **Units that end early.** Declare `incomplete_policy` (withhold_unit or common_interval) with a
   rationale. If dropout is related to the outcome (tumor-volume humane endpoints, deaths), any AUC
   comparison is biased; say so and suggest tumor-growth or time-to-event instead. The run flags groups
   whose withholding differs by more than 20 percentage points.
4. **Pairing.** For the same donor or culture under several treatments, declare
   `comparison.design=paired_t` and the pair policy/rationale in the
   [paired contract](references/input-and-model.md#paired-comparisons-0134).
   Pair by unit ID; never rename shared donors into independent IDs. The default
   independent-group Welch calculation remains available for different units.

## Run and report

Run `analyze`, `verify`, read `interpretation_facts.json`. Report the interval actually used, the
baseline rule, n per group with any withheld units, mean AUC with intervals, and the Welch differences
with the declared family. AUC summarises a window; it does not show when curves diverge, so show the curves.

## Evidence (0.13.2)

Independent base-R trapezoid and `t.test` agree within 1.8e-13 on four datasets. Calibration (2000 per
row): Welch coverage 0.949 with unequal variances, Holm family-wise error 0.039, simultaneous coverage
0.960. With outcome-related dropout under withhold_unit, coverage fell to 0.821 (descriptive), which is
why the flag exists. Live-agent misuse scenarios (interval shopping, dropout) passed. See
[the 0.13.2 record](../../validation/RELEASE_0.13.2.md).
