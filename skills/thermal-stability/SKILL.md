---
name: thermal-stability
description: Fit apparent melting temperatures (Tm) from nanoDSF ratio, intrinsic fluorescence, extrinsic-dye DSF or CD thermal unfolding curves with a declared number of two-state transitions (e.g. IgG CH2, Fab, CH3), sloping baselines, profile-F Tm intervals, derivative inflections for comparison with instrument software, replicate summaries and delta-Tm between formulations or variants. Use for developability and formulation screening; not DSC heat capacity, aggregation onset (Tagg) from scattering, or isothermal stability.
---

# Thermal stability (apparent Tm)

Read the [runtime procedure](../agentic-prism/references/runtime.md) and the
[input and model contract](references/input-and-model.md). Implemented in 0.13.1.

## Before fitting

1. **Technique and readout.** nanoDSF ratio, intrinsic fluorescence, extrinsic
   dye (SYPRO-type) or CD are fitted as unfolding curves. DSC (excess heat
   capacity) needs a peak model and is refused; light scattering measures
   aggregation onset (Tagg), not unfolding, and is refused. Dye DSF can show a
   post-transition decrease from dye quenching or aggregation; restrict
   `model.fit_range_c` before fitting rather than after a poor fit.
2. **Number of transitions** (`model.transitions`, 1–3) comes from the molecule
   and prior data (an IgG1 often shows CH2, Fab and CH3), recorded in
   `transitions_source`. Do not add transitions to absorb residual structure in
   this curve; if residuals show unexplained structure, report it.
3. **Conditions.** Declare scan rate, buffer and reversibility. Antibody
   unfolding is usually irreversible and scan-rate dependent: Tm is an apparent
   value and dH is a shape parameter. Compare only curves measured under the
   same conditions.
4. **Replicates.** Capillaries filled from one mix are technical repeats, not
   independent replicates. Declare `replicates.independent_unit: replicate_id`
   only for separate preparations or runs; comparisons need it.

## Run and report

Run `analyze`, `verify`, then read `interpretation_facts.json`. Report each
transition's apparent Tm with its profile-F interval, the derivative inflection
alongside (descriptive), replicate means with t intervals, and delta-Tm with
the Welch interval and Holm-adjusted p. Withheld transitions (baseline not
covered, overlapping transitions, amplitude near noise, open interval) are
reported as withheld, not as numbers. Assigning a transition to a domain needs
separate evidence (fragments, mutants or literature).

## Evidence (0.13.1)

Fits agree with base-R `nls(algorithm = "plinear")` on the same model (7
fixture curves; SSE within 2.4e-13 relative, Tm within 1.4e-7 °C). Each curve is
refitted with one fewer transition; unless the declared number wins (F test,
p < 0.01), every Tm of that curve is withheld (`declared_transitions_not_supported_by_data`),
and a transition next to an unresolved one is withheld too.

Two registered calibrations are kept. The first missed coverage for
multi-transition curves (0.914–0.928) and let overlapped transitions through;
after the structure test, neighbour gate and multistart profiles, the second
covered 0.943–0.963 except the noisy third transition (0.932, a retained miss).
0.9% of overlapped second-transition estimates still passed every gate and were
mostly wrong: treat Tm values of transitions closer than about two widths as
unsupported. The gates withhold more often (about 72% of noisy three-transition
curves reportable). The published worked-example gate is unmet; live-agent misuse
scenarios passed ([live-agent results](../../validation/agent-scenarios/0.13.1/RESULTS.md)).
See [the 0.13.1 record](../../validation/RELEASE_0.13.1.md).
