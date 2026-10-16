---
name: growth-rate-inhibition
description: Growth-rate inhibition (GR) metrics for cell proliferation and viability dose-response experiments - GR values from time-zero, treated and untreated counts, GR50, GEC50, GRinf, GRmax and GR_AOC per curve following Hafner et al. (2016) and the authors' reference code, with summaries across independent experiments. Use when comparing drug sensitivity across cell lines or conditions whose division rates or assay lengths differ; not for receptor binding or functional EC50 without cell growth (dose-response).
---

# Growth-rate inhibition (GR) metrics

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: growth_rate_inhibition`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Time-zero counts are required.** GR needs the cell count (or a readout
   proportional to it, such as nuclei or ATP calibrated to cell number) at the
   time of treatment, plus untreated controls at the endpoint. Without time-zero
   counts GR cannot be computed; do not invent x0 from seeding density unless it
   was measured, and say that relative viability IC50 then depends on division rate.
2. **Controls must grow.** Curves whose untreated control grew less than half a
   doubling are withheld (the GR denominator would amplify count noise); with less
   than one doubling GR is noisy (flagged). Report the control fold growth.
3. **One value per concentration.** Average technical wells by a declared rule
   first. Concentrations must be positive; untreated wells are the control
   column, not concentration 0.
4. **Independent experiments.** Summaries (geometric GR50, GRmax, GR_AOC with
   intervals) need at least two independent experiments (separate seeding days).
   Replicate wells or plates seeded together are not experiments.

## Run and report

Per curve report the fit type (sigmoid or flat), GEC50 and h_GR (potency and
steepness), GRinf (efficacy at high dose), GR50 only when it is reached inside the
tested range, GRmax and GR_AOC. GR near 0 means cytostatic; negative GR means cell
death (−1 = all cells killed). Do not extrapolate GR50 beyond the tested range.
GRmax follows the reference definition (minimum observed GR at the two highest
doses) and is biased low under noise; treat its interval as approximate. Respect
`must_mention` and failing items in `interpretation_facts.json`.

## Evidence (0.14.0)

See [the 0.14.0 record](../../validation/RELEASE_0.14.0.md): the authors' published
toy example reproduced (112 curves), the registered calibration with its GRmax
coverage miss, and the division-rate supplement.
