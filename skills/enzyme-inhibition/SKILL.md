---
name: enzyme-inhibition
description: Estimate Ki and describe the inhibition mechanism from enzyme initial rates measured across substrate and inhibitor concentrations - global fit of competitive, noncompetitive, uncompetitive or mixed (alpha) models with profile intervals, a declared mechanism, tight-binding and range gates, and summaries across independent experiments; plus Morrison tight-binding fits of inhibitor titrations when Ki is near the enzyme concentration. Use for reversible inhibitors of a single-substrate enzyme; not for ordinary IC50 curves at one substrate concentration (dose-response) or plain Michaelis-Menten without inhibitor (nonlinear-models).
---

# Enzyme inhibition: Ki and mechanism

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Read the [input contract](references/input-and-model.md). Use
`analysis_type: enzyme_inhibition`. Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before analysis

1. **Initial rates.** Each rate must come from the linear early phase (typically
   under 10% substrate conversion, no product inhibition). Ask how rates were
   obtained if it is unclear. Time-dependent or irreversible inhibition (rates that
   fall with preincubation time) needs a different analysis; do not fit it here.
2. **Declare the mechanism before the results.** Use prior knowledge (structure,
   binding site, earlier studies); declare `mixed` when it is not known. The runtime
   also fits the other models to describe compatibility, but never use that
   comparison to pick the reported mechanism after seeing which fits best. A
   declared mechanism that the data reject against the mixed model withholds Ki.
3. **Tight binding.** Ask for the total (active-site) enzyme concentration. If Ki
   comes out below about 10 times that concentration, free inhibitor is depleted and
   the standard equations are invalid: Ki is withheld. Use the Morrison option
   below instead. When the enzyme concentration is unknown, say tight binding cannot
   be excluded.
4. **Design.** Substrate should span roughly 0.2-5 Km and inhibitor up to several
   times Ki, with an uninhibited series in every experiment. Independent experiments
   (separate enzyme preparations or days) are the unit for summaries; replicate wells
   are not experiments.
5. **IC50 is not Ki.** An IC50 at one substrate concentration depends on that
   concentration and the mechanism; convert only with the matching equation and its
   assumptions, and say so.

## Tight binding: Morrison (0.14.3)

When the inhibitor is potent relative to the enzyme concentration (Ki_app within
about 100-fold of total enzyme, or the standard fit withheld Ki for tight binding),
add a `tight_binding` block and analyse an inhibitor titration at one substrate
concentration per experiment (see the contract).
- `enzyme: fixed` treats the declared total enzyme as fully active; use it only
  when active-site titration supports that. Otherwise `enzyme: estimate` fits active
  enzyme; Ki was biased when inactive enzyme was fixed as active in calibration.
- Ki needs the mechanism and, for competitive or uncompetitive inhibitors, Km with
  its source (Km is treated as known). For `mixed`, only Ki_app at that substrate
  concentration is reported. One substrate concentration cannot establish the
  mechanism; it must be declared from other evidence.
- Ki_app far below the enzyme concentration (titration regime) is poorly
  determined; recommend a lower enzyme concentration rather than quoting it.
- Do not report the IC50 as a potency in a tight-binding assay: it depends on the
  enzyme concentration (IC50 = Ki_app + Et/2).

## Run and report

Report the declared mechanism and why, Ki (and Km, Vmax; alpha for mixed) with
profile intervals per experiment and the geometric mean across experiments, the
model comparison as description only, and every diagnostic or withheld item. An
open alpha interval means the second inhibition component was not detected, not
that it is absent. Respect `must_mention` and failing items in
`interpretation_facts.json`.

## Evidence

See [the 0.14.2 record](../../validation/RELEASE_0.14.2.md): GraphPad's published
competitive example, a base-R `nls`/profile oracle and the registered calibration.
Morrison: [the 0.14.3 record](../../validation/RELEASE_0.14.3.md) (base-R oracle with a
failed first attempt kept, registered calibration, unmet published-example gate).
