---
name: first-in-human-dose
description: Calculate first-in-human starting-dose inputs from declared nonclinical data - NOAEL to human equivalent dose and maximum recommended starting dose (FDA 2005 body-surface-area table, body-weight exponent or mg/kg for large proteins), MABEL / pharmacologically active dose from receptor occupancy or effect concentrations, and allometric projection of human antibody PK for candidate doses - with sensitivity tables. Use when a project needs the numbers behind a FIH dose rationale; it never chooses or approves the dose.
---

# First-in-human dose calculations

For a new analysis, follow the shared [optional discussion and data-assessment
entry](../agentic-prism/references/analysis-entry.md). Preserve a choice already
made; do not repeat the menu on handoff.

Use `analysis_type: fih_dose` (0.15.2), a design-only analysis: every input is a
declared number with a source, and there is no data file. Read the
[input contract](references/input-and-model.md). Locate the runtime with the
[ordered procedure](../agentic-prism/references/runtime.md).

## Before calculating

1. **Modality decides the scaling.** Small molecules: body-surface-area divisors
   (FDA 2005 Table 1) by default. Proteins over 100 kDa given parenterally
   (antibodies): mg/kg, not body surface area; the runtime refuses BSA for them.
2. **Species.** Use the most sensitive species unless a named species is
   scientifically more appropriate (for an antibody, often the only
   pharmacologically relevant species); say why. Never pick a species because it
   gives a higher dose.
3. **Safety factor.** The guidance default is 10. Anything lower needs a stated,
   sourced reason; higher may be needed for steep dose-response or novel targets.
4. **MABEL** for agonists, immune activators or high-risk targets: declare each
   pharmacology target (receptor occupancy from a human KD, or an effect
   concentration such as an EC20 from human cells) with its source, and the human
   central volume. A KD in nM needs the molecular weight. Report every target;
   the lowest gives the MABEL.
5. **Allometry** needs animal PK parameters per kg in the stated units (convert
   first), the animal weight and exponents with a source (cynomolgus to human
   IgG: CL and Q 0.85, volumes 1.0 are common). Do not change exponents to reach
   a preferred dose; the sensitivity table shows their effect.

## Run and report

Report each route to a number with its inputs: HED per species and the selected
species, safety factor and MRSD (mg/kg and mg); each MABEL target with
concentration (ug/mL and nM) and dose; the MRSD/MABEL ratio; projected human Cmax,
AUC, half-life, time above target concentrations and RO for candidate doses; and
the sensitivity tables. Never call any value the FIH dose or approved dose: the
starting dose is an integrated decision under FDA (2005) and EMA (2017) guidance.
Respect `must_mention` and failing items in `interpretation_facts.json`.

## Evidence

See [the 0.15.2 record](../../validation/RELEASE_0.15.2.md): the published
NOAEL-to-HED examples reproduced, the base-R oracle, the unmet MABEL and
allometry worked-example gates, and the live-agent scenarios.
