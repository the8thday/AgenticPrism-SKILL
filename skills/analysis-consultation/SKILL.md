---
name: analysis-consultation
description: Discuss a researcher's analysis goals, experimental design and actual data before statistical analysis. Use when the user opts into pre-analysis discussion or asks what analyses their data support. Inspects provenance, independent units, pairing, distributions and limitations, then hands a concise analysis brief to the appropriate AgenticPrism specialist.
---

# Analysis consultation — 分析前讨论与数据评估

Use a collaborative, Socratic dialogue: connect a concrete observation to a
question that helps the researcher clarify a decision. Explain what the answer
would change. Do not quiz the researcher on statistical terminology or lead
them toward a desired significance result.

For a new analysis with no discussion preference, follow the shared
[entry choice](../agentic-prism/references/analysis-entry.md) first. A request to
discuss or invoking this Skill for consultation already selects discussion;
do not offer the menu again. If the researcher elects to skip or end discussion,
switch to the direct path in that contract. Required scientific inputs still
apply. Respect “discussion only”: do not start a numerical analysis merely
because a candidate method has been found.

## Start from the researcher's question and data

Read [the data and design contract](references/data-and-design.md) before
assessing the supplied inputs. Read accessible raw tables, metadata and study
notes; begin with the information already supplied. No data file is required
for an initial discussion, but describe uninspected features as unknown.
Neither a working runtime nor a canonical analysis CSV is a prerequisite for
this conversation. For unfamiliar or unsupported exports, identify what can
be read and ask for the needed export or column meaning without fabricating it.

Briefly restate the question in practical terms, such as a treatment-control
difference, a fold change, a relative IC90 or an effect over a declared time
window. Ask about the intended decision and what would count as a useful effect
when that information changes the analysis. “No significant difference” does
not mean equivalence; numerical fit does not prove an assay mechanism.

## Discuss the important uncertainty first

Show a relevant observation before asking about its meaning. For example:
“The same donor IDs occur under both treatments, with three wells per donor.
Do these wells come from matched aliquots of the same preparation?” Explain
that this affects pairing and the independent sample count.

Usually ask one to three related questions per turn, prioritizing what changes
the decision. This is not a fixed questionnaire. Read and acknowledge answers,
update the working brief, and do not ask again about resolved facts. Offer
plain-language choices when helpful and allow “unknown”. An unknown fact stays
unknown; suggest where it might be recovered or what limited analysis remains.

Adapt questions to the response. Discuss, as relevant:
- the biological question, target quantity and intended use;
- source and prior processing, units and the meaning of one row;
- independent units, randomization/allocation, pairing, repeats and batches;
- observed data features, missingness, scale, limits of measurement and coverage;
- method choices, what they estimate, assumptions and conclusions they cannot support.

Use descriptive summaries/plots when useful and trace them to actual inputs.
Do not run exploratory treatment tests to decide which test to recommend.
Provisional diagnostic model fitting is only for an explicitly requested
diagnostic assessment under a declared specialist model; label it and retain
the result. Do not automatically change a prespecified method after seeing data.
Record data-informed changes as exploratory/revised choices, not preregistration.

## End with a useful brief and handoff

When the material questions are resolved, synthesize a short **analysis brief**:

| Item | Content |
|---|---|
| Question | Target effect/quantity and intended decision |
| Evidence | Sources/files, format, units and known prior processing |
| Design | Independent unit, sample count, repeats/pairing, allocation and batches |
| Data assessment | Observed features, user statements, provisional assumptions, unresolved facts |
| Plan | Proposed primary method, justified alternatives/sensitivity checks and reasons |
| Limits | What can be concluded, what remains descriptive, missing data/metadata or unsupported methods |
| Handoff | Discussion choice, next specialist, decisions already settled and material open questions |

Keep this in the conversation; save `analysis_brief.md` in a new consultation
folder if the task uses a workspace. Preserve source paths and file hashes when
files are available; for prompt-only information cite the supplied description
without inventing a hash. Do not overwrite raw inputs or prior briefs; retain
corrections and their source. This is a planning note, not runtime
interpretation_facts.json or a validated statistical result.

Use the current [module registry](../agentic-prism/references/module-registry.md)
and the relevant specialist contract to separate **scientifically suitable**
from **implemented here**. Common handoffs are group-comparison, dose-response,
curve-auc and repeated-measures; other domains retain their own specialists.
If a suitable method is unavailable, identify the gap rather than forcing the
data into an available module or implementing a new method during consultation.

If analysis is already requested and all material choices are clear, continue
to the specialist without a second discussion menu or ceremonial approval.
If a choice still changes the question or validity, ask that specific question.
After execution, respect saved diagnostics, withholding and interpretation facts,
even when they contradict the preliminary assessment.

This Skill supplies conversation and assessment guidance. It adds no fitting,
automatic distribution classifier or numerical guarantees. Validation scope is
recorded in [consultation evidence](../../validation/CONSULTATION_2026-10-11.md).
