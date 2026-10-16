# Optional discussion before a new analysis

Apply this entry contract to a **new scientific analysis**, whether entered through
agentic-prism or a specialist directly. It is an agent conversation step, not a
CLI option or a new statistical method. Do not insert it into installation,
doctor, import-only conversion, verification, re-rendering, explanation of an
existing result, or development of the collection.

## Offer the choice once

If the user has not already chosen, begin with a short choice in their language:

> 分析前，是否先讨论一下研究目的、实验设计和数据适合的统计方法？
>
> 1. **先讨论，再分析**：结合数据，逐步澄清目标和关键设计问题。
> 2. **直接分析**：按你提供的描述开展分析，并保留必要的数据适用性检查。

Use the host's question/choice UI if available; otherwise present numbered
options in chat and wait for the answer. Do not claim the Markdown creates
buttons. Either choice is valid; no answer or elapsed time is not a choice.
While waiting, safe file inventory/inspection is fine; do not begin hypothesis
tests or model fitting before the choice is resolved.

An explicit request to discuss is already opt-in. “Skip discussion”, “直接分析”,
or equivalent is already opt-out: proceed without another menu. A bare request
to calculate an endpoint or run a named test does not by itself waive the choice.
Remember the choice for this study/task, its continuation and specialist
handoffs. Respect a standing preference declared by the user; do not re-ask
merely because another Skill is loaded. A user may change mode at any time.

- **Discuss:** follow [analysis-consultation](../../analysis-consultation/SKILL.md).
  It returns a brief and routes to the existing specialist; it does not replace
  that specialist's applicability contract.
- **Direct:** use the supplied description and actual data. Carry out the short
  assessment below and continue with the specialist when material inputs are
  available. Do not impose the full consultation questionnaire, a new approval
  step, or an extra “confirm execution” after an already authorized analysis.

## Agent assessment in both modes

Before fitting/testing, explain concisely what the data and description support:
the target quantity, independent unit and design; relevant observed features;
the method and its reason; and material limits or unresolved facts. A short
paragraph normally suffices in direct mode. Report evidence and practical
consequences, not hidden reasoning or a catalogue of every possible test.

Read available files before claiming to have checked them. Separate observations
(with file/sheet/column or row pointers), the researcher's statements and your
provisional assumptions. If only a description or summary is available, say so;
do not invent raw observations, distribution checks or figures. Discussion-only
work does not require installing or checking the numerical runtime first.

Check the relevant subset of source/preprocessing, format/units, independent vs
technical replicates, pairing/repeated measures, batch confounding, missingness,
zeros/BLQ, scale/variation and dose/time coverage. Follow the detailed
[data and design contract](../../analysis-consultation/references/data-and-design.md)
when interpreting these features. Do not choose a method solely from a normality
p value, select the smallest treatment p value, silently preprocess data, or
populate applicability flags from an example configuration.

Skipping discussion does not supply a missing experimental fact. If unknown or
contradictory information changes the estimand, independence, pairing, units,
scientific validity, or a required sourced acceptance limit, explain its
consequence and ask only the blocking question(s). Continue independent safe
inspection. When nothing essential is missing, proceed without further questions.

Keep discussion notes/briefs outside the numerical run directory. Specialist
outputs and verified interpretation_facts.json remain authoritative after
execution; an earlier proposed method or expectation cannot override them.
