---
name: drug-combination
description: Analyze two-agent inhibition matrices with named Bliss, Loewe, HSA and ZIP reference models, declared response normalization, independent-matrix uncertainty and model-disagreement diagnostics. Refuses synergy claims from unreplicated matrices or incompatible Loewe effect ranges.
---

# Drug combinations

Read [input contract](references/input-and-model.md) and the shared
[runtime](../agentic-prism/references/runtime.md). Establish the agents, doses,
cell line/condition, readout and normalization. All responses must share a
percent-inhibition scale. Technical wells are not independent experiments.

Use `analysis_type=drug_combination`. Complete matrices require single-agent
zero-dose axes and at least four positive doses per agent. Never pool different
conditions. Bliss assumes independent effects; HSA compares with the stronger
single agent; Loewe requires compatible invertible effect ranges; ZIP compares
conditional fitted response slices with a fitted independence reference.

Declare Loewe compatibility with a source before analysis. Unsupported
single-agent or conditional 4PL curves withhold dependent scores. Do not use
an extrapolated audit estimate as a reportable synergy result.

Equal-weight independent matrices support Student t intervals. An unreplicated
matrix yields descriptive scores only. Report all requested models and their
assumptions, including disagreement in sign. Never choose the favorable model
without disclosing others. A score is dose-grid- and assay-specific and does
not establish mechanism or clinical combination benefit.

Run `agentic-prism analyze --config CONFIG --output NEW_DIRECTORY`; use
`interpretation_facts.json`, failing items and retained calibration misses in
all narratives. See [release evidence](../../validation/RELEASE_0.13.0.md).
