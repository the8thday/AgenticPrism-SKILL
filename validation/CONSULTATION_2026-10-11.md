# Optional analysis consultation — 2026-10-11

Scope: local Skill/conversation update on the 0.13.4 numerical runtime, not a new
numerical release. Baseline commit: `ee85760612a8f03bd2d819a55dba257c76f376dd`.

The router and all direct specialists reference one shared entry contract. A
new analysis offers discussion or direct execution unless already chosen;
choices survive handoffs and continuation. Opt-in `analysis-consultation`
uses actual data and adaptive questions to produce a brief. Both paths include
a proportional assessment of design, estimand, observed features, suitable
methods and limits. Skipping the discussion does not supply missing scientific
facts. Verification, rendering and capability questions do not trigger the menu.

The menu is a host-agent conversational choice, with numbered chat fallback;
there is no runtime popup, CLI flag, new analysis type or new config field.
Conversation persistence relies on the host agent's available context, not a
new package-level state store. The brief is kept outside scientific run folders.

## Five evidence gates

| Gate | Evidence / applicability |
|---|---|
| Established implementation comparison | Not applicable to this change: no estimator or numerical method added. Exact runtime/representative-output comparisons recorded separately |
| Published worked example | UNMET for the consultation workflow; GraphPad/NIST links are conceptual guidance, not validation of this Skill |
| Calibration | Not applicable: no statistical decision threshold, test or interval introduced. No numerical calibration claim for agent judgment |
| Live-agent misuse / interaction | 6/6 scenarios and 8/8 conversation turns passed the [fixed criteria](consultation_20261011/CRITERIA.md), scored by the implementer from actual responses and artifacts; plain-chat choice only |
| Interpretation facts | Existing specialist facts and withholding stay authoritative; no new facts extractor or inference in the brief. Representative facts hashes are compared to the baseline |

## Checks

Behavioral prompts use hand-constructed synthetic files under
`fixtures/consultation_20261011/`, not biological data. Fresh independent contexts
receive the actual request and inputs without expected answers. Scoring is by
the implementer. The plain-chat fallback is tested; a particular host's graphical
question widget is not validated by these tests.

| Scenario | Observed result |
|---|---|
| New request without preference | Offered both modes and waited; after the user selected direct analysis, ran the declared paired test without another menu |
| Discussion with unclear replicate design | Inspected real rows, asked three material questions and explained their consequences; after clarification retained the complete treatment/plate-confounding limit and saved a brief |
| Direct specialist with explicit skip | Assessed the data, executed the specified paired analysis and verified report, configuration, results and facts without an interview or extra approval |
| Skip with a literal `<LLOQ` qualifier | Identified the affected record, rejected zero substitution for a log ratio and asked a narrow data question without forcing the full consultation |
| Continue an existing discussion | Preserved known choices/design, recommended an additive paired AU estimand and saved a sourced brief without numerical execution |
| Capability explanation only | Did not offer the menu or run an analysis; distinguished time-curve AUC from PK/NCA |

The two authorized runs each used six complete pairs and produced the specified
paired-t result and a report; scientific-artifact hashes were independently
rechecked. The four discussion/limited/capability-only workspaces contain no
numerical run manifest. These scoped observations do not establish performance
on every design or prove all possible agent behavior.

- [Behavioral results and artifact pointers](consultation_20261011/behavioral_results.json):
  six fresh contexts, two follow-up turns; no expected answers supplied.
- [Evidence layout](consultation_20261011/README.md),
  [initial prompts](consultation_20261011/scenarios.json),
  [follow-up prompts](consultation_20261011/followups.json), and
  [frozen input/Skill hashes](consultation_20261011/inputs.json).
- [Skill validation](consultation_20261011/skill_validation.json): all 28 folders
  passed skill-creator `quick_validate.py`.
- [Compatibility](consultation_20261011/compatibility.json): all 105 runtime source
  files are byte-identical to the baseline; three representative configurations
  (paired comparison, ICx, paired AUC) produced 28 byte-identical scientific
  artifacts, including interpretation facts. Reproduce with
  `scripts/validate_consultation_20261011.py --output <new-record.json>` under the
  collection `.venv/bin/python`.

Full numerical tests, old benchmarks, calibrations and package rebuilds were
not rerun for this instruction-only update. Source, runtime schemas, version
declarations, existing release records and saved runs remain unchanged.

The initial staged whitespace check flagged CRLF in the three new CSV fixtures.
Scoped CSV attributes preserve the tested bytes and recognize CRLF; the input
hashes and behavioral evidence were not rewritten to change line endings.
