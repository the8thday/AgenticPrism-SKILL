# Prism stored-data extraction — 2026-10-19

Local unreleased addition on development runtime 0.14.4. No version bump,
publication, new estimator or automated analysis mapping. Baseline:
`b99ed678a18c5bae8c6e8973563fc76c0977fc5f`.

## Implemented behavior

`inspect-prism --input FILE`, `import-prism --input FILE --output NEW_DIRECTORY`,
and `verify-prism --run DIRECTORY`. Runtime parsing uses Python's standard
library; the CLI loads numerical dependencies only for their existing commands.

Input values remain strings and retain their positions, including explicit blanks,
absent trailing values, excluded values, row titles, X and repeated subcolumns.
Native table display format and dataset storage format are separate fields.
For example, the native CV table's dataset declares SD storage and its first
mean/SD remains `39.5` / `7.77817409999999843`; no CV is calculated. High/low
distances remain distances. An implicit X series stays a start/interval definition;
no generated X values are inserted. XML alternative X/date representations are
separate, not silently substituted. Native exclusion/censoring flags are retained
as metadata, not treated as interpreted include/exclude decisions.

All native CSV members are inventoried, including saved results and unlinked
objects. All valid JSON metadata and member hashes are saved. XML Table/HugeTable
contents and other elements are preserved; opaque content remains recoverable
from a byte-identical original. Unsafe/malformed essential input fails without
creating an output directory. Recoverable per-object problems appear as issues.

## Evidence and limits

| Gate | Scope/status |
|---|---|
| Established implementation | Bounded structural checks against Python CSV/XML primitives and pinned real project content. This is not an independent native Prism GUI oracle; full equivalence UNMET. No upstream parser code copied. |
| Public worked example | 25 pinned public files: 1 native `.prism` and 24 `.pzfx`; 66 table objects, 7,571 stored fields (includes alternative X representations and saved results). Every extracted cell text and XML cell attribute checked against the stored source. This is data-extraction evidence, not reproduction of a published statistical analysis. |
| Calibration | Not applicable to deterministic extraction; no estimator or inferential claim added. Cross-version completeness is UNMET, not inferred from passing engineering tests. |
| Live-agent misuse | UNMET: no separate agent scenario run. Skill instructions and deterministic tests are not a behavioral-validation pass. |
| Interpretation facts | Extraction-only facts copy counts/limits/issues from the hashed import report; tested for source identity. No model interpretation is introduced. |

The only native project uses formatVersion `1-4-1`, created in Prism 10.4.0.527
and modified in 10.4.1.532 on Mac. Unknown native format versions retain raw
values while withholding unverified column roles. This is not a Windows/macOS
or Prism 10/11 compatibility matrix. The paired `.prism`/`.pzfx` example has
different save dates and is not claimed to be a lossless round trip.

Uninterpreted: native graphs/layouts/binary payloads and `.dt` semantics; XML
compressed Template; native exclusion/censoring semantics; unknown metadata.
No project rendering, analysis recomputation, `.pzf` binary support, official
MCP integration or automatic data conversion. A missing/unknown XML TableType
is unclassified, not silently promoted to experimental observations.

## Checks

- Focused suite: `tests/test_prism_import.py`, 42 passed on the final parser. It covers all
  public source hashes and stored values; CV versus SD, bounds versus distances,
  all 14 native input/7 result/8 unlinked objects, implicit series, repeated titles,
  XML exclusions/date/X errors/HugeTable, namespaces, subcolumn titles, missing
  versus blank versus literal NA, uneven lengths, future schemas, dimension
  mismatches, unsafe ZIP names, duplicate members/JSON keys, UTF-8/UTF-16 DTD
  refusal, size limits, output preservation, facts and hash tampering.
- First focused attempt: 39 passed, 2 failed due to inaccurate fixture expectations.
  The native dataset's format was actually `y_sd`, despite the table being `y_cv`;
  `table_display_format` was added to expose both. The XML exclusion fixture
  contained exactly `90.`, `60.`, `2.` (not just integer strings 90/60). Expectations
  were corrected from the original bytes; stored values were never changed.
- [Reproducible corpus and artifact checks](prism_import_2026-10-19/validation.json):
  `scripts/validate_prism_import.py --output NEW_DIRECTORY`. Demo outputs are
  retained for [native](prism_import_2026-10-19/demo_prism/import_report.json) and
  [XML](prism_import_2026-10-19/demo_pzfx/import_report.json), including snapshots,
  every table, metadata, interpretation facts and manifests.
- [Legacy check](prism_import_2026-10-19/legacy.json): 119 pre-existing Python
  implementation files (all except CLI) byte-identical to baseline. Six selected
  existing analysis CLI fixtures (equilibrium, paired comparison, dose response,
  ELISA, PCA and clustering) produced 61 byte-identical scientific artifacts.
  This is a scoped CLI regression check, not an exhaustive all-method rerun.
  No unchanged numerical benchmark or calibration was rerun.
- Both changed Skill folders passed skill-creator `quick_validate.py`. The first
  router description exceeded 1,024 characters after adding the new route; it
  was shortened, then validated successfully.
- Standard-library-only CLI smoke passed with `PYTHONPATH=src python -S -m
  agentic_prism.cli inspect-prism` on the native fixture (site packages disabled).
- A final conservative guard removes speculative dataset-level CV/high-low role
  aliases: those strings were observed as table display formats, not dataset
  storage formats. A synthetic unsupported dataset-format test now requires
  `unknown`, while preserving every stored value. Initial successful artifacts
  remain in `prism_import_2026-10-19_initial`; final artifacts are in
  `prism_import_2026-10-19`. All 43 exported demo table JSON files are byte-identical
  between these checks. Parser/test identity hashes distinguish the attempts.
- Full suite: `.venv/bin/python -m pytest -q` — **1,020 passed**, 8 existing
  `Pandas4Warning` warnings in `ada.py:209-210`, 457.48 seconds. This run collected
  the first 41 importer tests. The final conservative storage-format guard and
  one extra case were checked separately afterward: **42 importer tests passed**
  on the final parser (0.69 seconds). Current collection contains 1,021 tests;
  an all-1,021 final-snapshot rerun is not claimed. Existing numerical source
  remains unchanged and the final six-fixture CLI check passed again.
- New documentation links resolve; `git diff --check` passed. No install,
  publication or historical release-artifact modification was performed.

## Sources

Fixture provenance, commits, hashes and license notices:
[fixtures/prism_import](../fixtures/prism_import/README.md). Local research evidence:
[format research](../plans/prism-file-research-2026-10-19/README.md).
Data contract: [import-prism](../skills/import-prism/references/input-and-output.md).
The reader was written for this repository; it does not wrap R packages or copy
the unlicensed Python converter investigated during research.
