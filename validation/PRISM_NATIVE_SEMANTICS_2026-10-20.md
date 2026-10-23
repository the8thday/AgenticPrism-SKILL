# Native Prism structural compatibility and exclusion flags — 2026-10-20

Local unreleased extension of [the initial reader](PRISM_IMPORT_2026-10-19.md).
Runtime remains 0.14.4; extraction schema is now `prism-extraction/2`.
Earlier /1 artifacts and evidence are not rewritten; /1 manifests remain
verifiable. No new statistical estimator or automatic filtering is introduced.

## Changes and scope

- Scoped `1-2-0` profile: `text`, `y_single`, `y_replicates`, plus explicit X
  series. Other dataset storage formats in that version remain unknown.
  `1-4-1` keeps its previously checked storage roles.
- An X dataset with `format=series`, one `SeriesReplicate` and finite numeric
  `startValue`/`interval` is a virtual column in both profiles. Its metadata is
  preserved; it consumes no CSV field and generates no new values. Missing,
  boolean or malformed definitions do not silently shift subsequent Y columns.
- Native `DriverReplicate` flags now yield `excluded_by_row` at absolute CSV row
  positions. Common EXCLUDED with no range entry, and explicit range EXCLUDED,
  are decoded. Absence of EXCLUDED is a stored-bit observation, not proof of
  inclusion in a historical analysis. No source values are dropped or edited.
- Common EXCLUDED plus a range without it has ambiguous override precedence,
  so the affected cells remain null. Conflicting overlapping range flags also
  remain null. Invalid records/ranges make a column's flags unknown with an
  issue. Unfamiliar flags (including censoring) stay raw metadata.
- Absent CSV fields differ from stored empty strings. Sparse CSVs cannot expand
  into an unbounded rectangular flag matrix; the annotation limit retains all
  raw values while flag interpretation remains unknown. Reports/facts include
  marked, unmarked, unknown and unmapped stored-cell counts.

These are structural reading rules, not official native schema definitions.
The code does not mark every possible layout in a recognized format version as
supported. File-level provenance and raw metadata remain authoritative.

## Public evidence

All repository checks use pinned public source files or independently invented
small examples. No private project, filename, label, dataset ID, measurement,
private hash or private-output path is included in this record or its artifacts.
Local acceptance work on confidential inputs is kept outside the repository and
is not a redistributable fixture or a substitute for public/GUI evidence.

| Gate | Observed scope/status |
|---|---|
| Established implementation | Storage fidelity compared with CSV/XML primitives; synthetic tests exercise explicit native structures. Independent native Prism GUI/schema oracle remains UNMET. |
| Public worked example | Existing public `1-4-1` project and 24 XML files remain readable without changing stored text. `1-2-0` is covered in repository tests by invented fixtures, not claimed as another public native project or full-version benchmark. |
| Calibration | Not applicable: deterministic extraction, no estimator added. This does not establish all-version compatibility. |
| Live-agent misuse | UNMET: no independent live-agent scenario run. Deterministic tests and Skill wording are not a behavioral-validation pass. |
| Interpretation facts | Stored report counts and limits are copied into facts, with report hash. Exclusion summaries stay descriptive and never perform inference or filtering. |

## Checks

- Focused parser suite: **68 passed** (`tests/test_prism_import.py` plus
  `tests/test_prism_native_semantics.py`). The 26 additional tests use entirely
  invented values/IDs; they cover both version profiles, virtual X recognition,
  absolute row ranges, replicate ranges, default/common and per-row EXCLUDED,
  conflicting overlaps, missing/invalid flags, row offsets, absent versus blank
  fields, annotation limits, unknown versions/layouts, facts and old manifests.
- Reproducible artifact check:
  `.venv/bin/python scripts/validate_prism_import.py --output NEW_DIRECTORY`.
  [Current artifacts](prism_native_semantics_2026-10-20/validation.json) retain
  the same 25 public files, 66 tables and 7,571 stored fields. New metadata and
  flags intentionally change extraction JSON; the source data is unchanged.
- [Legacy check](prism_native_semantics_2026-10-20/legacy.json): all 119
  pre-existing non-CLI Python files remain byte-identical to pre-import baseline
  `b99ed678a18c5bae8c6e8973563fc76c0977fc5f`. Six selected existing analysis CLI
  fixtures produce the same 61 scientific artifacts. No unchanged numerical
  benchmark or calibration is rerun.
- Changed Skills pass skill-creator `quick_validate.py`; `git diff --check`
  passes.
- Final full suite: **1,047 passed**, 8 pre-existing `Pandas4Warning` warnings in
  `src/agentic_prism/ada.py:209-210`, 472.25 seconds. This full run uses the final
  parser and all 68 parser cases. No parser code changed after it started.
- Both old /1 and new /2 public demo manifests verify; all 43 demo tables retain
  identical raw content (native CSV strings; XML values/metadata). The reader CLI
  also works with site packages disabled (`python -S`).
- A local privacy check confirmed private source names, hashes and object IDs are
  absent from changed deliverables. The private file and its extraction remain
  outside the repository; no original input was modified or uploaded.

## Meaning of exclusion versus storage encoding

GraphPad explains that excluded cell values remain in a data table while being
omitted from analysis and graphing. It separately documents choosing which
datasets a particular analysis uses. These support keeping the stored marker
distinct from analysis selection. [Cell exclusion](https://www.graphpad.com/guides/prism/latest/user-guide/excluding_values.htm),
[dataset selection](https://www.graphpad.com/guides/prism/latest/user-guide/using_excluding_points.htm).

Those guides do **not** document `commonAttributes`/`cellAttributes` ZIP encoding
or establish replacement precedence. The reader explicitly leaves ambiguous
cases unresolved instead of inventing that precedence. Consult the
[updated contract](../skills/import-prism/references/input-and-output.md) for the
output representation and boundaries.
