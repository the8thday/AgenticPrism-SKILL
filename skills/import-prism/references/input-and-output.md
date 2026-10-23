# Prism extraction contract — prism-extraction/2

## Inputs and scope

One local GraphPad project: ZIP with `document.json` (`.prism`) or XML rooted at
`GraphPadPrismFile` (`.pzfx`). Inspect writes nothing. Import uses a byte snapshot,
creates a new directory and never overwrites one. Malformed essentials fail with
CLI exit 2; local recoverable problems are saved as issues. Exit 0 means the
operation completed with disclosed limits, not every semantic detail is known.

Limits: source/member 64 MiB, ZIP expanded 256 MiB, 10,000 archive members,
2,000,000 stored cells/CSV records and 2,000,000 flag-matrix positions per table,
10,000 columns. DTDs, encrypted ZIP members,
unsafe paths and duplicate archive names are refused. No `extractall`, formula
execution, XML external entity loading or inference. `.pzf` is unsupported.

## Artifacts

| File | Meaning |
|---|---|
| `source.prism` or `source.pzfx` | Byte-identical source snapshot; authoritative for opaque/uninterpreted data and original encoding |
| `project.json` | Source hash, versions, table file index; all valid native JSON metadata and complete member hash inventory, or XML non-table elements |
| `tables/table_NNNN.json` | One table/CSV's stored values, stable extraction ID, source ID/path, column descriptions and raw metadata |
| `import_report.json` | Readable JSON inventory, classifications, issues and limits; no numeric results recomputed |
| `interpretation_facts.json` | Extraction-only counts, issues and disclosures copied from the hashed report; no scientific inference |
| `import_manifest.json` | Runtime/schema/parser identity, original SHA-256, hashes of generated artifacts; verify with `verify-prism`, not analysis `verify` |

All coordinates are zero-based. Output IDs use table order, never titles as file
names; duplicate titles remain distinct. Manifest paths are relative. Hashes detect
changes relative to this local manifest; `verify-prism` also rejects files the
manifest does not record. They do not authenticate the external
sender or validate Prism's opaque `.dt` fields.

## Native `.prism`

`raw_rows` holds CSV field strings in source order. `""` is a stored empty field;
`[]` is a blank CSV record; a short record has absent trailing fields. No padding,
type coercion or trimming. `null` for the whole `raw_rows` means data unavailable,
with an issue. `source` identifies the CSV inside the preserved archive.

`columns[i]` describes physical CSV field i: X/Y/row-title kind, dataset ID,
dataset `declared_format`, `table_display_format`, subcolumn index, `stored_role`,
metadata pointer and full replicate metadata. A width mismatch clears positional
associations and retains descriptions in `unmapped_column_candidates`.
Dataset `cellAttributes`/`commonAttributes` and ranges remain verbatim.
`excluded_by_row[i]` is the stored EXCLUDED bit at absolute CSV row i: true,
false or null. It does not remove values. `firstRow`/`lastRow` metadata never
shifts the CSV or the flag coordinates. Row offsets remain descriptive.

The decoder supports `DriverReplicate` in the two scoped profiles. With no
per-row entry, common EXCLUDED applies; with an explicit EXCLUDED entry the
cell is marked excluded. No common EXCLUDED and a range without EXCLUDED means
not marked. Common EXCLUDED combined with a range lacking it is ambiguous
(replacement versus supplementary flags), so the result is null. Conflicting
overlapping ranges also give null. Single rows and inclusive `a~b` ranges are
supported; invalid ranges/records make that column's flags unknown, with an issue.
Other flags are preserved without inventing their meanings. Censoring is not decoded.
Recognized tokens are those observed in the pinned public project (`CENSORED`,
`SECTION_TITLE`, `NUMERIC`, `USER`, `FLOAT`, `INTEGER`, `POSITIVE`, `NOTCONV`,
`NOEDIT`) plus `EXCLUDED`. Any other token could be an unfamiliar exclusion
spelling: it raises `unrecognized_cell_attribute` and makes the cells it covers
(all rows for a common token) null unless they carry explicit `EXCLUDED`.
The `EXCLUDED` token itself has no public native example in this repository.

`exclusion_status` is `decoded`, `partial`, `invalid_metadata`, `unverified` or
`resource_limit`. Absent physical cells get null; a stored empty field still has
its stored flag status. Excessive rectangular annotation of a sparse/ragged CSV
is refused without losing raw data (`excluded_by_row=null` for the column).
`exclusion_summary` counts marked/not-marked/unknown mapped cells and unmapped
stored cells. The report and interpretation facts also save aggregate counts.
False means no EXCLUDED marker, **not** proven inclusion in any saved analysis.
No rule propagates an X-cell flag to Y cells or infers selected datasets.

| Stored role | Meaning; values remain strings |
|---|---|
| `replicate`, `value`, `text`, `row_title` | Stored value/label; replicate columns do not establish biological independence |
| `mean`, `sd`, `sem`, `n` | Stored summary components, not reconstructed observations |
| `upper_distance`, `lower_distance` | Distances stored alongside mean; not absolute limits |
| `unknown` | Unsupported/unverified format or version; do not infer a role |

For the pinned native sample, a table displayed as `y_cv` has datasets stored as
`y_sd`. Its first stored mean/SD is `39.5` / `7.77817409999999843`; both strings
remain unchanged. The matching XML CV table can store a CV instead. Never apply
one format's role rules to the other format.

Native `SeriesReplicate` metadata goes in `virtual_columns`; no X values are
generated and it consumes no CSV field in either supported layout. Recognition
requires an X dataset, `format=series`, one `SeriesReplicate`, and finite numeric
`startValue`/`interval`. Invalid or unknown-version series definitions do not
silently consume or remove a CSV column; inconsistent mappings are withheld. Result tables stay
`analysis_result` and retain their `result_links`; all valid analysis sheet
metadata, parameters and saved result JSON are in `project.json`. CSVs without
sheet links are `unlinked_csv`, not assumed observations. Native role
interpretation is withheld for an unknown or non-string `formatVersion`; raw
extraction can still succeed. Non-string sheet/table uids or result `dataSheet`
links are issues (`invalid_sheet_uid`, `missing_table_data`,
`invalid_result_link`), not crashes or guessed links.

## XML `.pzfx`

Each column has its own `values`: `{row, raw_text, excluded, attributes, xml}`.
Missing trailing rows are absent list entries, never empty strings or recycled
values. Explicit `<d/>` has `raw_text=""`. Text such as `NA`, `1,25`, `90.` and
date strings stays exactly as XML text (XML entity decoding/line normalization
applies; the snapshot preserves original bytes). No numeric inference.

`Excluded=1/true` becomes true, `0/false` or absent becomes false under the
documented XML attribute default; invalid values become null with an issue.
No row is removed. `Title` under a `Subcolumn` is metadata, not a data row.
Rich titles and table-level `SubColumnTitles` remain in `raw_xml` and separate
metadata. `XAdvancedColumn` is retained in `alternate_x_columns`, separately
from XColumn; date conversion or choosing a representation is a later explicit
step. `StartAt`/`Interval` stay metadata; no series is generated.

XML CV uses `cv_percent`; `upper-lower-limits` uses `upper_limit`/`lower_limit`;
`low-high` uses distances. Absent YFormat receives `value_unspecified` (no
assumed statistical role). Unknown versions/formats get `unknown`. `Table` and
`HugeTable` are both read; TableSequence governs order with a documented fallback.
The compressed `Template` (analyses, graphs, layouts) is preserved, not decoded.

## Compatibility boundaries

The public corpus is 1 native project and 24 XML projects, with source commits,
hashes and license notices in `fixtures/prism_import`. This is a bounded file
extraction test, not a complete version matrix or Prism GUI round-trip test.
All original bytes remain available even when a semantic mapping is unsupported.
No automated transformation to a specialist's analysis CSV is currently provided.


## Schema compatibility and evidence

New imports write `prism-extraction/2`; `verify-prism` accepts both /1 and /2
manifests without rewriting earlier outputs. This schema revision adds native
per-row exclusion states, summaries, scoped profile evidence and validated
virtual-series recognition. It does not change stored values or statistical
methods. Unknown versions never inherit a supported version's role rules.

The `1-2-0` profile is restricted to `text`, `y_single`, `y_replicates` and the
explicit X series structure; other storage layouts remain unknown even if the
version is recognized. `1-4-1` retains its previously checked summary roles.
Independent synthetic cases cover range/default/conflict/malformed behaviors;
private files are not distributed as validation fixtures.
[Follow-up evidence](../../../validation/PRISM_NATIVE_SEMANTICS_2026-10-20.md).

GraphPad describes excluded values as retained in a data table while omitted
from analysis and graphing, and distinguishes cell exclusion from choosing
which datasets an analysis uses ([cell exclusion](https://www.graphpad.com/guides/prism/latest/user-guide/excluding_values.htm),
[dataset selection](https://www.graphpad.com/guides/prism/latest/user-guide/using_excluding_points.htm)).
Those guides do not specify the native ZIP flag encoding; our range/bit decoder
is a bounded structural interpretation, not an official schema or GUI oracle.
