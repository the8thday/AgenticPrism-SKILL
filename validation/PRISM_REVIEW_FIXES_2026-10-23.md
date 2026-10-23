# Prism extraction review fixes — 2026-10-23

Local unreleased follow-up to [the initial reader](PRISM_IMPORT_2026-10-19.md)
and [the native compatibility record](PRISM_NATIVE_SEMANTICS_2026-10-20.md),
which remain unchanged. Runtime stays 0.14.4; schema stays `prism-extraction/2`.
No statistical method, fixture or historical artifact is changed. Reviewed commit:
`9a5954ad7628562a605ec617923f5aba44686d8d`.

## Defect found in review

A native cell-flag token that differs from `EXCLUDED` (for example `EXCLUDE`)
was decoded as "not marked" with status `decoded` and no issue. The public
native project contains no excluded cell, so the `EXCLUDED` token itself is not
demonstrated by any public native file in this repository; only synthetic tests
and the undistributed local acceptance file support it. Neither earlier record
states which native version that private file used. On 2026-10-23 the user
identified it: a `formatVersion 1-2-0` project (minFormatVersion `1-1-0`) created
and last saved by Prism 10.1.2 on Windows. Its datasets use only `y_single`,
`y_replicates`, `text` and one X `series`; its cell flags include `EXCLUDED`
(10 range records) and `SECTION_TITLE`, with common tokens `NUMERIC`, `USER` and
`NOEDIT`. It was read in place and is not distributed; no name, label, value or
hash from it is recorded. So `EXCLUDED` is observed in a real 1-2-0 file, but
there is still no public native example, and none for `1-4-1`. The fixed parser
raises no `unrecognized_cell_attribute` issue and no unknown flag state on it.

The pinned public project stores these tokens: `commonAttributes` with `NUMERIC`,
`USER`, `INTEGER`, `POSITIVE`, `NOTCONV`, `NOEDIT`; `cellAttributes` with
`SECTION_TITLE`, `CENSORED`, `FLOAT`, `POSITIVE`, `NUMERIC`. Its range rows lie
within each replicate's `firstRow`..`lastRow` (for example rows 20 and 24 with
firstRow 20, lastRow 26), consistent with the absolute-row reading already used.

## Changes

- Recognized flag tokens: the observed set above plus `EXCLUDED`. Any other
  token raises `unrecognized_cell_attribute`. Cells covered by it, or every row
  for an unrecognized common token, become null unless they carry explicit
  `EXCLUDED`; the column becomes `partial`. Explicit `EXCLUDED` behaviour and the
  existing ambiguity rules are unchanged.
- Range rows accept ASCII digits only (`\d` had accepted, for example, `١` as 1).
- A non-string `formatVersion` is an unknown version instead of a `TypeError`.
  Non-string sheet uids, table uids and result `dataSheet` links are issues
  (`invalid_sheet_uid`, `missing_table_data`, `invalid_result_link`); a table
  without a string uid no longer links to a literal `data/tables/None` CSV.
- `verify-prism` rejects files that the manifest does not record. All six
  committed demo extractions (/1 and /2) still verify.
- `scripts/validate_prism_import.py` replaces `assert` with explicit checks that
  also run under `python -O`.

## Checks

- Parser tests: **76 passed** (68 earlier plus 8 new cases for unrecognized
  tokens, recognized public tokens, non-ASCII digits, non-string identifiers and
  unrecorded files).
- [Reproducible artifacts](prism_review_fixes_2026-10-23/validation.json):
  `.venv/bin/python scripts/validate_prism_import.py --output validation/prism_review_fixes_2026-10-23`
  passed: 25 public files, 66 tables, 7,571 stored fields. `public_corpus.json` and
  `legacy.json` are byte-identical to the 2026-10-20 run; in both demo extractions
  every file except `import_manifest.json` (parser hash) is byte-identical. The
  six-fixture legacy CLI check again gave 61 byte-identical scientific artifacts.
- Full suite: **1,055 passed** (1,047 earlier plus 8 new), 8 existing `Pandas4Warning`
  warnings in `src/agentic_prism/ada.py:209-210`, 453.56 seconds, on the final parser
  ([checks.json](prism_review_fixes_2026-10-23/checks.json)).
- Both changed Skill folders pass `quick_validate.py`; `git diff --check` passes;
  standard-library-only `inspect-prism` (`python -S`) still works.

## Gates

Unchanged from the 2026-10-20 record: native Prism GUI/schema oracle UNMET,
live-agent misuse scenarios UNMET, cross-version completeness UNMET. A public
native file containing excluded cells would close the `EXCLUDED` token gap.
