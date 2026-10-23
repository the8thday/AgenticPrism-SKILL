"""Read-only, loss-aware Prism extraction. No fitting or numeric conversion.

The original bytes are authoritative. Parsed values are strings; descriptors
are evidence-bounded annotations, never instructions to transform the values.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import stat
import xml.etree.ElementTree as ET
import zipfile

from . import __version__

SCHEMA_VERSION = "prism-extraction/2"
MAX_SOURCE_BYTES = 64 * 1024 * 1024
MAX_EXPANDED_BYTES = 256 * 1024 * 1024
MAX_MEMBERS = 10000
MAX_CELLS = 2_000_000
MAX_COLUMNS = 10000

# Profiles are deliberately scoped to observed storage layouts, not all formats
# a particular Prism release might write. Display-format names are not storage
# contracts. Public evidence and independent synthetic cases are in validation/.
NATIVE_ROLES = {
    "y_single": ["value"], "text": ["text"],
    "y_sd": ["mean", "sd"], "y_sd_n": ["mean", "sd", "n"],
    "y_se": ["mean", "sem"], "y_se_n": ["mean", "sem", "n"],
    "y_plus_minus": ["mean", "upper_distance", "lower_distance"],
}
NATIVE_PROFILES = {
    "1-4-1": {"formats": frozenset(NATIVE_ROLES) | {"y_replicates"},
              "evidence": "pinned_native_1-4-1"},
    "1-2-0": {"formats": frozenset({"text", "y_single", "y_replicates"}),
              "evidence": "structural_profile_1-2-0"},
}
# Cell flags observed in the pinned public native project, plus EXCLUDED. Any
# other token could be an unfamiliar exclusion spelling, so it can never yield
# a "not marked" state; an explicit EXCLUDED still marks the cell.
KNOWN_CELL_FLAGS = frozenset({"EXCLUDED", "CENSORED", "SECTION_TITLE", "NUMERIC", "USER",
                              "FLOAT", "INTEGER", "POSITIVE", "NOTCONV", "NOEDIT"})
XML_ROLES = {
    "SD": ["mean", "sd"], "SDN": ["mean", "sd", "n"],
    "SE": ["mean", "sem"], "SEN": ["mean", "sem", "n"],
    "CV": ["mean", "cv_percent"], "CVN": ["mean", "cv_percent", "n"],
    "low-high": ["mean", "upper_distance", "lower_distance"],
    "upper-lower-limits": ["mean", "upper_limit", "lower_limit"],
    "text": ["text"],
}


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _issue(project, code, source, detail):
    project["issues"].append({"code": code, "source": source, "detail": detail})


def _json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def invalid(value):
        raise ValueError(f"Invalid JSON number: {value}")

    return json.loads(data, object_pairs_hook=pairs, parse_constant=invalid)


def _read_source(path):
    with Path(path).open("rb") as stream:
        data = stream.read(MAX_SOURCE_BYTES + 1)
    if len(data) > MAX_SOURCE_BYTES:
        raise ValueError("Source exceeds the 64 MiB read limit")
    return data


def _range(value):
    """Small, inclusive native ranges; never expand an unbounded range."""
    match = re.fullmatch(r"([0-9]+)(?:~([0-9]+))?", str(value))
    if not match:
        raise ValueError(f"Unsupported range: {value!r}")
    start = int(match[1]); end = int(match[2] or match[1])
    if end < start or end >= MAX_CELLS:
        raise ValueError("Invalid or excessive range")
    return start, end


def _replicates(dataset):
    if "replicates" in dataset and "replicate ranges" in dataset:
        raise ValueError("Ambiguous replicate metadata")
    if "replicates" in dataset:
        return dataset["replicates"]
    result = {}
    for item in dataset.get("replicate ranges", []):
        start, end = _range(item["range"])
        if end >= MAX_COLUMNS:
            raise ValueError("Excessive subcolumn count")
        for idx in range(start, end + 1):
            if idx in result:
                raise ValueError("Overlapping replicate ranges")
            result[idx] = item["replicate"]
    if result and set(result) != set(range(max(result) + 1)):
        raise ValueError("Non-contiguous replicate metadata")
    return [result[i] for i in range(len(result))]


def _native_columns(project, table, metadata, source, profile):
    columns = []
    virtual = []
    refs = []
    if table.get("rowTitlesDataSet"):
        refs.append(("row_title", table["rowTitlesDataSet"]))
    if table.get("xDataSet"):
        refs.append(("x", table["xDataSet"]))
    refs.extend(("y", uid) for uid in table.get("dataSets", []))
    for kind, uid in refs:
        name = f"data/sets/{uid}.json"
        dataset = metadata.get(name)
        if not isinstance(dataset, dict):
            raise ValueError(f"Missing or invalid dataset {uid}")
        reps = _replicates(dataset)
        fmt = dataset.get("format")
        if (profile and kind == "x" and fmt == "series" and len(reps) == 1
                and isinstance(reps[0], dict) and reps[0].get("@class") == "SeriesReplicate"
                and all(type(reps[0].get(k)) in (int, float) and math.isfinite(reps[0][k])
                        for k in ("startValue", "interval"))):
            virtual.append({"kind": kind, "dataset_id": uid, "metadata_source": name,
                            "metadata": dataset, "values": "not_materialized",
                            "role_evidence": profile["evidence"]})
            _issue(project, "implicit_series", name, "Series definition preserved separately; no generated X values or shifted Y columns")
            continue
        if not reps:
            raise ValueError(f"Missing replicate metadata for {uid}")
        roles = (["replicate"] * len(reps) if fmt == "y_replicates"
                 else NATIVE_ROLES.get(fmt, []))
        roles_known = bool(profile and fmt in profile["formats"] and len(roles) == len(reps)
                           and all(isinstance(r, dict) and r.get("@class") == "DriverReplicate" for r in reps))
        if not roles_known:
            roles = ["unknown"] * len(reps)
            _issue(project, "unverified_column_roles", name, f"format={fmt}; subcolumns={len(reps)}")
        for si, rep in enumerate(reps):
            if not isinstance(rep, dict):
                raise ValueError("Invalid replicate description")
            columns.append({"index": len(columns), "kind": kind,
                            "dataset_id": uid, "title": dataset.get("title"),
                            "subcolumn_index": si, "declared_format": fmt,
                            "table_display_format": table.get("dataFormat") if kind == "y" else table.get("xDataSetFormat") if kind == "x" else None,
                            "stored_role": "row_title" if kind == "row_title" and roles_known else roles[si],
                            "role_evidence": profile["evidence"] if roles_known else "unverified",
                            "metadata_source": name, "replicate_metadata": rep,
                            "exclusion_status": "not_interpreted"})
    return columns, virtual


def _flags(value):
    if not isinstance(value, list) or any(not isinstance(x, str) for x in value):
        raise ValueError("Flags must be a list of strings")
    return set(value)


def _native_exclusions(project, columns, rows, profile, source):
    """Decode the stored EXCLUDED bit, not an analysis-selection decision.

    Only non-conflicting explicit bits and unambiguous defaults are decoded.
    Whether a range replaces or supplements common flags is not assumed:
    common EXCLUDED plus a range without EXCLUDED therefore remains unknown.
    Ranges are absolute zero-based CSV row positions, not firstRow offsets.
    Unrecognized flag tokens make otherwise unmarked cells unknown.
    """
    if rows is not None and len(rows) * len(columns) > MAX_CELLS:
        _issue(project, "exclusion_annotation_limit", source,
               "Rectangular flag annotation exceeds limit; raw ragged values preserved")
        for column in columns:
            column.update(exclusion_status="resource_limit", excluded_by_row=None,
                          exclusion_evidence="unverified",
                          exclusion_scope="stored_cell_flag_only_not_analysis_selection")
        return
    for column in columns:
        rep = column["replicate_metadata"]
        nrows = len(rows) if rows is not None else 0
        states = [None] * nrows
        status = "unverified"
        if profile and rep.get("@class") == "DriverReplicate" and rows is not None:
            try:
                common_flags = _flags(rep.get("commonAttributes", []))
                common = "EXCLUDED" in common_flags
                unrecognized = common_flags - KNOWN_CELL_FLAGS
                common_unknown = bool(unrecognized) and not common
                overrides = rep.get("cellAttributes", [])
                if not isinstance(overrides, list):
                    raise ValueError("cellAttributes must be a list")
                # Bits: 1 range without EXCLUDED, 2 range with EXCLUDED,
                # 4 range without EXCLUDED but with an unrecognized token.
                mask = bytearray(nrows)
                visited = 0
                for item in overrides:
                    if not isinstance(item, dict) or set(item) != {"rows", "attributes"}:
                        raise ValueError("Unrecognized cell attribute record")
                    start, end = _range(item["rows"])
                    flags = _flags(item["attributes"])
                    if end >= nrows:
                        raise ValueError("Cell attribute range exceeds stored row count")
                    visited += end - start + 1
                    if visited > MAX_CELLS:
                        raise ValueError("Cell attribute expansion exceeds limit")
                    unrecognized |= flags - KNOWN_CELL_FLAGS
                    bit = 2 if "EXCLUDED" in flags else 5 if flags - KNOWN_CELL_FLAGS else 1
                    for row in range(start, end + 1):
                        mask[row] |= bit
                ambiguous = unknown = False
                states = []
                for m in mask:
                    if m == 0:
                        state = True if common else None if common_unknown else False
                        unknown |= state is None
                    elif m == 2:
                        state = True
                    elif m & 2 or common:
                        state = None; ambiguous = True
                    elif m & 4 or common_unknown:
                        state = None; unknown = True
                    else:
                        state = False
                    states.append(state)
                status = "partial" if ambiguous or unknown else "decoded"
                if ambiguous:
                    _issue(project, "ambiguous_exclusion_flags", source,
                           f"Column {column['index']}: conflicting/common-vs-range flags retained as unknown")
                if unrecognized:
                    _issue(project, "unrecognized_cell_attribute", source,
                           f"Column {column['index']}: {sorted(unrecognized)}; unmarked cells carrying them are unknown")
            except (ValueError, TypeError, KeyError) as exc:
                status = "invalid_metadata"
                _issue(project, "invalid_exclusion_metadata", source,
                       f"Column {column['index']}: {exc}")
        for row in range(nrows):
            if column["index"] >= len(rows[row]):
                states[row] = None  # Absent field, distinct from a stored empty field.
        column.update(exclusion_status=status, excluded_by_row=states,
                      exclusion_evidence="stored_EXCLUDED_marker" if status in ("decoded", "partial") else "unverified",
                      exclusion_scope="stored_cell_flag_only_not_analysis_selection")


def _exclusion_summary(columns, rows):
    marked = sum(v is True for c in columns for v in (c.get("excluded_by_row") or []))
    not_marked = sum(v is False for c in columns for v in (c.get("excluded_by_row") or []))
    mapped = sum(min(len(r), len(columns)) for r in rows) if rows is not None else 0
    return {
        "marked_excluded": marked,
        "not_marked_excluded": not_marked,
        "unknown_stored_cells": mapped - marked - not_marked,
        "unmapped_stored_cells": sum(max(0, len(r) - len(columns)) for r in rows) if rows is not None else None,
        "scope": "mapped_stored_cells_only; false does not establish use in an analysis",
    }


def _native(project, data):
    metadata = {}; blobs = {}; members = []
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        infos = archive.infolist()
        if len(infos) > MAX_MEMBERS or sum(i.file_size for i in infos) > MAX_EXPANDED_BYTES:
            raise ValueError("Archive exceeds member or expanded-size limit")
        seen = set()
        for info in infos:
            name = info.filename
            parts = PurePosixPath(name).parts
            if (name in seen or not parts or name.startswith("/") or "\\" in name
                    or ".." in parts or ":" in parts[0]
                    or stat.S_ISLNK(info.external_attr >> 16) or info.flag_bits & 1):
                raise ValueError(f"Unsafe, duplicate or encrypted archive member: {name}")
            seen.add(name)
            if info.file_size > MAX_SOURCE_BYTES:
                raise ValueError("Archive member exceeds 64 MiB limit")
            raw = archive.read(info)  # ZIP CRC checked; no extraction to paths.
            entry = {"name": name, "bytes": len(raw), "sha256": _sha(raw),
                     "handling": "directory" if info.is_dir() else "preserved_in_original"}
            if name.endswith(".json"):
                try:
                    metadata[name] = _json(raw)
                    entry["handling"] = "json_metadata"
                except (ValueError, UnicodeError) as exc:
                    _issue(project, "invalid_json", name, str(exc))
            if name.endswith(".csv"):
                blobs[name] = raw
            members.append(entry)
    doc = metadata.get("document.json")
    if not isinstance(doc, dict) or not isinstance(doc.get("sheets"), dict):
        raise ValueError("Not a supported Prism ZIP: document.json/sheets missing")
    project.update(format="prism", metadata=metadata, members=members,
                   version={k: doc.get(k) for k in ("formatVersion", "minFormatVersion", "minPrismVersion", "createdBy", "modifiedBy")})
    format_version = doc.get("formatVersion")
    profile = NATIVE_PROFILES.get(format_version) if isinstance(format_version, str) else None
    project["native_profile"] = {"version": format_version,
                                  "status": "scoped" if profile else "unknown",
                                  "storage_formats": sorted(profile["formats"]) if profile else []}
    if not profile:
        _issue(project, "unverified_native_version", "document.json", "Raw extraction only; native column roles withheld for this formatVersion")
    data_ids = doc["sheets"].get("data", [])
    if not isinstance(data_ids, list) or any(not isinstance(s, str) for s in data_ids):
        raise ValueError("Invalid data sheet sequence")
    if len(set(data_ids)) != len(data_ids):
        _issue(project, "duplicate_sheet_reference", "document.json", "Duplicate data sheet references retained in metadata")
    result_links = {}
    project["analyses"] = []
    for name, obj in metadata.items():
        if name.startswith("analyses/") and name.endswith("/sheet.json") and isinstance(obj, dict):
            project["analyses"].append({"source": name, "metadata": obj})
        if "/result_sheets/" in name and isinstance(obj, dict) and obj.get("dataSheet"):
            if isinstance(obj["dataSheet"], str):
                result_links.setdefault(obj["dataSheet"], []).append({"source": name, "metadata": obj})
            else:
                _issue(project, "invalid_result_link", name, "dataSheet is not a string; link not used")
    sheet_names = [f"data/sheets/{s}/sheet.json" for s in dict.fromkeys(data_ids)]
    sheet_names += sorted(n for n in metadata if n.startswith("data/sheets/") and n.endswith("/sheet.json") and n not in sheet_names)
    parsed_csv = {}
    cells = 0
    for name, raw in blobs.items():
        try:
            rows = []
            for row in csv.reader(io.StringIO(raw.decode("utf-8-sig"), newline=""), strict=True):
                cells += max(1, len(row))
                if cells > MAX_CELLS or len(row) > MAX_COLUMNS:
                    raise ValueError("CSV exceeds cell/column limit")
                rows.append(row)  # [] is a blank record; do not discard it.
            parsed_csv[name] = rows
        except (UnicodeError, csv.Error) as exc:
            _issue(project, "unreadable_csv", name, str(exc))
    linked = set()
    for name in sheet_names:
        sheet = metadata.get(name)
        if not isinstance(sheet, dict) or not isinstance(sheet.get("table"), dict):
            _issue(project, "missing_sheet", name, "Sheet or table metadata unavailable")
            continue
        table = sheet["table"]; sid = sheet.get("uid")
        if not isinstance(sid, str):
            _issue(project, "invalid_sheet_uid", name, "Sheet uid is not a string; not classified by uid")
            sid = None
        table_uid = table.get("uid")
        csv_name = f"data/tables/{table_uid}/data.csv" if isinstance(table_uid, str) else None
        linked.add(csv_name)
        rows = parsed_csv.get(csv_name)
        if rows is None:
            _issue(project, "missing_table_data", name, csv_name or "Table uid is not a string")
        try:
            columns, virtual = _native_columns(project, table, metadata, name, profile)
        except (ValueError, TypeError, KeyError) as exc:
            columns = []
            virtual = []
            _issue(project, "unmapped_columns", name, str(exc))
        storage = metadata.get(f"data/tables/{table.get('uid')}/content.json", {})
        if not isinstance(storage, dict):
            storage = {}
        width = max((len(r) for r in rows or []), default=0)
        if rows is not None and any(len(r) != width for r in rows):
            _issue(project, "ragged_csv", csv_name, "Unequal record widths/blank records preserved; absent fields are not padded")
        if rows is not None and (storage.get("numberOfRows") != len(rows)
                                or storage.get("numberOfColumns") != width):
            _issue(project, "dimension_mismatch", csv_name, "Observed CSV dimensions differ from content.json; no padding or truncation applied")
        if width != len(columns):
            _issue(project, "column_mapping_mismatch", csv_name, "CSV width and metadata columns differ; all column associations withheld")
            candidate = columns; columns = []
        else:
            candidate = []
        _native_exclusions(project, columns, rows, profile, csv_name)
        role = "input_data" if sid in data_ids else "analysis_result" if sid in result_links else "unclassified"
        project["tables"].append({"source_id": sid, "title": sheet.get("title"), "kind": role,
                                  "source": csv_name, "metadata_source": name,
                                  "declared_format": table.get("dataFormat"),
                                  "table_metadata": table, "storage_metadata": storage,
                                  "result_links": result_links.get(sid, []),
                                  "columns": columns, "unmapped_column_candidates": candidate,
                                  "virtual_columns": virtual,
                                  "exclusion_summary": _exclusion_summary(columns, rows),
                                  "raw_rows": rows, "observed_rows": len(rows) if rows is not None else None,
                                  "observed_columns": width})
    for name, rows in parsed_csv.items():
        if name not in linked:
            project["tables"].append({"source_id": name, "title": None, "kind": "unlinked_csv",
                                      "source": name, "columns": [], "raw_rows": rows,
                                      "exclusion_summary": _exclusion_summary([], rows),
                                      "observed_rows": len(rows), "observed_columns": max(map(len, rows), default=0)})
            _issue(project, "unlinked_csv", name, "Extracted, but not identified as input observations")
    for member in members:
        if member["name"] in parsed_csv:
            member["handling"] = "csv_extracted"
    project["limitations"] += [
        "Native roles use scoped 1-4-1 and 1-2-0 profiles; unsupported versions/layouts remain unverified, not a cross-version guarantee.",
        "Native EXCLUDED markers are decoded without removing values; ambiguous or unsupported flags stay null. False is not proof of analysis inclusion. Censoring remains uninterpreted.",
        "Native CV/CVN second values remain stored SD; high/low values remain stored distances in the verified fixture. No conversion applied.",
        "Graphs, layouts, binary payloads and .dt semantics are not decoded. SHA-256 hashes are ours, not a Prism .dt validation.",
        "Saved analysis results are historical file contents, not recomputed or scientifically validated results.",
    ]


def _tag(node):
    return node.tag.rsplit("}", 1)[-1]


def _children(node, name):
    return [e for e in node if _tag(e) == name]


def _text(node):
    return "" if node is None else "".join(node.itertext())


def _title(node):
    return _text(next(iter(_children(node, "Title")), None))


class _NoDTD(ET.TreeBuilder):
    def doctype(self, name, pubid, system):
        raise ValueError("DTD/entity declarations are not accepted")


def _xml(project, data):
    root = ET.fromstring(data, parser=ET.XMLParser(target=_NoDTD()))
    if _tag(root) != "GraphPadPrismFile":
        raise ValueError("Not a GraphPadPrismFile XML document")
    project.update(format="pzfx", version={"root_attributes": dict(root.attrib),
                   "created": [dict(e.attrib) for e in root.iter() if _tag(e) in ("OriginalVersion", "MostRecentVersion")]},
                   metadata={"root_attributes": dict(root.attrib)}, members=[], analyses=[])
    xml_version = root.get("PrismXMLVersion")
    known_version = xml_version == "5.00"
    if not known_version:
        _issue(project, "unverified_xml_version", "/GraphPadPrismFile", str(xml_version))
    table_nodes = [(i, e) for i, e in enumerate(root) if _tag(e) in ("Table", "HugeTable")]
    seq = [ref.get("ID") for s in _children(root, "TableSequence") for ref in _children(s, "Ref")]
    ids = [e.get("ID") for _, e in table_nodes]
    if len(set(ids)) != len(ids) or None in ids:
        _issue(project, "ambiguous_table_ids", "/GraphPadPrismFile", "Missing/duplicate IDs retained; output IDs use source positions")
    if len(set(seq)) != len(seq) or any(s not in ids for s in seq):
        _issue(project, "invalid_table_sequence", "/TableSequence", "Unresolved/duplicate references preserved")
    order = {key: i for i, key in enumerate(dict.fromkeys(seq))}
    if not seq and table_nodes:
        _issue(project, "document_order_fallback", "/TableSequence", "No sequence; document order used")
    table_nodes.sort(key=lambda pair: (order.get(pair[1].get("ID"), len(order)), pair[0]))
    cells = 0
    for position, table in table_nodes:
        source = f"/GraphPadPrismFile/{_tag(table)}[@document_child={position}]"
        columns = []; extra = []
        fmt = table.get("YFormat")
        for ci, column in enumerate(table):
            kind = _tag(column)
            if kind not in ("RowTitlesColumn", "XColumn", "YColumn", "XAdvancedColumn"):
                continue
            subs = _children(column, "Subcolumn")
            for child in column:
                if _tag(child) not in ("Title", "Subcolumn"):
                    _issue(project, "uninterpreted_column_element", source, _tag(child))
            if kind == "RowTitlesColumn":
                roles = ["row_title"] * len(subs)
            elif kind in ("XColumn", "XAdvancedColumn"):
                xf = table.get("XFormat")
                roles = ["x", "x_error"] if xf == "error" else ["x"] if xf in ("numbers", "date", "time", "text") else []
            else:
                roles = (["value_unspecified"] * len(subs) if fmt is None else
                         ["replicate"] * len(subs) if fmt == "replicates" else XML_ROLES.get(fmt, []))
            if not known_version or len(roles) != len(subs):
                roles = ["unknown"] * len(subs)
                _issue(project, "unverified_column_roles", source, f"{kind}: {table.get('XFormat') if kind.startswith('X') else fmt}")
            for si, sub in enumerate(subs):
                values = []
                for child in sub:
                    if _tag(child) not in ("Title", "d"):
                        _issue(project, "uninterpreted_subcolumn_element", source, _tag(child))
                for ri, node in enumerate(_children(sub, "d")):
                    cells += 1
                    if cells > MAX_CELLS:
                        raise ValueError("XML exceeds cell limit")
                    flag = node.get("Excluded")
                    excluded = {None: False, "0": False, "false": False, "1": True, "true": True}.get(flag)
                    if excluded is None:
                        _issue(project, "unknown_exclusion_flag", source, str(flag))
                    values.append({"row": ri, "raw_text": _text(node), "excluded": excluded,
                                   "attributes": dict(node.attrib), "xml": ET.tostring(node, encoding="unicode")})
                record = {"index": len(extra) if kind == "XAdvancedColumn" else len(columns),
                          "kind": kind, "title": _title(column), "source_column_index": ci,
                          "subcolumn_index": si, "subcolumn_title": _title(sub),
                          "column_attributes": dict(column.attrib), "subcolumn_attributes": dict(sub.attrib),
                          "stored_role": roles[si], "declared_format": table.get("XFormat") if kind.startswith("X") else fmt,
                          "role_evidence": "xml_format_mapping" if roles[si] != "unknown" else "unverified",
                          "values": values}
                (extra if kind == "XAdvancedColumn" else columns).append(record)
            if len(columns) + len(extra) > MAX_COLUMNS:
                raise ValueError("XML exceeds column limit")
        if extra:
            _issue(project, "alternate_x_representation", source, "XAdvancedColumn preserved separately; no date conversion or substitution")
        # Preserve the entire table, including notes, namespaces, rich labels and
        # unknown children. The original XML is the byte-exact representation.
        table_type = table.get("TableType")
        kind = ("analysis_result" if table_type == "Result" else "input_data"
                if table_type in ("XY", "OneWay", "TwoWay", "Contingency", "Survival", "PartsOfWhole") else "unclassified")
        if kind == "unclassified":
            _issue(project, "unclassified_table", source, f"TableType={table_type}; not assumed to be input observations")
        project["tables"].append({"source_id": table.get("ID"), "title": _title(table),
                                  "kind": kind,
                                  "source": source, "table_metadata": dict(table.attrib),
                                  "declared_format": fmt, "columns": columns,
                                  "alternate_x_columns": extra,
                                  "subcolumn_titles_xml": [ET.tostring(e, encoding="unicode") for e in _children(table, "SubColumnTitles")],
                                  "observed_rows": max((len(c["values"]) for c in columns + extra), default=0),
                                  "observed_columns": len(columns), "raw_xml": ET.tostring(table, encoding="unicode")})
        known = {"Title", "RowTitlesColumn", "XColumn", "XAdvancedColumn", "YColumn", "SubColumnTitles"}
        for child in table:
            if _tag(child) not in known:
                _issue(project, "uninterpreted_table_element", source, _tag(child))
    project["metadata"]["non_table_elements"] = [
        {"tag": _tag(e), "document_child": i, "xml": ET.tostring(e, encoding="unicode")}
        for i, e in enumerate(root) if _tag(e) not in ("Table", "HugeTable")]
    if any(_tag(e) == "Template" for e in root):
        _issue(project, "opaque_template", "/Template", "Compressed analysis/graph/layout payload preserved, not decoded")
    project["limitations"] += [
        "XML d values are text, never automatically converted to numbers, dates, NA or missing values.",
        "Explicit empty d is an empty string; an absent row in a subcolumn is absent, never recycled or filled.",
        "Template payloads and unknown XML metadata are preserved but not semantically decoded.",
        "Table/Subcolumn titles and excluded values are retained; exclusion flags are not applied to remove observations.",
    ]


def _parse(data, filename):
    project = {"schema_version": SCHEMA_VERSION, "source": {"name": filename, "bytes": len(data), "sha256": _sha(data)},
               "tables": [], "issues": [], "limitations": [
                   "Extraction does not identify independent units, infer units or validate an analysis.",
                   "Values are stored representations. No CV conversion, bound reconstruction, normalization, fitting or exclusion is performed."]}
    if zipfile.is_zipfile(io.BytesIO(data)):
        _native(project, data)
    else:
        _xml(project, data)
    for i, table in enumerate(project["tables"]):
        table["id"] = f"table_{i + 1:04d}"
    return project


def read_prism(path):
    """Return extracted project data without writing files (signature detection)."""
    return _parse(_read_source(path), Path(path).name)


def summarize(project):
    """Agent-facing inventory, explicitly separate from full table extraction."""
    return {"schema_version": SCHEMA_VERSION, "format": project["format"], "source": project["source"],
            "version": project["version"], "status": "extracted_with_limits",
            "tables": [{k: t.get(k) for k in ("id", "source_id", "title", "kind", "declared_format", "observed_rows", "observed_columns", "source", "exclusion_summary")} for t in project["tables"]],
            "counts": {kind: sum(t["kind"] == kind for t in project["tables"])
                       for kind in ("input_data", "analysis_result", "unclassified", "unlinked_csv")},
            "native_exclusions": {key: sum(t.get("exclusion_summary", {}).get(key, 0) or 0 for t in project["tables"])
                                  for key in ("marked_excluded", "not_marked_excluded", "unknown_stored_cells", "unmapped_stored_cells")}
                                  if project["format"] == "prism" else None,
            "issues": project["issues"], "limitations": project["limitations"]}


def _dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def import_prism(path, output):
    """Write a new, auditable extraction directory; never overwrite a directory."""
    data = _read_source(path)
    project = _parse(data, Path(path).name)
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    original = "source." + project["format"]
    (out / original).write_bytes(data)
    (out / "tables").mkdir()
    inventory = []
    for table in project["tables"]:
        destination = f"tables/{table['id']}.json"
        _dump(out / destination, table)
        inventory.append({"id": table["id"], "file": destination, "source": table["source"]})
    _dump(out / "project.json", {k: v for k, v in project.items() if k != "tables"} | {"tables": inventory})
    report = summarize(project)
    _dump(out / "import_report.json", report)
    _dump(out / "interpretation_facts.json", {
        "schema_version": SCHEMA_VERSION, "purpose": "extraction_only",
        "source_sha256": {"import_report.json": _sha((out / "import_report.json").read_bytes())},
        "primary": {"counts": report["counts"], "status": report["status"], "native_exclusions": report["native_exclusions"]},
        "must_mention": report["limitations"], "issues": report["issues"],
        "scientific_inference": "not_performed"})
    _dump(out / "import_manifest.json", {
        "schema_version": SCHEMA_VERSION, "runtime_version": __version__, "parser_sha256": _sha(Path(__file__).read_bytes()),
        "original": original, "original_sha256": _sha(data),
        "artifacts_sha256": {p.relative_to(out).as_posix(): _sha(p.read_bytes())
                             for p in sorted(out.rglob("*")) if p.is_file()}})
    return out


def verify_prism_import(directory):
    """Verify our recorded extraction hashes (not scientific validity)."""
    out = Path(directory)
    manifest = _json((out / "import_manifest.json").read_bytes())
    if manifest.get("schema_version") not in ("prism-extraction/1", SCHEMA_VERSION):
        raise ValueError("Unknown extraction manifest schema")
    for name, digest in manifest["artifacts_sha256"].items():
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or "\\" in name:
            raise ValueError("Unsafe manifest path")
        actual = out / name
        if not actual.resolve().is_relative_to(out.resolve()) or actual.is_symlink():
            raise ValueError("Manifest artifact escapes extraction directory")
        if _sha(actual.read_bytes()) != digest:
            raise ValueError(f"Artifact hash mismatch: {name}")
    if manifest["artifacts_sha256"].get(manifest["original"]) != manifest["original_sha256"]:
        raise ValueError("Source hash mismatch")
    present = {p.relative_to(out).as_posix() for p in out.rglob("*")
               if p.is_file() or p.is_symlink()} - {"import_manifest.json"}
    extra = sorted(present - set(manifest["artifacts_sha256"]))
    if extra:
        raise ValueError(f"Unrecorded artifact: {extra[0]}")
    return manifest
