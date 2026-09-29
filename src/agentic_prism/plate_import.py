"""Adapter from plate-reader grid CSVs plus plate-map layers to a canonical long table.

A grid CSV has a header row whose first cell is a free label followed by column
numbers 1..12 (96-well) or 1..24 (384-well), then one row per plate row A..H or
A..P. The importer never fills, averages, blank-subtracts or reinterprets values:
every output cell comes from exactly one declared layer, plate constant or
global constant, and wells whose key layer is blank are skipped and counted.
"""
from pathlib import Path
import csv
import json
import shutil
import string
import numpy as np
import pandas as pd
from .workflow import dump, sha

SHAPES = {(8, 12): 96, (16, 24): 384}
TARGETS = {
    "hts_qc": {"key": "role", "columns": ["plate_id", "well_id", "compound_id", "independent_unit_id", "row", "column", "role", "value"], "optional": [], "template": {"analysis_type": "hts_qc", "input": "data.csv", "plate_shape": None, "readout": None, "response_unit": None, "direction": None, "independent_wells": False, "randomized_layout": False, "majority_inactive": False}},
    "elisa_quantification": {
        "key": "role",
        "columns": ["plate_id", "well_id", "role", "sample_id", "concentration", "concentration_unit",
                    "dilution_factor", "response", "response_unit"],
        "optional": ["exclude", "exclusion_reason"],
        "template": {"analysis_type": "elisa_quantification", "input": "data.csv", "direction": None,
                     "assay": {"analyte": "", "response_definition": "", "matrix": "", "rationale": ""}}},
    "dose_response_4pl": {
        "key": "curve_id",
        "columns": ["observation_id", "plate_id", "well_id", "sample_id", "experiment_id", "curve_id",
                    "concentration", "concentration_unit", "response", "response_unit"],
        "optional": ["replicate_id", "exclude", "exclusion_reason"],
        "template": {"analysis_type": "dose_response_4pl", "input": "data.csv",
                     "assay": {"endpoint": None, "direction": None, "tested_agent": "", "response_definition": "",
                               "relative_half_response_supported": None, "rationale": ""}}},
}
AUTOMATIC = {"plate_id", "well_id", "response", "observation_id", "row", "column", "value"}


def read_grid(path):
    """Return {well_id: text} for a 96- or 384-well grid; cells are stripped strings."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = [r for r in csv.reader(handle)]
    while rows and not any(cell.strip() for cell in rows[-1]):
        rows.pop()
    if not rows:
        raise ValueError(f"Empty plate grid: {path}")
    header = [c.strip() for c in rows[0][1:]]
    body = rows[1:]
    shape = (len(body), len(header))
    if shape not in SHAPES:
        raise ValueError(f"{path}: grid is {shape[0]}x{shape[1]}; expected 8x12 (96-well) or 16x24 (384-well)")
    if header != [str(i) for i in range(1, shape[1] + 1)]:
        raise ValueError(f"{path}: header must be column numbers 1..{shape[1]}")
    cells = {}
    for letter, row in zip(string.ascii_uppercase, body):
        if row[0].strip().upper() != letter or len(row) != shape[1] + 1:
            raise ValueError(f"{path}: expected row {letter} with {shape[1]} values")
        for column, value in enumerate(row[1:], 1):
            cells[f"{letter}{column}"] = value.strip()
    return cells, SHAPES[shape]


def import_plates(manifest_path, output):
    src, out = Path(manifest_path).resolve(), Path(output).resolve()
    if out.exists():
        raise FileExistsError("Choose a new import directory")
    raw = json.loads(src.read_text())
    if set(raw) - {"format", "target", "source", "constants", "plates"} or raw.get("format") != "plate_grid_csv":
        raise ValueError("Unsupported plate manifest; format must be plate_grid_csv")
    if raw.get("target") not in TARGETS or not isinstance(raw.get("source"), str) or not raw["source"].strip():
        raise ValueError(f"Declare target ({sorted(TARGETS)}) and a source description")
    spec = TARGETS[raw["target"]]
    allowed = set(spec["columns"] + spec["optional"]) - AUTOMATIC
    constants = raw.get("constants", {})
    if not isinstance(constants, dict) or set(constants) - allowed:
        raise ValueError(f"Unsupported constants; allowed {sorted(allowed)}")
    if not isinstance(raw.get("plates"), list) or not raw["plates"]:
        raise ValueError("No plates declared")
    rows, files, plate_records, seen = [], {}, [], set()

    def register(path, role):
        key = str(path)
        files.setdefault(key, {"original_path": key, "sha256": sha(path), "roles": []})
        files[key]["roles"].append(role)

    for plate in raw["plates"]:
        if not isinstance(plate, dict) or set(plate) - {"plate_id", "readout", "layers", "constants"}:
            raise ValueError("Each plate needs plate_id, readout, layers and optional constants")
        pid = plate.get("plate_id")
        if not isinstance(pid, str) or not pid.strip() or pid in seen:
            raise ValueError("plate_id must be a unique non-empty string")
        seen.add(pid)
        layers, local = plate.get("layers", {}), plate.get("constants", {})
        if not isinstance(layers, dict) or not isinstance(local, dict) or set(layers) - allowed or set(local) - allowed:
            raise ValueError(f"Plate {pid}: unsupported layer or constant; allowed {sorted(allowed)}")
        if spec["key"] not in layers:
            raise ValueError(f"Plate {pid}: the {spec['key']} layer defines which wells are used and must be a grid")
        for column in spec["columns"] + spec["optional"]:
            sources = [name for name, found in (("layer", column in layers), ("plate constant", column in local),
                                                ("global constant", column in constants)) if found]
            if len(sources) > 1:
                raise ValueError(f"Plate {pid}: {column} is defined by more than one source {sources}")
            if column not in AUTOMATIC and column in spec["columns"] and not sources:
                raise ValueError(f"Plate {pid}: {column} needs a layer or constant")
        readout_path = (src.parent / plate["readout"]).resolve()
        readout, size = read_grid(readout_path)
        register(readout_path, f"{pid}:readout")
        grids = {}
        for column, name in layers.items():
            path = (src.parent / name).resolve()
            grids[column], layer_size = read_grid(path)
            register(path, f"{pid}:{column}")
            if layer_size != size:
                raise ValueError(f"Plate {pid}: layer {column} is {layer_size}-well but readout is {size}-well")
        used = [w for w in readout if grids[spec["key"]][w]]
        if not used:
            raise ValueError(f"Plate {pid}: no wells assigned in the {spec['key']} layer")
        for well in used:
            text = readout[well]
            try:
                value = float(text)
            except ValueError:
                raise ValueError(f"Plate {pid} well {well}: readout {text!r} is not numeric; resolve saturated or "
                                 "missing reads explicitly in the source before import") from None
            if not np.isfinite(value):
                raise ValueError(f"Plate {pid} well {well}: nonfinite readout")
            record = {"plate_id": pid, "well_id": well, "response": text}
            if raw["target"] == "hts_qc":
                record.update(row=string.ascii_uppercase.index(well[0])+1, column=int(well[1:]), value=text)
            if "observation_id" in spec["columns"]:
                record["observation_id"] = f"{pid}-{well}"
            for column in spec["columns"] + spec["optional"]:
                if column in grids:
                    record[column] = grids[column][well]
                elif column in local:
                    record[column] = str(local[column])
                elif column in constants:
                    record[column] = str(constants[column])
            rows.append(record)
        plate_records.append({"plate_id": pid, "plate_size": size, "wells_used": len(used),
                              "wells_skipped_blank_key": size - len(used),
                              "wells_skipped_with_readout": sum(1 for w in readout if readout[w] and w not in used)})
    present = [c for c in spec["columns"] + spec["optional"] if any(c in r for r in rows)]
    out.mkdir(parents=True)
    (out / "source_files").mkdir()
    for i, record in enumerate(files.values()):
        dest = out / "source_files" / f"{i:03d}_{Path(record['original_path']).name}"
        shutil.copyfile(record["original_path"], dest)
        record["snapshot_path"] = str(dest.relative_to(out))
    pd.DataFrame(rows, columns=present).to_csv(out / "data.csv", index=False)
    dump(out / "import_manifest.json", {"adapter": "plate_grid_csv_v1", "target": raw["target"], "source": raw["source"],
          "configuration": raw, "input_files": list(files.values()), "plates": plate_records,
          "data_sha256": sha(out / "data.csv"),
          "preprocessing": "None: no blank subtraction, averaging, outlier removal or unit conversion. Readout text copied "
                           "as-is; wells with a blank key layer skipped and counted."})
    dump(out / "config.template.json", {**spec["template"], "source": raw["source"]})
    return out
