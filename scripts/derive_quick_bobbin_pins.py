#!/usr/bin/env python3
"""Derive src/data/quick_bobbin_pins.json from the MAS catalogue (ABT #1220).

Bobbin::create_quick_bobbin can synthesise a two-row THT pinout for a core that has no
catalogue former. Every number it uses comes from this table, and every number in the table
comes from the catalogue records in MAS/data/bobbins.ndjson that state their pin geometry
(pitch AND rowDistance), joined to their core shape in MAS/data/core_shapes.ndjson -- except
the 2.54 mm grid itself, which is the IEC 62317 grid named in ABT #1220 and which the same
records confirm (the evidence count is written next to it).

Usage:  python3 scripts/derive_quick_bobbin_pins.py [--check]
        (run from the MKF root; --check exits 1 when the committed JSON differs)

Core dimensions are read the way MKF's CorePiece processors read them
(src/constructive_models/CorePiece.cpp), so the statistics refer to exactly the quantities
create_quick_bobbin measures at runtime:
  width  (core.get_width(), the pitch axis X)   = A
  height (core.get_height(), the column axis Y) = 2 B for the two-piece sets used here
  depth  (core.get_depth(), Z)                  = E for rm, C + max(0, K) for efd, C otherwise
  window height (processed winding window)      = 2 D
"""

import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BOBBINS = ROOT / "MAS" / "data" / "bobbins.ndjson"
SHAPES = ROOT / "MAS" / "data" / "core_shapes.ndjson"
OUTPUT = ROOT / "src" / "data" / "quick_bobbin_pins.json"

GRID = 0.00254          # IEC 62317 grid, 0.1 in
HALF_GRID = GRID / 2     # 1.27 mm: 3.81 mm (1.5 grid) is a catalogue pitch, so pitches snap to half-grid steps
SNAP_TOLERANCE = 0.0003  # a scraped metric pitch (5.0, 3.75, 3.8, 2.5) within 0.3 mm of a half-grid step is that step
MINIMUM_FOOTPRINTS_PER_BAND = 2

# Core shape families pooled into one construction class. Only families that have at least one
# catalogue former with pin geometry appear; create_quick_bobbin throws for any other family.
CLASSES = {
    "E-type": ["e", "ec", "efd", "er", "etd", "planarER"],
    "PQ": ["pq"],
    "RM": ["rm"],
    "EP": ["ep"],
}


def nominal(value):
    if isinstance(value, (int, float)):
        return float(value)
    if "nominal" in value:
        return float(value["nominal"])
    if "minimum" in value and "maximum" in value:
        return (float(value["minimum"]) + float(value["maximum"])) / 2
    if "maximum" in value:
        return float(value["maximum"])
    return float(value["minimum"])


def core_dimensions(shape):
    d = {k: nominal(v) for k, v in shape["dimensions"].items()}
    family = shape["family"]
    if family == "rm":
        depth = d["E"]
    elif family == "efd":
        depth = d["C"] + max(0.0, d.get("K", 0.0))
    else:
        depth = d["C"]
    return {"width": d["A"], "height": 2 * d["B"], "depth": depth, "windowHeight": 2 * d["D"]}


def snap(pitch):
    steps = round(pitch / HALF_GRID)
    snapped = steps * HALF_GRID
    return round(snapped, 6), abs(snapped - pitch) <= SNAP_TOLERANCE


def row_layout(pinout):
    """Pins and pitch per row, or None when the record does not say how its pins split."""
    pitch = pinout["pitch"]
    if pinout.get("numberPinsPerRow"):
        per_row = pinout["numberPinsPerRow"]
    else:
        rows = pinout.get("numberRows", 2)
        if pinout["numberPins"] % rows:
            return None
        per_row = [pinout["numberPins"] // rows] * rows
    pitches = pitch if isinstance(pitch, list) else [pitch] * len(per_row)
    if len(pitches) != len(per_row):
        return None
    return list(zip(per_row, pitches))


def row_span(pins, pitch, central):
    if pins < 2:
        return 0.0
    span = (pins - 1) * pitch
    if central is not None and pins % 2 == 0:
        span += central - pitch
    return span


def median(values):
    return round(statistics.median(values), 6)


def summary(values):
    return {"footprints": len(values), "minimum": round(min(values), 6), "median": median(values), "maximum": round(max(values), 6)}


def main():
    shapes = {}
    for line in SHAPES.read_text().splitlines():
        if line.strip():
            shape = json.loads(line)
            shapes[shape["name"]] = shape
            for alias in shape.get("aliases") or []:
                shapes.setdefault(alias, shape)
    family_class = {family: name for name, families in CLASSES.items() for family in families}

    skipped = Counter()
    pitch_samples = defaultdict(list)            # class -> [(width, snapped pitch, record index)]
    margins = defaultdict(list)                  # class -> [edge margin]
    clearances = defaultdict(list)               # (class, orientation) -> [clearance]
    pin_diameters, pin_lengths, rail_standoffs = [], [], []
    rail_outer_walls, rail_inner_walls, rail_end_walls = [], [], []
    grid_hits = 0
    grid_total = 0
    used = 0
    footprints = set()
    rails = set()
    for index, line in enumerate(BOBBINS.read_text().splitlines()):
        if not line.strip():
            continue
        record = json.loads(line)
        functional = record.get("functionalDescription", {})
        pinout = functional.get("pinout") or {}
        if pinout.get("pitch") is None or pinout.get("rowDistance") is None:
            continue
        if functional.get("base") is not None:
            # ABT #1173: a toroid base (family t) is not a former; its pins stand on the base, not on a rail.
            skipped["toroid base (family t), not a former"] += 1
            continue
        shape = shapes.get(functional.get("shape"))
        if shape is None:
            skipped["shape not in core_shapes.ndjson"] += 1
            continue
        klass = family_class.get(shape["family"])
        if klass is None:
            skipped["core family " + shape["family"] + " not in a class"] += 1
            continue
        orientation = functional.get("orientation")
        if orientation not in ("vertical", "horizontal"):
            skipped["no orientation"] += 1
            continue
        layout = row_layout(pinout)
        if layout is None:
            skipped["pin split over the rows not stated"] += 1
            continue
        used += 1
        # The catalogue repeats one footprint under several part numbers (chamber variants, pin
        # materials, duplicated scrapes); a footprint is one piece of evidence, however often listed.
        geometry = {k: v for k, v in pinout.items() if k != "pinDescription"}
        footprint = (shape["name"], orientation, json.dumps(geometry, sort_keys=True))
        if footprint in footprints:
            skipped["repeat of a footprint already counted"] += 1
            continue
        footprints.add(footprint)
        core = core_dimensions(shape)
        for _, pitch in layout:
            snapped, on_grid = snap(pitch)
            grid_total += 1
            grid_hits += on_grid
            pitch_samples[klass].append((core["width"], snapped, index))
        central = pinout.get("centralPitch")
        span = max(row_span(pins, pitch, central) for pins, pitch in layout)
        margins[klass].append((core["width"] - span) / 2)
        row_distance = pinout["rowDistance"]
        if orientation == "vertical":
            clearances[(klass, orientation)].append((row_distance - core["depth"]) / 2)
        else:
            clearances[(klass, orientation)].append((row_distance - core["windowHeight"]) / 2)
        description = pinout.get("pinDescription")
        if description and len(description.get("dimensions", [])) >= 3:
            pin_diameters.append(description["dimensions"][0])
            pin_lengths.append(description["dimensions"][2])
        dims = functional.get("dimensions", {})
        # The one verified pin-rail datum (Bobbin::get_pin_rail_distance, ABT #1207 corrected by
        # ABT #1249): vertical PQ, rail underside at H1/2 + H3 below the column centre; its standoff
        # beyond the core is that minus half the core height. The same drawings give the rail block
        # (a, b, b1): its walls around the pin row, measured from the pin surface.
        rail_labels = ("H1", "H3", "a", "b", "b1")
        rail = (shape["name"],) + tuple(json.dumps(dims.get(label)) for label in rail_labels)
        if (shape["family"] == "pq" and orientation == "vertical" and all(label in dims for label in rail_labels)
                and description and rail not in rails):
            rails.add(rail)
            half_row = row_distance / 2
            pin_radius = description["dimensions"][0] / 2
            rail_standoffs.append(nominal(dims["H1"]) / 2 + nominal(dims["H3"]) - core["height"] / 2)
            rail_outer_walls.append(nominal(dims["b"]) / 2 - half_row - pin_radius)
            rail_inner_walls.append(half_row - (nominal(dims["b"]) / 2 - nominal(dims["b1"])) - pin_radius)
            rail_end_walls.append(nominal(dims["a"]) / 2 - span / 2 - pin_radius)

    classes = {}
    for klass, families in CLASSES.items():
        samples = sorted(pitch_samples[klass])
        by_width = defaultdict(list)
        for width, pitch, index in samples:
            by_width[round(width, 6)].append((pitch, index))
        widths = sorted(by_width)
        # Mode pitch per catalogue core width; a tie goes to the smaller pitch (more pins).
        points = []
        for width in widths:
            counts = Counter(p for p, _ in by_width[width])
            best = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
            points.append((width, best, len({i for _, i in by_width[width]})))
        # Merge neighbouring widths with the same pitch into bands; boundaries half-way between.
        bands = []
        for width, pitch, records in points:
            if bands and bands[-1]["catalogueMode"] == pitch:
                bands[-1]["coreWidths"].append(width)
                bands[-1]["footprints"] += records
            else:
                bands.append({"catalogueMode": pitch, "coreWidths": [width], "footprints": records})
        result = []
        for position, band in enumerate(bands):
            lower = 0.0 if position == 0 else round((bands[position - 1]["coreWidths"][-1] + band["coreWidths"][0]) / 2, 6)
            upper = None if position == len(bands) - 1 else round((band["coreWidths"][-1] + bands[position + 1]["coreWidths"][0]) / 2, 6)
            entry = {"minimumCoreWidth": lower, "maximumCoreWidth": upper,
                     "catalogueCoreWidths": band["coreWidths"], "footprints": band["footprints"]}
            if band["footprints"] >= MINIMUM_FOOTPRINTS_PER_BAND:
                entry["pitch"] = band["catalogueMode"]
                entry["source"] = "catalogue mode"
            else:
                entry["pitch"] = GRID
                entry["catalogueMode"] = band["catalogueMode"]
                entry["source"] = ("standard grid: fewer than " + str(MINIMUM_FOOTPRINTS_PER_BAND) +
                                   " distinct catalogue footprints in this band, so the IEC 62317 2.54 mm grid pitch is used, not a catalogue fit")
            result.append(entry)
        entry = {"coreShapeFamilies": families, "pitchBands": result,
                 "edgeMargin": summary(margins[klass])}
        for orientation in ("vertical", "horizontal"):
            values = clearances.get((klass, orientation))
            entry[orientation + "RowClearance"] = summary(values) if values else None
        classes[klass] = entry

    table = {
        "_comment": ("Pin synthesis rules for Bobbin::create_quick_bobbin (ABT #1220). GENERATED by "
                     "scripts/derive_quick_bobbin_pins.py from MAS/data/bobbins.ndjson + core_shapes.ndjson; "
                     "do not edit by hand. Lengths in metres. MKF uses the median of every summary."),
        "_definitions": {
            "grid": "IEC 62317 2.54 mm grid; rowDistance is rounded UP to a multiple of it. Evidence: gridEvidence.",
            "pitchBands": ("pitch by core width A (core.get_width()): catalogue pitches snapped to 1.27 mm steps, mode per "
                           "catalogue core width (ties -> smaller pitch), equal neighbours merged, boundaries half-way between "
                           "catalogue widths"),
            "edgeMargin": "(core width - longest row span) / 2: distance from the core's end face to the outermost pin centre",
            "verticalRowClearance": "rowDistance / 2 - core depth / 2: rows straddle the core depth (Z)",
            "horizontalRowClearance": "rowDistance / 2 - core window height / 2: one row per end flange (Y)",
            "pinDiameter": "pinDescription.dimensions[0] of every distinct footprint that carries one",
            "pinLength": "pinDescription.dimensions[2] of every distinct footprint that carries one",
            "railStandoff": ("H1/2 + H3 - core height/2 of the vertical PQ records, one value per distinct (shape, H1, H3, a, b, b1) drawing: the pin rail's underside (the pin standoff) beyond the core's "
                             "outer face along the pin direction. H3 (bottom flange outer face -> pin standoff) was measured at scale "
                             "on the Miles-Platts PQ0010..PQ0080 drawings (ABT #1249; ABT #1207's c - H1/2 included the tab height). "
                             "It is the ONLY rail datum in the catalogue; no horizontal record locates its rail, so MKF applies the "
                             "same standoff beyond the core depth face for horizontal pins, and says so."),
            "railOuterWall": ("b/2 - rowDistance/2 - pin diameter/2 of the same records: rail plastic between a pin's surface and the "
                              "rail's outer edge, across the row. Bobbin::get_pin_rails sizes a quick bobbin's rail with it."),
            "railInnerWall": ("rowDistance/2 - (b/2 - b1) - pin diameter/2 of the same records: rail plastic between a pin's surface "
                              "and the rail's inner edge (the side facing the core), across the row."),
            "railEndWall": ("a/2 - longest row span/2 - pin diameter/2 of the same records: rail plastic beyond the outermost pin's "
                            "surface, along the row. A quick bobbin's rail is one block per row (its pinout has no central gap); "
                            "the Miles-Platts rails' centre gap a1 is not carried over. Rail height: from the bottom flange's outer "
                            "face down to the pin standoff, as on the drawings (no separate legs)."),
        },
        "grid": GRID,
        "gridEvidence": {"rowPitches": grid_total, "onHalfGridWithin": SNAP_TOLERANCE, "onHalfGrid": grid_hits},
        "minimumPinsPerRow": 2,
        "recordsWithPinGeometry": used,
        "distinctFootprints": len(footprints),
        "recordsSkipped": dict(sorted(skipped.items())),
        "classes": classes,
        "pinDiameter": summary(pin_diameters),
        "pinLength": summary(pin_lengths),
        "railStandoff": summary(rail_standoffs),
        "railOuterWall": summary(rail_outer_walls),
        "railInnerWall": summary(rail_inner_walls),
        "railEndWall": summary(rail_end_walls),
    }
    text = json.dumps(table, indent=4) + "\n"
    if "--check" in sys.argv:
        if OUTPUT.read_text() != text:
            print(OUTPUT, "differs from the derivation", file=sys.stderr)
            return 1
        print(OUTPUT, "reproduces")
        return 0
    OUTPUT.write_text(text)
    print("wrote", OUTPUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
