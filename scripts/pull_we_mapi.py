#!/usr/bin/env python3
"""ABT #357 phase 1: pull the Wuerth WE-MAPI catalogue from RedExpert's public
JSON API and produce DRAFT MAS 'magnetic' stubs carrying
manufacturerInfo.datasheetInfo (subtype "inductor").

Re-runnable: cached pulls are reused unless --force is given.

Outputs (all in OUT_DIR):
  product_list_4.json            raw catalogue pull (all power inductors)
  tc_measurements_48_4.json      raw temperature-id list for chart type 48 (L vs I)
  tc_values_48_4_<id>.json       raw L(I) curves per temperature id
  we_mapi_datasheet_records.ndjson   one MAS magnetic stub per WE-MAPI part
  we_mapi_unmapped_model_params.ndjson  RedExpert equivalent-circuit + warming
      coefficients per part (NO slot in the MAS schema yet -- kept OUT of the
      records on purpose; see magneticDatasheetInfo.model oneOf = chipBead only)

Usage: python3 pull_we_mapi.py [--force]
"""

import json
import re
import subprocess
import sys
import urllib.request
from datetime import date
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent
MAS_SCHEMAS = Path("/home/alf/OpenMagnetics/MAS/schemas")
PEAS_SCHEMAS = Path("/home/alf/OpenMagnetics/PEAS/schemas")

BASE = "https://redexpert.we-online.com/redexpert"
UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/151.0.0.0 Safari/537.36")

# Temperatures to pull L(I) curves at (Celsius). Ids resolved from the
# measurements endpoint at runtime.
CURVE_TEMPERATURES_C = [20, 100, -40]
MAX_POINTS_PER_TEMPERATURE = 15

MANUFACTURER = "Würth Elektronik"
FAMILY = "WE-MAPI"
MODULE_ID = 4          # RedExpert module 4 = Power Inductors
CHART_TYPE = 48        # chart type 48 = L vs I


# ----------------------------------------------------------------------------
# Pull helpers
# ----------------------------------------------------------------------------

def fetch(url: str, dest: Path, force: bool) -> dict | list:
    """GET url -> dest (cached), strip control chars, return parsed JSON."""
    if force or not dest.exists() or dest.stat().st_size == 0:
        req = urllib.request.Request(url, headers={
            "User-Agent": UA,
            "Accept": "application/json",
            "Accept-Encoding": "gzip",
        })
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = resp.read()
        if body[:2] == b"\x1f\x8b":  # gzip-sniff
            import gzip
            body = gzip.decompress(body)
        dest.write_bytes(body)
        print(f"pulled {url} -> {dest.name} ({len(body)} bytes)")
    raw = dest.read_bytes().decode("utf-8", "replace")
    raw = re.sub(r"[\x00-\x1f]", "", raw)  # raw control chars break json.loads
    return json.loads(raw)


# ----------------------------------------------------------------------------
# Conversion
# ----------------------------------------------------------------------------

def require(part: dict, key: str):
    """No silent fallbacks: a missing required source field is an error."""
    value = part.get(key)
    if value is None or value == "":
        raise ValueError(f"part {part.get('Order_Code')}: source field {key!r} "
                         "is missing/empty in the RedExpert record")
    return value


def downsample(points: list, limit: int) -> list:
    """Keep first+last, thin the middle evenly to <= limit points."""
    if len(points) <= limit:
        return points
    step = (len(points) - 1) / (limit - 1)
    indices = sorted({round(i * step) for i in range(limit)})
    return [points[i] for i in indices]


def convert_part(part: dict, curves_by_temperature: dict) -> tuple[dict, dict]:
    """RedExpert catalogue record -> (MAS magnetic stub, unmapped sidecar)."""
    order_code = require(part, "Order_Code")
    inductance = require(part, "Inductance")            # already Henries
    tol_minus = require(part, "TolMin")                 # fraction, e.g. 0.2
    tol_plus = require(part, "TolMax")

    electrical = {
        "subtype": "inductor",
        "inductance": {
            "nominal": inductance,
            "minimum": inductance * (1.0 - tol_minus),
            "maximum": inductance * (1.0 + tol_plus),
        },
        "dcResistance": {
            "nominal": require(part, "Resistance_DC_TYP"),   # Ohm, typical
            "maximum": require(part, "Resistance_DC_MAX"),   # Ohm, max
        },
        "ratedCurrents": [require(part, "Rated_Current")],
        "ratedCurrentPoints": [{
            "current": require(part, "Rated_Current"),
            # Rated_Current_Temperature is the deltaT criterion (K), 40 for MAPI
            "temperatureRise": require(part, "Rated_Current_Temperature"),
        }],
        "saturationCurrents": [{
            # Saturation_Current_DROP is a fraction (0.2 -> 20 %)
            "percentInductanceDrop": require(part, "Saturation_Current_DROP") * 100.0,
            "current": require(part, "Saturation_Current"),
        }],
        "selfResonantFrequency": require(part, "Self_Resonan_Frequency"),  # Hz
    }

    inductance_points = []
    reference_curve = curves_by_temperature.get(CURVE_TEMPERATURES_C[0], {}).get(part["ID"])
    for temperature in CURVE_TEMPERATURES_C:
        curve = curves_by_temperature.get(temperature, {}).get(part["ID"])
        if curve is None:
            continue  # coverage gap reported separately, not fatal
        if temperature != CURVE_TEMPERATURES_C[0] and curve == reference_curve:
            # RedExpert replicates the room-temperature measurement under every
            # temperature id for parts measured at a single temperature. Do not
            # fabricate temperature dependence -- keep only the 20 C points.
            continue
        for current_a, inductance_uh in downsample(curve, MAX_POINTS_PER_TEMPERATURE):
            inductance_points.append({
                "inductance": inductance_uh * 1e-6,  # RedExpert curves are uH
                "current": current_a,
                "temperature": temperature,
            })
    if inductance_points:
        electrical["inductancePoints"] = inductance_points

    datasheet_info = {
        "part": {
            "partNumber": order_code,
            "family": FAMILY,
            "caseCode": require(part, "Size"),
            "material": require(part, "Core_Material"),
            "numberOfWindings": 1,
            "shielded": require(part, "Shielding_Type") == "Shielded",
            "automotive": part.get("Is_AECQ_Text") == "1",
        },
        "electrical": [electrical],
        "thermal": {
            "operatingTemperature": {
                "minimum": require(part, "Temperature_MIN"),
                "maximum": require(part, "Temperature_OPT"),
            },
        },
        "mechanical": {
            "length": {"nominal": require(part, "Size_Length"),
                       "maximum": require(part, "Size_Length_Max")},
            "width": {"nominal": require(part, "Size_Width"),
                      "maximum": require(part, "Size_Width_Max")},
            "height": {"nominal": require(part, "Size_Height"),
                       "maximum": require(part, "Size_Height_Max")},
            "mounting": "smt",
            "assemblyType": require(part, "Assembling_Technology"),
        },
        "provenance": [{
            "source": "manufacturerDatabase",
            "sourceName": "Würth RedExpert JSON API",
            "sourceUrl": f"{BASE}/product/list/{MODULE_ID}",
            "retrievedDate": date.today().isoformat(),
        }],
    }

    record = {
        "manufacturerInfo": {
            "name": MANUFACTURER,
            "reference": order_code,
            "orderCode": order_code,
            "family": FAMILY,
            "status": "production",
            # canonical WE pattern, HEAD-verified (application/pdf)
            "datasheetUrl": ("https://www.we-online.com/components/products/"
                             f"datasheet/{order_code}.pdf"),
            "datasheetInfo": datasheet_info,
        },
    }

    # RedExpert fields with NO slot in the current MAS magneticDatasheetInfo
    # schema (model oneOf only defines the chip-bead variant). Kept out of the
    # record -- never emit schema-invalid objects -- and parked in a sidecar
    # for a future inductor model variant.
    sidecar = {
        "orderCode": order_code,
        "productId": part["ID"],
        "equivalentCircuit": {  # parallel tank fitted by WE
            "parallelResistance": part.get("Rp"),      # Ohm
            "parallelCapacitance": part.get("Cp"),     # F
            "parallelInductance": part.get("Lp"),      # H
            "r1": part.get("R1"),                      # Ohm
        },
        "warming": {  # WE self-heating model: deltaT = factor * P^exp (loss-based)
            "warmingFactor": part.get("Warming_Factor"),
            "warmingExp": part.get("Warming_Exp"),
            "rv": part.get("Rv"),
            "rt": part.get("Rt"),
        },
        "voltagePeak": part.get("Voltage_Peak"),       # V, rated peak voltage
        "footprint": part.get("Footprint"),            # m^2 pad area
    }
    return record, sidecar


# ----------------------------------------------------------------------------
# Validation (Draft 2020-12, refs resolved from local MAS + PEAS checkouts)
# ----------------------------------------------------------------------------

def build_registry():
    from referencing import Registry, Resource
    from referencing.jsonschema import DRAFT202012

    resources = []
    for root, base_uri in [(MAS_SCHEMAS, "https://psma.com/mas/"),
                           (PEAS_SCHEMAS, "https://psma.com/peas/")]:
        for path in root.rglob("*.json"):
            contents = json.loads(path.read_text())
            uri = base_uri + str(path.relative_to(root))
            resources.append((uri, Resource(contents=contents,
                                            specification=DRAFT202012)))
            declared = contents.get("$id")
            if declared and declared != uri:
                resources.append((declared, Resource(contents=contents,
                                                     specification=DRAFT202012)))
    return Registry().with_resources(resources)


def build_validators():
    import jsonschema

    registry = build_registry()
    magnetic_schema = json.loads((MAS_SCHEMAS / "magnetic.json").read_text())
    full = jsonschema.Draft202012Validator(magnetic_schema, registry=registry)
    datasheet_info_only = jsonschema.Draft202012Validator(
        {"$ref": "https://psma.com/mas/magnetic.json#/$defs/magneticDatasheetInfo"},
        registry=registry)
    return full, datasheet_info_only


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> int:
    force = "--force" in sys.argv

    catalogue = fetch(f"{BASE}/product/list/{MODULE_ID}",
                      OUT_DIR / f"product_list_{MODULE_ID}.json", force)["Data"]
    mapi = [p for p in catalogue if p.get("Series") == FAMILY]
    print(f"catalogue: {len(catalogue)} power inductors, {len(mapi)} {FAMILY}")

    measurements = fetch(f"{BASE}/tc/measurements/{CHART_TYPE}/{MODULE_ID}",
                         OUT_DIR / f"tc_measurements_{CHART_TYPE}_{MODULE_ID}.json",
                         force)
    temperature_ids = {m["Values"][0]: m["ID"] for m in measurements}

    curves_by_temperature = {}
    for temperature in CURVE_TEMPERATURES_C:
        if temperature not in temperature_ids:
            raise ValueError(f"no RedExpert measurement id for {temperature} C; "
                             f"available: {sorted(temperature_ids)}")
        meas_id = temperature_ids[temperature]
        curves = fetch(f"{BASE}/tc/valuesxms/{CHART_TYPE}/{MODULE_ID}/{meas_id}",
                       OUT_DIR / f"tc_values_{CHART_TYPE}_{MODULE_ID}_{meas_id}.json",
                       force)
        curves_by_temperature[temperature] = {c["ID"]: c["Values"] for c in curves}

    full_validator, dsi_validator = build_validators()

    records_path = OUT_DIR / "we_mapi_datasheet_records.ndjson"
    sidecar_path = OUT_DIR / "we_mapi_unmapped_model_params.ndjson"
    n_full_ok = n_dsi_ok = 0
    missing_curves = {t: [] for t in CURVE_TEMPERATURES_C}
    errors = []

    with records_path.open("w") as rec_f, sidecar_path.open("w") as side_f:
        for part in sorted(mapi, key=lambda p: p["Order_Code"]):
            record, sidecar = convert_part(part, curves_by_temperature)
            for temperature in CURVE_TEMPERATURES_C:
                if part["ID"] not in curves_by_temperature[temperature]:
                    missing_curves[temperature].append(part["Order_Code"])

            full_errors = list(full_validator.iter_errors(record))
            if not full_errors:
                n_full_ok += 1
            else:
                errors.append((part["Order_Code"], "full",
                               [e.message for e in full_errors[:3]]))
            dsi = record["manufacturerInfo"]["datasheetInfo"]
            dsi_errors = list(dsi_validator.iter_errors(dsi))
            if not dsi_errors:
                n_dsi_ok += 1
            else:
                errors.append((part["Order_Code"], "datasheetInfo",
                               [e.message for e in dsi_errors[:3]]))

            rec_f.write(json.dumps(record, ensure_ascii=False) + "\n")
            side_f.write(json.dumps(sidecar, ensure_ascii=False) + "\n")

    print(f"wrote {records_path} ({len(mapi)} records)")
    print(f"wrote {sidecar_path} (equivalent circuit + warming, unmapped)")
    print(f"validation: {n_full_ok}/{len(mapi)} pass as full MAS 'magnetic', "
          f"{n_dsi_ok}/{len(mapi)} datasheetInfo pass against its $def")
    for temperature, codes in missing_curves.items():
        if codes:
            print(f"missing L(I) curve at {temperature} C: {len(codes)} parts: "
                  f"{codes}")
    for order_code, kind, messages in errors[:20]:
        print(f"INVALID [{kind}] {order_code}: {messages}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
