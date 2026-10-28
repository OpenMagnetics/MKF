#!/usr/bin/env python3
"""Apply the web's legacy-default-voltage migration to a stored MAS test fixture.

The web app's old default excitation carried a -20.5 / +70.5 V primary voltage (mean +25 V,
against its own processed data: offset 0, 100 V peak to peak, duty 0.5), and the ABT #1670
reflection bug gave secondaries 0 / -100 V. MKF refuses both as not volt-second balanced
(ABT #1585, #1718). The web replaces exactly these waveforms when it loads a saved design
(WebFrontend src/services/loadMasIntoApp.js, migrateLegacyDefaultVoltages, ABT #1719). This
script is a line-for-line port of that function, so a fixture captured from the web can be
given the voltages the web itself would send today:

- a winding whose voltage is the legacy -20.5 / +70.5 V rectangle (as stored, or as the
  engine's equidistant samples of it) gets the web's default voltage object
  (Defaults.defaultOperatingPointExcitation.voltage: -50 / +50 V on the same switching
  instants, processed {dutyCycle 0.5, peakToPeak 100, offset 0, label rectangular});
- a secondary whose voltage is the legacy 0 / -100 V (or that divided by the turns ratio)
  gets {waveform: primary data / turns ratio, primary time}; its stale processed data and
  harmonics are dropped, as the web drops them (MKF recomputes both);
- the old two-switch-forward example secondary gets its reset interval ended where the
  volt-seconds cancel.

Only the "voltage" values of the excitations are rewritten, spliced into the original text,
so every other byte of the file stays as it was.

Usage:  python3 scripts/migrate_legacy_default_voltages.py <input.json> <output.json>
"""

import json
import sys

LEGACY_DEFAULT_PRIMARY_VOLTAGE = {"data": [-20.5, 70.5, 70.5, -20.5, -20.5], "timeFractions": [0, 0, 0.5, 0.5, 1]}
LEGACY_ABT1670_SECONDARY_VOLTAGE = {"data": [0, 0, -100, -100, 0], "timeFractions": [0, 0.5, 0.5, 1, 1]}
LEGACY_FORWARD_EXAMPLE_SECONDARY_VOLTAGE = {
    "data": [0, 7.291333333333334, 7.291333333333334, -6.135333333333334, -6.135333333333334, 0, 0],
    "timeFractions": [0, 0, 0.45, 0.45, 0.855, 0.855, 1],
    "balancedResetEndFraction": 0.9847876,
}
# WebSharedComponents/assets/js/defaults.js, defaultOperatingPointExcitation.voltage
DEFAULT_VOLTAGE = {
    "waveform": {"data": [-50, 50, 50, -50, -50], "time": [0, 0, 0.000005, 0.000005, 0.00001]},
    "processed": {"dutyCycle": 0.5, "peakToPeak": 100, "offset": 0, "label": "rectangular"},
}


def same_samples(actual, expected, tolerance):
    return isinstance(actual, list) and len(actual) == len(expected) \
        and all(abs(a - e) <= tolerance for a, e in zip(actual, expected))


def is_legacy_voltage(voltage, legacy, scale=1):
    waveform = (voltage or {}).get("waveform")
    time = (waveform or {}).get("time")
    period = time[-1] if isinstance(time, list) and time else None
    return period is not None and period > 0 \
        and same_samples(waveform.get("data"), [v * scale for v in legacy["data"]], 1e-9 * abs(scale)) \
        and same_samples(time, [f * period for f in legacy["timeFractions"]], period * 1e-9)


def is_sampled_two_level_voltage(voltage, levels):
    data = ((voltage or {}).get("waveform") or {}).get("data")
    if not isinstance(data, list) or len(data) < 8:
        return False
    tolerance = 1e-9 * max(abs(levels[0]), abs(levels[1]))
    on = lambda value, level: abs(value - level) <= tolerance
    return all(on(v, levels[0]) or on(v, levels[1]) for v in data) \
        and any(on(v, levels[0]) for v in data) and any(on(v, levels[1]) for v in data)


def resolve_dimension(value, description):
    if isinstance(value, (int, float)):
        return value
    if value.get("nominal") is not None:
        return value["nominal"]
    if value.get("minimum") is not None and value.get("maximum") is not None:
        return (value["minimum"] + value["maximum"]) / 2
    if value.get("maximum") is not None:
        return value["maximum"]
    if value.get("minimum") is not None:
        return value["minimum"]
    raise ValueError(f"{description} has no nominal, minimum or maximum")


def migrate(mas):
    """Returns {(operatingPointIndex, windingIndex): new voltage} for every replaced voltage."""
    replaced = {}
    turns_ratios = mas["inputs"]["designRequirements"].get("turnsRatios") or []
    for op_index, operating_point in enumerate(mas["inputs"].get("operatingPoints") or []):
        excitations = operating_point.get("excitationsPerWinding")
        if not isinstance(excitations, list):
            continue
        for w, excitation in enumerate(excitations):
            voltage = excitation.get("voltage")
            if is_legacy_voltage(voltage, LEGACY_FORWARD_EXAMPLE_SECONDARY_VOLTAGE):
                time = list(voltage["waveform"]["time"])
                time[4] = LEGACY_FORWARD_EXAMPLE_SECONDARY_VOLTAGE["balancedResetEndFraction"] * time[-1]
                time[5] = time[4]
                excitation["voltage"] = {"waveform": {"data": list(voltage["waveform"]["data"]), "time": time}}
                replaced[(op_index, w)] = excitation["voltage"]
                continue
            corrected = json.loads(json.dumps(DEFAULT_VOLTAGE))
            if is_legacy_voltage(voltage, LEGACY_DEFAULT_PRIMARY_VOLTAGE):
                corrected["waveform"]["time"] = list(voltage["waveform"]["time"])
            elif is_sampled_two_level_voltage(voltage, [-20.5, 70.5]) and (excitation.get("frequency") or 0) > 0:
                period = 1 / excitation["frequency"]
                corrected["waveform"]["time"] = [f * period for f in LEGACY_DEFAULT_PRIMARY_VOLTAGE["timeFractions"]]
            else:
                continue
            excitation["voltage"] = corrected
            replaced[(op_index, w)] = corrected
        primary_voltage = excitations[0].get("voltage") if excitations else None
        for w in range(1, len(excitations)):
            excitation = excitations[w]
            voltage = excitation.get("voltage")
            if not (voltage or {}).get("waveform") or not (primary_voltage or {}).get("waveform"):
                continue
            is_unscaled = is_legacy_voltage(voltage, LEGACY_ABT1670_SECONDARY_VOLTAGE) \
                or is_sampled_two_level_voltage(voltage, [0, -100])
            is_scaled = False
            if not is_unscaled and len(turns_ratios) > w - 1 and turns_ratios[w - 1] is not None:
                scale = 1 / resolve_dimension(turns_ratios[w - 1], f"Turns ratio {w}")
                is_scaled = is_legacy_voltage(voltage, LEGACY_ABT1670_SECONDARY_VOLTAGE, scale) \
                    or is_sampled_two_level_voltage(voltage, [0, -100 * scale])
            if not is_unscaled and not is_scaled:
                continue
            if len(turns_ratios) <= w - 1:
                raise ValueError(f"Turns ratio {w} is missing")
            turns_ratio = resolve_dimension(turns_ratios[w - 1], f"Turns ratio {w}")
            excitation["voltage"] = {"waveform": {
                "data": [v / turns_ratio for v in primary_voltage["waveform"]["data"]],
                "time": list(primary_voltage["waveform"]["time"]),
            }}
            replaced[(op_index, w)] = excitation["voltage"]
    return replaced


# --- splice: locate the byte span of each replaced "voltage" value in the original text ---

def _skip_ws(text, i):
    while text[i] in " \t\r\n":
        i += 1
    return i


def _value_span(text, i, path, wanted, spans):
    """Parses the JSON value at i; records spans of values whose path is in wanted. Returns end."""
    i = _skip_ws(text, i)
    if path in wanted:
        start = i
    c = text[i]
    if c == "{":
        i = _skip_ws(text, i + 1)
        if text[i] == "}":
            end = i + 1
        else:
            while True:
                key, i = json.decoder.scanstring(text, _skip_ws(text, i) + 1)
                i = _skip_ws(text, i)
                assert text[i] == ":"
                i = _value_span(text, i + 1, path + (key,), wanted, spans)
                i = _skip_ws(text, i)
                if text[i] == ",":
                    i += 1
                    continue
                assert text[i] == "}"
                end = i + 1
                break
    elif c == "[":
        i = _skip_ws(text, i + 1)
        if text[i] == "]":
            end = i + 1
        else:
            index = 0
            while True:
                i = _value_span(text, i, path + (index,), wanted, spans)
                i = _skip_ws(text, i)
                index += 1
                if text[i] == ",":
                    i += 1
                    continue
                assert text[i] == "]"
                end = i + 1
                break
    else:
        _, end = json.JSONDecoder().raw_decode(text, i)
    if path in wanted:
        spans[path] = (start, end)
    return end


def _js_numbers(value):
    """Integral floats as JavaScript's JSON.stringify writes them (-50, not -50.0)."""
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, list):
        return [_js_numbers(v) for v in value]
    if isinstance(value, dict):
        return {k: _js_numbers(v) for k, v in value.items()}
    return value


def main():
    source, target = sys.argv[1], sys.argv[2]
    text = open(source, encoding="utf-8").read()
    mas = json.loads(text)
    replaced = migrate(mas)
    if not replaced:
        raise SystemExit(f"{source}: no legacy default voltage found; nothing to migrate")
    wanted = {("inputs", "operatingPoints", op, "excitationsPerWinding", w, "voltage"): v for (op, w), v in replaced.items()}
    spans = {}
    _value_span(text, 0, (), wanted, spans)
    assert set(spans) == set(wanted)
    # indent unit of the file: the indentation of its first key
    first_line = text.split("\n")[1]
    unit = len(first_line) - len(first_line.lstrip(" "))
    out = text
    for path in sorted(spans, key=lambda p: spans[p][0], reverse=True):
        start, end = spans[path]
        line_start = out.rfind("\n", 0, start) + 1
        base = out[line_start:start]
        base = " " * (len(base) - len(base.lstrip(" ")))
        rendered = json.dumps(_js_numbers(wanted[path]), indent=unit, ensure_ascii=False).replace("\n", "\n" + base)
        out = out[:start] + rendered + out[end:]
    # the spliced document must parse to exactly the migrated document
    assert json.loads(out) == mas
    open(target, "w", encoding="utf-8").write(out)
    for (op, w), v in sorted(replaced.items()):
        print(f"{target}: operating point {op}, winding {w}: voltage -> {v['waveform']['data']}")


if __name__ == "__main__":
    main()
