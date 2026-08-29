# PLECS Export Example — Design

**Date:** 2026-04-06
**Status:** Approved (brainstorming complete)
**Related:**
- Existing PLECS exporter: `src/processors/CircuitSimulatorPlecs.cpp` (#57)
- Reference pattern (Simba parallel exporter): `tests/TestCircuitSimulatorInterface.cpp:79-104`

## Goal

Provide a small, video-demo-friendly CLI binary that turns a MAS `Magnetic` JSON
file into a `.plecs` snippet, suitable for copy-paste into a larger PLECS file.

The example also establishes a new top-level `examples/` convention for MKF
(currently MKF only builds the library and `MKF_tests`).

## Non-goals

- A general-purpose CLI for the rest of the MKF API.
- Support for full MAS files with operating points (only the `Magnetic` JSON shape
  used by the existing Simba exporter test fixtures).
- A new test suite — the underlying PLECS exporter is already tested via
  `tests/TestCircuitSimulatorPlecs.cpp`.
- New sample JSON files — existing fixtures in `tests/testData/` are reused.

## CLI

Positional only, no flag-parsing library:

```
plecs_export_example <input.json> <output.plecs> <subcircuit|symbol> [frequency_Hz] [temperature_C]
```

Behavior:

| Mode         | Required positional args                                            | Underlying call |
|--------------|---------------------------------------------------------------------|-----------------|
| `subcircuit` | `input`, `output`, `subcircuit`, `frequency_Hz`, `temperature_C`    | `CircuitSimulatorExporterPlecsModel::export_magnetic_as_subcircuit(magnetic, freq, temp)` |
| `symbol`     | `input`, `output`, `symbol`                                         | `CircuitSimulatorExporterPlecsModel::export_magnetic_as_symbol(magnetic)` |

Validation rules:

- `argc < 4` → print usage to `stderr`, exit 1.
- `mode` not in `{subcircuit, symbol}` → error, exit 1.
- `mode == subcircuit` and `argc < 6` → error "subcircuit requires frequency and temperature", exit 1.
- `frequency_Hz` / `temperature_C` not parseable as `double` → error, exit 1.
- Input file open failure → error, exit 1.
- JSON parse failure → propagate `nlohmann::json::exception` message to `stderr`, exit 1.
- Output file open failure → error, exit 1.

On success: prints `Wrote <output_path>\n` to `stdout` and exits 0.

## Layout (new files)

```
examples/
├── CMakeLists.txt          # plecs_export subdirectory hook
└── plecs_export/
    ├── CMakeLists.txt      # add_executable(plecs_export_example main.cpp); target_link_libraries(... MKF)
    └── main.cpp            # ~50 lines
```

Root `CMakeLists.txt` change:

```cmake
add_subdirectory(examples)
```

If MKF already defines an option pattern like `MKF_BUILD_TESTS`, mirror it as
`MKF_BUILD_EXAMPLES` (default `ON`) and gate the `add_subdirectory` on it.
Otherwise, unconditional.

## main.cpp data flow

1. Parse positional args, validate per the rules above.
2. Open input: `std::ifstream in(input_path);`
3. Parse JSON: `auto j = nlohmann::json::parse(in);`
4. Construct magnetic: `OpenMagnetics::Magnetic magnetic(j);`
   (This one-arg constructor is the pattern used by the Simba exporter tests
   at `tests/TestCircuitSimulatorInterface.cpp:82`.)
5. Construct exporter: `OpenMagnetics::CircuitSimulatorExporterPlecsModel exporter;`
6. Call the matching `export_magnetic_as_*`:
   - `subcircuit` → `exporter.export_magnetic_as_subcircuit(magnetic, freq, temp);`
   - `symbol`     → `exporter.export_magnetic_as_symbol(magnetic);`
   The returned `std::string` is the full content to write.
7. Write to output: `std::ofstream out(output_path); out << content;`
8. Print success line, return 0.

### Why the example writes the file (not the exporter)

Both `CircuitSimulatorExporterPlecsModel::export_magnetic_as_subcircuit` and
`export_magnetic_as_symbol` declare their `filePathOrFile` parameter as
`[[maybe_unused]]` (verified at `src/processors/CircuitSimulatorPlecs.cpp:622`
and `:668`). The functions return the generated string but do **not** write it
to disk. The example is therefore responsible for the `std::ofstream`.

## Sample inputs (reused, not added)

The existing Simba exporter fixtures in `tests/testData/` are top-level
`Magnetic` JSONs and parse cleanly through the same one-arg constructor.
Recommended for the video demo:

- `tests/testData/test_circuitsimulatorexporter_simba_json_87.json` (basic E-core)
- `tests/testData/test_circuitsimulatorexporter_simba_json_ur_74.json` (U/R core)
- `tests/testData/test_circuitsimulatorexporter_simba_json_toroidal_core_100.json` (toroidal)
- `tests/testData/test_circuitsimulatorexporter_simba_json_ep_core_113.json` (EP core)

## Demo commands for the video

```bash
# Subcircuit (with curve fitting at 100 kHz, 25 °C)
./build/plecs_export_example \
    tests/testData/test_circuitsimulatorexporter_simba_json_87.json \
    /tmp/inductor.plecs subcircuit 100e3 25

# Symbol only
./build/plecs_export_example \
    tests/testData/test_circuitsimulatorexporter_simba_json_87.json \
    /tmp/inductor_symbol.plecs symbol
```

Open the resulting `.plecs` file in PLECS, or copy-paste its contents into a
larger `.plecs` schematic.

## Error handling summary

All errors go to `stderr` with a one-line message and a non-zero exit code.
No stack traces, no recovery. The example is short enough that failure modes
are obvious from context.

## Risks and assumptions

- **`OpenMagnetics::Magnetic(json)` constructor signature.** Verified to exist
  and accept a `nlohmann::json` at `tests/TestCircuitSimulatorInterface.cpp:82`.
  Implementation should `#include` whatever header that test includes.
- **CMake target linking.** The example links against the `MKF` library target.
  Need to confirm the exact target name in the root `CMakeLists.txt` during
  implementation (likely `MKF` based on `MKF_tests`).
- **`MKF_BUILD_EXAMPLES` option.** If MKF has no existing `MKF_BUILD_*` options,
  the example builds unconditionally. Confirm during implementation.

## Out of scope (YAGNI)

- `--help` / `-h` text beyond the usage line
- Long-form flags (`--input`, `--output`, …)
- A flag-parsing library
- Reading frequency/temperature from a MAS operating point in the JSON
- Generating fresh sample JSON files
- A new test case for the example binary
- Cross-platform packaging
