# PLECS Export Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small CLI binary `plecs_export_example` that turns a MAS `Magnetic` JSON file into a `.plecs` snippet (subcircuit or symbol), and establish a top-level `examples/` directory in MKF.

**Architecture:** Single `main.cpp` (~60 lines) under a new `examples/plecs_export/` directory. It loads a JSON file, constructs `OpenMagnetics::Magnetic` via its `nlohmann::json`-accepting constructor, calls the existing `CircuitSimulatorExporterPlecsModel::export_magnetic_as_subcircuit` or `export_magnetic_as_symbol`, and writes the returned string to disk. A new `examples/CMakeLists.txt` registers the target; the root `CMakeLists.txt` is gated on a new `INCLUDE_MKF_EXAMPLES` option mirroring the existing `INCLUDE_MKF_TESTS` pattern.

**Tech Stack:** C++23, CMake, `nlohmann::json` (already in MKF), `MKF` shared/static library target, no flag-parsing library.

**Spec:** `docs/superpowers/specs/2026-04-06-plecs-export-example-design.md`

---

## File Structure

| File | Action | Purpose |
|------|--------|---------|
| `examples/CMakeLists.txt` | Create | Adds the `plecs_export` subdirectory |
| `examples/plecs_export/CMakeLists.txt` | Create | Defines the `plecs_export_example` executable target |
| `examples/plecs_export/main.cpp` | Create | The CLI itself |
| `CMakeLists.txt` | Modify (root, around line 490) | Add `INCLUDE_MKF_EXAMPLES` option + `add_subdirectory(examples)` gate |

No source files under `src/` are modified. No test files are added.

---

## Task 1: Create the `examples/plecs_export/main.cpp` CLI

**Files:**
- Create: `examples/plecs_export/main.cpp`

- [ ] **Step 1: Create the file with the full CLI**

Write `examples/plecs_export/main.cpp` with this exact content:

```cpp
// PLECS export example for MKF.
//
// Usage:
//   plecs_export_example <input.json> <output.plecs> <subcircuit|symbol> [frequency_Hz] [temperature_C]
//
// Reads a MAS Magnetic JSON file and writes a .plecs snippet that can be
// copy-pasted into a larger PLECS schematic file.

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include <nlohmann/json.hpp>

#include "processors/CircuitSimulatorInterface.h"

namespace {

void print_usage(std::ostream& os) {
    os << "Usage: plecs_export_example <input.json> <output.plecs> "
       << "<subcircuit|symbol> [frequency_Hz] [temperature_C]\n"
       << "  subcircuit mode requires frequency_Hz and temperature_C.\n"
       << "  symbol mode ignores frequency_Hz and temperature_C.\n";
}

int fail(const std::string& message) {
    std::cerr << "error: " << message << "\n";
    print_usage(std::cerr);
    return 1;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        return fail("not enough arguments");
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    const std::string mode = argv[3];

    if (mode != "subcircuit" && mode != "symbol") {
        return fail("mode must be 'subcircuit' or 'symbol' (got '" + mode + "')");
    }

    double frequency = 0.0;
    double temperature = 0.0;
    if (mode == "subcircuit") {
        if (argc < 6) {
            return fail("subcircuit mode requires frequency_Hz and temperature_C");
        }
        try {
            frequency = std::stod(argv[4]);
            temperature = std::stod(argv[5]);
        } catch (const std::exception& e) {
            return fail(std::string{"could not parse frequency/temperature: "} + e.what());
        }
    }

    std::ifstream input_file(input_path);
    if (!input_file.is_open()) {
        return fail("could not open input file: " + input_path);
    }

    nlohmann::json j;
    try {
        input_file >> j;
    } catch (const std::exception& e) {
        return fail(std::string{"failed to parse JSON: "} + e.what());
    }

    std::string content;
    try {
        OpenMagnetics::Magnetic magnetic(j);
        OpenMagnetics::CircuitSimulatorExporterPlecsModel exporter;
        if (mode == "subcircuit") {
            content = exporter.export_magnetic_as_subcircuit(magnetic, frequency, temperature);
        } else {
            content = exporter.export_magnetic_as_symbol(magnetic);
        }
    } catch (const std::exception& e) {
        return fail(std::string{"PLECS export failed: "} + e.what());
    }

    std::ofstream output_file(output_path);
    if (!output_file.is_open()) {
        return fail("could not open output file for writing: " + output_path);
    }
    output_file << content;
    if (!output_file) {
        return fail("failed while writing output file: " + output_path);
    }

    std::cout << "Wrote " << output_path << "\n";
    return 0;
}
```

- [ ] **Step 2: Commit**

```bash
git add examples/plecs_export/main.cpp
git commit -m "feat(examples): add plecs_export_example CLI source

Small CLI that loads a MAS Magnetic JSON and emits a .plecs snippet
(subcircuit or symbol) using CircuitSimulatorExporterPlecsModel.
Build wiring follows in subsequent commits.

Refs: #57"
```

---

## Task 2: Add per-target CMakeLists for the example

**Files:**
- Create: `examples/plecs_export/CMakeLists.txt`

- [ ] **Step 1: Create the per-target CMakeLists**

Write `examples/plecs_export/CMakeLists.txt` with this exact content:

```cmake
add_executable(plecs_export_example main.cpp)
target_link_libraries(plecs_export_example PRIVATE MKF)
# MKF's target_include_directories already exposes src/ and src/processors/
# as PUBLIC, so plecs_export_example inherits the include paths it needs.
set_target_properties(plecs_export_example PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
```

- [ ] **Step 2: Commit**

```bash
git add examples/plecs_export/CMakeLists.txt
git commit -m "build(examples): add CMake target for plecs_export_example

Defines the plecs_export_example executable, links against MKF, and
places the binary at \${CMAKE_BINARY_DIR} so the demo command line
matches ./build/plecs_export_example.

Refs: #57"
```

---

## Task 3: Add the `examples/` umbrella CMakeLists

**Files:**
- Create: `examples/CMakeLists.txt`

- [ ] **Step 1: Create the umbrella CMakeLists**

Write `examples/CMakeLists.txt` with this exact content:

```cmake
# MKF examples — small standalone binaries that demonstrate library usage.
# Add new examples by creating a subdirectory with its own CMakeLists.txt
# and listing it here.

add_subdirectory(plecs_export)
```

- [ ] **Step 2: Commit**

```bash
git add examples/CMakeLists.txt
git commit -m "build(examples): add umbrella CMakeLists for examples/

Aggregates example subdirectories. Currently only plecs_export.

Refs: #57"
```

---

## Task 4: Wire `examples/` into the root build

**Files:**
- Modify: `CMakeLists.txt` (root, around line 490 — the slot just after the commented-out `magnetic_adviser` block and before `if(INCLUDE_MKF_TESTS)`)

- [ ] **Step 1: Add the option and the gated `add_subdirectory`**

Find this region in the root `CMakeLists.txt` (around lines 490-493):

```cmake
# add_executable(magnetic_adviser "src/advisers/MagneticAdviser.cpp")
# target_link_libraries(magnetic_adviser MKF)

if(INCLUDE_MKF_TESTS)
```

Replace it with:

```cmake
# add_executable(magnetic_adviser "src/advisers/MagneticAdviser.cpp")
# target_link_libraries(magnetic_adviser MKF)

option(INCLUDE_MKF_EXAMPLES "Build small standalone example binaries under examples/" ON)
if(INCLUDE_MKF_EXAMPLES)
    add_subdirectory(examples)
endif()

if(INCLUDE_MKF_TESTS)
```

(Mirrors the existing `INCLUDE_MKF_TESTS` pattern. Default `ON` so the example builds by default; downstream users who don't want it set `-DINCLUDE_MKF_EXAMPLES=OFF`.)

- [ ] **Step 2: Configure CMake to verify the option and target are recognized**

Run from repo root:

```bash
cmake -S . -B build -G Ninja
```

Expected: configuration succeeds. The output should mention the new target `plecs_export_example` somewhere (or at least no error about `examples/CMakeLists.txt`). If the build dir already exists from prior work, this is a re-configure and should still succeed.

If configuration fails on the new lines, fix the syntax in `CMakeLists.txt` and re-run.

- [ ] **Step 3: Build only the new target**

```bash
cmake --build build --target plecs_export_example
```

Expected: a clean build of `main.cpp` plus any MKF dependencies, and `build/plecs_export_example` exists afterwards.

Verify the binary is present:

```bash
ls -l build/plecs_export_example
```

Expected: file exists and is executable.

Note: per `memory/feedback_no_full_compilation.md`, full MKF compilation is known to crash WSL. This task only builds the single new target, which compiles `main.cpp` and links against the already-built `MKF` library if it exists. If `MKF` itself has not yet been built in this checkout, the build will need to compile MKF as a dependency — in that case, stop and report this back rather than letting WSL die.

- [ ] **Step 4: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: gate examples/ subdirectory on INCLUDE_MKF_EXAMPLES

Mirrors the INCLUDE_MKF_TESTS option pattern. Defaults to ON so the
plecs_export_example binary builds by default; downstream consumers
can disable with -DINCLUDE_MKF_EXAMPLES=OFF.

Refs: #57"
```

---

## Task 5: Smoke-test the binary against an existing fixture

**Files:**
- None modified. This task only runs the binary.

- [ ] **Step 1: Run subcircuit mode against the existing Simba fixture**

```bash
./build/plecs_export_example \
    tests/testData/test_circuitsimulatorexporter_simba_json_87.json \
    /tmp/plecs_export_smoke_subcircuit.plecs \
    subcircuit 100e3 25
```

Expected stdout: `Wrote /tmp/plecs_export_smoke_subcircuit.plecs`
Expected exit code: 0
Expected file: `/tmp/plecs_export_smoke_subcircuit.plecs` exists and is non-empty.

Verify:

```bash
test -s /tmp/plecs_export_smoke_subcircuit.plecs && echo "non-empty OK"
head -3 /tmp/plecs_export_smoke_subcircuit.plecs
```

Expected: `non-empty OK` printed, and the first few lines look like a PLECS file (likely starting with `Plecs {` or a `%%` comment header from `assemble_plecs_file`).

- [ ] **Step 2: Run symbol mode against the same fixture**

```bash
./build/plecs_export_example \
    tests/testData/test_circuitsimulatorexporter_simba_json_87.json \
    /tmp/plecs_export_smoke_symbol.plecs \
    symbol
```

Expected stdout: `Wrote /tmp/plecs_export_smoke_symbol.plecs`
Expected exit code: 0
Expected file: `/tmp/plecs_export_smoke_symbol.plecs` exists and is non-empty.

Verify:

```bash
test -s /tmp/plecs_export_smoke_symbol.plecs && echo "non-empty OK"
```

- [ ] **Step 3: Run argument validation paths**

```bash
./build/plecs_export_example
echo "exit=$?"
```
Expected: error message + usage to stderr, exit code 1.

```bash
./build/plecs_export_example a.json b.plecs banana
echo "exit=$?"
```
Expected: `error: mode must be 'subcircuit' or 'symbol' (got 'banana')`, exit code 1.

```bash
./build/plecs_export_example a.json b.plecs subcircuit
echo "exit=$?"
```
Expected: `error: subcircuit mode requires frequency_Hz and temperature_C`, exit code 1.

```bash
./build/plecs_export_example /no/such/file.json /tmp/x.plecs symbol
echo "exit=$?"
```
Expected: `error: could not open input file: /no/such/file.json`, exit code 1.

- [ ] **Step 4: Clean up smoke-test artifacts**

```bash
rm -f /tmp/plecs_export_smoke_subcircuit.plecs /tmp/plecs_export_smoke_symbol.plecs
```

(Nothing to commit in this task — it's a manual verification task. If any step failed, fix the issue in the appropriate prior task and re-run from Task 4 step 2.)

---

## Self-Review Notes

**Spec coverage check:**

| Spec section | Covered by |
|--------------|------------|
| Layout (`examples/plecs_export/main.cpp`, two CMakeLists) | Tasks 1, 2, 3 |
| CLI positional shape | Task 1 step 1 (argv parsing) |
| `subcircuit` mode → `export_magnetic_as_subcircuit(magnetic, freq, temp)` | Task 1 step 1 |
| `symbol` mode → `export_magnetic_as_symbol(magnetic)` | Task 1 step 1 |
| Validation rules (argc, mode, freq/temp parse, file open, JSON parse, output write) | Task 1 step 1 + Task 5 step 3 |
| Example writes file itself (because `filePathOrFile` is `[[maybe_unused]]`) | Task 1 step 1 (`std::ofstream output_file`) |
| `OpenMagnetics::Magnetic(json)` constructor reuse | Task 1 step 1 |
| `INCLUDE_MKF_EXAMPLES` option mirroring `INCLUDE_MKF_TESTS` | Task 4 |
| `add_subdirectory(examples)` from root | Task 4 |
| Reuse existing `tests/testData/test_circuitsimulatorexporter_simba_json_87.json` | Task 5 |
| No new test, no new sample files | Confirmed by absence in plan |

**Type / signature consistency:**
- `OpenMagnetics::CircuitSimulatorExporterPlecsModel::export_magnetic_as_subcircuit(Magnetic, double frequency, double temperature)` — matches signature at `src/processors/CircuitSimulatorInterface.h:573` (later positional args use defaults).
- `OpenMagnetics::CircuitSimulatorExporterPlecsModel::export_magnetic_as_symbol(Magnetic)` — matches signature at `src/processors/CircuitSimulatorInterface.h:572`.
- `OpenMagnetics::Magnetic(nlohmann::json)` — verified in use at `tests/TestCircuitSimulatorInterface.cpp:82`.
- Include path `processors/CircuitSimulatorInterface.h` — `src/processors/` is on `MKF`'s PUBLIC include path (root `CMakeLists.txt:485`), so this resolves once we link against `MKF`.

**No placeholders:** all code shown in full, all commands shown in full, all expected outputs shown.
