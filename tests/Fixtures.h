#pragma once

#include "constructive_models/Mas.h"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace OpenMagneticsTesting::fixtures {

// =====================================================================
// Test fixture registry.
//
// Tests that hold MAS / magnetic / converter JSON fixtures inline as
// huge raw-string literals make their .cpp files slow to compile,
// bloat MKF_tests' .rodata, and hide what's actually fixture-data
// versus test logic. This registry loads those fixtures from on-disk
// NDJSON files once per process and serves them by name.
//
// On-disk layout: tests/fixtures/*.ndjson, one fixture per line, each
// line a JSON object with shape:
//     {"name": "<unique fixture name>", "data": { ... }}
//
// The first call to `get_*` does a single walk of tests/fixtures/,
// validates structure (every line is an object with a string "name"
// and a "data" field, no duplicate names), and populates an in-memory
// map. Subsequent calls are constant-time lookups. A missing fixture
// throws — there is no silent fallback.
//
// Why NDJSON, not a single JSON array:
//   * one-fixture-per-line keeps git diffs readable when a fixture
//     changes,
//   * file-level chunking by topology / test-area is natural,
//   * append-only growth doesn't reflow the whole document,
//   * partial-parse failures point at a specific line number, not
//     "somewhere in the array".
// =====================================================================

// Return the raw JSON of a named fixture. Throws std::runtime_error if
// the fixture is not registered (with a list of close matches if any).
nlohmann::json get_json(const std::string& name);

// Return a fully-loaded MAS bundle from a fixture whose "data" payload
// matches the shape that TestingUtils::mas_loader accepts (top-level
// keys "inputs", "magnetic", "outputs"). Throws if either the fixture
// name is unknown or the payload doesn't load as a MAS bundle.
OpenMagnetics::Mas get_mas(const std::string& name);

// Return every registered fixture name, sorted. Useful for diagnostics
// and to validate (in a startup smoke test) that the registry actually
// found the fixtures the test suite expects.
std::vector<std::string> list_names();

// Force-clear and re-load the registry. Only intended for tests of the
// registry itself; regular tests should not call this.
void reload_for_testing();

}  // namespace OpenMagneticsTesting::fixtures
