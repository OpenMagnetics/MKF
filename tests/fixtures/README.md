# Test fixtures registry

External NDJSON-backed fixtures loaded once per test process by
`OpenMagneticsTesting::fixtures::get_json(...)` / `get_mas(...)`.
See `tests/Fixtures.h` for the API and `tests/Fixtures.cpp` for the
loader implementation.

## Why this exists

The MKF test binary is ~1 GB. Most of that is .rodata from huge inline
JSON literals — Magnetic, Inputs, Mas — pasted as `R"({...})"` raw
strings inside .cpp files. Each such literal:

  * is parsed at compile time (slow on heavy test TUs),
  * is embedded in the binary (bloats disk / debug info),
  * obscures the actual test logic by overwhelming it with data.

Moving the big literals into NDJSON keeps them as data, not code.

## On-disk layout

```
tests/fixtures/
    README.md                  # this file
    <topic>.ndjson             # one fixture per line, name + payload
```

Each line of an `.ndjson` file is a single JSON object:

```json
{"name": "painter-cmc-toroidal-pq2625", "data": { ... full Magnetic JSON ... }}
```

Blank lines and lines starting with `#` or `//` are skipped.

## Inventory of extraction candidates

(Sweep done 2026-05-23 of all `tests/*.cpp` raw-string JSON literals.)

| File | LOC | Long raw strings | Total embedded lines | Notes |
|---|---:|---:|---:|---|
| `TestPainter.cpp` | 5728 | 3 | 1829 | **prime target** — Magnetic JSON, biggest single 1050 lines |
| `TestCircuitSimulatorInterface.cpp` | 5243 | 1 | 2269 | **NOT a JSON fixture** — SPICE deck templates, different refactor |
| `TestCoreLosses.cpp` | 3767 | 1 | 210 | candidate |
| `TestTopologyBuck.cpp` | 738 | 1 | 148 | candidate |
| `TestMagneticAdviser.cpp` | 5459 | 1 | 112 | candidate |
| `TestTopologyFlyback.cpp` | 1883 | 1 | 111 | candidate |
| `TestCoilAdviser.cpp` | 2568 | 1 | 91 | candidate |

Realistic extractable JSON across the test suite: **~2500 lines**.

Note: most `TestTopology*.cpp` files don't use raw-string literals —
they construct fixtures procedurally via `llcJson["operatingPoints"] =
json::array();` etc. That's a separate, lower-priority refactor: those
fixtures are smaller per-test-case and the test logic is intertwined
with the fixture construction. Procedural fixtures can be moved by
serialising the constructed JSON to a fixture file and loading it
back, but only when the procedural construction is itself repetitive
(e.g. many test cases sharing the same shape with small tweaks).

## Adding a fixture

1. Run the test that uses the literal you want to extract, with the
   literal `json::parse(...)` output dumped to a file. Or hand-author
   the JSON.
2. Wrap it: `{"name": "<unique-kebab-case-name>", "data": <your json>}`.
3. Append the line to an `.ndjson` file under `tests/fixtures/`.
4. Replace the inline literal in the .cpp with:
   ```cpp
   auto magneticJson = OpenMagneticsTesting::fixtures::get_json("<unique-kebab-case-name>");
   ```
   Or use `get_mas(...)` if the fixture is a full MAS bundle (with
   `inputs`/`magnetic`/`outputs` keys).
5. Rebuild — fixture names that collide across files cause a startup
   error, not a silent overwrite. Missing fixtures throw at lookup
   time with close-match suggestions.

## Conventions

- **Names are unique across the whole `tests/fixtures/` tree.** The
  loader's duplicate-name detector fails loud at startup; this is
  intentional so renaming via "add a copy with the new name, leave
  the old one for now" can't silently happen.
- **Names are kebab-case.** Easier to grep for.
- **Group fixtures by topic in files, not by test name.** A fixture
  may be used by multiple tests — file it under what it *is*
  (`painter-bobbins`), not under who first used it.
- **No comments inside the JSON** (NDJSON parses strict JSON per
  line). Annotate via the surrounding `.ndjson` file with
  `#`-prefixed lines between fixtures.
