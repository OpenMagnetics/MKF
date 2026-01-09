# MKF Architecture Documentation

## Overview

MKF (Magnetics Framework) is a C++ library for modeling, simulating, and optimizing magnetic components such as inductors and transformers. It provides both physical models (for calculating electromagnetic properties) and constructive models (for describing the physical structure of components).

## Project Structure

```
MKF/
├── src/                      # Source code
│   ├── advisers/             # High-level component selection algorithms
│   ├── constructive_models/  # Physical structure representations
│   ├── converter_models/     # Power converter topology models
│   ├── physical_models/      # Electromagnetic calculations
│   ├── processors/           # Data processing utilities
│   ├── support/              # Utility classes and infrastructure
│   └── *.h                   # Main interface headers
├── tests/                    # Unit and integration tests
├── MAS/                      # Magnetic Assembly Schema (data definitions)
├── output/                   # Generated test outputs
├── scripts/                  # Build and utility scripts
└── docs/                     # Documentation
```

## Core Concepts

### MAS (Magnetic Assembly Schema)

MAS is a JSON Schema that defines the data structures for magnetic components. The MKF library uses [quicktype](https://quicktype.io/) to generate C++ classes from these schemas, providing type-safe access to component data.

Key schema files in `MAS/schemas/`:
- `MAS.json` - Root schema combining all component types
- `magnetic.json` - Magnetic component definitions
- `inputs.json` - Operating conditions and requirements
- `outputs.json` - Calculation results
- `utils.json` - Shared utility types

### Module Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Application Layer                         │
│  (Advisers: CoreAdviser, CoilAdviser, WireAdviser, etc.)         │
├─────────────────────────────────────────────────────────────────┤
│                         Physical Models                           │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────────┐ │
│  │ CoreLosses   │ │ WindingLosses│ │ MagnetizingInductance    │ │
│  ├──────────────┤ ├──────────────┤ ├──────────────────────────┤ │
│  │ LeakageInd.  │ │ StrayCapac.  │ │ Reluctance               │ │
│  └──────────────┘ └──────────────┘ └──────────────────────────┘ │
├─────────────────────────────────────────────────────────────────┤
│                       Constructive Models                         │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────────┐ │
│  │    Core      │ │    Coil      │ │       Magnetic           │ │
│  ├──────────────┤ ├──────────────┤ ├──────────────────────────┤ │
│  │    Wire      │ │   Bobbin     │ │    Insulation            │ │
│  └──────────────┘ └──────────────┘ └──────────────────────────┘ │
├─────────────────────────────────────────────────────────────────┤
│                         Support Layer                             │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────────┐ │
│  │  Databases   │ │   Settings   │ │     Exceptions           │ │
│  ├──────────────┤ ├──────────────┤ ├──────────────────────────┤ │
│  │   Logger     │ │   Painter    │ │        Utils             │ │
│  └──────────────┘ └──────────────┘ └──────────────────────────┘ │
├─────────────────────────────────────────────────────────────────┤
│                     Data Layer (MAS Schema)                       │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                    MAS.hpp (generated)                    │   │
│  └──────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

## Key Components

### Constructive Models (`src/constructive_models/`)

These classes represent the physical structure of magnetic components:

| Class | Purpose |
|-------|---------|
| `Core` | Magnetic core with shape, material, and gap information |
| `Coil` | Winding configuration with turns, layers, and wire selection |
| `Wire` | Conductor specifications (round, litz, rectangular, foil) |
| `Bobbin` | Mechanical framework for coil placement |
| `Magnetic` | Complete magnetic component (core + coil) |
| `Insulation` | Electrical isolation materials and configurations |

### Physical Models (`src/physical_models/`)

These classes calculate electromagnetic properties:

| Class | Purpose |
|-------|---------|
| `CoreLosses` | Core power dissipation (Steinmetz, iGSE, etc.) |
| `WindingLosses` | Conductor losses (DC, skin effect, proximity) |
| `MagnetizingInductance` | Primary inductance calculation |
| `LeakageInductance` | Inter-winding leakage inductance |
| `StrayCapacitance` | Parasitic capacitance between windings |
| `Reluctance` | Magnetic circuit reluctance calculations |
| `MagneticField` | Field distribution analysis |
| `Permeability` | Material permeability modeling |

### Advisers (`src/advisers/`)

High-level algorithms for component selection and optimization:

| Class | Purpose |
|-------|---------|
| `MagneticAdviser` | Complete magnetic design optimization |
| `CoreAdviser` | Core shape and material selection |
| `CoilAdviser` | Winding strategy optimization |
| `WireAdviser` | Conductor selection for efficiency |

### Support Infrastructure (`src/support/`)

| Class | Purpose |
|-------|---------|
| `DatabaseManager` | Singleton managing component databases |
| `Settings` | Global configuration singleton |
| `Exceptions` | Typed exception hierarchy with error codes |
| `Logger` | Configurable logging infrastructure |
| `Painter` | Visualization output (SVG, etc.) |
| `Utils` | Mathematical and utility functions |

## Data Flow

### Typical Design Flow

```
1. Input Specification
   └── DesignRequirements (voltage, current, frequency, etc.)
   
2. Component Selection (Advisers)
   ├── CoreAdviser → Select core shape & material
   ├── WireAdviser → Select conductor type
   └── CoilAdviser → Determine winding configuration

3. Construction (Constructive Models)
   ├── Core::process_data() → Calculate core geometry
   ├── Core::process_gap() → Set gap dimensions
   ├── Coil::wind() → Generate winding layout
   └── Magnetic::validate() → Verify assembly

4. Analysis (Physical Models)
   ├── CoreLosses::calculate() → Core power loss
   ├── WindingLosses::calculate() → Copper losses
   ├── MagnetizingInductance::calculate() → Inductance
   └── StrayCapacitance::calculate() → Parasitic effects

5. Output
   └── Outputs (losses, temperatures, efficiency)
```

## Database System

MKF uses embedded JSON databases for component data:

| Database | Contents |
|----------|----------|
| `cores.ndjson` | Pre-built core configurations |
| `core_shapes.ndjson` | Core geometry definitions |
| `core_materials.ndjson` | Magnetic material properties |
| `wires.ndjson` | Wire specifications |
| `wire_materials.ndjson` | Conductor material properties |
| `bobbins.ndjson` | Bobbin definitions |
| `insulation_materials.ndjson` | Insulation properties |

Access via `DatabaseManager`:
```cpp
auto& dbm = DatabaseManager::getInstance();
auto core = dbm.findCore("E 42/21/15");
```

## Exception Handling

The library uses a typed exception hierarchy (see `src/support/Exceptions.h`):

```cpp
try {
    core.process_gap();
} catch (const GapNotProcessedException& e) {
    // Handle missing gap data
} catch (const CoreNotProcessedException& e) {
    // Handle unprocessed core
} catch (const OpenMagneticsException& e) {
    // Handle any library exception
    std::cerr << "Error " << e.code() << ": " << e.what() << std::endl;
}
```

## Configuration

Global settings via `Settings` singleton:
```cpp
auto settings = Settings::GetInstance();
settings->set_use_toroidal_cores(false);
settings->set_painter_mode(PainterModes::SVG);
```

## Building

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DINCLUDE_MKF_TESTS=ON ..
cmake --build .
ctest --output-on-failure
```

## Testing

Tests use Catch2 framework:
```bash
./MKF_tests                     # Run all tests
./MKF_tests [Utils]             # Run tests with tag
./MKF_tests "CeilFloat"         # Run specific test
./MKF_tests --list-tests        # List available tests
```

## Dependencies

| Library | Purpose |
|---------|---------|
| nlohmann/json | JSON parsing |
| magic-enum | Enum reflection |
| Catch2 | Testing framework |
| matplot++ | Visualization |
| rapidfuzz | String matching |
| cmrc | Resource embedding |
| levmar | Levenberg-Marquardt optimization |
| spline | Interpolation |
| svg | SVG generation |

## Contributing

1. Follow the `.clang-format` style
2. Add tests for new functionality
3. Document public APIs with Doxygen comments
4. Use the exception hierarchy for error handling
5. Use the Logger for diagnostic output

## License

See [LICENSE.md](LICENSE.md) for license information.
