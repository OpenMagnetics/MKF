# MKF - Magnetic Kernel Framework

## Project Overview
C++ library for magnetic component calculations, converter modeling, and circuit simulation.

## Architecture
- **Language:** C++17
- **Build System:** CMake
- **Test Framework:** Catch2
- **JSON:** nlohmann/json
- **Optional:** ngspice for circuit simulation

## Key Components

### Converter Models (`src/converter_models/`)
- `PushPull.cpp` - Push-pull converter with center-tapped transformer
- `Flyback.cpp` - Flyback converter
- `Forward.cpp` - Forward converters (single-switch, two-switch, active clamp)
- `BuckBoost.cpp` - Buck/Boost converters
- `LLC.cpp`, `DAB.cpp` - Resonant converters

### Core Classes
- `OperatingPoint` - Operating point calculations
- `Magnetic` - Core, winding, coil data
- `Inputs` - Design inputs and requirements
- `NgspiceRunner` - SPICE simulation interface

### Utilities
- `Painter` - Waveform visualization (outputs SVG)
- `Utils` - Helper functions
- `TestingUtils` - Test utilities

## Directory Structure
```
src/
├── converter_models/     # Converter topology implementations
├── models/              # Magnetic component models
├── processors/          # ngspice runner
├── support/             # Utils, Painter, TestingUtils
└── tests/               # Unit tests
```

## Building

### Standard Build
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Run Tests
```bash
cd build
ctest --output-on-failure
```

### Generate Test Waveforms
Tests output SVG files to `tests/output/` for visual inspection.

## Converter Model Implementation

### Adding a New Converter
1. Create class in `src/converter_models/`
2. Inherit from base converter class
3. Implement required methods:
   - `process()` - Calculate operating points
   - `generate_ngspice_circuit()` - Create SPICE netlist
   - `simulate_and_extract_operating_points()` - Run simulation

### ngspice Circuit Design Tips
- **Coupling:** Use pairwise K statements for coupled inductors
- **Snubbers:** Add RC/RCD snubbers across switches for stability
- **Damping:** Include series resistance (10-100mΩ) with inductors
- **Tolerances:** RELTOL=0.001 for cleaner waveforms
- **Ringing:** If oscillations occur, add snubbers or reduce coupling coefficient

### Common Issues
- **Flux imbalance:** Ensure all windings on same core are mutually coupled
- **Non-periodic results:** Check all windings are properly coupled
- **Convergence errors:** Relax tolerances or add snubber resistors

## Testing

### Create New Test
```cpp
TEST_CASE("Test_MyConverter", "[converter-model][my-converter]") {
    json inputsJson;
    // ... setup inputs
    
    MyConverter converter(inputsJson);
    auto inputs = converter.process();
    
    // Validate results
    REQUIRE(inputs.get_operating_points().size() > 0);
}
```

### Visual Debugging
Use Painter to save waveforms:
```cpp
Painter painter("output.svg", false, true);
painter.paint_waveform(waveform);
painter.export_svg();
```

## Integration with WebLibMKF
This library is compiled to WebAssembly by WebLibMKF project. Changes here require rebuilding WebLibMKF to update the web frontend.

## Dependencies
- MAS (Magnetic Application Standard)
- nlohmann/json
- magic_enum
- Catch2 (testing)
- ngspice (optional, for simulation)