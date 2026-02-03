# Quick Start

This guide walks you through your first magnetic component simulation with MKF.

## Basic Concepts

MKF uses a layered architecture:

```
┌─────────────────────────────────────────────┐
│           Advisers (Design Optimization)     │
│   CoreAdviser, CoilAdviser, WireAdviser     │
├─────────────────────────────────────────────┤
│         Physical Models (Calculations)       │
│  CoreLosses, WindingLosses, Inductance...   │
├─────────────────────────────────────────────┤
│      Constructive Models (Components)        │
│       Magnetic, Core, Coil, Wire            │
├─────────────────────────────────────────────┤
│            Data Layer (MAS Schema)           │
│         JSON schemas, DatabaseManager        │
└─────────────────────────────────────────────┘
```

## Example 1: Calculate Core Reluctance

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Load a core from the database
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Create a reluctance model (Zhang method)
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        OpenMagnetics::ReluctanceModels::ZHANG
    );

    // Calculate reluctance
    auto result = reluctanceModel->get_core_reluctance(core);

    std::cout << "Core: E 42/21/15" << std::endl;
    std::cout << "Total reluctance: " << result.get_core_reluctance() << " H^-1" << std::endl;
    std::cout << "Ungapped reluctance: " << result.get_ungapped_core_reluctance() << " H^-1" << std::endl;

    return 0;
}
```

## Example 2: Calculate Core Losses

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Load core with material
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Create operating point excitation
    MAS::OperatingPointExcitation excitation;

    // Set up magnetic flux density waveform
    MAS::SignalDescriptor magneticFluxDensity;
    MAS::Processed processed;
    processed.set_peak(0.1);  // 100 mT peak
    processed.set_offset(0.0);
    magneticFluxDensity.set_processed(processed);
    excitation.set_magnetic_flux_density(magneticFluxDensity);

    // Set frequency
    excitation.set_frequency(100000);  // 100 kHz

    // Calculate losses using iGSE model
    OpenMagnetics::CoreLosses coreLosses;
    auto losses = coreLosses.calculate_core_losses(core, excitation, 25.0);

    std::cout << "Core losses at 100kHz, 100mT:" << std::endl;
    std::cout << "  Total: " << losses.get_core_losses() << " W" << std::endl;
    std::cout << "  Volumetric: " << losses.get_volumetric_losses().value() << " W/m³" << std::endl;

    return 0;
}
```

## Example 3: Use the Database Manager

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Get the database manager singleton
    auto& db = OpenMagnetics::DatabaseManager::getInstance();

    // List available core shapes
    auto coreShapes = db.get_core_shape_names();
    std::cout << "Available core shapes: " << coreShapes.size() << std::endl;

    // Find cores by shape family
    auto eCores = db.get_cores_by_shape_family("E");
    std::cout << "E-cores available: " << eCores.size() << std::endl;

    // Find materials by type
    auto ferrites = db.get_materials_by_type("ferrite");
    std::cout << "Ferrite materials: " << ferrites.size() << std::endl;

    // Load a specific wire
    auto wire = OpenMagnetics::find_wire_by_name("AWG 20");
    std::cout << "AWG 20 diameter: " << wire.get_conducting_diameter().value() << " m" << std::endl;

    return 0;
}
```

## Example 4: Configure Settings

```cpp
#include "OpenMagnetics.h"

int main() {
    // Access the global settings singleton
    auto& settings = OpenMagnetics::Settings::GetInstance();

    // Configure core adviser behavior
    settings.set_core_adviser_include_stacks(true);
    settings.set_core_adviser_include_distributed_gaps(true);

    // Configure wire adviser
    settings.set_wire_adviser_include_litz(true);
    settings.set_wire_adviser_include_foil(false);

    // Set default reluctance model
    settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::ZHANG);

    // Enable verbose output for debugging
    settings.set_verbose(true);

    return 0;
}
```

## Python Quick Start

```python
import PyMKF

# Load a core
core = PyMKF.find_core_by_name("E 42/21/15")

# Calculate reluctance
result = PyMKF.calculate_reluctance(core, model="ZHANG")
print(f"Reluctance: {result['core_reluctance']} H^-1")

# Calculate core losses
losses = PyMKF.calculate_core_losses(
    core=core,
    frequency=100000,  # Hz
    peak_flux_density=0.1,  # T
    temperature=25  # C
)
print(f"Core losses: {losses['total']} W")
```

## Next Steps

- Explore [Physical Models](../models/index.md) for detailed model documentation
- Learn about [Design Advisers](../algorithms/core-adviser.md) for optimization
- Check the [Examples](examples.md) page for more use cases
- Review [Settings Reference](../settings/reference.md) for configuration options
