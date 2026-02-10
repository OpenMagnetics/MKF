# OpenMagnetics Documentation

Welcome to the documentation for **MKF (Magnetics Knowledge Foundation)** - the C++ simulation engine powering OpenMagnetics.

## What is OpenMagnetics?

OpenMagnetics is an open-source framework for modeling magnetic components such as inductors and transformers. MKF provides:

- **Physical Models**: Accurate electromagnetic calculations for core losses, winding losses, inductance, magnetic field distribution, and more
- **Constructive Models**: Physical component representations including cores, coils, wires, and bobbins
- **Design Advisers**: Optimization algorithms to help you find the best core, wire, and coil configurations for your requirements
- **Multi-language Support**: Core engine in C++ with bindings for Python (PyMKF) and .NET (MKFNet)

## Key Features

### Multiple Competing Models

MKF implements multiple scientific models for each electromagnetic phenomenon, allowing you to choose the most appropriate model for your application:

| Phenomenon | Available Models |
|------------|-----------------|
| Air Gap Reluctance | Zhang, Muehlethaler, Partridge, Stenglein, Balakrishnan, Classic |
| Core Losses | Steinmetz, iGSE, MSE, NSE, Roshen, Albach, Barg, Proprietary |
| Winding Skin Effect | Dowell, Albach, Kutkut, Ferreira, and more |
| Proximity Effect | Ferreira, Lammeraner, Albach, Dowell |

### Design Optimization

The adviser system helps you explore the design space efficiently:

- **CoreAdviser**: Find optimal core shapes and materials for your power and frequency requirements
- **CoilAdviser**: Optimize winding configurations for minimum losses
- **WireAdviser**: Select the best wire type and dimensions
- **MagneticAdviser**: Full magnetic component optimization

### Extensive Database

MKF includes embedded databases of:

- Commercial magnetic cores from major manufacturers
- Core materials (ferrites, powder cores)
- Wire specifications (round, litz, rectangular, foil, planar)
- Bobbin designs

## Quick Example

```cpp
#include "OpenMagnetics.h"

// Create a core from the database
auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

// Calculate reluctance
auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
    OpenMagnetics::ReluctanceModels::ZHANG
);
auto result = reluctanceModel->get_core_reluctance(core);

std::cout << "Core reluctance: " << result.get_core_reluctance() << " H^-1" << std::endl;
```

## Getting Started

1. **[Installation](getting-started/installation.md)** - Build MKF from source
2. **[Quick Start](getting-started/quick-start.md)** - Your first magnetic simulation
3. **[Examples](getting-started/examples.md)** - Common use cases and patterns

## Documentation Sections

- **[Physical Models](models/index.md)** - Detailed documentation of electromagnetic models with scientific references
- **[Settings Reference](settings/reference.md)** - Configuration options for customizing behavior
- **[Algorithms](algorithms/core-adviser.md)** - How the design advisers work
- **[API Reference](api/index.md)** - C++ API documentation (Doxygen)

## Contributing

MKF is open source and welcomes contributions. Visit the [GitHub repository](https://github.com/OpenMagnetics/MKF) to:

- Report issues
- Submit pull requests
- Request features

## License

MKF is released under the MIT License.
