# API Reference

This section provides detailed C++ API documentation for MKF.

## Doxygen Documentation

The full API documentation is generated using Doxygen and is available in the `docs/api/html/` directory after building.

To generate the documentation locally:

```bash
# From the repository root
doxygen Doxyfile
```

Then open `docs/api/html/index.html` in your browser.

## Key Classes

### Core Classes

| Class | Description |
|-------|-------------|
| `OpenMagnetics::Core` | Represents a magnetic core with shape, material, and gap information |
| `OpenMagnetics::Coil` | Represents winding configuration with turns, layers, and wire info |
| `OpenMagnetics::Wire` | Represents wire specifications (round, litz, rectangular, foil) |
| `OpenMagnetics::Magnetic` | Complete magnetic assembly combining core and coil |
| `OpenMagnetics::Bobbin` | Bobbin/former specifications for wound components |

### Physical Model Classes

| Class | Description |
|-------|-------------|
| `OpenMagnetics::CoreLossesModel` | Base class for core loss calculations |
| `OpenMagnetics::ReluctanceModel` | Base class for reluctance calculations |
| `OpenMagnetics::WindingLosses` | Winding loss calculations (DC, skin, proximity) |
| `OpenMagnetics::MagneticField` | Magnetic field strength calculations |
| `OpenMagnetics::Inductance` | Inductance calculations |

### Adviser Classes

| Class | Description |
|-------|-------------|
| `OpenMagnetics::CoreAdviser` | Core selection optimization |
| `OpenMagnetics::CoilAdviser` | Winding configuration optimization |
| `OpenMagnetics::WireAdviser` | Wire selection optimization |
| `OpenMagnetics::MagneticAdviser` | Complete magnetic design optimization |

### Support Classes

| Class | Description |
|-------|-------------|
| `OpenMagnetics::DatabaseManager` | Singleton for database access |
| `OpenMagnetics::Settings` | Global configuration singleton |
| `OpenMagnetics::Logger` | Logging facility |

## Namespaces

### OpenMagnetics

Main namespace containing all MKF classes and functions.

```cpp
namespace OpenMagnetics {
    // Core functionality
    class Core;
    class Coil;
    class Wire;
    class Magnetic;

    // Physical models
    class CoreLossesModel;
    class ReluctanceModel;
    // ...

    // Support functions
    Core find_core_by_name(const std::string& name);
    Wire find_wire_by_name(const std::string& name);
    // ...
}
```

### MAS

Auto-generated namespace from JSON schemas containing data structures.

```cpp
namespace MAS {
    // Data structures
    struct OperatingPoint;
    struct OperatingPointExcitation;
    struct DesignRequirements;
    struct SignalDescriptor;
    // ...

    // Enumerations
    enum class WireType;
    enum class GapType;
    enum class CoreShapeFamily;
    // ...
}
```

## Common Patterns

### Factory Pattern

Models use factories for instantiation:

```cpp
// Create specific model
auto model = OpenMagnetics::ReluctanceModel::factory(
    OpenMagnetics::ReluctanceModels::ZHANG
);

// Create default model
auto model = OpenMagnetics::ReluctanceModel::factory();
```

### Singleton Access

```cpp
// Database access
auto& db = OpenMagnetics::DatabaseManager::getInstance();

// Settings access
auto& settings = OpenMagnetics::Settings::GetInstance();
```

### JSON Serialization

All MAS structures support JSON serialization:

```cpp
// Serialize to JSON
json j;
to_json(j, myObject);
std::string jsonString = j.dump();

// Deserialize from JSON
from_json(j, myObject);
```

## Error Handling

MKF uses typed exceptions defined in `Exceptions.h`:

```cpp
try {
    auto core = OpenMagnetics::find_core_by_name("NonExistent");
} catch (const OpenMagnetics::CoreNotFoundException& e) {
    std::cerr << "Core not found: " << e.what() << std::endl;
} catch (const OpenMagnetics::DatabaseException& e) {
    std::cerr << "Database error: " << e.what() << std::endl;
}
```

Common exception types:
- `CoreNotFoundException`
- `MaterialException`
- `ModelNotAvailableException`
- `CalculationException`
- `GapException`
- `CoilNotProcessedException`

## Thread Safety

- `DatabaseManager` is thread-safe for read operations
- `Settings` modifications should be done before multi-threaded calculations
- Model instances are not thread-safe; create separate instances per thread
