# Physical Models Overview

MKF implements multiple scientific models for electromagnetic calculations. Each phenomenon has several competing implementations, allowing you to choose the most appropriate model for your application.

## Model Selection Philosophy

MKF follows a **multi-model approach**:

1. **Multiple implementations** for each physical phenomenon
2. **Documented accuracy** through error metrics from validation tests
3. **Scientific references** linking to original papers
4. **Default recommendations** based on comprehensive testing

## Available Physical Models

### Core-Related Models

| Category | Models | Default |
|----------|--------|---------|
| [Core Losses](core-losses.md) | Steinmetz, iGSE, MSE, NSE, Roshen, Albach, Barg, Proprietary | iGSE |
| [Reluctance](reluctance.md) | Zhang, Muehlethaler, Partridge, Stenglein, Balakrishnan, Classic | Zhang |
| [Thermal](thermal.md) | Kazimierczuk, Maniktala, TDK, Dixon, Amidon | Maniktala |

### Winding-Related Models

| Category | Models | Default |
|----------|--------|---------|
| [Skin Effect](winding-losses.md#skin-effect) | Dowell, Wojda, Albach, Payne, Kutkut, Ferreira | Wire-dependent |
| [Proximity Effect](winding-losses.md#proximity-effect) | Rossmanith, Wang, Ferreira, Lammeraner, Albach, Dowell | Wire-dependent |

### Field Models

| Category | Models | Default |
|----------|--------|---------|
| [Magnetic Field](magnetic-field.md) | Binns-Lawrenson, Lammeraner, Dowell, Wang, Albach | Binns-Lawrenson |
| [Fringing Field](magnetic-field.md#fringing-effect) | Roshen, Albach | Roshen |

### Other Models

| Category | Models | Default |
|----------|--------|---------|
| [Capacitance](capacitance.md) | Koch, Albach, Duerdoth, Massarini | Albach |
| [Inductance](inductance.md) | Based on reluctance model | - |

## Using Models Programmatically

### Factory Pattern

All models use a factory pattern for instantiation:

```cpp
// Create a specific model
auto model = OpenMagnetics::ReluctanceModel::factory(
    OpenMagnetics::ReluctanceModels::ZHANG
);

// Or use default settings
auto model = OpenMagnetics::ReluctanceModel::factory();
```

### Querying Model Information

Each model class provides static methods for introspection:

```cpp
// Get descriptions of all reluctance models
auto info = OpenMagnetics::ReluctanceModel::get_models_information();
for (const auto& [name, description] : info) {
    std::cout << name << ": " << description << std::endl;
}

// Get error metrics
auto errors = OpenMagnetics::ReluctanceModel::get_models_errors();

// Get external references (papers, books)
auto links = OpenMagnetics::ReluctanceModel::get_models_external_links();
```

### Configuring Default Models

Use the Settings singleton to configure defaults:

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();

// Set default reluctance model
settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::MUEHLETHALER);

// Set default skin effect model
settings.set_winding_skin_effect_losses_model(
    OpenMagnetics::WindingSkinEffectLossesModels::ALBACH
);
```

## Model Accuracy

Each model page includes:

- **Mean error**: Average deviation from measured values
- **Standard deviation**: Variability of predictions
- **Validation dataset**: Description of test conditions

Error metrics are computed against published experimental data and internal test suites.

## Adding Custom Models

MKF's model architecture supports extension. To add a custom model:

1. Inherit from the appropriate base class (e.g., `ReluctanceModel`)
2. Implement the required virtual methods
3. Register in the factory function

See the source code for examples of model implementations.
