# Coretemperature

*Auto-generated from `src/physical_models/CoreTemperature.h`*

## Overview

MKF provides multiple implementations for coretemperature calculations. Each model is based on published scientific research and has been validated against experimental data.

## Available Models

### Maniktala

Based on "witching Power Supplies A - Z, 2nd edition" by Sanjaya Maniktala

**Validation Error:** 24.8% mean deviation

**Reference:** [https://www.goodreads.com/book/show/12042906-switching-power-supplies-a-z](https://www.goodreads.com/book/show/12042906-switching-power-supplies-a-z)

### Kazimierczuk

Based on "High-Frequency Magnetic Components 2nd Edition" by Marian Kazimierczuk

**Validation Error:** 25.8% mean deviation

**Reference:** [https://www.goodreads.com/book/show/18861402-high-frequency-magnetic-components](https://www.goodreads.com/book/show/18861402-high-frequency-magnetic-components)

### TDK

Based on "Ferrites and Accessories" by TDK

**Validation Error:** 51.5% mean deviation

**Reference:** [https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/](https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/)

### Dixon

Based on "Design of Flyback Transformers and Filter Inductors" by Lloyd H. Dixon

**Validation Error:** 24.6% mean deviation

**Reference:** [https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086](https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086)

### Amidon

Based on "Iron Powder Core Loss Characteristics" by Amidon

**Validation Error:** 25.4% mean deviation

**Reference:** [https://www.amidoncorp.com/product_images/specifications/1-38.pdf](https://www.amidoncorp.com/product_images/specifications/1-38.pdf)

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|
| Maniktala | 24.8% | [Link](https://www.goodreads.com/book/show/12042906-switching-power-supplies-a-z) |
| Kazimierczuk | 25.8% | [Link](https://www.goodreads.com/book/show/18861402-high-frequency-magnetic-components) |
| TDK | 51.5% | [Link](https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/) |
| Dixon | 24.6% | [Link](https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086) |
| Amidon | 25.4% | [Link](https://www.amidoncorp.com/product_images/specifications/1-38.pdf) |

## Usage

```cpp
#include "physical_models/CoreTemperature.h"

// Create a specific model
auto model = OpenMagnetics::CoreTemperatureModel::factory(
    OpenMagnetics::CoreTemperatureModels::MANIKTALA
);

// Or use the default model
auto model = OpenMagnetics::CoreTemperatureModel::factory();
```

## Model Selection

To set a default model globally:

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_coretemperature_model(OpenMagnetics::CoreTemperatureModels::...);
```