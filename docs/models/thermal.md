# Coretemperature

*Auto-generated from `CoreTemperature.h` with comprehensive documentation*

Thermal management is critical for magnetic component reliability. MKF provides models for:
- Core temperature rise from core losses
- Winding hot-spot temperature
- Thermal resistance estimation

Temperature affects:
- Core losses (Steinmetz parameters are temperature-dependent)
- Winding resistance (copper has +0.39%/°C temperature coefficient)
- Saturation flux density (decreases with temperature)
- Component lifetime (every 10°C reduces life by ~50%)

## Available Models

### Maniktala

Maniktala's empirical formula for temperature rise:

$$\Delta T = \left(\frac{P_{total}}{A_{surface}}\right)^{0.833}$$

Where:
- $P_{total}$ is total power loss (W)
- $A_{surface}$ is effective cooling surface area (cm²)
- Result is temperature rise in °C

**Simple and practical** for initial estimates and design iteration.

**Reference:** Maniktala, S. "Switching Power Supplies A-Z." Newnes, 2012.

**Validation Error:** 24.8% mean deviation

**Reference:** [Maniktala, S. 'Switching Power Supplies A-Z.' Newnes, 2012](https://www.elsevier.com/books/switching-power-supplies-a-z/maniktala/978-0-12-386533-5)

### Kazimierczuk

Kazimierczuk's thermal model separates core and winding contributions:

$$\Delta T_{core} = R_{th,core} \cdot P_{core}$$
$$\Delta T_{winding} = R_{th,winding} \cdot P_{winding}$$

Thermal resistances are calculated from geometry and material properties.
More accurate than empirical formulas but requires more parameters.

**Reference:** Kazimierczuk, M.K. "High-Frequency Magnetic Components." Wiley, 2014.

**Validation Error:** 25.8% mean deviation

**Reference:** [Kazimierczuk, M.K. 'High-Frequency Magnetic Components.' Wiley, 2014](https://www.wiley.com/en-us/High+Frequency+Magnetic+Components%2C+2nd+Edition-p-9781118717790)

### TDK

TDK's thermal models are empirically derived from extensive testing of their
ferrite cores. Provides manufacturer-specific accuracy for TDK materials.

**Validation Error:** 51.5% mean deviation

**Reference:** [TDK Ferrites and Accessories Application Notes](https://www.tdk-electronics.tdk.com/download/531372/4730039e8720de4e83e3e32614ca3efc/ferrites-magnetic-design-tool.pdf)

### Dixon

Dixon's thermal resistance approximation for pot cores and similar shapes:

$$R_{th} \approx \frac{50}{A_{surface}^{0.7}}$$ (°C/W)

Where $A_{surface}$ is in cm².

**Validation Error:** 24.6% mean deviation

**Reference:** [Dixon, L. 'Magnetics Design for Switching Power Supplies.' Texas Instruments](https://www.ti.com/lit/ml/slup123/slup123.pdf)

### Amidon

Amidon's thermal data is empirically derived for their iron powder and ferrite cores.
Useful when designing with Amidon (Micrometals) materials.

**Validation Error:** 25.4% mean deviation

**Reference:** [Amidon Corp. Technical Specifications](https://www.amidoncorp.com/specs/)

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|
| Maniktala | 24.8% | [Link](https://www.goodreads.com/book/show/12042906-switching-power-supplies-a-z) |
| Kazimierczuk | 25.8% | [Link](https://www.goodreads.com/book/show/18861402-high-frequency-magnetic-components) |
| TDK | 51.5% | [Link](https://www.tdk-electronics.tdk.com/download/531536/badc7640e8213233c951b4540e3745e2/) |
| Dixon | 24.6% | [Link](https://www.ti.com/lit/ml/slup076/slup076.pdf?ts=1679429443086) |
| Amidon | 25.4% | [Link](https://www.amidoncorp.com/product_images/specifications/1-38.pdf) |

## Thermal Model Selection

| Application | Recommended | Notes |
|-------------|-------------|-------|
| Quick estimates | Maniktala | Simple, conservative |
| Detailed analysis | Kazimierczuk | Separates loss sources |
| TDK cores | TDK | Manufacturer data |
| Worst-case | Use lowest thermal resistance estimate | Safety margin |

**Thermal Resistance Models:**
The thermal resistance between core/winding and ambient depends on:
- Cooling method (natural convection, forced air, heat sink)
- Surface treatment (painted, bare, potted)
- Mounting orientation (vertical, horizontal)
- Ambient temperature and altitude

**Design Margin:** Add 10-15°C margin for component reliability.

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

## Configuring Default Model

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_coretemperature_model(OpenMagnetics::CoreTemperatureModels::...);
```