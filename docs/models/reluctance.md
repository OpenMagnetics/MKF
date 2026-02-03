# Reluctance

*Auto-generated from `Reluctance.h` with comprehensive documentation*

The air gap in a magnetic core is where energy is stored in an inductor. Accurately calculating
gap reluctance is essential for predicting inductance and ensuring designs meet specifications.

The magnetic field doesn't pass straight through the gap—it **fringes outward**, effectively
increasing the gap's cross-sectional area. Different models account for this fringing in
different ways, with significant impact on inductance predictions.

$$R_{gap} = \frac{l_g}{\mu_0 \cdot A_{eff}}$$

Where $A_{eff} > A_{core}$ due to fringing, and the various models differ in how they
compute this effective area.

## Available Models

### Zhang

Zhang's model (2020) treats total reluctance as parallel combination of internal and fringing:

$$R_g = R_{in} \parallel R_{fr}$$

The internal reluctance handles direct flux, while fringing reluctance accounts for
flux spreading around gap edges:

$$R_{fr} = \frac{\pi}{\mu_0 \cdot C \cdot \ln\left(\frac{2h + l_g}{l_g}\right)}$$

Where $C$ is the perimeter of core cross-section and $h$ is winding window height.

**Advantages:**
- Excellent accuracy for typical power magnetics
- Simple analytical formulation
- Default model in MKF

**Validated error:** ~3.6% mean deviation from measurements.

**Validation Error:** 11.6% mean deviation

**Reference:** [Zhang et al. 'Improved Calculation Method for Inductance Value of the Air-Gap Inductor.' CIYCEE, 2020](https://ieeexplore.ieee.org/document/9332553)

### Muehlethaler

Mühlethaler's 3D approach divides the fringing field into distinct geometric regions,
each with its own reluctance contribution:

$$R_{total} = \left(\sum_i \frac{1}{R_i}\right)^{-1}$$

Regions include:
- Direct gap path
- Corner fringing (4 corners)
- Edge fringing (4 edges)
- Top/bottom fringing

**Best for:**
- Planar cores with complex geometries
- Very large gaps (relative to core dimensions)
- Cores with multiple distributed gaps
- EI, EE, and planar E cores

**Validated error:** ~1.4% mean deviation (best accuracy).

**Validation Error:** 11.1% mean deviation

**Reference:** [Mühlethaler et al. 'A Novel Approach for 3D Air Gap Reluctance Calculations.' ECCE Asia, 2011](https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf)

### Partridge

Partridge's model uses empirical correction factors for fringing:

$$R_g = \frac{l_g}{\mu_0 \cdot A_c \cdot F}$$

Where $F$ is a fringing factor that depends on gap-to-core dimension ratios.
Simple and fast for quick estimates.

**Validation Error:** 12.4% mean deviation

**Reference:** [Partridge, 'Analysis of magnetic circuits with air gaps.' AIEE Trans., 1937](https://ieeexplore.ieee.org/document/1085)

### Effective Area

Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.

**Validation Error:** 17.5% mean deviation

**Reference:** [https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33)

### Effective Length

Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.

**Validation Error:** 17.5% mean deviation

**Reference:** [https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33)

### Stenglein

Stenglein's model provides improved accuracy for gapped ferrite cores by including
detailed geometric factors and validated against extensive FEA simulations.
Particularly accurate for E-cores and pot cores.

**Validation Error:** 14.3% mean deviation

**Reference:** [Stenglein et al. 'The Reluctance of Large Air Gaps in Ferrite Cores.' APEC, 2018](https://ieeexplore.ieee.org/document/8423500)

### Balakrishnan

Balakrishnan's model uses a permeance-based approach with geometric correction factors.
Well-suited for planar magnetics and integrated magnetics applications.

**Validation Error:** 13.7% mean deviation

**Reference:** [Balakrishnan et al. 'Air-gap reluctance and inductance calculations.' IEEE PESC, 1997](https://ieeexplore.ieee.org/document/602562)

### Classic

The classic reluctance formula assumes uniform magnetic field (no fringing):

$$R_g = \frac{l_g}{\mu_0 \cdot A_c}$$

Where:
- $R_g$ is gap reluctance (H⁻¹)
- $l_g$ is gap length (m)
- $A_c$ is core cross-sectional area (m²)
- $\mu_0$ is permeability of free space

**Limitation:** Ignores fringing entirely, leading to overestimated reluctance
and underestimated inductance. Only useful for rough estimates or very small gaps.

**Validation Error:** 28.4% mean deviation

**Reference:** [Standard magnetic circuit theory](https://en.wikipedia.org/wiki/Magnetic_reluctance)

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|
| Zhang | 11.6% | [Link](https://ieeexplore.ieee.org/document/9332553) |
| Muehlethaler | 11.1% | [Link](https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/) |
| Partridge | 12.4% | [Link](https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook) |
| Effective Area | 17.5% | [Link](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33) |
| Effective Length | 17.5% | [Link](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33) |
| Stenglein | 14.3% | [Link](https://ieeexplore.ieee.org/document/7695271/) |
| Balakrishnan | 13.7% | [Link](https://ieeexplore.ieee.org/document/602560) |
| Classic | 28.4% | [Link](https://en.wikipedia.org/wiki/Magnetic_reluctance) |

## Model Selection Guide

| Core Type | Recommended Model | Notes |
|-----------|-------------------|-------|
| E-cores (standard) | Zhang | Good balance of accuracy and simplicity |
| Planar E-cores | Mühlethaler | Handles wide, flat geometries |
| Pot cores | Zhang or Stenglein | Enclosed geometry |
| Toroids | Classic | No air gap, or distributed gap |
| RM cores | Zhang | Similar to E-cores |
| PQ cores | Zhang or Mühlethaler | Choose based on gap size |
| Very large gaps | Mühlethaler | >20% of center leg width |

**Default Recommendation:** Use **Zhang** for most applications. Switch to **Mühlethaler**
for planar cores or when maximum accuracy is required.

**Remember:** Fringing always *reduces* reluctance (increases inductance). If measured
inductance exceeds calculations, fringing is likely underestimated.

## Usage

```cpp
#include "physical_models/Reluctance.h"

// Create a specific model
auto model = OpenMagnetics::ReluctanceModel::factory(
    OpenMagnetics::ReluctanceModels::ZHANG
);

// Or use the default model
auto model = OpenMagnetics::ReluctanceModel::factory();
```

## Configuring Default Model

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::...);
```