# Reluctance

*Auto-generated from `src/physical_models/Reluctance.h`*

## Overview

MKF provides multiple implementations for reluctance calculations. Each model is based on published scientific research and has been validated against experimental data.

## Available Models

### Zhang

Based on "Improved Calculation Method for Inductance Value of the Air-Gap Inductor" by Xinsheng Zhang.

**Validation Error:** 11.6% mean deviation

**Reference:** [https://ieeexplore.ieee.org/document/9332553](https://ieeexplore.ieee.org/document/9332553)

### Muehlethaler

Based on "A Novel Approach for 3D Air Gap Reluctance Calculations" by Jonas Mühlethaler.

**Validation Error:** 11.1% mean deviation

**Reference:** [https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/](https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/)

**Related Article:** [/musings/10_gap_reluctance_and_muehlethaler_method](/musings/10_gap_reluctance_and_muehlethaler_method)

### Partridge

Based on the method described in page 8-11 from "Transformer and Inductor Design Handbook Fourth Edition" by Colonel Wm. T. McLyman.

**Validation Error:** 12.4% mean deviation

**Reference:** [https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook](https://www.goodreads.com/book/show/30187347-transformer-and-inductor-design-handbook)

### Effective Area

Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.

**Validation Error:** 17.5% mean deviation

**Reference:** [https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33)

### Effective Length

Based on the method described in page 60 from "High-Frequency Magnetic Components, Second Edition" by Marian Kazimierczuk.

**Validation Error:** 17.5% mean deviation

**Reference:** [https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33](https://www.goodreads.com/book/show/18227470-high-frequency-magnetic-components?ref=nav_sb_ss_1_33)

### Stenglein

Based on "The Reluctance of Large Air Gaps in Ferrite Cores" by Erika Stenglein.

**Validation Error:** 14.3% mean deviation

**Reference:** [https://ieeexplore.ieee.org/document/7695271/](https://ieeexplore.ieee.org/document/7695271/)

**Related Article:** [/musings/11_inductance_variables_and_stenglein_method](/musings/11_inductance_variables_and_stenglein_method)

### Balakrishnan

Based on "Air-gap reluctance and inductance calculations for magnetic circuits using a Schwarz-Christoffel transformation" by A. Balakrishnan.

**Validation Error:** 13.7% mean deviation

**Reference:** [https://ieeexplore.ieee.org/document/602560](https://ieeexplore.ieee.org/document/602560)

### Classic

Based on the reluctance of a uniform magnetic circuit

**Validation Error:** 28.4% mean deviation

**Reference:** [https://en.wikipedia.org/wiki/Magnetic_reluctance](https://en.wikipedia.org/wiki/Magnetic_reluctance)

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

## Model Selection

To set a default model globally:

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::...);
```