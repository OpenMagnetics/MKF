# Corelosses

*Auto-generated from `src/physical_models/CoreLosses.h`*

## Overview

MKF provides multiple implementations for corelosses calculations. Each model is based on published scientific research and has been validated against experimental data.

## Available Models

### Steinmetz

Based on "On the law of hysteresis" by Charles Proteus Steinmetz

**Validation Error:** 39.4% mean deviation

**Reference:** [https://ieeexplore.ieee.org/document/1457110](https://ieeexplore.ieee.org/document/1457110)

### iGSE

Based on "Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters" by Charles R. Sullivan

**Validation Error:** 35.8% mean deviation

**Reference:** [http://inductor.thayerschool.org/papers/IGSE.pdf](http://inductor.thayerschool.org/papers/IGSE.pdf)

### Barg

Based on "Core Loss Calculation of Symmetric Trapezoidal Magnetic Flux Density Waveform" by Sobhi Barg

**Validation Error:** 37.4% mean deviation

**Reference:** [https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf](https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf)

### Roshen

Based on "Ferrite Core Loss for Power Magnetic Components Design" and "A Practical, Accurate and Very General Core Loss Model for Nonsinusoidal Waveforms" by Waseem Roshen

**Validation Error:** 48.8% mean deviation

**Reference:** [https://ieeexplore.ieee.org/document/4052433](https://ieeexplore.ieee.org/document/4052433)

**Related Article:** [/musings/4_roshen_method_core_losses](/musings/4_roshen_method_core_losses)

### Albach

Based on "Calculating Core Losses in Transformers for Arbitrary Magnetizing Currents A Comparison of Different Approaches" by Manfred Albach

**Validation Error:** 35.7% mean deviation

**Reference:** [https://ieeexplore.ieee.org/iel3/3925/11364/00548774.pdf](https://ieeexplore.ieee.org/iel3/3925/11364/00548774.pdf)

### NSE

Based on "Measurement and Loss Model of Ferrites with Non-sinusoidal Waveforms" by Alex Van den Bossche

**Validation Error:** 35.8% mean deviation

**Reference:** [http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf](http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf)

### MSE

Based on "Calculation of Losses in Ferro- and Ferrimagnetic Materials Based on the Modified Steinmetz Equation" by Jürgen Reinert

**Validation Error:** 35.7% mean deviation

**Reference:** [https://www.mikrocontroller.net/attachment/129490/Modified_Steinmetz.pdf](https://www.mikrocontroller.net/attachment/129490/Modified_Steinmetz.pdf)

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|
| Steinmetz | 39.4% | [Link](https://ieeexplore.ieee.org/document/1457110) |
| iGSE | 35.8% | [Link](http://inductor.thayerschool.org/papers/IGSE.pdf) |
| Barg | 37.4% | [Link](https://miun.diva-portal.org/smash/get/diva2:1622559/FULLTEXT01.pdf) |
| Roshen | 48.8% | [Link](https://ieeexplore.ieee.org/document/4052433) |
| Albach | 35.7% | [Link](https://ieeexplore.ieee.org/iel3/3925/11364/00548774.pdf) |
| NSE | 35.8% | [Link](http://web.eecs.utk.edu/~dcostine/ECE482/Spring2015/materials/magnetics/NSE.pdf) |
| MSE | 35.7% | [Link](https://www.mikrocontroller.net/attachment/129490/Modified_Steinmetz.pdf) |

## Usage

```cpp
#include "physical_models/CoreLosses.h"

// Create a specific model
auto model = OpenMagnetics::CoreLossesModel::factory(
    OpenMagnetics::CoreLossesModels::STEINMETZ
);

// Or use the default model
auto model = OpenMagnetics::CoreLossesModel::factory();
```

## Model Selection

To set a default model globally:

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_corelosses_model(OpenMagnetics::CoreLossesModels::...);
```