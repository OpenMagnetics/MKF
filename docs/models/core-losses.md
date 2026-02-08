# Corelosses

*Auto-generated from `CoreLosses.h` with comprehensive documentation*

Core losses in magnetic materials arise from two primary mechanisms: **hysteresis losses** (energy dissipated
in reorienting magnetic domains) and **eddy current losses** (currents induced in the conductive material).
Accurate prediction of these losses is critical for thermal management and efficiency optimization in
transformers, inductors, and other magnetic components.

MKF implements multiple core loss models, each with different trade-offs between accuracy and complexity.
The choice of model depends on the waveform type, available material data, and required precision.

## Available Models

### Steinmetz

The original Steinmetz equation (1892) empirically relates core losses to frequency and flux density:

$$P_v = k \cdot f^\alpha \cdot \hat{B}^\beta$$

Where:
- $P_v$ is volumetric power loss (W/m³)
- $f$ is frequency (Hz)
- $\hat{B}$ is peak flux density (T)
- $k$, $\alpha$, $\beta$ are material-specific Steinmetz coefficients

**Limitations:** Only valid for sinusoidal excitation. Using this model with non-sinusoidal
waveforms (triangular, trapezoidal, PWM) will produce significant errors.

**Validation Error:** 39.4% mean deviation

**Reference:** [Steinmetz, C.P. 'On the law of hysteresis.' AIEE Trans., 1892](https://ieeexplore.ieee.org/document/4768063)

### iGSE

The **improved Generalized Steinmetz Equation** (iGSE) extends Steinmetz to arbitrary waveforms:

$$P_v = \frac{1}{T} \int_0^T k_i \left| \frac{dB}{dt} \right|^\alpha (\Delta B)^{\beta - \alpha} dt$$

Where $k_i = \frac{k}{(2\pi)^{\alpha-1} \int_0^{2\pi} |\cos\theta|^\alpha 2^{\beta-\alpha} d\theta}$

**Key insight:** Core losses depend on the rate of change of flux density (dB/dt), not just
peak values. This allows accurate predictions for:
- Triangular waveforms (inductors in CCM)
- Trapezoidal waveforms (half-bridge converters)
- Arbitrary PWM waveforms

**Advantage:** Uses the same Steinmetz parameters (k, α, β) as the original equation,
requiring no additional material characterization.

**Validation Error:** 35.8% mean deviation

**Reference:** [Li, Abdallah, Sullivan. 'Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms.' IEEE COMPEL, 2001](http://inductor.thayerschool.org/papers/IGSE.pdf)

### Barg

Barg's model extends loss separation for modern ferrite materials with improved
temperature and frequency dependence modeling. It separates static and dynamic
hysteresis contributions for better accuracy at high frequencies.

**Validation Error:** 37.4% mean deviation

**Reference:** [Barg et al. 'A Dynamic Hysteresis Loss Model.' IEEE Trans. Power Electronics, 2021](https://ieeexplore.ieee.org/document/9530729)

### Roshen

Roshen's physics-based model separates hysteresis and eddy current contributions:

$$P_v = P_{hyst} + P_{eddy} + P_{excess}$$

Each component is derived from material physics:
- **Hysteresis:** From B-H loop integration
- **Eddy currents:** From material resistivity and geometry
- **Excess losses:** Anomalous losses from domain wall motion

**Advantages:**
- Better accuracy at extreme operating conditions
- Handles DC bias effects
- Temperature dependence modeled physically

**Disadvantage:** Requires more material parameters than Steinmetz-based methods.

**Validation Error:** 48.8% mean deviation

**Reference:** [Roshen, W. 'Ferrite Core Loss for Power Magnetic Components Design.' IEEE Trans. Magnetics, 1991](https://ieeexplore.ieee.org/document/4052433)

### Albach

Albach's model combines loss separation with improved high-frequency accuracy:

$$P_v = C_h f B^\gamma + C_e (fB)^2 / \rho + C_a (fB)^{1.5}$$

Where:
- $C_h$, $\gamma$ are hysteresis parameters
- $C_e$ is the classical eddy current coefficient
- $C_a$ accounts for anomalous (excess) losses
- $\rho$ is material resistivity

This provides good accuracy across a wide frequency range.

**Validation Error:** 35.7% mean deviation

**Reference:** [Albach et al. 'Calculating core losses in transformers.' IEEE Trans. Magnetics, 2011](https://ieeexplore.ieee.org/document/5680877)

### NSE

The **Natural Steinmetz Extension** (NSE) uses direct integration without equivalent frequency:

$$P_v = \left(\frac{\Delta B}{2}\right)^{\beta-\alpha} \frac{k}{T} \int_0^T \left|\frac{dB}{dt}\right|^\alpha dt$$

This provides a more physically meaningful interpretation by relating losses directly
to the actual flux rate of change over the entire period.

**Validation Error:** 35.8% mean deviation

**Reference:** [Van den Bossche et al. 'Evaluation of the Steinmetz Parameters.' IEEE Trans. Magnetics, 2008](https://ieeexplore.ieee.org/document/4544745)

### MSE

The **Modified Steinmetz Equation** (MSE) uses an equivalent frequency concept:

$$P_v = k \cdot f_{eq}^{\alpha-1} \cdot f_r \cdot \hat{B}^\beta$$

Where $f_{eq} = \frac{2}{\Delta B^2 \pi^2} \int_0^T \left(\frac{dB}{dt}\right)^2 dt$

The equivalent frequency captures the effect of non-sinusoidal waveforms by computing
what sinusoidal frequency would produce the same losses.

**Validation Error:** 35.7% mean deviation

**Reference:** [Venkatachalam et al. 'Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms.' IEEE COMPEL, 2002](https://www.psma.com/sites/default/files/uploads/files/Limits%20of%20Accuracy%20of%20Steinmetz%20Equation-Handout.pdf)

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

## Model Selection Guide

| Waveform Type | Recommended Model | Notes |
|---------------|-------------------|-------|
| Sinusoidal | Steinmetz | Original equation, most manufacturer data available |
| Triangular | iGSE | Best for inductor applications (CCM) |
| Trapezoidal | iGSE or MSE | Good for half-bridge topologies |
| Arbitrary PWM | iGSE | Most general non-sinusoidal model |
| DC bias present | Roshen | Physics-based handles bias effects |
| Wide temp range | Roshen or Albach | Better temperature extrapolation |

**Default Recommendation:** Use **iGSE** for most power electronics applications.
It handles non-sinusoidal waveforms well and uses readily available Steinmetz parameters.

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

## Configuring Default Model

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
// settings.set_corelosses_model(OpenMagnetics::CoreLossesModels::...);
```