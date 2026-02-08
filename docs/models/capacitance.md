# Stray Capacitance

*Comprehensive documentation for parasitic capacitance models*

Parasitic (stray) capacitance in magnetic components affects:
- High-frequency behavior and resonances
- EMI performance (common-mode noise paths)
- Switching transients in power converters
- Self-resonant frequency of inductors

MKF models several capacitance components:
- **Turn-to-turn capacitance:** Between adjacent turns in a layer
- **Layer-to-layer capacitance:** Between winding layers
- **Winding-to-core capacitance:** Between windings and grounded core
- **Primary-to-secondary capacitance:** Coupling between windings (transformers)

## Capacitance Components

### Turn-to-Turn Capacitance

Capacitance between adjacent turns in the same layer:

$$C_{tt} = \varepsilon_0 \varepsilon_r \frac{l_{turn}}{d_{tt}} \cdot k_{fringe}$$

Where:
- $l_{turn}$ is the mean turn length
- $d_{tt}$ is the turn-to-turn spacing
- $k_{fringe}$ accounts for fringing fields

### Layer-to-Layer Capacitance

For windings with multiple layers:

$$C_{ll} = \varepsilon_0 \varepsilon_r \frac{A_{overlap}}{d_{ll}}$$

Where:
- $A_{overlap}$ is the overlapping area between layers
- $d_{ll}$ is the layer-to-layer spacing (insulation thickness)

### Winding-to-Core Capacitance

$$C_{wc} = \varepsilon_0 \varepsilon_r \frac{A_{winding}}{d_{wc}}$$

This capacitance provides a path for common-mode noise in isolated converters.

## Available Models

### Koch

Koch's model uses parallel-plate approximation with fringing corrections:

$$C = \varepsilon_0 \varepsilon_r \frac{A}{d} \cdot k_{fringe}$$

Where $k_{fringe}$ accounts for field spreading at conductor edges.
Good for turn-to-turn and layer-to-layer capacitance.

**Reference:** [Koch et al. 'Stray Capacitance in Inductors.' IEEE Trans. Power Electronics, 2020](https://ieeexplore.ieee.org/document/9206060)

### Massarini

Massarini's analytical model for winding capacitance includes:
- Multi-layer effects
- Wire geometry (round, rectangular)
- Insulation thickness

Provides separate formulas for different capacitance components.

**Reference:** [Massarini, Kazimierczuk. 'Self-capacitance of inductors.' IEEE Trans. Power Electronics, 1997](https://ieeexplore.ieee.org/document/668030)

### Albach

Albach's comprehensive capacitance model accounts for:
- 2D field distribution
- Non-uniform winding spacing
- Core proximity effects

Provides good accuracy across various winding configurations.

**Reference:** [Albach et al. 'Calculating stray capacitance.' IEEE Trans. Magnetics, 2011](https://ieeexplore.ieee.org/document/5680877)

### Duerdoth

Duerdoth's coil self-capacitance model uses energy-based approach:

$$C_{self} = \frac{2 W_E}{V^2}$$

Where $W_E$ is the electric field energy stored between turns.

**Reference:** [Duerdoth, W.T. 'Equivalent Capacitance of Transformer Windings.' Wireless Engineer, 1946](https://ieeexplore.ieee.org/document/7866793)

## Self-Resonant Frequency

The self-resonant frequency (SRF) of an inductor is where capacitive and inductive
reactances cancel:

$$f_{SRF} = \frac{1}{2\pi\sqrt{L \cdot C_{stray}}}$$

Above the SRF, the inductor behaves as a capacitor. For effective filtering,
operate well below the SRF (typically f < SRF/10).

## Model Comparison

| Model | Best For | Complexity |
|-------|----------|------------|
| Koch | Turn-to-turn, layer-to-layer | Medium |
| Massarini | Complete self-capacitance | Medium |
| Albach | General purpose | High |
| Duerdoth | Quick estimates | Low |

**Default Model:** Albach provides good balance of accuracy and generality.

## Capacitance Minimization Strategies

| Technique | Effect | Trade-off |
|-----------|--------|-----------|
| Increase turn spacing | Reduces turn-to-turn C | Larger winding window |
| Section winding | Reduces layer-to-layer C | More complex construction |
| Shield windings | Controls coupling paths | Additional losses |
| Bank winding | Reduces distributed C | Higher leakage inductance |
| Thicker insulation | Reduces all C | Thermal resistance |
| Progressive winding | Reduces layer-to-layer C | Complex winding pattern |

## EMI Considerations

### Common-Mode Capacitance

In isolated converters, capacitance between primary and secondary windings ($C_{ps}$)
and winding-to-core capacitance ($C_{wc}$) create common-mode noise paths:

$$I_{CM} = C_{ps} \cdot \frac{dV}{dt}$$

To minimize common-mode noise:
- Use shield windings (Faraday shields)
- Increase primary-secondary spacing
- Consider winding orientation

### Y-Capacitors

External Y-capacitors can shunt common-mode currents but are limited by safety
regulations (typically 4.7nF max for Class II equipment).

## Usage

```cpp
#include "physical_models/StrayCapacitance.h"

// Configure capacitance model
auto& settings = OpenMagnetics::Settings::GetInstance();
settings.set_stray_capacitance_model(
    OpenMagnetics::StrayCapacitanceModels::ALBACH
);

// Calculate capacitance
auto capacitance = OpenMagnetics::StrayCapacitance();
auto C_stray = capacitance.calculate_stray_capacitance(magnetic, operatingPoint);

// Get individual components
auto C_tt = capacitance.get_turn_to_turn_capacitance();
auto C_ll = capacitance.get_layer_to_layer_capacitance();
auto C_wc = capacitance.get_winding_to_core_capacitance();
```

## Design Guidelines

### For EMI-Sensitive Applications
- Minimize primary-secondary capacitance
- Consider adding Faraday shields
- Place windings symmetrically

### For High-Frequency Inductors
- Keep SRF at least 10x operating frequency
- Use single-layer winding if possible
- Consider air-core for very high frequencies

### For Resonant Converters
- Capacitance may be a design parameter
- Account for capacitance in resonant tank design
- Verify with impedance analyzer measurements
