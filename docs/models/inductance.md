# Inductance

*Comprehensive documentation for inductance calculation models*

MKF calculates inductance using magnetic circuit analysis with reluctance-based methods.
The fundamental relationship is:

$$L = \frac{N^2}{\mathcal{R}_{total}}$$

Where $N$ is the number of turns and $\mathcal{R}_{total}$ is the total magnetic circuit
reluctance (sum of core and gap reluctances).

For coupled inductors and transformers, MKF uses **matrix theory** to compute the full
inductance matrix including self and mutual inductances.

## Magnetizing Inductance

Primary (magnetizing) inductance is calculated from magnetic circuit analysis:

$$L_m = \frac{N^2}{\mathcal{R}_{core} + \mathcal{R}_{gap}}$$

Where:
- **Core reluctance:** $\mathcal{R}_{core} = \frac{l_e}{\mu_0 \mu_r A_e}$
- **Gap reluctance:** Uses the selected reluctance model (Zhang, Muhlethaler, etc.)

The effective permeability with an air gap is:

$$\mu_{eff} = \frac{\mu_r}{1 + \mu_r \cdot l_g / l_e}$$

For gapped cores, the air gap dominates the reluctance, making inductance relatively
insensitive to core permeability variations.

## Inductance Matrix for Multi-Winding Components

For transformers with multiple windings, MKF computes the full inductance matrix:

$$[L] = [N]^T [\mathcal{P}] [N]$$

Where:
- $[N]$ is the turns matrix (diagonal for standard transformers)
- $[\mathcal{P}]$ is the permeance matrix

The inductance matrix elements are:
- **Diagonal elements ($L_{ii}$):** Self-inductance of winding $i$
- **Off-diagonal elements ($L_{ij}$):** Mutual inductance between windings $i$ and $j$

## Leakage Inductance

Leakage inductance represents flux that doesn't link all windings. It's critical for:
- Transformer voltage regulation
- Resonant converter design (LLC, phase-shifted full-bridge)
- EMI and voltage spikes

MKF calculates leakage inductance using energy-based methods:

$$L_{leak} = \frac{2 W_{leak}}{I^2} = \frac{1}{I^2} \int_{volume} \frac{B^2}{\mu_0} dV$$

The magnetic field distribution is computed using the selected magnetic field model
(Binns-Lawrenson, Dowell, etc.), then integrated over the winding window volume.

**Reference:** [Spreen, J.H. "Electrical Terminal Representation of Conductor Loss in Transformers." IEEE Trans. Power Electronics, 1990](https://ieeexplore.ieee.org/document/63145)

### Leakage Inductance Factors

Leakage inductance depends on:

| Factor | Effect | Design Implication |
|--------|--------|-------------------|
| Winding spacing | Linear increase | Minimize gaps between windings |
| Number of layers | Quadratic increase | Interleave to reduce |
| Winding width | Inverse | Use full bobbin width |
| Window height | Linear decrease | Taller windows help |
| Interleaving | Significant reduction | P-S-P-S beats P-P-S-S |

### Rogowski Factor

The Rogowski factor corrects for non-uniform field distribution at winding ends:

$$k_R = 1 - \frac{1 - e^{-\pi h/w}}{\pi h/w}$$

Where $h$ is winding height and $w$ is winding width.

## Configuration Settings

Key settings affecting inductance calculations:

| Setting | Effect |
|---------|--------|
| `reluctance_model` | Determines gap fringing accuracy |
| `magnetic_field_strength_model` | Affects leakage calculation |
| `magnetic_field_include_fringing` | Include gap fringing in field calc |
| `magnetizing_inductance_include_air_inductance` | Include air-core inductance |
| `leakage_inductance_grid_auto_scaling` | Auto-adjust grid density |
| `leakage_inductance_grid_precision_level_planar` | Grid precision for planar |
| `leakage_inductance_grid_precision_level_wound` | Grid precision for wound |

**For resonant converters:** Accurate leakage inductance is critical. Use Muhlethaler
reluctance model and Binns-Lawrenson field model with fine grid resolution.

## Usage

```cpp
#include "physical_models/Inductance.h"
#include "physical_models/LeakageInductance.h"

// Create magnetic component
auto magnetic = OpenMagnetics::Magnetic(...);

// Calculate magnetizing inductance
auto inductance = OpenMagnetics::MagnetizingInductance();
auto L_m = inductance.calculate_inductance(magnetic, operatingPoint);

// Calculate leakage inductance
auto leakage = OpenMagnetics::LeakageInductance();
auto L_leak = leakage.calculate_leakage_inductance(magnetic, operatingPoint);

// Get full inductance matrix for multi-winding transformer
auto L_matrix = inductance.calculate_inductance_matrix(magnetic, operatingPoint);
```

## Design Guidelines

### For Inductors
- Gap length determines inductance (assuming high-permeability core)
- Use Zhang or Muhlethaler model for accurate gap fringing
- Check saturation: $B_{peak} = \frac{L \cdot I_{peak}}{N \cdot A_e}$

### For Transformers
- Magnetizing inductance affects duty cycle in flyback
- Leakage inductance causes voltage spikes and affects regulation
- Interleave windings to reduce leakage (but increases capacitance)

### For Resonant Converters
- Leakage inductance is a design parameter, not a parasitic
- LLC: Leakage acts as series resonant inductor
- PSFB: Leakage enables ZVS transitions
- Use accurate models and verify with measurements
