# Inductance

MKF calculates inductance using magnetic circuit analysis with reluctance-based methods.
The fundamental relationship is:

$$L = \frac{N^2}{\mathcal{R}_{total}}$$

Where $N$ is the number of turns and $\mathcal{R}_{total}$ is the total magnetic circuit
reluctance (sum of core and gap reluctances).

For coupled inductors and transformers, MKF uses **matrix theory** to compute the full
inductance matrix including self and mutual inductances.

## Available Models

### Reluctance-based

Primary inductance calculation from magnetic circuit analysis:

$$L = \frac{N^2}{\mathcal{R}_{core} + \mathcal{R}_{gap}}$$

Core reluctance: $\mathcal{R}_{core} = \frac{l_e}{\mu_0 \mu_r A_e}$

Gap reluctance: Uses selected reluctance model (Zhang, Mühlethaler, etc.)

For multi-winding components, the inductance matrix is:

$$[L] = [N]^T [\mathcal{P}] [N]$$

Where $[\mathcal{P}]$ is the permeance matrix and $[N]$ is the turns matrix.

## Leakage Inductance

Leakage inductance represents flux that doesn't link all windings. It's critical for:
- Transformer voltage regulation
- Resonant converter design (LLC, phase-shifted full-bridge)
- EMI and voltage spikes

MKF calculates leakage inductance using energy-based methods:

$$L_{leak} = \frac{2 W_{leak}}{I^2} = \frac{1}{I^2} \int_{volume} \frac{B^2}{\mu_0} dV$$

The magnetic field distribution is computed using the selected magnetic field model
(Binns-Lawrenson, Dowell, etc.), then integrated over the winding window volume.

**Reference:** Spreen, J.H. "Electrical Terminal Representation of Conductor Loss in
Transformers." IEEE Trans. Power Electronics, 1990.
[IEEE](https://ieeexplore.ieee.org/document/63145)

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|

## Inductance Calculation Settings

Key settings affecting inductance calculations:

| Setting | Effect |
|---------|--------|
| `reluctance_model` | Determines gap fringing accuracy |
| `magnetic_field_strength_model` | Affects leakage calculation |
| `magnetic_field_include_fringing` | Include gap fringing in field calc |
| `leakage_inductance_grid_auto_scaling` | Auto-adjust grid density |

**For resonant converters:** Accurate leakage inductance is critical. Use Mühlethaler
reluctance model and Binns-Lawrenson field model with fine grid resolution.
