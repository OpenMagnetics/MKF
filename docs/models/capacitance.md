# Capacitance

Parasitic (stray) capacitance in magnetic components affects:
- High-frequency behavior and resonances
- EMI performance (common-mode noise paths)
- Switching transients in power converters
- Self-resonant frequency of inductors

MKF models several capacitance components:
- **Turn-to-turn capacitance:** Between adjacent turns
- **Layer-to-layer capacitance:** Between winding layers
- **Winding-to-core capacitance:** Between windings and grounded core
- **Primary-to-secondary capacitance:** Coupling between windings

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

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|

## Capacitance Minimization Strategies

| Technique | Effect | Trade-off |
|-----------|--------|-----------|
| Increase turn spacing | Reduces turn-to-turn C | Larger winding window |
| Section winding | Reduces layer-to-layer C | More complex construction |
| Shield windings | Controls coupling paths | Additional losses |
| Bank winding | Reduces distributed C | Higher leakage inductance |
| Thicker insulation | Reduces all C | Thermal resistance |

**Default Model:** Albach provides good balance of accuracy and generality.
