# Winding Losses

*Comprehensive documentation for winding loss models*

Winding losses in magnetic components arise from two frequency-dependent mechanisms beyond
DC resistance:

1. **Skin Effect:** AC current concentrates near the conductor surface, reducing effective
   cross-sectional area and increasing resistance.

2. **Proximity Effect:** Magnetic fields from adjacent conductors induce eddy currents,
   causing additional losses that can dominate at high frequencies.

MKF provides multiple models for each effect, with the optimal choice depending on
wire type (round, litz, foil, rectangular) and frequency range.

## Skin Effect Models

The skin depth $\delta$ determines where current flows in a conductor:

$$\delta = \sqrt{\frac{\rho}{\pi f \mu}}$$

For frequencies where $\delta$ is comparable to or smaller than the conductor radius,
skin effect significantly increases AC resistance.

### Dowell

Dowell's classic 1D analysis for layered windings:

$$F_R = \frac{\xi}{2} \left[ \frac{\sinh(\xi) + \sin(\xi)}{\cosh(\xi) - \cos(\xi)} \right]$$

Where $\xi = h/\delta$ is the ratio of conductor height to skin depth.

Originally developed for foil windings, widely used for all wire types despite
being a 1D approximation.

**Reference:** [Dowell, P.L. 'Effects of eddy currents in transformer windings.' Proc. IEE, 1966](https://ieeexplore.ieee.org/document/5256365)

### Wojda

Wojda's extension improves accuracy for round conductors by using proper 2D field
solutions. Recommended for round wire windings at moderate frequencies.

**Reference:** [Wojda, R. 'Winding resistance and winding power loss of high-frequency power inductors.' IEEE Trans. Industry Appl., 2012](https://ieeexplore.ieee.org/document/6205273)

### Albach

Albach's comprehensive model includes:
- Full 2D field solution for round wires
- Improved accuracy at high frequencies
- Better handling of dense winding packing

Particularly accurate for litz wire bundles.

**Reference:** [Albach et al. 'Calculating core losses in transformers.' IEEE Trans. Magnetics, 2011](https://ieeexplore.ieee.org/document/5680877)

### Payne

Payne's model provides empirically-validated corrections for practical winding
configurations, including non-ideal layer stacking and partial layers.

**Reference:** [Payne, A. 'Skin Effect, Proximity Effect and the Resistance of Conductors.' 2014](http://g3rbj.co.uk/wp-content/uploads/2014/05/Skin-Effect-revised.pdf)

### Kutkut

Kutkut's model for litz wire accounts for bundle-level and strand-level skin effects
separately. Essential for accurate litz wire loss predictions at high frequencies.

**Reference:** [Kutkut et al. 'AC resistance of planar power inductors.' IEEE Trans. Power Electronics, 1996](https://ieeexplore.ieee.org/document/498992)

### Ferreira

Ferreira's analytical approach provides closed-form expressions suitable for
optimization algorithms. Good for initial design iterations.

**Reference:** [Ferreira, J.A. 'Improved analytical modeling of conductive losses.' IEEE Trans. Power Electronics, 1994](https://ieeexplore.ieee.org/document/80877)

## Proximity Effect Models

Proximity effect losses arise from eddy currents induced by external magnetic fields.
For multi-layer windings, these losses can far exceed skin effect losses:

$$P_{prox} \propto H^2 \cdot f^2$$

Where $H$ is the magnetic field from adjacent conductors.

### Rossmanith

Rossmanith's model uses orthogonal field decomposition for improved accuracy
in complex winding arrangements. Handles both radial and axial field components.

**Reference:** [Rossmanith et al. 'Improved Proximity-Effect Loss Calculation.' IEEE Trans. Industry Appl., 2009](https://ieeexplore.ieee.org/document/4912922)

### Lammeraner

Lammeraner's classical analysis provides the foundation for proximity effect
calculations in transformer windings.

**Reference:** [Lammeraner, J. 'Eddy Currents.' Iliffe Books, 1966](https://link.springer.com/book/10.1007/978-94-011-7095-6)

### Dowell (Proximity)

Dowell's 1D model extended to proximity effect. The loss factor increases
dramatically with number of layers:

$$F_R = \frac{\xi}{2} \left[ ... + \frac{2(m^2-1)}{3}(...) \right]$$

Where $m$ is the number of layers. For 10 layers, proximity losses can be
50x higher than skin effect alone.

**Reference:** [Dowell, P.L. 'Effects of eddy currents in transformer windings.' Proc. IEE, 1966](https://ieeexplore.ieee.org/document/5256365)

### Albach (Proximity)

Albach's 2D field solution for proximity effect in round conductors.
Provides better accuracy than 1D models, especially for sparse windings.

**Reference:** [Albach et al. 'Calculating core losses in transformers.' IEEE Trans. Magnetics, 2011](https://ieeexplore.ieee.org/document/5680877)

### Ferreira (Proximity)

Ferreira's proximity effect model complements his skin effect model,
providing consistent analytical framework for winding loss optimization.

**Reference:** [Ferreira, J.A. 'Improved analytical modeling of conductive losses.' IEEE Trans. Power Electronics, 1994](https://ieeexplore.ieee.org/document/80877)

## Model Selection by Wire Type

| Wire Type | Skin Model | Proximity Model | Notes |
|-----------|------------|-----------------|-------|
| Round solid | Wojda or Albach | Rossmanith | 2D models preferred |
| Litz wire | Kutkut or Albach | Albach | Strand-level effects |
| Foil | Dowell | Dowell | 1D model appropriate |
| Rectangular | Dowell or Payne | Dowell | Layer-based analysis |
| Planar PCB | Dowell | Dowell | Similar to foil |

## Design Tips

- Use litz wire when $d/\delta > 1.5$ (wire diameter to skin depth ratio)
- Interleave primary and secondary to reduce proximity effect
- Fewer layers with thicker wire often beats many thin layers
- Consider foil for very high currents at lower frequencies

## Usage

```cpp
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingProximityEffectLosses.h"

// Configure skin effect model
auto& settings = OpenMagnetics::Settings::GetInstance();
settings.set_winding_skin_effect_losses_model(
    OpenMagnetics::WindingSkinEffectLossesModels::ALBACH
);

// Configure proximity effect model
settings.set_winding_proximity_effect_losses_model(
    OpenMagnetics::WindingProximityEffectLossesModels::ROSSMANITH
);

// Or use factory for specific models
auto skinModel = OpenMagnetics::WindingSkinEffectLossesModel::factory(
    OpenMagnetics::WindingSkinEffectLossesModels::KUTKUT
);
```
