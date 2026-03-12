# Magnetic Field

Accurate magnetic field distribution calculations are essential for:
- Winding loss predictions (proximity effect depends on local field)
- Leakage inductance calculations
- Thermal analysis (loss distribution)
- Saturation checking (local flux density)

MKF implements several magnetic field models, each with different trade-offs between
accuracy and computational cost.

## Available Models

### Binns-Lawrenson

The Binns-Lawrenson method uses Fourier series expansion to solve the magnetic field
in the winding window:

$$H(x,y) = \sum_{n=1}^{\infty} \sum_{m=1}^{\infty} H_{nm} \sin(n\pi x/a) \sin(m\pi y/b)$$

**Advantages:**
- Handles arbitrary conductor positions
- Good accuracy for typical transformer geometries
- Computationally efficient

**Default model in MKF** for most calculations.

**Reference:** [Binns, Lawrenson. 'Analysis and Computation of Electric and Magnetic Field Problems.' Pergamon Press](https://onlinelibrary.wiley.com/doi/book/10.1002/0470020423)

### Dowell

Dowell's 1D field model assumes uniform current sheets:

$$H(y) = \frac{N \cdot I}{w} \cdot \frac{y}{h}$$

Where $w$ is winding width and $h$ is window height.

**Limitations:** 1D approximation, best for full-width layers.
**Advantage:** Very fast computation.

**Reference:** [Dowell, P.L. 'Effects of eddy currents in transformer windings.' Proc. IEE, 1966](https://ieeexplore.ieee.org/document/5256365)

### Wang

Wang's model provides improved 2D solutions for planar magnetic structures,
accounting for end effects and non-ideal layer positioning.

**Reference:** [Wang et al. 'Analysis of Planar E-Core Transformer.' IEEE Trans. Magnetics, 2009](https://ieeexplore.ieee.org/document/4798971)

### Albach

Albach's 2D analytical solution for round conductors provides accurate field
distributions even with partial layers and non-uniform spacing.

**Reference:** [Albach et al. 'Calculating core losses in transformers.' IEEE Trans. Magnetics, 2011](https://ieeexplore.ieee.org/document/5680877)

## Fringing Field Models

Gap fringing affects the field distribution near air gaps, causing:
- Increased losses in conductors near the gap
- Modified flux density distribution
- Potential localized saturation

### Roshen Fringing Model

Roshen's analytical model for fringing field around air gaps:

$$H_{fringe}(r) = H_g \cdot \frac{2}{\pi} \arctan\left(\frac{l_g}{2r}\right)$$

Where $r$ is the distance from gap edge.

**Reference:** Roshen, W. "Fringing Field Formulas and Winding Loss Due to Fringing Field."
IEEE Trans. Magnetics, 2007. [IEEE](https://ieeexplore.ieee.org/document/4276863)

### Albach Fringing Model

Albach's 2D field solution includes fringing effects in the overall field calculation,
providing smooth transitions between gap and non-gap regions.

**Design Tip:** Keep conductors at least 2-3 gap lengths away from the air gap to
minimize fringing-induced losses.

## Model Comparison

| Model | Error | Reference |
|-------|-------|-----------|

## Model Selection Guide

| Application | Recommended Model | Notes |
|-------------|-------------------|-------|
| Standard transformers | Binns-Lawrenson | Best general accuracy |
| Quick estimates | Dowell | Fast 1D approximation |
| Planar magnetics | Wang or Albach | Handle wide, flat windings |
| Near air gaps | Include fringing models | Critical for loss accuracy |

**Grid Resolution Settings:**
- `magnetic_field_number_points_x`: X-axis grid points (default: 20)
- `magnetic_field_number_points_y`: Y-axis grid points (default: 20)

Increase for more accurate leakage inductance; decrease for faster computation.
