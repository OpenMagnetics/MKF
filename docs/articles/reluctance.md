# Air Gap Reluctance Models: From Classic to Zhang

*How fringing fields affect your inductance calculations*

*Published: 2026-02-03*

---

The air gap in a magnetic core is where energy is stored in an inductor. Accurately calculating the gap reluctance is essential for predicting inductance and ensuring your design meets specifications.

However, the magnetic field doesn't simply pass straight through the gap - it fringes outward, effectively increasing the gap's cross-sectional area. Different models account for this fringing in different ways.

## Classic Model

The classic reluctance formula assumes a uniform magnetic field:

$$R_g = \frac{l_g}{\mu_0 \cdot A_c}$$

Where:
- $R_g$ is the gap reluctance (H⁻¹)
- $l_g$ is the gap length (m)
- $A_c$ is the core cross-sectional area (m²)
- $\mu_0$ is the permeability of free space

This model ignores fringing entirely, leading to overestimated reluctance and underestimated inductance. It's useful only for rough estimates or very small gaps.

## Zhang Model

Xinsheng Zhang's model (2020) treats the total reluctance as a parallel combination of internal and fringing components:

$$R_g = R_{in} \parallel R_{fr}$$

The internal reluctance handles the direct flux path, while the fringing reluctance accounts for flux that spreads around the gap edges:

$$R_{fr} = \frac{\pi}{\mu_0 \cdot C \cdot \ln\left(\frac{2h + l_g}{l_g}\right)}$$

Where $C$ is the perimeter of the core cross-section and $h$ is the winding window height. This model provides excellent accuracy for typical power magnetics applications and is the default in OpenMagnetics.

## Mühlethaler Model

Jonas Mühlethaler's 3D approach divides the fringing field into distinct regions, each with its own reluctance contribution. This is particularly valuable for:

- Planar cores with complex geometries
- Very large gaps (relative to core dimensions)
- Cores with multiple gaps

The model requires more geometric information but provides superior accuracy for challenging designs.

## Practical Recommendations

For typical E-cores and similar shapes with moderate gaps, **Zhang** provides the best accuracy with minimal complexity. It's the recommended default.

Use **Mühlethaler** when designing with planar cores, unusual geometries, or very large gaps where standard models show significant error.

The **Classic** model is useful for quick estimates or when you want to understand the baseline reluctance before fringing effects.

Remember: fringing always reduces reluctance (increases inductance). If your measured inductance exceeds calculations, fringing is likely the cause.

---

## References

- Standard magnetic circuit theory. [Wikipedia](https://en.wikipedia.org/wiki/Magnetic_reluctance)
- Zhang, X. et al. "Improved Calculation Method for Inductance Value of the Air-Gap Inductor." CIYCEE, 2020. [IEEE](https://ieeexplore.ieee.org/document/9332553)
- Mühlethaler, J. "A Novel Approach for 3D Air Gap Reluctance Calculations." ECCE Asia, 2011. [Paper](https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf)

---

*This article was generated from the [OpenMagnetics](https://github.com/OpenMagnetics/MKF) documentation. OpenMagnetics is an open-source framework for magnetic component design.*