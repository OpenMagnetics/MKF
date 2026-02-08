# Understanding Core Loss Models in Magnetic Design

*A deep dive into Steinmetz, iGSE, and beyond*

*Published: 2026-02-03*

---

When designing magnetic components like transformers and inductors, accurately predicting core losses is crucial for thermal management and efficiency optimization. Core losses arise from hysteresis and eddy currents in the magnetic material, and their magnitude depends on the operating frequency, flux density, and material properties.

In this article, we'll explore the various core loss models implemented in OpenMagnetics, their theoretical foundations, and when to use each one.

## Steinmetz Equation

The Steinmetz equation, first proposed by Charles Proteus Steinmetz in 1892, is the foundation of most core loss calculations:

$$P_v = k \cdot f^\alpha \cdot B^\beta$$

Where:
- $P_v$ is the volumetric power loss (W/m³)
- $f$ is the frequency (Hz)
- $B$ is the peak magnetic flux density (T)
- $k$, $\alpha$, $\beta$ are material-specific coefficients

The beauty of Steinmetz lies in its simplicity - manufacturers provide the k, α, β parameters for their materials, making loss estimation straightforward. However, it assumes sinusoidal excitation, which limits its applicability to power electronics with non-sinusoidal waveforms.

## Improved Generalized Steinmetz Equation (iGSE)

The iGSE, developed by Charles Sullivan at Dartmouth, extends Steinmetz to handle arbitrary waveforms:

$$P_v = \frac{1}{T} \int_0^T k_i \left| \frac{dB}{dt} \right|^\alpha (\Delta B)^{\beta - \alpha} dt$$

The key insight is that core losses depend on the rate of change of flux density (dB/dt), not just its peak value. This allows accurate predictions for:
- Triangular waveforms (typical in inductors)
- Trapezoidal waveforms (common in half-bridge converters)
- Arbitrary PWM waveforms

iGSE uses the same Steinmetz parameters (k, α, β) but applies them to instantaneous dB/dt, making it practical since no additional material characterization is needed.

## Roshen Model

Waseem Roshen's model takes a more physics-based approach, separating hysteresis and eddy current losses:

$$P_v = P_{hyst} + P_{eddy} + P_{excess}$$

The hysteresis component is derived from the material's B-H loop, while eddy current losses are calculated from the material's resistivity and geometry. This separation provides better accuracy when:
- Operating far from the Steinmetz characterization conditions
- Using materials with strong temperature dependence
- Designing with DC bias

However, Roshen requires more material parameters than Steinmetz-based methods.

## Choosing the Right Model

For most designs, **iGSE** offers the best balance of accuracy and practicality. It handles non-sinusoidal waveforms well and uses readily available Steinmetz parameters.

Use **Steinmetz** when you have purely sinusoidal excitation or need a quick estimate.

Consider **Roshen** when you need physics-based accuracy, especially with significant DC bias or extreme operating conditions.

OpenMagnetics automatically selects the most appropriate model based on available material data, but you can override this to compare predictions or match specific validation data.

---

## References

- Steinmetz, C.P. "On the law of hysteresis." AIEE Transactions, 1892. [IEEE](https://ieeexplore.ieee.org/document/1457110)
- Li, J., Abdallah, T., Sullivan, C.R. "Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters." [Paper](http://inductor.thayerschool.org/papers/IGSE.pdf)
- Roshen, W. "Ferrite Core Loss for Power Magnetic Components Design." IEEE Trans. Magnetics, 1991. [IEEE](https://ieeexplore.ieee.org/document/4052433)

---

*This article was generated from the [OpenMagnetics](https://github.com/OpenMagnetics/MKF) documentation. OpenMagnetics is an open-source framework for magnetic component design.*