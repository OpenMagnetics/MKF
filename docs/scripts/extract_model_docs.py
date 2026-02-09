#!/usr/bin/env python3
"""
Extract model documentation from C++ source files and generate Markdown documentation.

This script parses the physical model header files to extract:
- Model descriptions from get_models_information()
- Model errors from get_models_errors()
- External links from get_models_external_links()
- Internal links from get_models_internal_links()

It also includes comprehensive hardcoded documentation for models
where the embedded documentation is insufficient.
"""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class ModelInfo:
    """Information about a single model implementation."""
    name: str
    description: str
    error: Optional[float] = None
    external_link: Optional[str] = None
    internal_link: Optional[str] = None


@dataclass
class ModelCategory:
    """A category of models (e.g., CoreLosses, Reluctance)."""
    name: str
    class_name: str
    models: List[ModelInfo]
    source_file: str


# Comprehensive documentation for model categories
# This supplements the auto-extracted content with scientific context

# Constants for duplicated string literals
REF_ALBACH_2011 = "https://ieeexplore.ieee.org/document/5680877"
REF_TEXT_ALBACH = "Albach et al. 'Calculating core losses in transformers.' IEEE Trans. Magnetics, 2011"
REF_DOWELL_1966 = "https://ieeexplore.ieee.org/document/5256365"
REF_TEXT_DOWELL = "Dowell, P.L. 'Effects of eddy currents in transformer windings.' Proc. IEE, 1966"
LINK_WINDING_LOSSES = "winding-losses.md"
LINK_INDUCTANCE = "inductance.md"
LINK_THERMAL = "thermal.md"
CODE_BLOCK_CPP = "```cpp"

DETAILED_DOCS = {
    "core-losses.md": {
        "intro": """
Core losses in magnetic materials arise from two primary mechanisms: **hysteresis losses** (energy dissipated
in reorienting magnetic domains) and **eddy current losses** (currents induced in the conductive material).
Accurate prediction of these losses is critical for thermal management and efficiency optimization in
transformers, inductors, and other magnetic components.

MKF implements multiple core loss models, each with different trade-offs between accuracy and complexity.
The choice of model depends on the waveform type, available material data, and required precision.
""",
        "models": {
            "Steinmetz": {
                "description": """
The original Steinmetz equation (1892) empirically relates core losses to frequency and flux density:

$$P_v = k \\cdot f^\\alpha \\cdot \\hat{B}^\\beta$$

Where:
- $P_v$ is volumetric power loss (W/m³)
- $f$ is frequency (Hz)
- $\\hat{B}$ is peak flux density (T)
- $k$, $\\alpha$, $\\beta$ are material-specific Steinmetz coefficients

**Limitations:** Only valid for sinusoidal excitation. Using this model with non-sinusoidal
waveforms (triangular, trapezoidal, PWM) will produce significant errors.
""",
                "reference": "https://ieeexplore.ieee.org/document/4768063",
                "reference_text": "Steinmetz, C.P. 'On the law of hysteresis.' AIEE Trans., 1892"
            },
            "iGSE": {
                "description": """
The **improved Generalized Steinmetz Equation** (iGSE) extends Steinmetz to arbitrary waveforms:

$$P_v = \\frac{1}{T} \\int_0^T k_i \\left| \\frac{dB}{dt} \\right|^\\alpha (\\Delta B)^{\\beta - \\alpha} dt$$

Where $k_i = \\frac{k}{(2\\pi)^{\\alpha-1} \\int_0^{2\\pi} |\\cos\\theta|^\\alpha 2^{\\beta-\\alpha} d\\theta}$

**Key insight:** Core losses depend on the rate of change of flux density (dB/dt), not just
peak values. This allows accurate predictions for:
- Triangular waveforms (inductors in CCM)
- Trapezoidal waveforms (half-bridge converters)
- Arbitrary PWM waveforms

**Advantage:** Uses the same Steinmetz parameters (k, α, β) as the original equation,
requiring no additional material characterization.
""",
                "reference": "https://inductor.thayerschool.org/papers/IGSE.pdf",
                "reference_text": "Li, Abdallah, Sullivan. 'Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms.' IEEE COMPEL, 2001"
            },
            "MSE": {
                "description": """
The **Modified Steinmetz Equation** (MSE) uses an equivalent frequency concept:

$$P_v = k \\cdot f_{eq}^{\\alpha-1} \\cdot f_r \\cdot \\hat{B}^\\beta$$

Where $f_{eq} = \\frac{2}{\\Delta B^2 \\pi^2} \\int_0^T \\left(\\frac{dB}{dt}\\right)^2 dt$

The equivalent frequency captures the effect of non-sinusoidal waveforms by computing
what sinusoidal frequency would produce the same losses.
""",
                "reference": "https://www.psma.com/sites/default/files/uploads/files/Limits%20of%20Accuracy%20of%20Steinmetz%20Equation-Handout.pdf",
                "reference_text": "Venkatachalam et al. 'Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms.' IEEE COMPEL, 2002"
            },
            "NSE": {
                "description": """
The **Natural Steinmetz Extension** (NSE) uses direct integration without equivalent frequency:

$$P_v = \\left(\\frac{\\Delta B}{2}\\right)^{\\beta-\\alpha} \\frac{k}{T} \\int_0^T \\left|\\frac{dB}{dt}\\right|^\\alpha dt$$

This provides a more physically meaningful interpretation by relating losses directly
to the actual flux rate of change over the entire period.
""",
                "reference": "https://ieeexplore.ieee.org/document/4544745",
                "reference_text": "Van den Bossche et al. 'Evaluation of the Steinmetz Parameters.' IEEE Trans. Magnetics, 2008"
            },
            "Roshen": {
                "description": """
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
""",
                "reference": "https://ieeexplore.ieee.org/document/4052433",
                "reference_text": "Roshen, W. 'Ferrite Core Loss for Power Magnetic Components Design.' IEEE Trans. Magnetics, 1991"
            },
            "Albach": {
                "description": """
Albach's model combines loss separation with improved high-frequency accuracy:

$$P_v = C_h f B^\\gamma + C_e (fB)^2 / \\rho + C_a (fB)^{1.5}$$

Where:
- $C_h$, $\\gamma$ are hysteresis parameters
- $C_e$ is the classical eddy current coefficient
- $C_a$ accounts for anomalous (excess) losses
- $\\rho$ is material resistivity

This provides good accuracy across a wide frequency range.
""",
                "reference": REF_ALBACH_2011,
                "reference_text": REF_TEXT_ALBACH
            },
            "Barg": {
                "description": """
Barg's model extends loss separation for modern ferrite materials with improved
temperature and frequency dependence modeling. It separates static and dynamic
hysteresis contributions for better accuracy at high frequencies.
""",
                "reference": "https://ieeexplore.ieee.org/document/9530729",
                "reference_text": "Barg et al. 'A Dynamic Hysteresis Loss Model.' IEEE Trans. Power Electronics, 2021"
            },
            "Proprietary": {
                "description": """
Manufacturer-specific loss models that use detailed material characterization data
provided by core manufacturers. These models may include:
- Multi-temperature Steinmetz coefficients
- DC bias correction factors
- Frequency-dependent parameter adjustment

Contact the core manufacturer for specific model parameters.
""",
                "reference": "",
                "reference_text": "Manufacturer-specific data sheets"
            }
        },
        "selection_guide": """
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
"""
    },

    "reluctance.md": {
        "intro": """
The air gap in a magnetic core is where energy is stored in an inductor. Accurately calculating
gap reluctance is essential for predicting inductance and ensuring designs meet specifications.

The magnetic field doesn't pass straight through the gap—it **fringes outward**, effectively
increasing the gap's cross-sectional area. Different models account for this fringing in
different ways, with significant impact on inductance predictions.

$$R_{gap} = \\frac{l_g}{\\mu_0 \\cdot A_{eff}}$$

Where $A_{eff} > A_{core}$ due to fringing, and the various models differ in how they
compute this effective area.
""",
        "models": {
            "Classic": {
                "description": """
The classic reluctance formula assumes uniform magnetic field (no fringing):

$$R_g = \\frac{l_g}{\\mu_0 \\cdot A_c}$$

Where:
- $R_g$ is gap reluctance (H⁻¹)
- $l_g$ is gap length (m)
- $A_c$ is core cross-sectional area (m²)
- $\\mu_0$ is permeability of free space

**Limitation:** Ignores fringing entirely, leading to overestimated reluctance
and underestimated inductance. Only useful for rough estimates or very small gaps.
""",
                "reference": "https://en.wikipedia.org/wiki/Magnetic_reluctance",
                "reference_text": "Standard magnetic circuit theory"
            },
            "Zhang": {
                "description": """
Zhang's model (2020) treats total reluctance as parallel combination of internal and fringing:

$$R_g = R_{in} \\parallel R_{fr}$$

The internal reluctance handles direct flux, while fringing reluctance accounts for
flux spreading around gap edges:

$$R_{fr} = \\frac{\\pi}{\\mu_0 \\cdot C \\cdot \\ln\\left(\\frac{2h + l_g}{l_g}\\right)}$$

Where $C$ is the perimeter of core cross-section and $h$ is winding window height.

**Advantages:**
- Excellent accuracy for typical power magnetics
- Simple analytical formulation
- Default model in MKF

**Validated error:** ~3.6% mean deviation from measurements.
""",
                "reference": "https://ieeexplore.ieee.org/document/9332553",
                "reference_text": "Zhang et al. 'Improved Calculation Method for Inductance Value of the Air-Gap Inductor.' CIYCEE, 2020"
            },
            "Muehlethaler": {
                "description": """
Mühlethaler's 3D approach divides the fringing field into distinct geometric regions,
each with its own reluctance contribution:

$$R_{total} = \\left(\\sum_i \\frac{1}{R_i}\\right)^{-1}$$

Regions include:
- Direct gap path
- Corner fringing (4 corners)
- Edge fringing (4 edges)
- Top/bottom fringing

**Best for:**
- Planar cores with complex geometries
- Very large gaps (relative to core dimensions)
- Cores with multiple distributed gaps
- EI, EE, and planar E cores

**Validated error:** ~1.4% mean deviation (best accuracy).
""",
                "reference": "https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf",
                "reference_text": "Mühlethaler et al. 'A Novel Approach for 3D Air Gap Reluctance Calculations.' ECCE Asia, 2011"
            },
            "Partridge": {
                "description": """
Partridge's model uses empirical correction factors for fringing:

$$R_g = \\frac{l_g}{\\mu_0 \\cdot A_c \\cdot F}$$

Where $F$ is a fringing factor that depends on gap-to-core dimension ratios.
Simple and fast for quick estimates.
""",
                "reference": "https://ieeexplore.ieee.org/document/1085",
                "reference_text": "Partridge, 'Analysis of magnetic circuits with air gaps.' AIEE Trans., 1937"
            },
            "Stenglein": {
                "description": """
Stenglein's model provides improved accuracy for gapped ferrite cores by including
detailed geometric factors and validated against extensive FEA simulations.
Particularly accurate for E-cores and pot cores.
""",
                "reference": "https://ieeexplore.ieee.org/document/8423500",
                "reference_text": "Stenglein et al. 'The Reluctance of Large Air Gaps in Ferrite Cores.' APEC, 2018"
            },
            "Balakrishnan": {
                "description": """
Balakrishnan's model uses a permeance-based approach with geometric correction factors.
Well-suited for planar magnetics and integrated magnetics applications.
""",
                "reference": "https://ieeexplore.ieee.org/document/602562",
                "reference_text": "Balakrishnan et al. 'Air-gap reluctance and inductance calculations.' IEEE PESC, 1997"
            }
        },
        "selection_guide": """
## Model Selection Guide

| Core Type | Recommended Model | Notes |
|-----------|-------------------|-------|
| E-cores (standard) | Zhang | Good balance of accuracy and simplicity |
| Planar E-cores | Mühlethaler | Handles wide, flat geometries |
| Pot cores | Zhang or Stenglein | Enclosed geometry |
| Toroids | Classic | No air gap, or distributed gap |
| RM cores | Zhang | Similar to E-cores |
| PQ cores | Zhang or Mühlethaler | Choose based on gap size |
| Very large gaps | Mühlethaler | >20% of center leg width |

**Default Recommendation:** Use **Zhang** for most applications. Switch to **Mühlethaler**
for planar cores or when maximum accuracy is required.

**Remember:** Fringing always *reduces* reluctance (increases inductance). If measured
inductance exceeds calculations, fringing is likely underestimated.
"""
    },

    "winding-losses.md": {
        "intro": """
Winding losses in magnetic components arise from two frequency-dependent mechanisms beyond
DC resistance:

1. **Skin Effect:** AC current concentrates near the conductor surface, reducing effective
   cross-sectional area and increasing resistance.

2. **Proximity Effect:** Magnetic fields from adjacent conductors induce eddy currents,
   causing additional losses that can dominate at high frequencies.

MKF provides multiple models for each effect, with the optimal choice depending on
wire type (round, litz, foil, rectangular) and frequency range.
""",
        "skin_effect_intro": """
## Skin Effect Models

The skin depth $\\delta$ determines where current flows in a conductor:

$$\\delta = \\sqrt{\\frac{\\rho}{\\pi f \\mu}}$$

For frequencies where $\\delta$ is comparable to or smaller than the conductor radius,
skin effect significantly increases AC resistance.
""",
        "models": {
            "Dowell (Skin)": {
                "description": """
Dowell's classic 1D analysis for layered windings:

$$F_R = \\frac{\\xi}{2} \\left[ \\frac{\\sinh(\\xi) + \\sin(\\xi)}{\\cosh(\\xi) - \\cos(\\xi)} \\right]$$

Where $\\xi = h/\\delta$ is the ratio of conductor height to skin depth.

Originally developed for foil windings, widely used for all wire types despite
being a 1D approximation.
""",
                "reference": REF_DOWELL_1966,
                "reference_text": REF_TEXT_DOWELL
            },
            "Wojda (Skin)": {
                "description": """
Wojda's extension improves accuracy for round conductors by using proper 2D field
solutions. Recommended for round wire windings at moderate frequencies.
""",
                "reference": "https://ieeexplore.ieee.org/document/6205273",
                "reference_text": "Wojda, R. 'Winding resistance and winding power loss of high-frequency power inductors.' IEEE Trans. Industry Appl., 2012"
            },
            "Albach (Skin)": {
                "description": """
Albach's comprehensive model includes:
- Full 2D field solution for round wires
- Improved accuracy at high frequencies
- Better handling of dense winding packing

Particularly accurate for litz wire bundles.
""",
                "reference": REF_ALBACH_2011,
                "reference_text": REF_TEXT_ALBACH
            },
            "Payne (Skin)": {
                "description": """
Payne's model provides empirically-validated corrections for practical winding
configurations, including non-ideal layer stacking and partial layers.
""",
                "reference": "https://g3rbj.co.uk/wp-content/uploads/2014/05/Skin-Effect-revised.pdf",
                "reference_text": "Payne, A. 'Skin Effect, Proximity Effect and the Resistance of Conductors.' 2014"
            },
            "Kutkut (Skin)": {
                "description": """
Kutkut's model for litz wire accounts for bundle-level and strand-level skin effects
separately. Essential for accurate litz wire loss predictions at high frequencies.
""",
                "reference": "https://ieeexplore.ieee.org/document/498992",
                "reference_text": "Kutkut et al. 'AC resistance of planar power inductors.' IEEE Trans. Power Electronics, 1996"
            },
            "Ferreira (Skin)": {
                "description": """
Ferreira's analytical approach provides closed-form expressions suitable for
optimization algorithms. Good for initial design iterations.
""",
                "reference": "https://ieeexplore.ieee.org/document/80877",
                "reference_text": "Ferreira, J.A. 'Improved analytical modeling of conductive losses.' IEEE Trans. Power Electronics, 1994"
            }
        },
        "proximity_effect_intro": """
## Proximity Effect Models

Proximity effect losses arise from eddy currents induced by external magnetic fields.
For multi-layer windings, these losses can far exceed skin effect losses:

$$P_{prox} \\propto H^2 \\cdot f^2$$

Where $H$ is the magnetic field from adjacent conductors.
""",
        "proximity_models": {
            "Rossmanith": {
                "description": """
Rossmanith's model uses orthogonal field decomposition for improved accuracy
in complex winding arrangements. Handles both radial and axial field components.
""",
                "reference": "https://ieeexplore.ieee.org/document/4912922",
                "reference_text": "Rossmanith et al. 'Improved Proximity-Effect Loss Calculation.' IEEE Trans. Industry Appl., 2009"
            },
            "Lammeraner": {
                "description": """
Lammeraner's classical analysis provides the foundation for proximity effect
calculations in transformer windings.
""",
                "reference": "https://link.springer.com/book/10.1007/978-94-011-7095-6",
                "reference_text": "Lammeraner, J. 'Eddy Currents.' Iliffe Books, 1966"
            },
            "Dowell (Proximity)": {
                "description": """
Dowell's 1D model extended to proximity effect. The loss factor increases
dramatically with number of layers:

$$F_R = \\frac{\\xi}{2} \\left[ ... + \\frac{2(m^2-1)}{3}(...) \\right]$$

Where $m$ is the number of layers. For 10 layers, proximity losses can be
50x higher than skin effect alone.
""",
                "reference": REF_DOWELL_1966,
                "reference_text": REF_TEXT_DOWELL
            },
            "Albach (Proximity)": {
                "description": """
Albach's 2D field solution for proximity effect in round conductors.
Provides better accuracy than 1D models, especially for sparse windings.
""",
                "reference": REF_ALBACH_2011,
                "reference_text": REF_TEXT_ALBACH
            },
            "Ferreira (Proximity)": {
                "description": """
Ferreira's proximity effect model complements his skin effect model,
providing consistent analytical framework for winding loss optimization.
""",
                "reference": "https://ieeexplore.ieee.org/document/80877",
                "reference_text": "Ferreira, J.A. 'Improved analytical modeling of conductive losses.' IEEE Trans. Power Electronics, 1994"
            }
        },
        "selection_guide": """
## Model Selection by Wire Type

| Wire Type | Skin Model | Proximity Model | Notes |
|-----------|------------|-----------------|-------|
| Round solid | Wojda or Albach | Rossmanith | 2D models preferred |
| Litz wire | Kutkut or Albach | Albach | Strand-level effects |
| Foil | Dowell | Dowell | 1D model appropriate |
| Rectangular | Dowell or Payne | Dowell | Layer-based analysis |
| Planar PCB | Dowell | Dowell | Similar to foil |

**Design Tips:**
- Use litz wire when $d/\\delta > 1.5$ (wire diameter to skin depth ratio)
- Interleave primary and secondary to reduce proximity effect
- Fewer layers with thicker wire often beats many thin layers
- Consider foil for very high currents at lower frequencies
"""
    },

    "inductance.md": {
        "intro": """
MKF calculates inductance using magnetic circuit analysis with reluctance-based methods.
The fundamental relationship is:

$$L = \\frac{N^2}{\\mathcal{R}_{total}}$$

Where $N$ is the number of turns and $\\mathcal{R}_{total}$ is the total magnetic circuit
reluctance (sum of core and gap reluctances).

For coupled inductors and transformers, MKF uses **matrix theory** to compute the full
inductance matrix including self and mutual inductances.
""",
        "models": {
            "Reluctance-based": {
                "description": """
Primary inductance calculation from magnetic circuit analysis:

$$L = \\frac{N^2}{\\mathcal{R}_{core} + \\mathcal{R}_{gap}}$$

Core reluctance: $\\mathcal{R}_{core} = \\frac{l_e}{\\mu_0 \\mu_r A_e}$

Gap reluctance: Uses selected reluctance model (Zhang, Mühlethaler, etc.)

For multi-winding components, the inductance matrix is:

$$[L] = [N]^T [\\mathcal{P}] [N]$$

Where $[\\mathcal{P}]$ is the permeance matrix and $[N]$ is the turns matrix.
""",
                "reference": "",
                "reference_text": "Classical magnetic circuit theory"
            }
        },
        "leakage_section": """
## Leakage Inductance

Leakage inductance represents flux that doesn't link all windings. It's critical for:
- Transformer voltage regulation
- Resonant converter design (LLC, phase-shifted full-bridge)
- EMI and voltage spikes

MKF calculates leakage inductance using energy-based methods:

$$L_{leak} = \\frac{2 W_{leak}}{I^2} = \\frac{1}{I^2} \\int_{volume} \\frac{B^2}{\\mu_0} dV$$

The magnetic field distribution is computed using the selected magnetic field model
(Binns-Lawrenson, Dowell, etc.), then integrated over the winding window volume.

**Reference:** Spreen, J.H. "Electrical Terminal Representation of Conductor Loss in
Transformers." IEEE Trans. Power Electronics, 1990.
[IEEE](https://ieeexplore.ieee.org/document/63145)
""",
        "selection_guide": """
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
"""
    },

    "magnetic-field.md": {
        "intro": """
Accurate magnetic field distribution calculations are essential for:
- Winding loss predictions (proximity effect depends on local field)
- Leakage inductance calculations
- Thermal analysis (loss distribution)
- Saturation checking (local flux density)

MKF implements several magnetic field models, each with different trade-offs between
accuracy and computational cost.
""",
        "models": {
            "Binns-Lawrenson": {
                "description": """
The Binns-Lawrenson method uses Fourier series expansion to solve the magnetic field
in the winding window:

$$H(x,y) = \\sum_{n=1}^{\\infty} \\sum_{m=1}^{\\infty} H_{nm} \\sin(n\\pi x/a) \\sin(m\\pi y/b)$$

**Advantages:**
- Handles arbitrary conductor positions
- Good accuracy for typical transformer geometries
- Computationally efficient

**Default model in MKF** for most calculations.
""",
                "reference": "https://onlinelibrary.wiley.com/doi/book/10.1002/0470020423",
                "reference_text": "Binns, Lawrenson. 'Analysis and Computation of Electric and Magnetic Field Problems.' Pergamon Press"
            },
            "Dowell": {
                "description": """
Dowell's 1D field model assumes uniform current sheets:

$$H(y) = \\frac{N \\cdot I}{w} \\cdot \\frac{y}{h}$$

Where $w$ is winding width and $h$ is window height.

**Limitations:** 1D approximation, best for full-width layers.
**Advantage:** Very fast computation.
""",
                "reference": REF_DOWELL_1966,
                "reference_text": REF_TEXT_DOWELL
            },
            "Wang": {
                "description": """
Wang's model provides improved 2D solutions for planar magnetic structures,
accounting for end effects and non-ideal layer positioning.
""",
                "reference": "https://ieeexplore.ieee.org/document/4798971",
                "reference_text": "Wang et al. 'Analysis of Planar E-Core Transformer.' IEEE Trans. Magnetics, 2009"
            },
            "Albach": {
                "description": """
Albach's 2D analytical solution for round conductors provides accurate field
distributions even with partial layers and non-uniform spacing.
""",
                "reference": REF_ALBACH_2011,
                "reference_text": REF_TEXT_ALBACH
            }
        },
        "fringing_section": """
## Fringing Field Models

Gap fringing affects the field distribution near air gaps, causing:
- Increased losses in conductors near the gap
- Modified flux density distribution
- Potential localized saturation

### Roshen Fringing Model

Roshen's analytical model for fringing field around air gaps:

$$H_{fringe}(r) = H_g \\cdot \\frac{2}{\\pi} \\arctan\\left(\\frac{l_g}{2r}\\right)$$

Where $r$ is the distance from gap edge.

**Reference:** Roshen, W. "Fringing Field Formulas and Winding Loss Due to Fringing Field."
IEEE Trans. Magnetics, 2007. [IEEE](https://ieeexplore.ieee.org/document/4276863)

### Albach Fringing Model

Albach's 2D field solution includes fringing effects in the overall field calculation,
providing smooth transitions between gap and non-gap regions.

**Design Tip:** Keep conductors at least 2-3 gap lengths away from the air gap to
minimize fringing-induced losses.
""",
        "selection_guide": """
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
"""
    },

    "capacitance.md": {
        "intro": """
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
""",
        "models": {
            "Koch": {
                "description": """
Koch's model uses parallel-plate approximation with fringing corrections:

$$C = \\varepsilon_0 \\varepsilon_r \\frac{A}{d} \\cdot k_{fringe}$$

Where $k_{fringe}$ accounts for field spreading at conductor edges.
Good for turn-to-turn and layer-to-layer capacitance.
""",
                "reference": "https://ieeexplore.ieee.org/document/9206060",
                "reference_text": "Koch et al. 'Stray Capacitance in Inductors.' IEEE Trans. Power Electronics, 2020"
            },
            "Massarini": {
                "description": """
Massarini's analytical model for winding capacitance includes:
- Multi-layer effects
- Wire geometry (round, rectangular)
- Insulation thickness

Provides separate formulas for different capacitance components.
""",
                "reference": "https://ieeexplore.ieee.org/document/668030",
                "reference_text": "Massarini, Kazimierczuk. 'Self-capacitance of inductors.' IEEE Trans. Power Electronics, 1997"
            },
            "Albach": {
                "description": """
Albach's comprehensive capacitance model accounts for:
- 2D field distribution
- Non-uniform winding spacing
- Core proximity effects

Provides good accuracy across various winding configurations.
""",
                "reference": "https://ieeexplore.ieee.org/document/5680877",
                "reference_text": "Albach et al. 'Calculating stray capacitance.' IEEE Trans. Magnetics, 2011"
            },
            "Duerdoth": {
                "description": """
Duerdoth's coil self-capacitance model uses energy-based approach:

$$C_{self} = \\frac{2 W_E}{V^2}$$

Where $W_E$ is the electric field energy stored between turns.
""",
                "reference": "https://ieeexplore.ieee.org/document/7866793",
                "reference_text": "Duerdoth, W.T. 'Equivalent Capacitance of Transformer Windings.' Wireless Engineer, 1946"
            }
        },
        "selection_guide": """
## Capacitance Minimization Strategies

| Technique | Effect | Trade-off |
|-----------|--------|-----------|
| Increase turn spacing | Reduces turn-to-turn C | Larger winding window |
| Section winding | Reduces layer-to-layer C | More complex construction |
| Shield windings | Controls coupling paths | Additional losses |
| Bank winding | Reduces distributed C | Higher leakage inductance |
| Thicker insulation | Reduces all C | Thermal resistance |

**Default Model:** Albach provides good balance of accuracy and generality.
"""
    },

    "thermal.md": {
        "intro": """
Thermal management is critical for magnetic component reliability. MKF provides models for:
- Core temperature rise from core losses
- Winding hot-spot temperature
- Thermal resistance estimation

Temperature affects:
- Core losses (Steinmetz parameters are temperature-dependent)
- Winding resistance (copper has +0.39%/°C temperature coefficient)
- Saturation flux density (decreases with temperature)
- Component lifetime (every 10°C reduces life by ~50%)
""",
        "models": {
            "Maniktala": {
                "description": """
Maniktala's empirical formula for temperature rise:

$$\\Delta T = \\left(\\frac{P_{total}}{A_{surface}}\\right)^{0.833}$$

Where:
- $P_{total}$ is total power loss (W)
- $A_{surface}$ is effective cooling surface area (cm²)
- Result is temperature rise in °C

**Simple and practical** for initial estimates and design iteration.

**Reference:** Maniktala, S. "Switching Power Supplies A-Z." Newnes, 2012.
""",
                "reference": "https://www.elsevier.com/books/switching-power-supplies-a-z/maniktala/978-0-12-386533-5",
                "reference_text": "Maniktala, S. 'Switching Power Supplies A-Z.' Newnes, 2012"
            },
            "Kazimierczuk": {
                "description": """
Kazimierczuk's thermal model separates core and winding contributions:

$$\\Delta T_{core} = R_{th,core} \\cdot P_{core}$$
$$\\Delta T_{winding} = R_{th,winding} \\cdot P_{winding}$$

Thermal resistances are calculated from geometry and material properties.
More accurate than empirical formulas but requires more parameters.

**Reference:** Kazimierczuk, M.K. "High-Frequency Magnetic Components." Wiley, 2014.
""",
                "reference": "https://www.wiley.com/en-us/High+Frequency+Magnetic+Components%2C+2nd+Edition-p-9781118717790",
                "reference_text": "Kazimierczuk, M.K. 'High-Frequency Magnetic Components.' Wiley, 2014"
            },
            "TDK": {
                "description": """
TDK's thermal models are empirically derived from extensive testing of their
ferrite cores. Provides manufacturer-specific accuracy for TDK materials.
""",
                "reference": "https://www.tdk-electronics.tdk.com/download/531372/4730039e8720de4e83e3e32614ca3efc/ferrites-magnetic-design-tool.pdf",
                "reference_text": "TDK Ferrites and Accessories Application Notes"
            },
            "Dixon": {
                "description": """
Dixon's thermal resistance approximation for pot cores and similar shapes:

$$R_{th} \\approx \\frac{50}{A_{surface}^{0.7}}$$ (°C/W)

Where $A_{surface}$ is in cm².
""",
                "reference": "https://www.ti.com/lit/ml/slup123/slup123.pdf",
                "reference_text": "Dixon, L. 'Magnetics Design for Switching Power Supplies.' Texas Instruments"
            },
            "Amidon": {
                "description": """
Amidon's thermal data is empirically derived for their iron powder and ferrite cores.
Useful when designing with Amidon (Micrometals) materials.
""",
                "reference": "https://www.amidoncorp.com/specs/",
                "reference_text": "Amidon Corp. Technical Specifications"
            }
        },
        "selection_guide": """
## Thermal Model Selection

| Application | Recommended | Notes |
|-------------|-------------|-------|
| Quick estimates | Maniktala | Simple, conservative |
| Detailed analysis | Kazimierczuk | Separates loss sources |
| TDK cores | TDK | Manufacturer data |
| Worst-case | Use lowest thermal resistance estimate | Safety margin |

**Thermal Resistance Models:**
The thermal resistance between core/winding and ambient depends on:
- Cooling method (natural convection, forced air, heat sink)
- Surface treatment (painted, bare, potted)
- Mounting orientation (vertical, horizontal)
- Ambient temperature and altitude

**Design Margin:** Add 10-15°C margin for component reliability.
"""
    }
}


def extract_map_entries(content: str, method_name: str) -> Dict[str, str]:
    """Extract key-value pairs from a static map-returning method."""
    # Find the method
    pattern = rf'static\s+std::map<std::string,\s*(?:std::string|double)>\s+{method_name}\s*\(\s*\)\s*\{{'
    match = re.search(pattern, content, re.MULTILINE)
    if not match:
        return {}

    # Find the method body
    start = match.end()
    brace_count = 1
    end = start
    while brace_count > 0 and end < len(content):
        if content[end] == '{':
            brace_count += 1
        elif content[end] == '}':
            brace_count -= 1
        end += 1

    method_body = content[start:end-1]

    # Extract entries: information["Key"] = R"(value)" or information["Key"] = value
    entries = {}

    # Pattern for R"(...)" strings
    raw_pattern = r'\w+\["([^"]+)"\]\s*=\s*R"\(([^)]*)\)"'
    for match in re.finditer(raw_pattern, method_body, re.DOTALL):
        key, value = match.groups()
        entries[key] = value.strip()

    # Pattern for regular strings
    str_pattern = r'\w+\["([^"]+)"\]\s*=\s*"([^"]*)"'
    for match in re.finditer(str_pattern, method_body):
        key, value = match.groups()
        if key not in entries:
            entries[key] = value

    # Pattern for numeric values
    num_pattern = r'\w+\["([^"]+)"\]\s*=\s*([\d.]+)'
    for match in re.finditer(num_pattern, method_body):
        key, value = match.groups()
        if key not in entries:
            entries[key] = value

    return entries


def extract_models_from_header(header_path: Path) -> Optional[ModelCategory]:
    """Extract model information from a header file."""
    content = header_path.read_text(encoding='utf-8', errors='ignore')

    # Check if this file has model information methods
    if 'get_models_information' not in content:
        return None

    # Determine the model category name from the file
    stem = header_path.stem
    category_name = stem.replace('_', ' ').title()

    # Find the main class name - using a non-backtracking regex pattern
    # Pattern explanation:
    #   class\s+(\w+Model) - matches 'class' followed by whitespace and the class name ending in 'Model'
    #   (?:\s*:\s*public\s+\w+)? - optionally matches inheritance (whitespace tolerant)
    #   \s*\{ - matches optional whitespace and opening brace
    # Using possessive-style quantifiers by avoiding nested optional groups that could backtrack
    class_pattern = r'class\s+(\w+Model)(?:\s*:\s*public\s+\w+)?\s*\{'
    class_match = re.search(class_pattern, content)
    class_name = class_match.group(1) if class_match else stem + "Model"

    # Extract information
    info_map = extract_map_entries(content, 'get_models_information')
    errors_map = extract_map_entries(content, 'get_models_errors')
    external_links_map = extract_map_entries(content, 'get_models_external_links')
    internal_links_map = extract_map_entries(content, 'get_models_internal_links')

    if not info_map:
        return None

    models = []
    for name, description in info_map.items():
        error = None
        if name in errors_map:
            try:
                error = float(errors_map[name])
            except ValueError:
                pass

        external_link = external_links_map.get(name)
        internal_link = internal_links_map.get(name)

        models.append(ModelInfo(
            name=name,
            description=description,
            error=error,
            external_link=external_link,
            internal_link=internal_link
        ))

    return ModelCategory(
        name=category_name,
        class_name=class_name,
        models=models,
        source_file=header_path.name
    )


def _generate_header(output_name: str, categories: List[ModelCategory],
                     detailed: dict) -> List[str]:
    """Generate the markdown header section."""
    lines = []

    # Determine main title from first category or detailed docs
    main_title = categories[0].name if categories else output_name.replace(
        ".md", "").replace("-", " ").title()

    lines.append(f"# {main_title}")
    lines.append("")

    # Add source attribution
    if categories:
        source_files = ", ".join(f"`{c.source_file}`" for c in categories)
        lines.append(f"*Auto-generated from {source_files} with comprehensive documentation*")
        lines.append("")

    # Add introduction from detailed docs
    if "intro" in detailed:
        lines.append(detailed["intro"].strip())
        lines.append("")

    return lines


def _generate_model_section(model: ModelInfo, detailed_models: dict) -> List[str]:
    """Generate markdown for a single model."""
    lines = []
    lines.append(f"### {model.name}")
    lines.append("")

    # Use detailed description if available, otherwise extracted
    if model.name in detailed_models:
        dm = detailed_models[model.name]
        lines.append(dm["description"].strip())
        lines.append("")

        if model.error is not None:
            lines.append(f"**Validation Error:** {model.error:.1%} mean deviation")
            lines.append("")

        # Use detailed reference if available
        if dm.get("reference"):
            lines.append(f"**Reference:** [{dm['reference_text']}]({dm['reference']})")
        elif model.external_link:
            lines.append(f"**Reference:** [{model.external_link}]({model.external_link})")
        lines.append("")
    else:
        # Fall back to extracted description
        lines.append(model.description)
        lines.append("")

        if model.error is not None:
            lines.append(f"**Validation Error:** {model.error:.1%} mean deviation")
            lines.append("")

        if model.external_link:
            lines.append(f"**Reference:** [{model.external_link}]({model.external_link})")
            lines.append("")

    return lines


def _generate_proximity_section(detailed: dict) -> List[str]:
    """Generate proximity effect section for winding losses."""
    lines = []
    lines.append(detailed["proximity_effect_intro"].strip())
    lines.append("")

    proximity_models = detailed.get("proximity_models", {})
    for name, pm in proximity_models.items():
        lines.append(f"### {name}")
        lines.append("")
        lines.append(pm["description"].strip())
        lines.append("")
        if pm.get("reference"):
            lines.append(f"**Reference:** [{pm['reference_text']}]({pm['reference']})")
            lines.append("")

    return lines


def _generate_category_content(cat: ModelCategory, detailed: dict,
                               show_category_header: bool) -> List[str]:
    """Generate content for a single category."""
    lines = []

    if show_category_header:
        lines.append(f"## {cat.name}")
        lines.append("")

    lines.append("## Available Models")
    lines.append("")

    # Process each model
    detailed_models = detailed.get("models", {})
    for model in cat.models:
        lines.extend(_generate_model_section(model, detailed_models))

    # Add skin effect intro for winding losses
    if "skin_effect_intro" in detailed:
        lines.append(detailed["skin_effect_intro"].strip())
        lines.append("")

    # Add proximity effect section for winding losses
    if "proximity_effect_intro" in detailed:
        lines.extend(_generate_proximity_section(detailed))

    # Add fringing section for magnetic field
    if "fringing_section" in detailed:
        lines.append(detailed["fringing_section"].strip())
        lines.append("")

    # Add leakage section for inductance
    if "leakage_section" in detailed:
        lines.append(detailed["leakage_section"].strip())
        lines.append("")

    return lines


def _generate_comparison_table(categories: List[ModelCategory]) -> List[str]:
    """Generate model comparison table."""
    lines = [
        "## Model Comparison",
        "",
        "| Model | Error | Reference |",
        "|-------|-------|-----------|",
    ]

    for cat in categories:
        for model in cat.models:
            error_str = f"{model.error:.1%}" if model.error is not None else "N/A"
            ref_str = f"[Link]({model.external_link})" if model.external_link else "N/A"
            lines.append(f"| {model.name} | {error_str} | {ref_str} |")

    lines.append("")
    return lines


def _generate_usage_section(categories: List[ModelCategory]) -> List[str]:
    """Generate usage code example section."""
    if not categories:
        return []

    cat = categories[0]
    lines = [
        "## Usage",
        "",
        CODE_BLOCK_CPP,
        f"#include \"physical_models/{cat.source_file}\"",
        "",
        "// Create a specific model",
        f"auto model = OpenMagnetics::{cat.class_name}::factory(",
    ]

    if cat.models:
        enum_name = cat.class_name.replace("Model", "Models")
        model_enum = cat.models[0].name.upper().replace(" ", "_").replace("-", "_")
        lines.append(f"    OpenMagnetics::{enum_name}::{model_enum}")

    setting_name = cat.name.lower().replace(" ", "_")
    lines.extend([
        ");",
        "",
        "// Or use the default model",
        f"auto model = OpenMagnetics::{cat.class_name}::factory();",
        "```",
        "",
        "## Configuring Default Model",
        "",
        CODE_BLOCK_CPP,
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        f"// settings.set_{setting_name}_model(OpenMagnetics::{cat.class_name.replace('Model', 'Models')}::...);",
        "```",
    ])

    return lines


def generate_enhanced_markdown(output_name: str, categories: List[ModelCategory]) -> str:
    """Generate enhanced Markdown using DETAILED_DOCS plus extracted data."""
    detailed = DETAILED_DOCS.get(output_name, {})
    lines = []

    # Generate header section
    lines.extend(_generate_header(output_name, categories, detailed))

    # Process each category
    show_category_header = len(categories) > 1
    for cat in categories:
        lines.extend(_generate_category_content(cat, detailed, show_category_header))

    # Add model comparison table
    lines.extend(_generate_comparison_table(categories))

    # Add selection guide
    if "selection_guide" in detailed:
        lines.append(detailed["selection_guide"].strip())
        lines.append("")

    # Add usage section
    lines.extend(_generate_usage_section(categories))

    # Add usage section
    if categories:
        cat = categories[0]
        lines.extend([
            "## Usage",
            "",
            CODE_BLOCK_CPP,
            f"#include \"physical_models/{cat.source_file}\"",
            "",
            "// Create a specific model",
            f"auto model = OpenMagnetics::{cat.class_name}::factory(",
        ])

        if cat.models:
            enum_name = cat.class_name.replace("Model", "Models")
            model_enum = cat.models[0].name.upper().replace(" ", "_").replace("-", "_")
            lines.append(f"    OpenMagnetics::{enum_name}::{model_enum}")

        lines.extend([
            ");",
            "",
            "// Or use the default model",
            f"auto model = OpenMagnetics::{cat.class_name}::factory();",
            "```",
            "",
            "## Configuring Default Model",
            "",
            CODE_BLOCK_CPP,
            "auto& settings = OpenMagnetics::Settings::GetInstance();",
        ])

        setting_name = cat.name.lower().replace(" ", "_")
        lines.extend([
            f"// settings.set_{setting_name}_model(OpenMagnetics::{cat.class_name.replace('Model', 'Models')}::...);",
            "```",
        ])

    return "\n".join(lines)


def generate_model_markdown(category: ModelCategory) -> str:
    """Generate basic Markdown documentation for a model category (fallback)."""
    lines = [
        f"# {category.name}",
        "",
        f"*Auto-generated from `src/physical_models/{category.source_file}`*",
        "",
        "## Overview",
        "",
        f"MKF provides multiple implementations for {category.name.lower()} calculations. "
        "Each model is based on published scientific research and has been validated against experimental data.",
        "",
        "## Available Models",
        "",
    ]

    for model in category.models:
        lines.append(f"### {model.name}")
        lines.append("")
        lines.append(model.description)
        lines.append("")

        if model.error is not None:
            lines.append(f"**Validation Error:** {model.error:.1%} mean deviation")
            lines.append("")

        if model.external_link:
            lines.append(f"**Reference:** [{model.external_link}]({model.external_link})")
            lines.append("")

        if model.internal_link and model.internal_link.strip():
            lines.append(f"**Related Article:** [{model.internal_link}]({model.internal_link})")
            lines.append("")

    # Add comparison table
    lines.extend([
        "## Model Comparison",
        "",
        "| Model | Error | Reference |",
        "|-------|-------|-----------|",
    ])

    for model in category.models:
        error_str = f"{model.error:.1%}" if model.error is not None else "N/A"
        ref_str = f"[Link]({model.external_link})" if model.external_link else "N/A"
        lines.append(f"| {model.name} | {error_str} | {ref_str} |")

    lines.extend([
        "",
        "## Usage",
        "",
        CODE_BLOCK_CPP,
        f"#include \"physical_models/{category.source_file}\"",
        "",
        "// Create a specific model",
        f"auto model = OpenMagnetics::{category.class_name}::factory(",
    ])

    # Add enum value based on first model
    if category.models:
        enum_name = category.class_name.replace("Model", "Models")
        model_enum = category.models[0].name.upper().replace(" ", "_").replace("-", "_")
        lines.append(f"    OpenMagnetics::{enum_name}::{model_enum}")

    lines.extend([
        ");",
        "",
        "// Or use the default model",
        f"auto model = OpenMagnetics::{category.class_name}::factory();",
        "```",
        "",
        "## Model Selection",
        "",
        "To set a default model globally:",
        "",
        CODE_BLOCK_CPP,
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
    ])

    # Add settings call if applicable
    setting_name = category.name.lower().replace(" ", "_")
    lines.extend([
        f"// settings.set_{setting_name}_model(OpenMagnetics::{category.class_name.replace('Model', 'Models')}::...);",
        "```",
    ])

    return "\n".join(lines)


def main():
    """Main entry point."""
    # Determine paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent  # docs/scripts -> docs -> project_root
    src_dir = project_root / "src" / "physical_models"
    docs_dir = project_root / "docs" / "models"

    if not src_dir.exists():
        print(f"Error: Source directory not found: {src_dir}")
        sys.exit(1)

    docs_dir.mkdir(parents=True, exist_ok=True)

    # File mapping for output names
    file_mapping = {
        "CoreLosses.h": "core-losses.md",
        "Reluctance.h": "reluctance.md",
        "WindingLosses.h": "winding-losses.md",
        "WindingSkinEffectLosses.h": "winding-losses.md",  # Append to winding-losses
        "WindingProximityEffectLosses.h": "winding-losses.md",  # Append to winding-losses
        "MagneticField.h": "magnetic-field.md",
        "CoreTemperature.h": "thermal.md",
        "ThermalResistance.h": "thermal.md",  # Append to thermal
        "StrayCapacitance.h": "capacitance.md",
        "Inductance.h": "inductance.md",
        "LeakageInductance.h": "inductance.md",  # Append to inductance
    }

    # Process each header file
    categories = {}
    for header_file in src_dir.glob("*.h"):
        if header_file.name not in file_mapping:
            continue

        category = extract_models_from_header(header_file)
        if category:
            output_name = file_mapping[header_file.name]
            if output_name not in categories:
                categories[output_name] = []
            categories[output_name].append(category)
            print(f"Extracted {len(category.models)} models from {header_file.name}")

    # Generate enhanced markdown files
    for output_name, cats in categories.items():
        output_path = docs_dir / output_name

        # Use enhanced generation with DETAILED_DOCS
        content = generate_enhanced_markdown(output_name, cats)

        output_path.write_text(content)
        print(f"Generated {output_path}")

    # Generate files from DETAILED_DOCS that weren't covered by extraction
    for output_name, detailed in DETAILED_DOCS.items():
        output_path = docs_dir / output_name
        if not output_path.exists():
            # Generate from detailed docs alone
            content = generate_enhanced_markdown(output_name, [])
            output_path.write_text(content)
            print(f"Generated {output_path} from detailed docs")

    print(f"\nModel documentation generated in {docs_dir}")


if __name__ == "__main__":
    main()
