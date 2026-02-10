#!/usr/bin/env python3
"""
Generate standalone Substack-ready markdown articles from model documentation.

This script creates self-contained articles suitable for publication on Substack
or other blog platforms, including full citations and explanations.
"""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass
from datetime import datetime


@dataclass
class ModelArticle:
    """Article content for a model category."""
    title: str
    subtitle: str
    introduction: str
    models: List[dict]
    conclusion: str


def get_article_content() -> Dict[str, ModelArticle]:
    """Define article content for each model category."""
    return {
        "core-losses": ModelArticle(
            title="Understanding Core Loss Models in Magnetic Design",
            subtitle="A deep dive into Steinmetz, iGSE, and beyond",
            introduction="""
When designing magnetic components like transformers and inductors, accurately predicting core losses is crucial for thermal management and efficiency optimization. Core losses arise from hysteresis and eddy currents in the magnetic material, and their magnitude depends on the operating frequency, flux density, and material properties.

In this article, we'll explore the various core loss models implemented in OpenMagnetics, their theoretical foundations, and when to use each one.
""".strip(),
            models=[
                {
                    "name": "Steinmetz Equation",
                    "key": "Steinmetz",
                    "content": """
The Steinmetz equation, first proposed by Charles Proteus Steinmetz in 1892, is the foundation of most core loss calculations:

$$P_v = k \\cdot f^\\alpha \\cdot B^\\beta$$

Where:
- $P_v$ is the volumetric power loss (W/m³)
- $f$ is the frequency (Hz)
- $B$ is the peak magnetic flux density (T)
- $k$, $\\alpha$, $\\beta$ are material-specific coefficients

The beauty of Steinmetz lies in its simplicity - manufacturers provide the k, α, β parameters for their materials, making loss estimation straightforward. However, it assumes sinusoidal excitation, which limits its applicability to power electronics with non-sinusoidal waveforms.
"""
                },
                {
                    "name": "Improved Generalized Steinmetz Equation (iGSE)",
                    "key": "iGSE",
                    "content": """
The iGSE, developed by Charles Sullivan at Dartmouth, extends Steinmetz to handle arbitrary waveforms:

$$P_v = \\frac{1}{T} \\int_0^T k_i \\left| \\frac{dB}{dt} \\right|^\\alpha (\\Delta B)^{\\beta - \\alpha} dt$$

The key insight is that core losses depend on the rate of change of flux density (dB/dt), not just its peak value. This allows accurate predictions for:
- Triangular waveforms (typical in inductors)
- Trapezoidal waveforms (common in half-bridge converters)
- Arbitrary PWM waveforms

iGSE uses the same Steinmetz parameters (k, α, β) but applies them to instantaneous dB/dt, making it practical since no additional material characterization is needed.
"""
                },
                {
                    "name": "Roshen Model",
                    "key": "Roshen",
                    "content": """
Waseem Roshen's model takes a more physics-based approach, separating hysteresis and eddy current losses:

$$P_v = P_{hyst} + P_{eddy} + P_{excess}$$

The hysteresis component is derived from the material's B-H loop, while eddy current losses are calculated from the material's resistivity and geometry. This separation provides better accuracy when:
- Operating far from the Steinmetz characterization conditions
- Using materials with strong temperature dependence
- Designing with DC bias

However, Roshen requires more material parameters than Steinmetz-based methods.
"""
                }
            ],
            conclusion="""
## Choosing the Right Model

For most designs, **iGSE** offers the best balance of accuracy and practicality. It handles non-sinusoidal waveforms well and uses readily available Steinmetz parameters.

Use **Steinmetz** when you have purely sinusoidal excitation or need a quick estimate.

Consider **Roshen** when you need physics-based accuracy, especially with significant DC bias or extreme operating conditions.

OpenMagnetics automatically selects the most appropriate model based on available material data, but you can override this to compare predictions or match specific validation data.
"""
        ),
        "reluctance": ModelArticle(
            title="Air Gap Reluctance Models: From Classic to Zhang",
            subtitle="How fringing fields affect your inductance calculations",
            introduction="""
The air gap in a magnetic core is where energy is stored in an inductor. Accurately calculating the gap reluctance is essential for predicting inductance and ensuring your design meets specifications.

However, the magnetic field doesn't simply pass straight through the gap - it fringes outward, effectively increasing the gap's cross-sectional area. Different models account for this fringing in different ways.
""".strip(),
            models=[
                {
                    "name": "Classic Model",
                    "key": "Classic",
                    "content": """
The classic reluctance formula assumes a uniform magnetic field:

$$R_g = \\frac{l_g}{\\mu_0 \\cdot A_c}$$

Where:
- $R_g$ is the gap reluctance (H⁻¹)
- $l_g$ is the gap length (m)
- $A_c$ is the core cross-sectional area (m²)
- $\\mu_0$ is the permeability of free space

This model ignores fringing entirely, leading to overestimated reluctance and underestimated inductance. It's useful only for rough estimates or very small gaps.
"""
                },
                {
                    "name": "Zhang Model",
                    "key": "Zhang",
                    "content": """
Xinsheng Zhang's model (2020) treats the total reluctance as a parallel combination of internal and fringing components:

$$R_g = R_{in} \\parallel R_{fr}$$

The internal reluctance handles the direct flux path, while the fringing reluctance accounts for flux that spreads around the gap edges:

$$R_{fr} = \\frac{\\pi}{\\mu_0 \\cdot C \\cdot \\ln\\left(\\frac{2h + l_g}{l_g}\\right)}$$

Where $C$ is the perimeter of the core cross-section and $h$ is the winding window height. This model provides excellent accuracy for typical power magnetics applications and is the default in OpenMagnetics.
"""
                },
                {
                    "name": "Mühlethaler Model",
                    "key": "Muehlethaler",
                    "content": """
Jonas Mühlethaler's 3D approach divides the fringing field into distinct regions, each with its own reluctance contribution. This is particularly valuable for:

- Planar cores with complex geometries
- Very large gaps (relative to core dimensions)
- Cores with multiple gaps

The model requires more geometric information but provides superior accuracy for challenging designs.
"""
                }
            ],
            conclusion="""
## Practical Recommendations

For typical E-cores and similar shapes with moderate gaps, **Zhang** provides the best accuracy with minimal complexity. It's the recommended default.

Use **Mühlethaler** when designing with planar cores, unusual geometries, or very large gaps where standard models show significant error.

The **Classic** model is useful for quick estimates or when you want to understand the baseline reluctance before fringing effects.

Remember: fringing always reduces reluctance (increases inductance). If your measured inductance exceeds calculations, fringing is likely the cause.
"""
        ),
    }


def generate_article(key: str, article: ModelArticle) -> str:
    """Generate a standalone markdown article."""
    today = datetime.now().strftime("%Y-%m-%d")

    lines = [
        f"# {article.title}",
        "",
        f"*{article.subtitle}*",
        "",
        f"*Published: {today}*",
        "",
        "---",
        "",
        article.introduction,
        "",
    ]

    for model in article.models:
        lines.extend([
            f"## {model['name']}",
            "",
            model['content'].strip(),
            "",
        ])

    lines.extend([
        article.conclusion.strip(),
        "",
        "---",
        "",
        "## References",
        "",
    ])

    # Add references based on model keys
    references = {
        "Steinmetz": "Steinmetz, C.P. \"On the law of hysteresis.\" AIEE Transactions, 1892. [IEEE](https://ieeexplore.ieee.org/document/1457110)",
        "iGSE": "Li, J., Abdallah, T., Sullivan, C.R. \"Accurate Prediction of Ferrite Core Loss with Nonsinusoidal Waveforms Using Only Steinmetz Parameters.\" [Paper](http://inductor.thayerschool.org/papers/IGSE.pdf)",
        "Roshen": "Roshen, W. \"Ferrite Core Loss for Power Magnetic Components Design.\" IEEE Trans. Magnetics, 1991. [IEEE](https://ieeexplore.ieee.org/document/4052433)",
        "Zhang": "Zhang, X. et al. \"Improved Calculation Method for Inductance Value of the Air-Gap Inductor.\" CIYCEE, 2020. [IEEE](https://ieeexplore.ieee.org/document/9332553)",
        "Muehlethaler": "Mühlethaler, J. \"A Novel Approach for 3D Air Gap Reluctance Calculations.\" ECCE Asia, 2011. [Paper](https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf)",
        "Classic": "Standard magnetic circuit theory. [Wikipedia](https://en.wikipedia.org/wiki/Magnetic_reluctance)",
    }

    for model in article.models:
        model_key = model.get('key', model['name'])
        if model_key in references:
            lines.append(f"- {references[model_key]}")

    lines.extend([
        "",
        "---",
        "",
        "*This article was generated from the [OpenMagnetics](https://github.com/OpenMagnetics/MKF) documentation. "
        "OpenMagnetics is an open-source framework for magnetic component design.*",
    ])

    return "\n".join(lines)


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent  # docs/scripts -> docs -> project_root
    output_dir = project_root / "docs" / "articles"

    output_dir.mkdir(parents=True, exist_ok=True)

    articles = get_article_content()

    for key, article in articles.items():
        content = generate_article(key, article)
        output_path = output_dir / f"{key}.md"
        output_path.write_text(content)
        print(f"Generated {output_path}")

    # Create index
    index_lines = [
        "# Articles",
        "",
        "Standalone articles suitable for Substack or blog publication.",
        "",
        "## Available Articles",
        "",
    ]

    for key, article in articles.items():
        index_lines.append(f"- [{article.title}]({key}.md)")

    index_lines.extend([
        "",
        "## Generating Articles",
        "",
        "Run the generation script:",
        "",
        "```bash",
        "python scripts/generate_substack.py",
        "```",
        "",
        "Articles are generated from model documentation with added context and references.",
    ])

    index_path = output_dir / "index.md"
    index_path.write_text("\n".join(index_lines))
    print(f"Generated {index_path}")


if __name__ == "__main__":
    main()
