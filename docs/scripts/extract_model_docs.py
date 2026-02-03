#!/usr/bin/env python3
"""
Extract model documentation from C++ source files and generate Markdown documentation.

This script parses the physical model header files to extract:
- Model descriptions from get_models_information()
- Model errors from get_models_errors()
- External links from get_models_external_links()
- Internal links from get_models_internal_links()
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

    # Find the main class name
    class_pattern = r'class\s+(\w+Model)\s*(?::\s*public\s+\w+)?\s*\{'
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


def generate_model_markdown(category: ModelCategory) -> str:
    """Generate Markdown documentation for a model category."""
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
        "```cpp",
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
        "```cpp",
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

    # Generate markdown files
    for output_name, cats in categories.items():
        output_path = docs_dir / output_name

        if len(cats) == 1:
            content = generate_model_markdown(cats[0])
        else:
            # Combine multiple categories into one file
            parts = []
            for cat in cats:
                parts.append(generate_model_markdown(cat))
            content = "\n\n---\n\n".join(parts)

        output_path.write_text(content)
        print(f"Generated {output_path}")

    # Create placeholder files for models without extractable info
    placeholders = {
        "inductance.md": "Inductance",
        "winding-losses.md": "Winding Losses",
        "thermal.md": "Thermal Models",
        "capacitance.md": "Stray Capacitance",
        "magnetic-field.md": "Magnetic Field",
    }

    for filename, title in placeholders.items():
        output_path = docs_dir / filename
        if not output_path.exists():
            content = f"# {title}\n\n*Documentation coming soon. See source files in `src/physical_models/` for details.*\n"
            output_path.write_text(content)
            print(f"Created placeholder {output_path}")

    print(f"\nModel documentation generated in {docs_dir}")


if __name__ == "__main__":
    main()
