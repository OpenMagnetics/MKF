#!/usr/bin/env python3
"""
Extract settings documentation from Settings.h and generate Markdown documentation.

This script parses the Settings class to extract:
- Setting names (from getter/setter methods)
- Default values (from member initialization)
- Types (bool, size_t, double, string, etc.)
- Descriptions (from comments if available)
"""

import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class Setting:
    """Information about a single setting."""
    name: str
    getter: str
    setter: str
    type: str
    default_value: Optional[str] = None
    description: Optional[str] = None
    category: str = "General"


def extract_member_default(content: str, member_name: str) -> Optional[str]:
    """Extract the default value for a member variable."""
    # Pattern: _memberName = value; or _memberName{value};
    pattern = rf'{member_name}\s*[={{]\s*([^;}}]+)[;}}]'
    match = re.search(pattern, content)
    if match:
        value = match.group(1).strip()
        # Clean up the value
        if value.startswith('"') and value.endswith('"'):
            return value  # Keep quotes for strings
        return value
    return None


def extract_type_from_getter(content: str, getter_name: str) -> str:
    """Extract the return type from a getter method."""
    pattern = rf'(\w+(?:<[^>]+>)?(?:\s+const)?)\s+{getter_name}\s*\('
    match = re.search(pattern, content)
    if match:
        return match.group(1).strip()
    return "unknown"


def categorize_setting(name: str) -> str:
    """Categorize a setting based on its name."""
    name_lower = name.lower()

    if 'painter' in name_lower or 'color' in name_lower:
        return "Visualization"
    elif 'coil' in name_lower or 'winding' in name_lower:
        return "Coil/Winding"
    elif 'core' in name_lower:
        return "Core"
    elif 'wire' in name_lower:
        return "Wire"
    elif 'adviser' in name_lower:
        return "Advisers"
    elif 'magnetic_field' in name_lower or 'leakage' in name_lower:
        return "Magnetic Field"
    elif 'harmonic' in name_lower or 'input' in name_lower:
        return "Input Processing"
    elif 'model' in name_lower or 'reluctance' in name_lower or 'losses' in name_lower:
        return "Physical Models"
    else:
        return "General"


def format_type(cpp_type: str) -> str:
    """Format C++ type for documentation."""
    type_map = {
        "bool": "`bool`",
        "size_t": "`size_t`",
        "double": "`double`",
        "int": "`int`",
        "std::string": "`string`",
        "PainterModes": "`PainterModes`",
        "std::optional<double>": "`double` (optional)",
        "std::optional<MagneticFieldStrengthModels>": "`MagneticFieldStrengthModels` (optional)",
        "std::vector<CoreLossesModels>": "`vector<CoreLossesModels>`",
        "MagneticFieldStrengthModels": "`MagneticFieldStrengthModels`",
        "MagneticFieldStrengthFringingEffectModels": "`MagneticFieldStrengthFringingEffectModels`",
        "ReluctanceModels": "`ReluctanceModels`",
        "CoreTemperatureModels": "`CoreTemperatureModels`",
        "CoreThermalResistanceModels": "`CoreThermalResistanceModels`",
        "WindingSkinEffectLossesModels": "`WindingSkinEffectLossesModels`",
        "WindingProximityEffectLossesModels": "`WindingProximityEffectLossesModels`",
        "StrayCapacitanceModels": "`StrayCapacitanceModels`",
    }
    return type_map.get(cpp_type, f"`{cpp_type}`")


def format_default(value: Optional[str], cpp_type: str) -> str:
    """Format default value for documentation."""
    if value is None:
        return "N/A"

    # Handle booleans
    if value == "true":
        return "`true`"
    elif value == "false":
        return "`false`"

    # Handle enums
    if "::" in value:
        return f"`{value.split('::')[-1]}`"

    # Handle strings
    if value.startswith('"'):
        return f"`{value}`"

    # Handle std::nullopt
    if value == "std::nullopt":
        return "`null`"

    return f"`{value}`"


def parse_settings_header(header_path: Path) -> List[Setting]:
    """Parse Settings.h and extract all settings."""
    content = header_path.read_text(encoding='utf-8', errors='ignore')

    settings = []

    # Find all getter methods: <type> get_<name>() const;
    getter_pattern = r'(\w+(?:<[^>]+>)?)\s+(get_\w+)\s*\(\s*\)\s*const\s*;'

    for match in re.finditer(getter_pattern, content):
        return_type = match.group(1)
        getter_name = match.group(2)

        # Derive setting name from getter
        setting_name = getter_name.replace("get_", "")

        # Find corresponding setter
        setter_name = f"set_{setting_name}"

        # Find member variable name
        member_name = f"_{setting_name}"

        # Extract default value
        default_value = extract_member_default(content, member_name)

        # Categorize
        category = categorize_setting(setting_name)

        settings.append(Setting(
            name=setting_name,
            getter=getter_name,
            setter=setter_name,
            type=return_type,
            default_value=default_value,
            category=category
        ))

    return settings


def generate_settings_markdown(settings: List[Setting]) -> str:
    """Generate Markdown documentation for settings."""
    lines = [
        "# Settings Reference",
        "",
        "*Auto-generated from `src/support/Settings.h`*",
        "",
        "## Overview",
        "",
        "MKF uses a singleton `Settings` class to configure global behavior. "
        "Access it via `Settings::GetInstance()`.",
        "",
        "```cpp",
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        "",
        "// Example: Enable verbose logging",
        "settings.set_verbose(true);",
        "",
        "// Example: Configure reluctance model",
        "settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::ZHANG);",
        "```",
        "",
    ]

    # Group settings by category
    categories: Dict[str, List[Setting]] = {}
    for setting in settings:
        if setting.category not in categories:
            categories[setting.category] = []
        categories[setting.category].append(setting)

    # Sort categories
    category_order = [
        "General",
        "Physical Models",
        "Core",
        "Coil/Winding",
        "Wire",
        "Advisers",
        "Magnetic Field",
        "Input Processing",
        "Visualization",
    ]

    sorted_categories = sorted(
        categories.keys(),
        key=lambda x: category_order.index(x) if x in category_order else 999
    )

    # Generate table of contents
    lines.append("## Categories")
    lines.append("")
    for category in sorted_categories:
        anchor = category.lower().replace("/", "").replace(" ", "-")
        lines.append(f"- [{category}](#{anchor})")
    lines.append("")

    # Generate each category
    for category in sorted_categories:
        anchor = category.lower().replace("/", "").replace(" ", "-")
        lines.append(f"## {category}")
        lines.append("")
        lines.append("| Setting | Type | Default | Description |")
        lines.append("|---------|------|---------|-------------|")

        for setting in sorted(categories[category], key=lambda s: s.name):
            name = f"`{setting.name}`"
            type_str = format_type(setting.type)
            default_str = format_default(setting.default_value, setting.type)
            desc = setting.description or get_description(setting.name)

            lines.append(f"| {name} | {type_str} | {default_str} | {desc} |")

        lines.append("")

    # Add usage examples
    lines.extend([
        "## Usage Examples",
        "",
        "### Configuring Core Adviser",
        "",
        "```cpp",
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        "",
        "// Include stacked core configurations",
        "settings.set_core_adviser_include_stacks(true);",
        "",
        "// Include distributed gap configurations",
        "settings.set_core_adviser_include_distributed_gaps(true);",
        "",
        "// Enable pruning for faster searches",
        "settings.set_core_adviser_enable_intermediate_pruning(true);",
        "settings.set_core_adviser_maximum_magnetics_after_filtering(500);",
        "```",
        "",
        "### Configuring Wire Adviser",
        "",
        "```cpp",
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        "",
        "// Enable litz wire",
        "settings.set_wire_adviser_include_litz(true);",
        "",
        "// Disable foil wire",
        "settings.set_wire_adviser_include_foil(false);",
        "",
        "// Allow rectangular in toroids",
        "settings.set_wire_adviser_allow_rectangular_in_toroidal_cores(true);",
        "```",
        "",
        "### Configuring Physical Models",
        "",
        "```cpp",
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        "",
        "// Set reluctance model",
        "settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::MUEHLETHALER);",
        "",
        "// Set magnetic field model",
        "settings.set_magnetic_field_strength_model(",
        "    OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON",
        ");",
        "",
        "// Set winding loss models",
        "settings.set_winding_skin_effect_losses_model(",
        "    OpenMagnetics::WindingSkinEffectLossesModels::ALBACH",
        ");",
        "```",
        "",
        "### Resetting to Defaults",
        "",
        "```cpp",
        "auto& settings = OpenMagnetics::Settings::GetInstance();",
        "settings.reset();  // Reset all settings to defaults",
        "```",
    ])

    return "\n".join(lines)


def get_description(name: str) -> str:
    """Get a human-readable description for a setting."""
    descriptions = {
        "verbose": "Enable verbose logging output",
        "use_toroidal_cores": "Include toroidal cores in searches",
        "use_concentric_cores": "Include concentric (E, U, etc.) cores in searches",
        "inputs_trim_harmonics": "Trim insignificant harmonics from input waveforms",
        "inputs_number_points_sampled_waveforms": "Number of points for waveform sampling",
        "magnetizing_inductance_include_air_inductance": "Include air-core inductance in calculations",
        "coil_allow_margin_tape": "Allow margin tape in coil designs",
        "coil_allow_insulated_wire": "Allow pre-insulated wire",
        "coil_fill_sections_with_margin_tape": "Fill empty sections with margin tape",
        "coil_wind_even_if_not_fit": "Attempt winding even if marginal fit",
        "coil_delimit_and_compact": "Delimit sections and compact windings",
        "coil_try_rewind": "Try rewinding if initial attempt fails",
        "coil_include_additional_coordinates": "Include detailed turn coordinates",
        "coil_equalize_margins": "Equalize margins on both sides",
        "coil_only_one_turn_per_layer_in_contiguous_rectangular": "Limit rectangular wire to one turn per layer",
        "coil_maximum_layers_planar": "Maximum layers for planar transformers",
        "use_only_cores_in_stock": "Only use cores marked as in stock",
        "painter_number_points_x": "Grid resolution (X) for visualization",
        "painter_number_points_y": "Grid resolution (Y) for visualization",
        "painter_mode": "Visualization mode (CONTOUR, etc.)",
        "painter_logarithmic_scale": "Use logarithmic scale in plots",
        "painter_include_fringing": "Include fringing field in visualizations",
        "painter_draw_spacer": "Draw spacer elements",
        "painter_simple_litz": "Use simplified litz wire rendering",
        "painter_advanced_litz": "Use advanced litz wire rendering",
        "magnetic_field_number_points_x": "Grid resolution (X) for field calculations",
        "magnetic_field_number_points_y": "Grid resolution (Y) for field calculations",
        "magnetic_field_include_fringing": "Include fringing in field calculations",
        "coil_mesher_inside_turns_factor": "Factor for turn mesh density",
        "leakage_inductance_grid_auto_scaling": "Auto-scale grid for leakage calculations",
        "leakage_inductance_grid_precision_level_planar": "Grid precision for planar inductors",
        "leakage_inductance_grid_precision_level_wound": "Grid precision for wound inductors",
        "coil_adviser_maximum_number_wires": "Max wire candidates per coil design",
        "core_adviser_include_stacks": "Include core stacking combinations",
        "core_adviser_include_distributed_gaps": "Include distributed gap configurations",
        "core_adviser_include_margin": "Include margin in core designs",
        "core_adviser_enable_intermediate_pruning": "Prune candidates between filter stages",
        "core_adviser_maximum_magnetics_after_filtering": "Max candidates after filtering",
        "wire_adviser_include_planar": "Include planar wire in searches",
        "wire_adviser_include_foil": "Include foil wire in searches",
        "wire_adviser_include_rectangular": "Include rectangular wire in searches",
        "wire_adviser_include_litz": "Include litz wire in searches",
        "wire_adviser_include_round": "Include round wire in searches",
        "wire_adviser_allow_rectangular_in_toroidal_cores": "Allow rectangular wire in toroids",
        "harmonic_amplitude_threshold_quick_mode": "Use quick mode for harmonic filtering",
        "harmonic_amplitude_threshold": "Threshold for harmonic filtering",
        "magnetic_field_strength_model": "Default magnetic field model",
        "magnetic_field_strength_fringing_effect_model": "Default fringing effect model",
        "leakage_inductance_magnetic_field_strength_model": "Model for leakage inductance",
        "reluctance_model": "Default reluctance model",
        "core_temperature_model": "Default core temperature model",
        "core_thermal_resistance_model": "Default thermal resistance model",
        "winding_skin_effect_losses_model": "Default skin effect model",
        "winding_proximity_effect_losses_model": "Default proximity effect model",
        "stray_capacitance_model": "Default capacitance model",
        "preferred_core_material_ferrite_manufacturer": "Preferred ferrite manufacturer",
        "preferred_core_material_powder_manufacturer": "Preferred powder core manufacturer",
        "core_cross_referencer_allow_different_core_material_type": "Allow cross-referencing different material types",
    }
    return descriptions.get(name, "")


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent.parent  # docs/scripts -> docs -> project_root
    settings_path = project_root / "src" / "support" / "Settings.h"
    output_path = project_root / "docs" / "settings" / "reference.md"

    if not settings_path.exists():
        print(f"Error: Settings.h not found: {settings_path}")
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Parse settings
    settings = parse_settings_header(settings_path)
    print(f"Extracted {len(settings)} settings")

    # Generate markdown
    content = generate_settings_markdown(settings)
    output_path.write_text(content)
    print(f"Generated {output_path}")


if __name__ == "__main__":
    main()
