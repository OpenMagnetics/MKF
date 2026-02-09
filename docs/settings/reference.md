# Settings Reference

*Auto-generated from `src/support/Settings.h`*

## Overview

MKF uses a singleton `Settings` class to configure global behavior. Access it via `Settings::GetInstance()`.

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();

// Example: Enable verbose logging
settings.set_verbose(true);

// Example: Configure reluctance model
settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::ZHANG);
```

## Categories

- [General](#general)
- [Physical Models](#physical-models)
- [Core](#core)
- [Coil/Winding](#coilwinding)
- [Wire](#wire)
- [Magnetic Field](#magnetic-field)
- [Input Processing](#input-processing)
- [Visualization](#visualization)

## General

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `magnetizing_inductance_include_air_inductance` | `bool` | N/A | Include air-core inductance in calculations |
| `verbose` | `bool` | `false` | Enable verbose logging output |

## Physical Models

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `reluctance_model` | `ReluctanceModels` | N/A | Default reluctance model |
| `stray_capacitance_model` | `StrayCapacitanceModels` | N/A | Default capacitance model |

## Core

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `core_adviser_enable_intermediate_pruning` | `bool` | N/A | Prune candidates between filter stages |
| `core_adviser_include_distributed_gaps` | `bool` | N/A | Include distributed gap configurations |
| `core_adviser_include_margin` | `bool` | N/A | Include margin in core designs |
| `core_adviser_include_stacks` | `bool` | N/A | Include core stacking combinations |
| `core_adviser_maximum_magnetics_after_filtering` | `size_t` | N/A | Max candidates after filtering |
| `core_cross_referencer_allow_different_core_material_type` | `bool` | N/A | Allow cross-referencing different material types |
| `core_losses_model_names` | `vector<CoreLossesModels>` | N/A |  |
| `core_temperature_model` | `CoreTemperatureModels` | N/A | Default core temperature model |
| `core_thermal_resistance_model` | `CoreThermalResistanceModels` | N/A | Default thermal resistance model |
| `preferred_core_material_ferrite_manufacturer` | `string` | N/A | Preferred ferrite manufacturer |
| `preferred_core_material_powder_manufacturer` | `string` | N/A | Preferred powder core manufacturer |
| `use_concentric_cores` | `bool` | N/A | Include concentric (E, U, etc.) cores in searches |
| `use_only_cores_in_stock` | `bool` | N/A | Only use cores marked as in stock |
| `use_toroidal_cores` | `bool` | N/A | Include toroidal cores in searches |
| `wire_adviser_allow_rectangular_in_toroidal_cores` | `bool` | N/A | Allow rectangular wire in toroids |

## Coil/Winding

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `coil_adviser_maximum_number_wires` | `size_t` | N/A | Max wire candidates per coil design |
| `coil_allow_insulated_wire` | `bool` | N/A | Allow pre-insulated wire |
| `coil_allow_margin_tape` | `bool` | N/A | Allow margin tape in coil designs |
| `coil_delimit_and_compact` | `bool` | N/A | Delimit sections and compact windings |
| `coil_equalize_margins` | `bool` | N/A | Equalize margins on both sides |
| `coil_fill_sections_with_margin_tape` | `bool` | N/A | Fill empty sections with margin tape |
| `coil_include_additional_coordinates` | `bool` | N/A | Include detailed turn coordinates |
| `coil_maximum_layers_planar` | `size_t` | N/A | Maximum layers for planar transformers |
| `coil_mesher_inside_turns_factor` | `double` | N/A | Factor for turn mesh density |
| `coil_only_one_turn_per_layer_in_contiguous_rectangular` | `bool` | N/A | Limit rectangular wire to one turn per layer |
| `coil_try_rewind` | `bool` | N/A | Try rewinding if initial attempt fails |
| `coil_wind_even_if_not_fit` | `bool` | N/A | Attempt winding even if marginal fit |
| `winding_proximity_effect_losses_model` | `WindingProximityEffectLossesModels` | N/A | Default proximity effect model |
| `winding_skin_effect_losses_model` | `WindingSkinEffectLossesModels` | N/A | Default skin effect model |

## Wire

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `wire_adviser_include_foil` | `bool` | N/A | Include foil wire in searches |
| `wire_adviser_include_litz` | `bool` | N/A | Include litz wire in searches |
| `wire_adviser_include_planar` | `bool` | N/A | Include planar wire in searches |
| `wire_adviser_include_rectangular` | `bool` | N/A | Include rectangular wire in searches |
| `wire_adviser_include_round` | `bool` | N/A | Include round wire in searches |

## Magnetic Field

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `leakage_inductance_grid_auto_scaling` | `bool` | N/A | Auto-scale grid for leakage calculations |
| `leakage_inductance_grid_precision_level_planar` | `double` | N/A | Grid precision for planar inductors |
| `leakage_inductance_grid_precision_level_wound` | `double` | N/A | Grid precision for wound inductors |
| `leakage_inductance_magnetic_field_strength_model` | `MagneticFieldStrengthModels` | N/A | Model for leakage inductance |
| `magnetic_field_include_fringing` | `bool` | N/A | Include fringing in field calculations |
| `magnetic_field_mirroring_dimension` | `int` | N/A |  |
| `magnetic_field_number_points_x` | `size_t` | N/A | Grid resolution (X) for field calculations |
| `magnetic_field_number_points_y` | `size_t` | N/A | Grid resolution (Y) for field calculations |
| `magnetic_field_strength_fringing_effect_model` | `MagneticFieldStrengthFringingEffectModels` | N/A | Default fringing effect model |
| `magnetic_field_strength_model` | `MagneticFieldStrengthModels` | N/A | Default magnetic field model |

## Input Processing

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `harmonic_amplitude_threshold` | `double` | N/A | Threshold for harmonic filtering |
| `harmonic_amplitude_threshold_quick_mode` | `bool` | N/A | Use quick mode for harmonic filtering |
| `inputs_number_points_sampled_waveforms` | `size_t` | N/A | Number of points for waveform sampling |
| `inputs_trim_harmonics` | `bool` | N/A | Trim insignificant harmonics from input waveforms |

## Visualization

| Setting | Type | Default | Description |
|---------|------|---------|-------------|
| `painter_advanced_litz` | `bool` | N/A | Use advanced litz wire rendering |
| `painter_cci_coordinates_path` | `string` | N/A |  |
| `painter_color_bobbin` | `string` | N/A |  |
| `painter_color_copper` | `string` | N/A |  |
| `painter_color_current_density` | `string` | N/A |  |
| `painter_color_enamel` | `string` | N/A |  |
| `painter_color_etfe` | `string` | N/A |  |
| `painter_color_fep` | `string` | N/A |  |
| `painter_color_ferrite` | `string` | N/A |  |
| `painter_color_fr4` | `string` | N/A |  |
| `painter_color_insulation` | `string` | N/A |  |
| `painter_color_lines` | `string` | N/A |  |
| `painter_color_magnetic_field_maximum` | `string` | N/A |  |
| `painter_color_magnetic_field_minimum` | `string` | N/A |  |
| `painter_color_margin` | `string` | N/A |  |
| `painter_color_pfa` | `string` | N/A |  |
| `painter_color_silk` | `string` | N/A |  |
| `painter_color_spacer` | `string` | N/A |  |
| `painter_color_tca` | `string` | N/A |  |
| `painter_color_text` | `string` | N/A |  |
| `painter_draw_spacer` | `bool` | N/A | Draw spacer elements |
| `painter_include_fringing` | `bool` | N/A | Include fringing field in visualizations |
| `painter_logarithmic_scale` | `bool` | N/A | Use logarithmic scale in plots |
| `painter_magnetic_field_strength_model` | `optional<MagneticFieldStrengthModels>` | N/A |  |
| `painter_maximum_value_colorbar` | `optional<double>` | N/A |  |
| `painter_minimum_value_colorbar` | `optional<double>` | N/A |  |
| `painter_mirroring_dimension` | `int` | N/A |  |
| `painter_mode` | `PainterModes` | N/A | Visualization mode (CONTOUR, etc.) |
| `painter_number_points_x` | `size_t` | N/A | Grid resolution (X) for visualization |
| `painter_number_points_y` | `size_t` | N/A | Grid resolution (Y) for visualization |
| `painter_simple_litz` | `bool` | N/A | Use simplified litz wire rendering |

## Usage Examples

### Configuring Core Adviser

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();

// Include stacked core configurations
settings.set_core_adviser_include_stacks(true);

// Include distributed gap configurations
settings.set_core_adviser_include_distributed_gaps(true);

// Enable pruning for faster searches
settings.set_core_adviser_enable_intermediate_pruning(true);
settings.set_core_adviser_maximum_magnetics_after_filtering(500);
```

### Configuring Wire Adviser

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();

// Enable litz wire
settings.set_wire_adviser_include_litz(true);

// Disable foil wire
settings.set_wire_adviser_include_foil(false);

// Allow rectangular in toroids
settings.set_wire_adviser_allow_rectangular_in_toroidal_cores(true);
```

### Configuring Physical Models

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();

// Set reluctance model
settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::MUEHLETHALER);

// Set magnetic field model
settings.set_magnetic_field_strength_model(
    OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON
);

// Set winding loss models
settings.set_winding_skin_effect_losses_model(
    OpenMagnetics::WindingSkinEffectLossesModels::ALBACH
);
```

### Resetting to Defaults

```cpp
auto& settings = OpenMagnetics::Settings::GetInstance();
settings.reset();  // Reset all settings to defaults
```