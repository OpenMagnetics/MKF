# Examples

This page provides practical examples for common magnetic component design tasks, organized from basic to advanced.

## Table of Contents

- [Creating a Core](#creating-a-core)
- [Winding a Coil](#winding-a-coil)
- [Calculating Core Losses](#calculating-core-losses)
- [Calculating Winding Losses](#calculating-winding-losses)
- [Using the Design Adviser](#using-the-design-adviser)
- [Designing a Buck Inductor](#designing-a-buck-inductor)
- [Designing a Flyback Transformer](#designing-a-flyback-transformer)
- [Comparing Models](#comparing-models)
- [Working with Waveforms](#working-with-waveforms)
- [Configuration Settings](#configuration-settings)
- [Exporting to JSON](#exporting-to-json)
- [Python Examples (PyMKF)](#python-examples-pymkf)

---

## Creating a Core

Create a core from scratch by specifying shape, material, and gapping.

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Find a core shape by name
    auto shape = OpenMagnetics::find_core_shape_by_name("E 42/21/15");

    // Find a core material by name
    auto material = OpenMagnetics::find_core_material_by_name("3C95");

    // Create core functional description
    MAS::CoreFunctionalDescription functionalDescription;
    functionalDescription.set_shape(shape);
    functionalDescription.set_material(material);
    functionalDescription.set_number_stacks(1);

    // Add gapping (1mm subtractive gap)
    std::vector<MAS::CoreGap> gapping;
    MAS::CoreGap gap;
    gap.set_type(MAS::GapType::SUBTRACTIVE);
    gap.set_length(0.001);  // 1 mm
    gapping.push_back(gap);
    functionalDescription.set_gapping(gapping);

    // Create and process the core
    OpenMagnetics::Core core(functionalDescription);
    core.process_data();
    core.process_gap();

    // Access calculated parameters
    auto effectiveParams = core.get_processed_description()->get_effective_parameters();
    std::cout << "Core: " << core.get_name().value() << std::endl;
    std::cout << "Effective area: " << effectiveParams.get_effective_area() << " m²" << std::endl;
    std::cout << "Effective length: " << effectiveParams.get_effective_length() << " m" << std::endl;
    std::cout << "Effective volume: " << effectiveParams.get_effective_volume() << " m³" << std::endl;

    return 0;
}
```

### Loading a Core from Database

```cpp
// Simpler: Load a pre-defined core directly
auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

// List available cores
auto& db = OpenMagnetics::DatabaseManager::getInstance();
auto coreNames = db.get_core_names();
std::cout << "Available cores: " << coreNames.size() << std::endl;

// Filter by shape family
auto eCores = db.get_cores_by_shape_family("E");
auto toroids = db.get_cores_by_shape_family("T");
```

---

## Winding a Coil

Wind a coil on a core with specified turns and wire.

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Load or create a core
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Define winding functional description
    std::vector<MAS::WindingFunctionalDescription> windingDescriptions;

    // Primary winding: 50 turns of round wire
    MAS::WindingFunctionalDescription primary;
    primary.set_name("Primary");
    primary.set_number_turns(50);
    primary.set_number_parallels(1);
    primary.set_wire("Round 0.5 - Grade 1");
    windingDescriptions.push_back(primary);

    // Secondary winding: 10 turns with 3 parallels
    MAS::WindingFunctionalDescription secondary;
    secondary.set_name("Secondary");
    secondary.set_number_turns(10);
    secondary.set_number_parallels(3);
    secondary.set_wire("Round 1.0 - Grade 1");
    windingDescriptions.push_back(secondary);

    // Create coil functional description
    MAS::CoilFunctionalDescription coilFunctional;
    coilFunctional.set_bobbin("E 42/21/15");
    coilFunctional.set_windings_description(windingDescriptions);

    // Wind the coil
    OpenMagnetics::CoilWrapper coilWrapper;
    auto coil = coilWrapper.wind(coilFunctional, core);

    // Check result
    if (coil.get_turns_description()) {
        std::cout << "Winding successful!" << std::endl;
        auto layers = coil.get_layers_description();
        std::cout << "Number of layers: " << layers.size() << std::endl;
    }

    return 0;
}
```

### Winding with Interleaving

```cpp
// Define interleaving pattern: Primary-Secondary-Primary
std::vector<size_t> pattern = {0, 1, 0};  // Winding indices

MAS::CoilFunctionalDescription coilFunctional;
coilFunctional.set_bobbin("E 42/21/15");
coilFunctional.set_windings_description(windingDescriptions);
coilFunctional.set_sections_description(pattern);

// This will create P-S-P interleaving for reduced leakage inductance
```

---

## Calculating Core Losses

Calculate power dissipation in the magnetic core.

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Load a core
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Define operating point excitation
    MAS::OperatingPointExcitation excitation;
    excitation.set_frequency(100000);  // 100 kHz

    // Set magnetic flux density (sinusoidal)
    MAS::SignalDescriptor magneticFluxDensity;
    MAS::Processed processed;
    processed.set_peak(0.1);       // 100 mT peak
    processed.set_peak_to_peak(0.2);  // 200 mT peak-to-peak
    processed.set_offset(0.0);     // No DC bias
    magneticFluxDensity.set_processed(processed);
    excitation.set_magnetic_flux_density(magneticFluxDensity);

    // Calculate losses using iGSE model (default)
    OpenMagnetics::CoreLosses coreLossesCalculator;
    auto result = coreLossesCalculator.calculate_core_losses(core, excitation, 25.0);

    std::cout << "Operating conditions:" << std::endl;
    std::cout << "  Frequency: 100 kHz" << std::endl;
    std::cout << "  Flux density: 100 mT peak" << std::endl;
    std::cout << "  Temperature: 25°C" << std::endl;
    std::cout << std::endl;
    std::cout << "Results:" << std::endl;
    std::cout << "  Total core losses: " << result.get_core_losses() << " W" << std::endl;
    std::cout << "  Volumetric losses: " << result.get_volumetric_losses().value() << " W/m³" << std::endl;

    return 0;
}
```

### Using Different Core Loss Models

```cpp
// Steinmetz (sinusoidal only)
auto steinmetzModel = OpenMagnetics::CoreLossesModel::factory(
    OpenMagnetics::CoreLossesModels::STEINMETZ
);

// iGSE (arbitrary waveforms)
auto igseModel = OpenMagnetics::CoreLossesModel::factory(
    OpenMagnetics::CoreLossesModels::IGSE
);

// Roshen (physics-based)
auto roshenModel = OpenMagnetics::CoreLossesModel::factory(
    OpenMagnetics::CoreLossesModels::ROSHEN
);

// Calculate with specific model
auto losses = igseModel->get_core_losses(core, excitation, temperature);
```

---

## Calculating Winding Losses

Calculate DC, skin effect, and proximity effect losses.

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Create a complete magnetic assembly
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // ... (set up coil as shown above)

    MAS::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    // Define operating point
    MAS::OperatingPoint operatingPoint;
    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(25);
    operatingPoint.set_conditions(conditions);

    // Define current waveform
    MAS::OperatingPointExcitation excitation;
    excitation.set_frequency(100000);

    MAS::SignalDescriptor current;
    MAS::Processed processed;
    processed.set_peak(5.0);
    processed.set_rms(3.5);
    current.set_processed(processed);
    excitation.set_current(current);

    operatingPoint.set_excitations_per_winding({excitation});

    // Calculate winding losses
    OpenMagnetics::WindingLosses windingLossesCalculator;
    auto result = windingLossesCalculator.calculate_losses(magnetic, operatingPoint, 25.0);

    // Print breakdown
    std::cout << "Winding Losses Breakdown:" << std::endl;
    for (size_t i = 0; i < result.get_winding_losses_per_winding().size(); ++i) {
        const auto& winding = result.get_winding_losses_per_winding()[i];
        std::cout << "  Winding " << i + 1 << ":" << std::endl;
        std::cout << "    DC (ohmic): " << winding.get_ohmic_losses() << " W" << std::endl;
        if (winding.get_skin_effect_losses()) {
            std::cout << "    Skin effect: " << winding.get_skin_effect_losses().value() << " W" << std::endl;
        }
        if (winding.get_proximity_effect_losses()) {
            std::cout << "    Proximity effect: " << winding.get_proximity_effect_losses().value() << " W" << std::endl;
        }
    }

    return 0;
}
```

---

## Using the Design Adviser

Let the adviser find optimal core, wire, and winding configurations.

```cpp
#include "OpenMagnetics.h"
#include <iostream>

int main() {
    // Define design requirements
    MAS::DesignRequirements requirements;

    // Inductance requirement
    MAS::DimensionWithTolerance inductance;
    inductance.set_minimum(95e-6);   // 95 µH minimum
    inductance.set_nominal(100e-6);  // 100 µH nominal
    inductance.set_maximum(110e-6);  // 110 µH maximum
    requirements.set_magnetizing_inductance(inductance);

    // For transformer: set turns ratios
    std::vector<MAS::DimensionWithTolerance> turnsRatios;
    MAS::DimensionWithTolerance ratio;
    ratio.set_nominal(5.0);  // 5:1 turns ratio
    turnsRatios.push_back(ratio);
    requirements.set_turns_ratios(turnsRatios);

    // Define operating point
    MAS::OperatingPoint operatingPoint;

    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(40);
    operatingPoint.set_conditions(conditions);

    // Primary excitation with waveforms
    MAS::OperatingPointExcitation primaryExcitation;
    primaryExcitation.set_frequency(100000);

    // Triangular current waveform
    MAS::SignalDescriptor current;
    MAS::Waveform currentWaveform;
    currentWaveform.set_data({0, 2.0, 0, -2.0, 0});
    currentWaveform.set_time({0, 2.5e-6, 5e-6, 7.5e-6, 10e-6});
    current.set_waveform(currentWaveform);
    primaryExcitation.set_current(current);

    // Square voltage waveform
    MAS::SignalDescriptor voltage;
    MAS::Waveform voltageWaveform;
    voltageWaveform.set_data({100, 100, -100, -100, 100});
    voltageWaveform.set_time({0, 5e-6, 5e-6, 10e-6, 10e-6});
    voltage.set_waveform(voltageWaveform);
    primaryExcitation.set_voltage(voltage);

    operatingPoint.set_excitations_per_winding({primaryExcitation});

    // Run the magnetic adviser
    OpenMagnetics::MagneticAdviser adviser;
    auto recommendations = adviser.get_advised_magnetic(requirements, operatingPoint, 10);

    // Display results
    std::cout << "Top " << recommendations.size() << " magnetic designs:" << std::endl;
    for (size_t i = 0; i < recommendations.size(); ++i) {
        auto& rec = recommendations[i];
        std::cout << "\n" << i + 1 << ". ";
        std::cout << rec.get_core().get_name().value_or("Unknown") << " - ";
        std::cout << rec.get_core().get_material_name() << std::endl;
        std::cout << "   Efficiency: " << rec.get_outputs().get_efficiency().value_or(0) * 100 << "%" << std::endl;
    }

    return 0;
}
```

---

## Designing a Buck Inductor

Complete workflow for a DC-DC buck converter inductor.

```cpp
#include "OpenMagnetics.h"
#include <iostream>
#include <cmath>

int main() {
    // Buck converter specifications
    double Vin = 24.0;           // Input voltage (V)
    double Vout = 5.0;           // Output voltage (V)
    double Iout = 3.0;           // Output current (A)
    double fsw = 200000;         // Switching frequency (Hz)
    double rippleRatio = 0.3;    // Current ripple ratio (30%)

    // Calculate required inductance
    double duty = Vout / Vin;
    double deltaI = Iout * rippleRatio;
    double L = (Vin - Vout) * duty / (fsw * deltaI);

    std::cout << "Buck Converter Design:" << std::endl;
    std::cout << "  Input: " << Vin << " V" << std::endl;
    std::cout << "  Output: " << Vout << " V @ " << Iout << " A" << std::endl;
    std::cout << "  Frequency: " << fsw / 1000 << " kHz" << std::endl;
    std::cout << "  Required inductance: " << L * 1e6 << " µH" << std::endl;
    std::cout << "  Peak current: " << Iout + deltaI / 2 << " A" << std::endl;
    std::cout << std::endl;

    // Set up design requirements
    MAS::DesignRequirements requirements;
    MAS::DimensionWithTolerance inductance;
    inductance.set_nominal(L);
    inductance.set_minimum(L * 0.9);
    requirements.set_magnetizing_inductance(inductance);

    // Operating point with triangular current
    MAS::OperatingPoint operatingPoint;
    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(40);
    operatingPoint.set_conditions(conditions);

    MAS::OperatingPointExcitation excitation;
    excitation.set_frequency(fsw);

    // Triangular current: DC + ripple
    double Ipeak = Iout + deltaI / 2;
    double Ivalley = Iout - deltaI / 2;
    double Ton = duty / fsw;
    double Toff = (1 - duty) / fsw;

    MAS::SignalDescriptor current;
    MAS::Waveform waveform;
    waveform.set_data({Ivalley, Ipeak, Ivalley});
    waveform.set_time({0, Ton, Ton + Toff});
    current.set_waveform(waveform);
    excitation.set_current(current);

    operatingPoint.set_excitations_per_winding({excitation});

    // Get recommendations
    OpenMagnetics::MagneticAdviser adviser;
    auto designs = adviser.get_advised_magnetic(requirements, operatingPoint, 5);

    std::cout << "Recommended inductors:" << std::endl;
    for (size_t i = 0; i < designs.size(); ++i) {
        auto& design = designs[i];
        std::cout << i + 1 << ". " << design.get_core().get_name().value_or("Unknown");
        std::cout << " (" << design.get_core().get_material_name() << ")" << std::endl;
    }

    return 0;
}
```

---

## Designing a Flyback Transformer

Complete flyback transformer design workflow.

```cpp
#include "OpenMagnetics.h"
#include <iostream>
#include <cmath>

int main() {
    // Flyback specifications
    double Vin_min = 90;         // Minimum input voltage (V)
    double Vin_max = 375;        // Maximum input voltage (V)
    double Vout = 12.0;          // Output voltage (V)
    double Iout = 2.0;           // Output current (A)
    double fsw = 100000;         // Switching frequency (Hz)
    double Dmax = 0.45;          // Maximum duty cycle
    double efficiency = 0.85;    // Target efficiency
    double Vdiode = 0.5;         // Output diode drop (V)

    // Calculate turns ratio and inductance
    double Vout_total = Vout + Vdiode;
    double n = (Vin_min * Dmax) / (Vout_total * (1 - Dmax));  // Turns ratio Np/Ns

    double Pout = Vout * Iout;
    double Pin = Pout / efficiency;

    // Primary inductance (for CCM at minimum input)
    double Ip_peak = 2 * Pin / (Vin_min * Dmax);
    double Lp = (Vin_min * Dmax) / (fsw * Ip_peak);

    std::cout << "Flyback Transformer Design:" << std::endl;
    std::cout << "  Input range: " << Vin_min << "-" << Vin_max << " V" << std::endl;
    std::cout << "  Output: " << Vout << " V @ " << Iout << " A" << std::endl;
    std::cout << "  Frequency: " << fsw / 1000 << " kHz" << std::endl;
    std::cout << "  Turns ratio (Np/Ns): " << n << ":1" << std::endl;
    std::cout << "  Primary inductance: " << Lp * 1e6 << " µH" << std::endl;
    std::cout << "  Peak primary current: " << Ip_peak << " A" << std::endl;
    std::cout << std::endl;

    // Design requirements
    MAS::DesignRequirements requirements;

    MAS::DimensionWithTolerance inductance;
    inductance.set_nominal(Lp);
    inductance.set_minimum(Lp * 0.9);
    requirements.set_magnetizing_inductance(inductance);

    std::vector<MAS::DimensionWithTolerance> turnsRatios;
    MAS::DimensionWithTolerance ratio;
    ratio.set_nominal(n);
    turnsRatios.push_back(ratio);
    requirements.set_turns_ratios(turnsRatios);

    requirements.set_isolation_voltage(3000);  // 3 kV isolation

    // Operating point at worst case (min input)
    MAS::OperatingPoint operatingPoint;
    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(40);
    operatingPoint.set_conditions(conditions);

    // Primary excitation
    MAS::OperatingPointExcitation primaryExcitation;
    primaryExcitation.set_frequency(fsw);

    // Triangular primary current
    MAS::SignalDescriptor primaryCurrent;
    MAS::Waveform primaryWaveform;
    double Ton = Dmax / fsw;
    double Tperiod = 1.0 / fsw;
    primaryWaveform.set_data({0, Ip_peak, 0, 0});
    primaryWaveform.set_time({0, Ton, Ton, Tperiod});
    primaryCurrent.set_waveform(primaryWaveform);
    primaryExcitation.set_current(primaryCurrent);

    // Secondary excitation
    MAS::OperatingPointExcitation secondaryExcitation;
    secondaryExcitation.set_frequency(fsw);

    double Is_peak = Ip_peak * n;
    MAS::SignalDescriptor secondaryCurrent;
    MAS::Waveform secondaryWaveform;
    secondaryWaveform.set_data({0, 0, Is_peak, 0});
    secondaryWaveform.set_time({0, Ton, Ton, Tperiod});
    secondaryCurrent.set_waveform(secondaryWaveform);
    secondaryExcitation.set_current(secondaryCurrent);

    operatingPoint.set_excitations_per_winding({primaryExcitation, secondaryExcitation});

    // Get recommendations
    OpenMagnetics::MagneticAdviser adviser;
    auto designs = adviser.get_advised_magnetic(requirements, operatingPoint, 5);

    std::cout << "Recommended flyback transformers:" << std::endl;
    for (size_t i = 0; i < designs.size(); ++i) {
        auto& design = designs[i];
        std::cout << i + 1 << ". " << design.get_core().get_name().value_or("Unknown");
        std::cout << " (" << design.get_core().get_material_name() << ")" << std::endl;
    }

    return 0;
}
```

---

## Comparing Models

Compare different physical models for the same calculation.

### Comparing Reluctance Models

```cpp
#include "OpenMagnetics.h"
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Add 1mm gap
    auto gapping = core.get_functional_description().get_gapping();
    MAS::CoreGap gap;
    gap.set_length(0.001);
    gap.set_type(MAS::GapType::SUBTRACTIVE);
    gapping.push_back(gap);
    core.get_functional_description().set_gapping(gapping);
    core.process_gap();

    std::vector<std::pair<OpenMagnetics::ReluctanceModels, std::string>> models = {
        {OpenMagnetics::ReluctanceModels::ZHANG, "Zhang"},
        {OpenMagnetics::ReluctanceModels::MUEHLETHALER, "Muehlethaler"},
        {OpenMagnetics::ReluctanceModels::PARTRIDGE, "Partridge"},
        {OpenMagnetics::ReluctanceModels::STENGLEIN, "Stenglein"},
        {OpenMagnetics::ReluctanceModels::BALAKRISHNAN, "Balakrishnan"},
        {OpenMagnetics::ReluctanceModels::CLASSIC, "Classic"}
    };

    // Get model error information
    auto errors = OpenMagnetics::ReluctanceModel::get_models_errors();

    std::cout << std::left << std::setw(15) << "Model"
              << std::setw(20) << "Reluctance (H⁻¹)"
              << std::setw(15) << "Fringing"
              << std::setw(15) << "Error" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    for (const auto& [modelType, name] : models) {
        auto model = OpenMagnetics::ReluctanceModel::factory(modelType);
        auto result = model->get_core_reluctance(core);

        std::cout << std::left << std::setw(15) << name
                  << std::setw(20) << std::scientific << result.get_core_reluctance()
                  << std::setw(15) << std::fixed << std::setprecision(2)
                  << result.get_maximum_fringing_factor()
                  << std::setw(15) << std::fixed << std::setprecision(1)
                  << errors[name] * 100 << "%" << std::endl;
    }

    return 0;
}
```

### Comparing Core Loss Models

```cpp
std::vector<std::pair<OpenMagnetics::CoreLossesModels, std::string>> models = {
    {OpenMagnetics::CoreLossesModels::STEINMETZ, "Steinmetz"},
    {OpenMagnetics::CoreLossesModels::IGSE, "iGSE"},
    {OpenMagnetics::CoreLossesModels::MSE, "MSE"},
    {OpenMagnetics::CoreLossesModels::ALBACH, "Albach"}
};

for (const auto& [modelType, name] : models) {
    auto model = OpenMagnetics::CoreLossesModel::factory(modelType);
    auto losses = model->get_core_losses(core, excitation, 25.0);
    std::cout << name << ": " << losses.get_core_losses() << " W" << std::endl;
}
```

---

## Working with Waveforms

Define custom excitation waveforms.

```cpp
#include "OpenMagnetics.h"

// Method 1: Using processed values (for simple waveforms)
MAS::SignalDescriptor signal1;
MAS::Processed processed;
processed.set_peak(0.1);
processed.set_peak_to_peak(0.2);
processed.set_rms(0.0707);
processed.set_offset(0.0);
signal1.set_processed(processed);

// Method 2: Using waveform data points
MAS::SignalDescriptor signal2;
MAS::Waveform waveform;
waveform.set_data({0, 0.1, 0, -0.1, 0});  // Values
waveform.set_time({0, 2.5e-6, 5e-6, 7.5e-6, 10e-6});  // Time points
signal2.set_waveform(waveform);

// Method 3: Using harmonics
MAS::SignalDescriptor signal3;
MAS::Harmonics harmonics;
std::vector<MAS::Harmonic> harmonicList;

MAS::Harmonic fundamental;
fundamental.set_frequency(100000);
fundamental.set_amplitude(0.1);
fundamental.set_phase(0);
harmonicList.push_back(fundamental);

MAS::Harmonic third;
third.set_frequency(300000);
third.set_amplitude(0.01);
third.set_phase(0);
harmonicList.push_back(third);

harmonics.set_amplitudes(harmonicList);
signal3.set_harmonics(harmonics);
```

---

## Configuration Settings

Customize library behavior through settings.

```cpp
#include "OpenMagnetics.h"

int main() {
    auto& settings = OpenMagnetics::Settings::GetInstance();

    // === Core Adviser Settings ===
    settings.set_core_adviser_include_stacks(true);       // Allow stacked cores
    settings.set_core_adviser_include_distributed_gaps(true);
    settings.set_core_adviser_enable_intermediate_pruning(true);
    settings.set_core_adviser_maximum_magnetics_after_filtering(500);

    // === Wire Adviser Settings ===
    settings.set_wire_adviser_include_litz(true);
    settings.set_wire_adviser_include_round(true);
    settings.set_wire_adviser_include_rectangular(true);
    settings.set_wire_adviser_include_foil(false);

    // === Coil Settings ===
    settings.set_coil_allow_margin_tape(true);
    settings.set_coil_allow_insulated_wire(true);
    settings.set_coil_wind_even_if_not_fit(false);
    settings.set_coil_try_rewind(true);

    // === Physical Model Defaults ===
    settings.set_reluctance_model(OpenMagnetics::ReluctanceModels::ZHANG);
    settings.set_magnetic_field_strength_model(
        OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON
    );
    settings.set_winding_skin_effect_losses_model(
        OpenMagnetics::WindingSkinEffectLossesModels::ALBACH
    );

    // === Database Preferences ===
    settings.set_use_only_cores_in_stock(false);
    settings.set_use_toroidal_cores(true);
    settings.set_use_concentric_cores(true);
    settings.set_preferred_core_material_ferrite_manufacturer("TDK");

    // === Visualization Settings ===
    settings.set_painter_number_points_x(50);
    settings.set_painter_number_points_y(100);
    settings.set_painter_include_fringing(true);

    // === Reset to defaults ===
    settings.reset();

    return 0;
}
```

---

## Exporting to JSON

Export results for analysis or integration with other tools.

```cpp
#include "OpenMagnetics.h"
#include <fstream>
#include <iostream>

int main() {
    // Perform calculations
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        OpenMagnetics::ReluctanceModels::ZHANG
    );
    auto reluctanceResult = reluctanceModel->get_core_reluctance(core);

    // Export core to JSON
    json coreJson;
    to_json(coreJson, core);

    std::ofstream coreFile("core_data.json");
    coreFile << coreJson.dump(2);
    coreFile.close();

    // Export calculation results
    json resultJson;
    to_json(resultJson, reluctanceResult);

    std::ofstream resultFile("reluctance_result.json");
    resultFile << resultJson.dump(2);
    resultFile.close();

    std::cout << "Exported to core_data.json and reluctance_result.json" << std::endl;

    // Import from JSON
    std::ifstream inputFile("core_data.json");
    json importedJson;
    inputFile >> importedJson;

    OpenMagnetics::Core importedCore;
    from_json(importedJson, importedCore);

    return 0;
}
```

---

## Python Examples (PyMKF)

MKF can be used from Python via PyMKF bindings.

### Installation

```bash
# From MKF repository root
pip install -e . -vvv
```

### Basic Usage

```python
import PyMKF

# Find a core
core = PyMKF.find_core_by_name("E 42/21/15")

# Calculate reluctance
result = PyMKF.calculate_reluctance(core, model="ZHANG")
print(f"Reluctance: {result['core_reluctance']:.2e} H⁻¹")
print(f"Fringing factor: {result['maximum_fringing_factor']:.2f}")
```

### Core Losses

```python
import PyMKF

# Define operating conditions
operating_point = {
    "frequency": 100000,
    "magneticFluxDensity": {
        "processed": {
            "peak": 0.1,
            "offset": 0
        }
    }
}

# Calculate losses
core = PyMKF.find_core_by_name("E 42/21/15")
losses = PyMKF.calculate_core_losses(core, operating_point, model="IGSE")
print(f"Core losses: {losses['coreLosses']:.3f} W")
```

### Design Adviser

```python
import PyMKF

# Define requirements
inputs = {
    "designRequirements": {
        "magnetizingInductance": {
            "nominal": 100e-6
        }
    },
    "operatingPoints": [{
        "conditions": {"ambientTemperature": 25},
        "excitationsPerWinding": [{
            "frequency": 100000,
            "current": {
                "processed": {"peak": 2.0, "rms": 1.4}
            }
        }]
    }]
}

# Get recommendations
magnetics = PyMKF.calculate_advised_magnetics(inputs, max_results=5)
for i, mag in enumerate(magnetics):
    print(f"{i+1}. {mag['core']['name']}")
```

### Configuration

```python
import PyMKF

# Get current settings
settings = PyMKF.get_settings()

# Modify settings
settings["wireAdviserIncludeLitz"] = True
settings["wireAdviserIncludeFoil"] = False
settings["coreAdviserIncludeStacks"] = True

# Apply settings
PyMKF.set_settings(settings)
```

---

## More Resources

- **[PyOpenMagnetics](https://github.com/OpenMagnetics/PyOpenMagnetics)** - High-level Python wrapper with additional features
- **[API Reference](../api/index.md)** - Complete C++ API documentation
- **[Physical Models](../models/index.md)** - Detailed model documentation
- **[Algorithm Flowcharts](../algorithms/core-adviser.md)** - How advisers work internally
