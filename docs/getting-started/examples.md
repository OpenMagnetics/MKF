# Examples

This page provides practical examples for common magnetic component design tasks.

## Designing a Flyback Transformer

```cpp
#include "OpenMagnetics.h"

int main() {
    // Define design requirements
    MAS::DesignRequirements requirements;
    requirements.set_magnetizing_inductance(100e-6);  // 100 µH
    requirements.set_maximum_voltage(400);            // 400 V input
    requirements.set_maximum_current(2.5);            // 2.5 A
    requirements.set_isolation_voltage(3000);         // 3 kV isolation

    // Define operating point
    MAS::OperatingPoint operatingPoint;
    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(40);
    operatingPoint.set_conditions(conditions);

    MAS::OperatingPointExcitation excitation;
    excitation.set_frequency(100000);  // 100 kHz
    operatingPoint.set_excitations_per_winding({excitation});

    // Use MagneticAdviser to find optimal design
    OpenMagnetics::MagneticAdviser adviser;
    auto designs = adviser.get_advised_magnetic(requirements, operatingPoint, 10);

    // Print top 3 designs
    for (size_t i = 0; i < std::min(designs.size(), size_t(3)); ++i) {
        auto& design = designs[i];
        std::cout << "Design " << i + 1 << ":" << std::endl;
        std::cout << "  Core: " << design.get_core().get_name() << std::endl;
        std::cout << "  Score: " << design.get_score() << std::endl;
    }

    return 0;
}
```

## Calculating Winding Losses with Proximity Effect

```cpp
#include "OpenMagnetics.h"

int main() {
    // Load a complete magnetic assembly
    auto magnetic = OpenMagnetics::find_magnetic_by_name("example_transformer");

    // Define operating point
    MAS::OperatingPoint operatingPoint;
    MAS::OperatingConditions conditions;
    conditions.set_ambient_temperature(25);
    operatingPoint.set_conditions(conditions);

    // Set up current waveform for primary winding
    MAS::OperatingPointExcitation primaryExcitation;
    primaryExcitation.set_frequency(100000);

    MAS::SignalDescriptor current;
    MAS::Processed processed;
    processed.set_peak(5.0);  // 5 A peak
    processed.set_rms(3.5);   // RMS value
    current.set_processed(processed);
    primaryExcitation.set_current(current);

    operatingPoint.set_excitations_per_winding({primaryExcitation});

    // Calculate winding losses including proximity effect
    OpenMagnetics::WindingLosses windingLosses;
    auto result = windingLosses.calculate_losses(magnetic, operatingPoint, 25.0);

    std::cout << "Winding Losses Breakdown:" << std::endl;
    for (const auto& element : result.get_winding_losses_per_winding()) {
        std::cout << "  Ohmic (DC): " << element.get_ohmic_losses() << " W" << std::endl;
        if (element.get_skin_effect_losses()) {
            std::cout << "  Skin effect: " << element.get_skin_effect_losses().value() << " W" << std::endl;
        }
        if (element.get_proximity_effect_losses()) {
            std::cout << "  Proximity effect: " << element.get_proximity_effect_losses().value() << " W" << std::endl;
        }
    }

    return 0;
}
```

## Comparing Reluctance Models

```cpp
#include "OpenMagnetics.h"
#include <iostream>
#include <iomanip>

int main() {
    // Load a gapped core
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Add a center gap
    auto gapping = core.get_functional_description().get_gapping();
    MAS::CoreGap gap;
    gap.set_length(0.001);  // 1mm gap
    gap.set_type(MAS::GapType::SUBTRACTIVE);
    gapping.push_back(gap);
    core.get_functional_description().set_gapping(gapping);
    core.process_gap();

    // Compare all reluctance models
    std::vector<OpenMagnetics::ReluctanceModels> models = {
        OpenMagnetics::ReluctanceModels::ZHANG,
        OpenMagnetics::ReluctanceModels::MUEHLETHALER,
        OpenMagnetics::ReluctanceModels::PARTRIDGE,
        OpenMagnetics::ReluctanceModels::STENGLEIN,
        OpenMagnetics::ReluctanceModels::BALAKRISHNAN,
        OpenMagnetics::ReluctanceModels::CLASSIC
    };

    std::cout << std::setw(20) << "Model" << std::setw(20) << "Reluctance (H^-1)"
              << std::setw(20) << "Fringing Factor" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (auto modelType : models) {
        auto model = OpenMagnetics::ReluctanceModel::factory(modelType);
        auto result = model->get_core_reluctance(core);

        std::string modelName;
        switch (modelType) {
            case OpenMagnetics::ReluctanceModels::ZHANG: modelName = "Zhang"; break;
            case OpenMagnetics::ReluctanceModels::MUEHLETHALER: modelName = "Muehlethaler"; break;
            case OpenMagnetics::ReluctanceModels::PARTRIDGE: modelName = "Partridge"; break;
            case OpenMagnetics::ReluctanceModels::STENGLEIN: modelName = "Stenglein"; break;
            case OpenMagnetics::ReluctanceModels::BALAKRISHNAN: modelName = "Balakrishnan"; break;
            case OpenMagnetics::ReluctanceModels::CLASSIC: modelName = "Classic"; break;
            default: modelName = "Unknown";
        }

        std::cout << std::setw(20) << modelName
                  << std::setw(20) << std::scientific << result.get_core_reluctance()
                  << std::setw(20) << std::fixed << result.get_maximum_fringing_factor()
                  << std::endl;
    }

    return 0;
}
```

## Finding Optimal Wire for a Given Frequency

```cpp
#include "OpenMagnetics.h"

int main() {
    // Define requirements
    double frequency = 200000;  // 200 kHz
    double rms_current = 3.0;   // 3 A RMS
    double temperature = 40;    // 40°C

    // Configure wire adviser
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.set_wire_adviser_include_litz(true);
    settings.set_wire_adviser_include_round(true);
    settings.set_wire_adviser_include_rectangular(false);
    settings.set_wire_adviser_include_foil(false);

    // Use WireAdviser
    OpenMagnetics::WireAdviser adviser;
    auto wires = adviser.get_advised_wire(frequency, rms_current, temperature, 20);

    std::cout << "Top wire recommendations for " << frequency/1000 << " kHz:" << std::endl;
    for (size_t i = 0; i < std::min(wires.size(), size_t(5)); ++i) {
        auto& wire = wires[i];
        std::cout << "  " << i+1 << ". " << wire.get_name();
        if (wire.get_type() == MAS::WireType::LITZ) {
            std::cout << " (Litz, " << wire.get_number_conductors().value() << " strands)";
        }
        std::cout << std::endl;
    }

    return 0;
}
```

## Temperature Rise Estimation

```cpp
#include "OpenMagnetics.h"

int main() {
    // Load core and calculate losses
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    // Set up excitation
    MAS::OperatingPointExcitation excitation;
    excitation.set_frequency(100000);

    MAS::SignalDescriptor magneticFluxDensity;
    MAS::Processed processed;
    processed.set_peak(0.15);  // 150 mT
    magneticFluxDensity.set_processed(processed);
    excitation.set_magnetic_flux_density(magneticFluxDensity);

    // Calculate core losses
    OpenMagnetics::CoreLosses coreLossesModel;
    auto coreLosses = coreLossesModel.calculate_core_losses(core, excitation, 25.0);
    double totalLosses = coreLosses.get_core_losses();

    // Estimate temperature rise
    auto temperatureModel = OpenMagnetics::CoreTemperatureModel::factory(
        OpenMagnetics::CoreTemperatureModels::MANIKTALA
    );
    double ambientTemp = 25.0;
    double coreTemp = temperatureModel->get_core_temperature(core, totalLosses, ambientTemp);

    std::cout << "Core losses: " << totalLosses << " W" << std::endl;
    std::cout << "Ambient temperature: " << ambientTemp << " C" << std::endl;
    std::cout << "Estimated core temperature: " << coreTemp << " C" << std::endl;
    std::cout << "Temperature rise: " << coreTemp - ambientTemp << " C" << std::endl;

    return 0;
}
```

## Exporting Results to JSON

```cpp
#include "OpenMagnetics.h"
#include <fstream>

int main() {
    // Perform calculations
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");

    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(
        OpenMagnetics::ReluctanceModels::ZHANG
    );
    auto result = reluctanceModel->get_core_reluctance(core);

    // Convert result to JSON
    json j;
    to_json(j, result);

    // Write to file
    std::ofstream file("reluctance_result.json");
    file << j.dump(2);  // Pretty print with 2-space indent
    file.close();

    std::cout << "Results exported to reluctance_result.json" << std::endl;

    return 0;
}
```

## Using Multiple Operating Points

```cpp
#include "OpenMagnetics.h"

int main() {
    auto core = OpenMagnetics::find_core_by_name("E 42/21/15");
    OpenMagnetics::CoreLosses coreLossesModel;

    // Test multiple operating conditions
    std::vector<double> frequencies = {50000, 100000, 200000, 500000};
    std::vector<double> flux_densities = {0.05, 0.10, 0.15, 0.20};

    std::cout << "Core Losses Map (W):" << std::endl;
    std::cout << "Freq\\B  ";
    for (double B : flux_densities) {
        std::cout << std::setw(10) << B << "T ";
    }
    std::cout << std::endl;

    for (double freq : frequencies) {
        std::cout << std::setw(6) << freq/1000 << "kHz ";
        for (double B : flux_densities) {
            MAS::OperatingPointExcitation excitation;
            excitation.set_frequency(freq);

            MAS::SignalDescriptor magneticFluxDensity;
            MAS::Processed processed;
            processed.set_peak(B);
            magneticFluxDensity.set_processed(processed);
            excitation.set_magnetic_flux_density(magneticFluxDensity);

            try {
                auto losses = coreLossesModel.calculate_core_losses(core, excitation, 25.0);
                std::cout << std::setw(10) << std::fixed << std::setprecision(3)
                          << losses.get_core_losses() << " ";
            } catch (...) {
                std::cout << std::setw(10) << "N/A" << " ";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
```
