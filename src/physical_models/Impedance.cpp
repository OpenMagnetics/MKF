#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/Reluctance.h"
#include "physical_models/StrayCapacitance.h"
#include "physical_models/MagnetizingInductance.h"
#include "Constants.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include <cmath>


namespace OpenMagnetics {

std::complex<double>  Impedance::calculate_impedance(MagneticWrapper magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_impedance(core, coil, frequency, temperature);
}

std::complex<double> Impedance::calculate_impedance(CoreWrapper core, CoilWrapper coil, double frequency, double temperature) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    double numberTurns = coil.get_functional_description()[0].get_number_turns();
    double reluctanceCoreUnityPermeability = reluctanceModel->get_core_reluctance(core, 1).get_core_reluctance();

    OpenMagnetics::ComplexPermeability complexPermeabilityObj;
    auto coreMaterial = core.resolve_material();
    auto [complexPermeabilityRealPart, complexPermeabilityImaginaryPart] = complexPermeabilityObj.get_complex_permeability(coreMaterial, frequency);

    auto angularFrequency = 2 * std::numbers::pi * frequency;
    double airCoredInductance = numberTurns * numberTurns / reluctanceCoreUnityPermeability;
    auto inductiveImpedance = angularFrequency * airCoredInductance * std::complex<double>(complexPermeabilityImaginaryPart, -complexPermeabilityRealPart);

    double capacitance;
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil);
    }
    else {
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance_among_windings(coil);
        auto windingsKey = std::make_pair(coil.get_functional_description()[0].get_name(), coil.get_functional_description()[0].get_name());
        capacitance = capacitanceMatrix[windingsKey];
    }

    auto capacitiveImpedance = std::complex<double>(0, 1.0 / (angularFrequency * capacitance));

    auto impedance = 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance);

    return impedance;
}

double Impedance::calculate_self_resonant_frequency(MagneticWrapper magnetic, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_self_resonant_frequency(core, coil, temperature);
}

double Impedance::calculate_self_resonant_frequency(CoreWrapper core, CoilWrapper coil, double temperature) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double capacitance;
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil);
    }
    else {
        if (!coil.get_turns_description()) {
            coil.wind();
        }
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance_among_windings(coil);
        auto windingsKey = std::make_pair(coil.get_functional_description()[0].get_name(), coil.get_functional_description()[0].get_name());

        capacitance = capacitanceMatrix[windingsKey];
    }

    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_cooling(std::nullopt);
    conditions.set_ambient_temperature(temperature);
    operatingPoint.set_conditions(conditions);
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil, nullptr).get_magnetizing_inductance().get_nominal().value();

    double selfResonantFrequency = 1.0 / (2 * std::numbers::pi * sqrt(magnetizingInductance * capacitance));

    return selfResonantFrequency;
}

} // namespace OpenMagnetics
