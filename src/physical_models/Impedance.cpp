#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/Reluctance.h"
#include "physical_models/StrayCapacitance.h"
#include "physical_models/MagnetizingInductance.h"
#include "constructive_models/NumberTurns.h"
#include "support/Utils.h"
#include <cmath>


namespace OpenMagnetics {

std::complex<double> Impedance::calculate_impedance(Magnetic magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_impedance(core, coil, frequency, temperature);
}

std::complex<double> Impedance::calculate_impedance(Core core, Coil coil, double frequency, double temperature) {
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
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil).get_capacitance_among_windings().value();
        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
    }

    auto capacitiveImpedance = std::complex<double>(0, 1.0 / (angularFrequency * capacitance));

    auto impedance = 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance);

    return impedance;
}

int64_t Impedance::calculate_minimum_number_turns(Magnetic magnetic, Inputs inputs) {
    NumberTurns numberTurns(1, inputs.get_design_requirements());
    if (!inputs.get_design_requirements().get_minimum_impedance()) {
        throw std::invalid_argument("Missing impedance requirement");
    }
    auto minimumImpedanceRequirement = inputs.get_design_requirements().get_minimum_impedance().value();
    auto temperature = inputs.get_maximum_temperature();
    auto timeout = _maximumNumberTurns;
    int64_t calculatedNumberTurns = -1;
    while (true) {
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        calculatedNumberTurns = numberTurnsCombination[0];
        magnetic.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(calculatedNumberTurns);
        auto selfResonantFrequency = calculate_self_resonant_frequency(magnetic);

        bool validDesign = true;
        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto requiredImpedanceMagnitude = impedanceAtFrequency.get_impedance().get_magnitude();
            auto impedanceMagnitude = abs(calculate_impedance(magnetic, frequency, temperature));
            if (impedanceMagnitude < requiredImpedanceMagnitude) {
                validDesign = false;
            }
            if (frequency > 0.25 * selfResonantFrequency) {  // hardcoded 20% of SRF
                validDesign = false;
                return -1;
            }
        }

        if (validDesign) {
            break;
        }

        if (timeout == 0) {
            return -1;
        }
        timeout--;
    }

    return calculatedNumberTurns;
}

double Impedance::calculate_q_factor(Magnetic magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_q_factor(core, coil, frequency, temperature);
}

double Impedance::calculate_q_factor(Core core, Coil coil, double frequency, double temperature) {
    auto impedance = calculate_impedance(core, coil, frequency, temperature);
    return impedance.imag() / impedance.real();
}

double Impedance::calculate_self_resonant_frequency(Magnetic magnetic, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_self_resonant_frequency(core, coil, temperature);
}

double Impedance::calculate_self_resonant_frequency(Core core, Coil coil, double temperature) {
    double capacitance;
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil);
    }
    else {
        if (!coil.get_turns_description()) {
            coil.wind();
        }
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil).get_capacitance_among_windings().value();

        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
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
