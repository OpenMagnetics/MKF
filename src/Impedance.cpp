#include "Impedance.h"
#include "ComplexPermeability.h"
#include "Reluctance.h"
#include "CoreLosses.h"
#include "Constants.h"
#include "Utils.h"
#include <cmath>


namespace OpenMagnetics {

std::complex<double>  Impedance::calculate_impedance(MagneticWrapper magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_impedance(core, coil, frequency, temperature);
}

double capacitance_turn_to_turn(double turnDiameter, double wireRadius, double centerSeparation) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double epsilon0 = Constants().vacuumPermittivity;
    double ctt = pow(std::numbers::pi, 2) * turnDiameter * epsilon0 / log(centerSeparation / (2 * wireRadius) + sqrt(pow(centerSeparation / (2 * wireRadius), 2) - 1));
    return ctt;
}

double capacitance_turn_to_shield(double turnDiameter, double wireRadius, double distance) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    double epsilon0 = Constants().vacuumPermittivity;
    double cts = 2 * pow(std::numbers::pi, 2) * turnDiameter * epsilon0 / log(distance / wireRadius + sqrt(pow(distance / wireRadius, 2) - 1));
    return cts;
}

double cab(double n, double ctt, double cts) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    if (n == 2) {
        return ctt + cts / 2;
    }
    else if (n == 3) {
        return ctt / 2 + cts / 2;
    }
    else {
        double cabValue = cab(n - 2, ctt, cts);
        return (cabValue * ctt / 2) / (cabValue  + ctt / 2) + cts / 2;
    }
}

double cas(double n, double ctt, double cts) {
    // According to https://sci-hub.st/https://ieeexplore.ieee.org/document/793378
    if (n == 1) {
        return cts;
    }
    else {
        double casValue = cas(n - 1, ctt, cts);
        return (casValue * ctt) / (casValue  + ctt) + cts;
    }
}

std::complex<double> Impedance::calculate_impedance(CoreWrapper core, CoilWrapper coil, double frequency, double temperature) {
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

    // We are assuming one layer for now, because i fogot to push to good code in Osma
    auto wireRadius = coil.resolve_wire(0).get_maximum_conducting_width() / 2;
    double distanceTurnsToCore = coil.resolve_bobbin().get_processed_description()->get_column_thickness() + coil.resolve_wire(0).get_maximum_outer_width() / 2;
    double turnDiameter = 2 * std::numbers::pi * (coil.resolve_bobbin().get_processed_description()->get_column_width().value() + wireRadius);
    double centerSeparation = coil.resolve_wire(0).get_maximum_outer_width();
    if (coil.get_turns_description()) {
        if (coil.get_turns_description().value().size() > 1) {
            auto firstTurn = coil.get_turns_description().value()[0];
            auto secondTurn = coil.get_turns_description().value()[1];
            centerSeparation = sqrt(pow(firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0], 2) + pow(firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1], 2));
        }
    }

    double ctt = capacitance_turn_to_turn(turnDiameter, wireRadius, centerSeparation);
    double cts = capacitance_turn_to_shield(turnDiameter, wireRadius, distanceTurnsToCore);
    double cabValue = cab(numberTurns, ctt, cts);
    double casValue = cas(numberTurns, ctt, cts);
    double C2 = 2 * cabValue * casValue / (4 * cabValue - casValue);

    auto capacitance = C2 * 2;
    if (coil.get_layers_description()) {
        capacitance *= coil.get_layers_by_winding_index(0).size();
    }
    auto capacitiveImpedance = std::complex<double>(0, 1.0 / (angularFrequency * capacitance));

    auto coreLossesModel = OpenMagnetics::CoreLossesModel::factory(CoreLossesModels::LOSS_FACTOR);

    auto seriesResistance = coreLossesModel->get_core_losses_series_resistance(core, frequency, temperature, abs(inductiveImpedance));
    auto resistiveImpedance = std::complex<double>(seriesResistance, 0);

    auto impedance = 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance + 1.0 / resistiveImpedance);

    return impedance;
}

} // namespace OpenMagnetics
