#include "Impedance.h"
#include "ComplexPermeability.h"
#include "Reluctance.h"
#include "Constants.h"
#include "Utils.h"
#include <cmath>


namespace OpenMagnetics {

std::complex<double>  Impedance::calculate_impedance(MagneticWrapper magnetic, double frequency) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_impedance(core, coil, frequency);
}

std::complex<double> Impedance::calculate_impedance(CoreWrapper core, CoilWrapper coil, double frequency) {
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

    // auto capacitance= 4e-11;
    // auto capacitiveImpedance = std::complex<double>(0, 1.0 / (angularFrequency * capacitance));
    // // auto resistiveImpedance = std::complex<double>(1e5, 0);

    // auto impedance = 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance);
    // // auto impedance = 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance + 1.0 / resistiveImpedance);
    // std::cout << "**********************************" << std::endl;
    // std::cout << "inductiveImpedance: " << inductiveImpedance << std::endl;
    // std::cout << "capacitiveImpedance: " << capacitiveImpedance << std::endl;
    // // std::cout << "resistiveImpedance: " << resistiveImpedance << std::endl;
    // std::cout << "impedance: " << impedance << std::endl;

    return inductiveImpedance;
}

} // namespace OpenMagnetics
