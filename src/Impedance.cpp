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
    auto impedance = angularFrequency * airCoredInductance * std::complex<double>(complexPermeabilityImaginaryPart, -complexPermeabilityRealPart);

    return impedance;
}

} // namespace OpenMagnetics
