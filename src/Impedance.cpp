#include "Impedance.h"
#include "ComplexPermeability.h"
#include "MagnetizingInductance.h"
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

    OpenMagnetics::ComplexPermeability complexPermeabilityObj;
    auto coreMaterial = core.resolve_material();
    auto [complexPermeabilityRealPart, complexPermeabilityImaginaryPart] = complexPermeabilityObj.get_complex_permeability(coreMaterial, frequency);
    MagnetizingInductance magnetizingInductance;
    auto airCoredSolenoidInductance = magnetizingInductance.calculate_inductance_air_solenoid(core, coil);

    auto complexPermeability = constants.vacuumPermeability * std::complex<double>(complexPermeabilityRealPart, -complexPermeabilityImaginaryPart);
    auto angularFrequency = 2 * std::numbers::pi * frequency;
    auto impedance = std::complex<double>(0, angularFrequency * airCoredSolenoidInductance) * complexPermeability;

    return impedance;
}

} // namespace OpenMagnetics
