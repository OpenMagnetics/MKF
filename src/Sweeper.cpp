#include "Sweeper.h"
#include "ComplexPermeability.h"
#include "Constants.h"
#include <cmath>


namespace OpenMagnetics {


Curve2D Sweeper::sweep_impedance_over_frequency(MagneticWrapper magnetic, double start, double stop, size_t numberElements, std::string title) {
    auto frequencies = linear_spaced_array(start, stop, numberElements);
    std::vector<double> impedances;
    for (auto frequency : frequencies) {
        auto impedance = abs(OpenMagnetics::Impedance().calculate_impedance(magnetic, frequency));
        impedances.push_back(impedance);
    }

    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return Curve2D(frequencies, impedances, title);
}

} // namespace OpenMagnetics
