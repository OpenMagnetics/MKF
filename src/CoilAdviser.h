#pragma once
#include "WireAdviser.h"
#include "CoilWrapper.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoilAdviser : public WireAdviser {
    public:

        std::vector<std::pair<MasWrapper, double>> get_advised_coil(MasWrapper mas, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_coil(std::vector<WireWrapper>* wires, MasWrapper mas, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_coil_for_pattern(std::vector<WireWrapper>* wires, MasWrapper mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> withstandVoltageForWires, size_t maximumNumberResults);
        std::vector<std::vector<WireSolidInsulationRequirements>> get_solid_insulation_requirements_for_wires(InputsWrapper inputs);

};


} // namespace OpenMagnetics