#pragma once
#include "WireAdviser.h"
#include "CoilWrapper.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoilAdviser : public WireAdviser {
    private:
        bool _allowMarginTape = true;
        bool _allowInsulatedWire = true;
    public:

        std::vector<MasWrapper> get_advised_coil(MasWrapper mas, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_coil(std::vector<WireWrapper>* wires, MasWrapper mas, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_coil_for_pattern(std::vector<WireWrapper>* wires, MasWrapper mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference);
        std::vector<std::vector<WireSolidInsulationRequirements>> get_solid_insulation_requirements_for_wires(InputsWrapper& inputs, std::vector<size_t> pattern, size_t repetitions);
        void set_allow_margin_tape(bool value) {
            _allowMarginTape = value;
        }
        void set_allow_insulated_wire(bool value) {
            _allowInsulatedWire = value;
        }

};


} // namespace OpenMagnetics