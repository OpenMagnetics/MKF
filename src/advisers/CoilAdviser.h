#pragma once
#include "advisers/WireAdviser.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Mas.h"
#include "support/Utils.h"
#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class CoilAdviser : public WireAdviser {
    private:
        bool _allowMarginTape = true;
        bool _allowInsulatedWire = true;
        OpenMagnetics::WireAdviser _wireAdviser;
        std::optional<WireStandard> _commonWireStandard = defaults.commonWireStandard;
    public:

        std::vector<Mas> get_advised_coil(Mas mas, size_t maximumNumberResults=1);
        std::vector<Mas> get_advised_coil(std::vector<Wire>* wires, Mas mas, size_t maximumNumberResults=1);
        std::vector<Section> get_advised_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions);
        std::vector<Section> get_advised_planar_sections(Mas mas, std::vector<size_t> pattern, size_t repetitions);
        std::vector<Mas> get_advised_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, std::vector<WireSolidInsulationRequirements> solidInsulationRequirementsForWires, size_t maximumNumberResults, std::string reference);
        std::vector<Mas> get_advised_planar_coil_for_pattern(std::vector<Wire>* wires, Mas mas, std::vector<size_t> pattern, size_t repetitions, size_t maximumNumberResults, std::string reference);
        void set_allow_margin_tape(bool value) {
            _allowMarginTape = value;
        }
        void set_allow_insulated_wire(bool value) {
            _allowInsulatedWire = value;
        }
        void set_common_wire_standard(std::optional<WireStandard> commonWireStandard) {
            _commonWireStandard = commonWireStandard;
        }
        std::optional<WireStandard> get_common_wire_standard() {
            return _commonWireStandard;
        }

};


} // namespace OpenMagnetics