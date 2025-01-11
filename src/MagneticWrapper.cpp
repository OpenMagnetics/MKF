#include <MAS.hpp>
#include "MagneticWrapper.h"

namespace OpenMagnetics {

BobbinWrapper MagneticWrapper::get_bobbin() {
    return get_mutable_coil().resolve_bobbin();
}

std::vector<WireWrapper> MagneticWrapper::get_wires() {
    return get_mutable_coil().get_wires();
}

std::vector<double> MagneticWrapper::get_turns_ratios() {
    return get_mutable_coil().get_turns_ratios();
}

WireWrapper MagneticWrapper::get_wire(size_t windingIndex) {
    return get_mutable_coil().resolve_wire(windingIndex);
}

std::string MagneticWrapper::get_reference() {
    if (get_manufacturer_info()) {
        if (get_manufacturer_info()->get_reference()) {
            return get_manufacturer_info()->get_reference().value();
        }
    }
    return "Custom component made with OpenMagnetic";
}

        
} // namespace OpenMagnetics
