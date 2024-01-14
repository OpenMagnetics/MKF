#include <MAS.hpp>
#include "MagneticWrapper.h"

namespace OpenMagnetics {

BobbinWrapper MagneticWrapper::get_bobbin() {
    return get_mutable_coil().resolve_bobbin();
}

std::vector<WireWrapper> MagneticWrapper::get_wires() {
    return get_mutable_coil().get_wires();
}

WireWrapper MagneticWrapper::get_wire(size_t windingIndex) {
    return get_mutable_coil().resolve_wire(windingIndex);
}

        
} // namespace OpenMagnetics
