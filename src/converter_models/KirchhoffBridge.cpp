#include "converter_models/KirchhoffBridge.h"

#include "KirchhoffApi.hpp"   // Kirchhoff::api — string/JSON facade only (no KH MAS types cross the boundary)

#include <stdexcept>

namespace OpenMagnetics {

namespace {
// KH's api returns "Exception: ..." on error (never throws across the .so boundary). Re-raise MKF-side.
nlohmann::json parse_or_throw(const std::string& r, const char* what) {
    if (r.rfind("Exception", 0) == 0)
        throw std::runtime_error(std::string("Kirchhoff ") + what + ": " + r);
    return nlohmann::json::parse(r);
}
}  // namespace

nlohmann::json design_converter_tas(const std::string& topology, const nlohmann::json& spec) {
    return parse_or_throw(Kirchhoff::api::design_tas(topology, spec.dump()), "design_tas");
}

Inputs design_requirements_from_tas(const nlohmann::json& tas) {
    nlohmann::json mainInputs = parse_or_throw(Kirchhoff::api::main_magnetic_inputs(tas.dump()),
                                               "main_magnetic_inputs");
    return mainInputs.get<Inputs>();
}

nlohmann::json simulate_tas(const nlohmann::json& tas, const nlohmann::json& fidelity) {
    // simulate_ngspice returns {"success":false,...} (not an "Exception" string) when it has no simulator,
    // so pass the raw JSON through rather than treating that as an error.
    return nlohmann::json::parse(Kirchhoff::api::simulate_ngspice(tas.dump(), fidelity.dump()));
}

}  // namespace OpenMagnetics
