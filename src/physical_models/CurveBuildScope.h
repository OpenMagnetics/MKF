#pragma once
// Internal to InitialPermeability.cpp / AmplitudePermeability.cpp: the loop-invariant state of ONE magnetisation-curve
// build (InitialPermeability::get_magnetisation_curve). The build evaluates the same const material at one temperature
// for ~2600 field strengths; nothing in that loop writes to the material, so what depends only on (material,
// temperature) is computed at its first use inside the build and read from here afterwards, by the evaluations on that
// very object at that very temperature (bit-equal). The guard in get_magnetisation_curve saves and restores it.
#include "MAS.hpp"
#include <array>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics::curve_build {

struct Scope {
    const MAS::CoreMaterial* material = nullptr;
    double temperature = 0;
    std::string interpolatorKey;
    // Roshen loop inputs at the build's temperature: saturation H, saturation B, coercive force, remanence.
    std::optional<std::array<double, 4>> loopInputs;
    // InitialPermeability::get_only_temperature_dependent_indexes of the material's points.
    std::optional<std::vector<size_t>> temperatureIndexes;
};

inline thread_local Scope scope;

inline bool same_temperature(double a, double b) {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// The scope's loop inputs, when this evaluation is the build's own (same object, bit-equal temperature).
inline std::optional<std::array<double, 4>>* loop_inputs_for(const MAS::CoreMaterial& material, double temperature) {
    if (scope.material == &material && same_temperature(scope.temperature, temperature)) {
        return &scope.loopInputs;
    }
    return nullptr;
}

} // namespace OpenMagnetics::curve_build
