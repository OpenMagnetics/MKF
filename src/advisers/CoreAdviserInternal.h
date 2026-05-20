#pragma once

// Internal helpers shared across the CoreAdviser translation units
// (CoreAdviser.cpp, CoreAdviserGapping.cpp, CoreAdviserMaterials.cpp,
// CoreAdviserPipeline.cpp, ...). Not part of the public API; only sibling
// .cpp files in src/advisers/ should include this header.
//
// Phase 5 (structural cleanup): these helpers were originally defined in
// an anonymous namespace inside CoreAdviser.cpp. They were promoted here
// once the file was split so the same compact log-format / sinusoidal-
// excitation builder / Steinmetz+Proprietary model-pair factory can be
// reused without duplication.

#include "constructive_models/Magnetic.h"
#include "physical_models/CoreLosses.h"
#include "support/Logger.h"
#include <memory>
#include <string>

namespace OpenMagnetics {

// Compact replacement for the verbatim
//   logEntry("There are " + std::to_string(c.size()) + " magnetics after the X filter.", "CoreAdviser");
// idioms that used to be scattered across the pipeline functions. Format:
//   "After <stage>: <n>"   /   "Pruned to <n> before <stage>"
inline void log_stage(const std::string& stage, size_t count) {
    logEntry("After " + stage + ": " + std::to_string(count), "CoreAdviser");
}
inline void log_pruned(const std::string& stage, size_t count) {
    logEntry("Pruned to " + std::to_string(count) + " before " + stage, "CoreAdviser");
}

// Sinusoidal-excitation builder used by both add_ferrite_materials_by_losses
// and add_ferrite_materials_by_impedance.
inline OperatingPointExcitation make_sinusoidal_excitation(double peak, double offset, double frequency) {
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(offset);
    processed.set_peak(peak);
    processed.set_peak_to_peak(2 * peak);
    magneticFluxDensity.set_processed(processed);
    OperatingPointExcitation excitation;
    excitation.set_magnetic_flux_density(magneticFluxDensity);
    excitation.set_frequency(frequency);
    return excitation;
}

struct CoreLossesModelPair {
    std::shared_ptr<CoreLossesModel> steinmetz;
    std::shared_ptr<CoreLossesModel> proprietary;
};
inline CoreLossesModelPair make_default_core_losses_model_pair() {
    return {
        CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}})),
        CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}))
    };
}

// Forward declarations for free helpers that still live in CoreAdviser.cpp
// but are called from sibling TUs (CoreAdviserPipeline.cpp, CoreAdviserFilters.cpp).
void add_initial_turns_by_inductance(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, Inputs inputs);
void add_alternative_materials(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, Inputs inputs);
std::vector<double> normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration);
void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring);

} // namespace OpenMagnetics
