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
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

namespace OpenMagnetics {

// Compact replacement for the verbatim
//   logEntry("There are " + std::to_string(c.size()) + " magnetics after the X filter.", "CoreAdviser");
// idioms that used to be scattered across the pipeline functions. Format:
//   "After <stage>: <n>"   /   "Pruned to <n> before <stage>"
//
// CORE_ADVISER_PROFILE=1 in the environment additionally streams the
// elapsed time between consecutive log_stage calls to stderr, for
// per-stage timing investigations (e.g. finding the bottleneck in slow
// converter wizards). Off by default — production callers and the
// WebFrontend WASM build see only the existing log entries.
// Shared cursor for the CORE_ADVISER_PROFILE stream, so that log_stage and
// log_probe below report deltas against the same previous mark rather than two
// independent ones (otherwise a probe between two stages would leave the stage
// delta double-counting the probe's own span).
inline std::chrono::steady_clock::time_point& profile_cursor() {
    thread_local auto prev = std::chrono::steady_clock::now();
    return prev;
}

inline void profile_mark(const std::string& label, size_t count) {
    if (const char* dbg = std::getenv("CORE_ADVISER_PROFILE"); dbg && dbg[0] == '1') {
        auto& prev = profile_cursor();
        auto now = std::chrono::steady_clock::now();
        auto dt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - prev).count();
        std::cerr << "[CA-PROFILE] " << label << ": " << dt_ms << " ms (count=" << count << ")\n";
        prev = now;
    }
}

inline void log_stage(const std::string& stage, size_t count) {
    logEntry("After " + stage + ": " + std::to_string(count), "CoreAdviser");
    profile_mark(stage, count);
}

// Profile-only probe for stretches that emit no stage log at all. ABT #859: the
// DMC run spends minutes both BEFORE the first logged stage (dataset build,
// input pre-processing) and AFTER the last one (per-candidate winding in
// post_process_and_cut), and neither span was visible to CORE_ADVISER_PROFILE
// because that stream is driven exclusively by log_stage. Unlike log_stage this
// writes nothing at all unless CORE_ADVISER_PROFILE=1, so it can sit on a hot
// path without changing what production logging emits.
inline void log_probe(const std::string& label, size_t count = 0) {
    profile_mark(label, count);
}
inline void log_pruned(const std::string& stage, size_t count) {
    logEntry("Pruned to " + std::to_string(count) + " before " + stage, "CoreAdviser");
}

// Sinusoidal-excitation builder used by both add_ferrite_materials_by_losses
// and add_ferrite_materials_by_impedance.
inline OperatingPointExcitation make_sinusoidal_excitation(double peak, double offset, double frequency) {
    SignalDescriptor magneticFluxDensity;
    ProcessedWaveform processed;
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
void add_initial_turns_by_inductance(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, const Inputs& inputs);
void add_alternative_materials(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, Inputs inputs);
std::vector<double> normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration);
void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring);

} // namespace OpenMagnetics
