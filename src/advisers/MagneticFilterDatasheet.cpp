#include "advisers/MagneticFilter.h"
#include "processors/Inputs.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace OpenMagnetics {

// ABT #19: gate catalogue parts by their OWN datasheet-published electrical
// limits, no-op for designed/custom magnetics that carry no datasheet. No
// analytical physics here — operating values come from `Inputs`, datasheet
// limits from the MAS model. See MagneticFilterDatasheetLimits doc comment in
// MagneticFilter.h for the full contract.

namespace {

// Headroom margin of a single (limit, operating) pair: 1.0 = limit untouched,
// 0.0 = operating at or above the limit. Clamped to [0,1] so the worst-case
// minimum across limits stays a meaningful [0,1] score.
double datasheet_headroom(double limit, double operating) {
    if (limit <= 0) {
        return 1.0;  // no meaningful limit published
    }
    return std::clamp((limit - operating) / limit, 0.0, 1.0);
}

// Pick the datasheet electrical entry to validate against. Prefer the entry
// whose numberTurns matches the candidate coil's turns; otherwise fall back to
// the most conservative entry (smallest published rated current) rather than
// guessing. Returns an index into `electricals`.
size_t select_electrical_entry(const std::vector<MagneticDatasheetElectrical>& electricals,
                               Magnetic* magnetic) {
    if (electricals.size() == 1) {
        return 0;
    }

    // Candidate coil turns. For a CMC/DMC every winding shares a turn count, so
    // the first winding is representative.
    std::optional<double> coilTurns;
    const auto& windings = magnetic->get_coil().get_functional_description();
    if (!windings.empty()) {
        coilTurns = static_cast<double>(windings[0].get_number_turns());
    }

    if (coilTurns) {
        for (size_t i = 0; i < electricals.size(); ++i) {
            auto numberTurns = electricals[i].get_number_turns();
            if (numberTurns && std::abs(numberTurns.value() - coilTurns.value()) < 0.5) {
                return i;
            }
        }
    }

    // No turn match (or turns unset on every entry): most conservative entry =
    // smallest published rated current. Entries without a rated current are
    // only used as a last resort (no entry publishes one).
    size_t mostConservative = 0;
    double smallestRated = std::numeric_limits<double>::max();
    bool anyRated = false;
    for (size_t i = 0; i < electricals.size(); ++i) {
        auto ratedCurrents = electricals[i].get_rated_currents();
        if (ratedCurrents && !ratedCurrents->empty()) {
            double smallest = *std::min_element(ratedCurrents->begin(), ratedCurrents->end());
            if (smallest < smallestRated) {
                smallestRated = smallest;
                mostConservative = i;
                anyRated = true;
            }
        }
    }
    return anyRated ? mostConservative : 0;
}

}  // namespace

std::pair<bool, double> MagneticFilterDatasheetLimits::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    // Pass-through for any magnetic without published datasheet electricals —
    // i.e. every designed/custom part. Neutral score so it is neither rejected
    // nor penalised relative to catalogue parts.
    auto manufacturerInfo = magnetic->get_manufacturer_info();
    if (!manufacturerInfo || !manufacturerInfo->get_datasheet_info()) {
        return {true, 1.0};
    }
    auto datasheetInfo = manufacturerInfo->get_datasheet_info().value();
    auto electricals = datasheetInfo.get_electrical();
    if (!electricals || electricals->empty()) {
        return {true, 1.0};
    }

    const auto& entry = electricals.value()[select_electrical_entry(electricals.value(), magnetic)];

    auto ratedCurrents = entry.get_rated_currents();
    auto saturationCurrentPeak = entry.get_saturation_current_peak();
    auto ratedVoltageAc = entry.get_rated_voltage_ac();
    auto ratedVoltageDc = entry.get_rated_voltage_dc();

    bool valid = true;
    double worstScore = 1.0;        // smallest headroom across all checked limits
    bool anyLimitChecked = false;

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        const auto& excitations = operatingPoint.get_excitations_per_winding();

        // saturationCurrentPeak gates the largest winding peak current.
        double maxWindingPeak = 0.0;
        bool haveAnyPeak = false;

        for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
            const auto& excitation = excitations[windingIndex];

            // --- rated current (RMS), per winding ----------------------------
            if (ratedCurrents && !ratedCurrents->empty() &&
                excitation.get_current() && excitation.get_current()->get_processed() &&
                excitation.get_current()->get_processed()->get_rms()) {
                // Single-entry array ⇒ same limit for every winding; multi-entry
                // ⇒ per-winding (windings beyond the published list are skipped,
                // not gated against a guessed limit).
                std::optional<double> limit;
                if (ratedCurrents->size() == 1) {
                    limit = ratedCurrents.value()[0];
                } else if (windingIndex < ratedCurrents->size()) {
                    limit = ratedCurrents.value()[windingIndex];
                }
                if (limit) {
                    double currentRms = excitation.get_current()->get_processed()->get_rms().value();
                    anyLimitChecked = true;
                    worstScore = std::min(worstScore, datasheet_headroom(limit.value(), currentRms));
                    if (currentRms > limit.value()) {
                        valid = false;
                    }
                }
            }

            // collect peak for the saturation-current check
            if (saturationCurrentPeak && excitation.get_current() &&
                excitation.get_current()->get_processed() &&
                excitation.get_current()->get_processed()->get_peak()) {
                maxWindingPeak = std::max(maxWindingPeak,
                                          excitation.get_current()->get_processed()->get_peak().value());
                haveAnyPeak = true;
            }

            // --- rated AC voltage (RMS), per winding -------------------------
            if (ratedVoltageAc && excitation.get_voltage() &&
                excitation.get_voltage()->get_processed() &&
                excitation.get_voltage()->get_processed()->get_rms()) {
                double voltageRms = excitation.get_voltage()->get_processed()->get_rms().value();
                anyLimitChecked = true;
                worstScore = std::min(worstScore, datasheet_headroom(ratedVoltageAc.value(), voltageRms));
                if (voltageRms > ratedVoltageAc.value()) {
                    valid = false;
                }
            }

            // --- rated DC voltage (offset), per winding ----------------------
            if (ratedVoltageDc && excitation.get_voltage() &&
                excitation.get_voltage()->get_processed()) {
                double voltageDc = std::abs(excitation.get_voltage()->get_processed()->get_offset());
                anyLimitChecked = true;
                worstScore = std::min(worstScore, datasheet_headroom(ratedVoltageDc.value(), voltageDc));
                if (voltageDc > ratedVoltageDc.value()) {
                    valid = false;
                }
            }
        }

        // --- saturation current peak (worst winding peak) --------------------
        if (saturationCurrentPeak && haveAnyPeak) {
            anyLimitChecked = true;
            worstScore = std::min(worstScore, datasheet_headroom(saturationCurrentPeak.value(), maxWindingPeak));
            if (maxWindingPeak > saturationCurrentPeak.value()) {
                valid = false;
            }
        }
    }

    // Datasheet present but no limit field we know how to check (e.g. only
    // `inductance`): nothing to gate on ⇒ neutral pass. This is the "only if it
    // exists" contract — a part is never rejected for a limit it doesn't publish.
    if (!anyLimitChecked) {
        return {true, 1.0};
    }

    return {valid, worstScore};
}

}  // namespace OpenMagnetics
