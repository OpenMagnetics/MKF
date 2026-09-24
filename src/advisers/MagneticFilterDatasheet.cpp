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
// limits. No analytical physics here — operating values come from `Inputs`,
// datasheet limits from the MAS model. See MagneticFilterDatasheetLimits doc
// comment in MagneticFilter.h for the full contract.

namespace {

// Pick the datasheet electrical entry to validate against. Prefer the entry
// whose numberTurns matches the candidate coil's turns; otherwise fall back to
// the most conservative entry (smallest published rated current) rather than
// guessing. Returns an index into `electricals`. A datasheet-only part has no
// coil to match turns against, so it goes straight to the conservative choice.
size_t select_electrical_entry(const std::vector<MagneticDatasheetElectrical>& electricals,
                               Magnetic* magnetic) {
    if (electricals.size() == 1) {
        return 0;
    }

    // Candidate coil turns. For a CMC/DMC every winding shares a turn count, so
    // the first winding is representative.
    std::optional<double> coilTurns;
    if (magnetic->has_coil()) {
        const auto& windings = magnetic->get_coil().get_functional_description();
        if (!windings.empty()) {
            coilTurns = static_cast<double>(windings[0].get_number_turns());
        }
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

// The rated currents the entry publishes: the plain array, or -- when only the
// ΔT-qualified table is given -- its most conservative (smallest) current, as a
// single limit shared by every winding.
std::optional<std::vector<double>> rated_currents(const MagneticDatasheetElectrical& entry) {
    auto ratedCurrents = entry.get_rated_currents();
    if (ratedCurrents && !ratedCurrents->empty()) {
        return ratedCurrents;
    }
    auto ratedCurrentPoints = entry.get_rated_current_points();
    if (ratedCurrentPoints && !ratedCurrentPoints->empty()) {
        double smallest = std::numeric_limits<double>::max();
        for (const auto& point : ratedCurrentPoints.value()) {
            smallest = std::min(smallest, point.get_current());
        }
        return std::vector<double>{smallest};
    }
    return std::nullopt;
}

// The saturation current to gate the peak against: the smallest the datasheet
// states, across the unqualified peak and every |dL/L| criterion of its table.
// The smallest belongs to the tightest roll-off criterion, so a part is never
// passed on a looser criterion than one it also publishes.
std::optional<double> saturation_current(const MagneticDatasheetElectrical& entry) {
    std::optional<double> smallest = entry.get_saturation_current_peak();
    auto saturationCurrents = entry.get_saturation_currents();
    if (saturationCurrents) {
        for (const auto& point : saturationCurrents.value()) {
            if (!smallest || point.get_current() < smallest.value()) {
                smallest = point.get_current();
            }
        }
    }
    return smallest;
}

std::optional<MagneticDatasheetElectrical> selected_entry(Magnetic* magnetic) {
    auto manufacturerInfo = magnetic->get_manufacturer_info();
    if (!manufacturerInfo || !manufacturerInfo->get_datasheet_info()) {
        return std::nullopt;
    }
    auto electricals = manufacturerInfo->get_datasheet_info()->get_electrical();
    if (!electricals || electricals->empty()) {
        return std::nullopt;
    }
    return electricals.value()[select_electrical_entry(electricals.value(), magnetic)];
}

}  // namespace

bool MagneticFilterDatasheetLimits::applies_to(Magnetic* magnetic) const {
    // Only a part whose datasheet publishes at least one limit this filter
    // checks: every designed magnetic, and a catalogue part that states only an
    // inductance, have nothing here to be judged on.
    auto entry = selected_entry(magnetic);
    if (!entry) {
        return false;
    }
    return rated_currents(entry.value()) || saturation_current(entry.value()) || entry->get_rated_voltage_ac() || entry->get_rated_voltage_dc();
}

std::pair<bool, double> MagneticFilterDatasheetLimits::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    if (!applies_to(magnetic)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "MagneticFilterDatasheetLimits: '" + magnetic->get_reference() + "' publishes no datasheet limit to check; the caller must consult applies_to() first");
    }
    const auto entry = selected_entry(magnetic).value();

    auto ratedCurrents = rated_currents(entry);
    auto saturationCurrent = saturation_current(entry);
    auto ratedVoltageAc = entry.get_rated_voltage_ac();
    auto ratedVoltageDc = entry.get_rated_voltage_dc();

    bool valid = true;
    // The largest fraction of any published limit the operating point uses: 0
    // is untouched, 1 is at the limit, above 1 is over it. Lower is better, the
    // convention every other filter follows, so the adviser's usual invert=true
    // ranks the most comfortable part first.
    double utilisation = 0.0;
    bool anyLimitChecked = false;
    auto check = [&](double limit, double operating) {
        if (limit <= 0) {
            return;  // no meaningful limit published
        }
        anyLimitChecked = true;
        utilisation = std::max(utilisation, operating / limit);
        if (operating > limit) {
            valid = false;
        }
    };

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        const auto& excitations = operatingPoint.get_excitations_per_winding();

        // The saturation current gates the largest winding peak current.
        double maxWindingPeak = 0.0;
        bool haveAnyPeak = false;

        for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
            const auto& excitation = excitations[windingIndex];

            // --- rated current (RMS), per winding ----------------------------
            if (ratedCurrents &&
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
                    check(limit.value(), excitation.get_current()->get_processed()->get_rms().value());
                }
            }

            // collect peak for the saturation-current check
            if (saturationCurrent && excitation.get_current() &&
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
                check(ratedVoltageAc.value(), excitation.get_voltage()->get_processed()->get_rms().value());
            }

            // --- rated DC voltage (offset), per winding ----------------------
            if (ratedVoltageDc && excitation.get_voltage() &&
                excitation.get_voltage()->get_processed()) {
                check(ratedVoltageDc.value(), std::abs(excitation.get_voltage()->get_processed()->get_offset()));
            }
        }

        // --- saturation current (worst winding peak) -------------------------
        if (saturationCurrent && haveAnyPeak) {
            check(saturationCurrent.value(), maxWindingPeak);
        }
    }

    // The part publishes limits, but the operating point carries none of the
    // quantities they limit (e.g. a rated current and no current RMS): there is
    // nothing to compare, which is the caller's input, not the part's merit.
    if (!anyLimitChecked) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "MagneticFilterDatasheetLimits: the operating point carries none of the processed currents or voltages that '" + magnetic->get_reference() + "' publishes limits for");
    }

    return {valid, utilisation};
}

}  // namespace OpenMagnetics
