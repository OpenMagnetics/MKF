#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/CoreTemperature.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace OpenMagnetics {

std::pair<bool, double> MagneticFilterDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    auto depth = core.get_depth();

    if (magnetic->get_coil().get_layers_description()) {
        auto coilDepth = magnetic->get_mutable_core().get_columns()[0].get_depth();
        auto layers = magnetic->get_coil().get_layers_description().value();
        for (auto layer : layers) {
            coilDepth += layer.get_dimensions()[0] * 2;
        }
        depth = std::max(depth, coilDepth);
    }

    double volume = core.get_width() * core.get_height() * depth;

    // Hard-reject planar-family cores when wiring is wound wire. PLANAR_E / PLANAR_EL /
    // PLANAR_ER / LP cores assume PCB/planar conductors; a wound coil cannot make
    // efficient use of the wide-shallow window and the resulting design is always
    // dominated by a properly-shaped core. (Note: winding-window width/height are not
    // populated for planar shapes, so we key off the shape family enum.)
    if (inputs && inputs->get_wiring_technology() == WiringTechnology::WOUND) {
        auto family = magnetic->get_mutable_core().get_shape_family();
        if (family == CoreShapeFamily::PLANAR_E ||
            family == CoreShapeFamily::PLANAR_EL ||
            family == CoreShapeFamily::PLANAR_ER ||
            family == CoreShapeFamily::LP) {
            return {false, volume};
        }
    }

    return {true, volume};
}

std::pair<bool, double> MagneticFilterTurnCount::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto coil = magnetic->get_coil();

    // Sum N across all windings. For a CMC the windings are equal so this
    // is 2 × N; for a single inductor it's just N.
    double totalTurns = 0;
    for (const auto& winding : coil.get_functional_description()) {
        totalTurns += static_cast<double>(winding.get_number_turns());
    }

    if (totalTurns <= 0) {
        // Turns not yet assigned — return a neutral, non-rejecting score.
        // Caller is responsible for ordering (run after add_initial_turns_*).
        return {true, 0.0};
    }

    // Manufacturability proxy = wire length used ≈ N × (core size). We must
    // multiply by a size proxy: scoring N alone is monotonically minimised
    // by the largest core in the catalogue (because for a fixed |Z| target
    // Z = µ_complex · ω · N² · Aₑ / lₑ, N ∝ 1/√(µ·Aₑ/lₑ), so bigger core
    // → fewer turns). That made the CMC adviser pick T 140-class toroids
    // for plain 230 V / 5 A / 1 kΩ@150 kHz wizard defaults — see
    // TestTopologyCmc::Test_Cmc_AdviserMustNotPickOversizedToroid_WizardDefaults.
    //
    // get_width() returns the OD for toroids and the A-dimension for
    // two-piece sets; both scale linearly with overall core size and are
    // a reasonable stand-in for mean turn length without forcing the
    // filter to instantiate a bobbin (which the cheap pre-filter pipeline
    // explicitly avoids). Falls back to raw N if width is unavailable.
    // Avoid get_width() here: it throws CoreNotProcessedException on
    // unprocessed cores, and exception throwing in the hot pre-filter loop
    // (thousands of cores × multiple passes) costs seconds. Read the
    // processed description directly instead.
    auto& core = magnetic->get_mutable_core();
    auto processed = core.get_processed_description();
    if (!processed) {
        return {true, totalTurns};
    }
    double widthProxy = processed->get_width();
    if (!std::isfinite(widthProxy) || widthProxy <= 0.0) {
        return {true, totalTurns};
    }

    return {true, totalTurns * widthProxy};
}

std::pair<bool, double> MagneticFilterTurnsRatios::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_turns_ratios().size() > 0) {
        auto magneticTurnsRatios = magnetic->get_turns_ratios();
        if (magneticTurnsRatios.size() != inputs->get_design_requirements().get_turns_ratios().size()) {
            return {false, 0.0};
        }
        for (size_t i = 0; i < inputs->get_design_requirements().get_turns_ratios().size(); ++i) {
            auto turnsRatioRequirement = inputs->get_design_requirements().get_turns_ratios()[i];
            // TODO: Try all combinations of turns ratios, not just the default order
            if (!check_requirement(turnsRatioRequirement, magneticTurnsRatios[i])) {
                return {false, 0.0};
            }
            scoring += abs(resolve_dimensional_values(turnsRatioRequirement) - resolve_dimensional_values(magneticTurnsRatios[i]));
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterMaximumDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_maximum_dimensions()) {
        auto maximumDimensions = inputs->get_design_requirements().get_maximum_dimensions().value();
        auto magneticDimensions = magnetic->get_maximum_dimensions();
        scoring = sqrt(pow(maximumDimensions.get_width().value() - magneticDimensions[0], 2) + pow(maximumDimensions.get_height().value() - magneticDimensions[1], 2)+ pow(maximumDimensions.get_depth().value() - magneticDimensions[2], 2));
        if (!magnetic->fits(maximumDimensions, true)) {
            valid = false;
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double volume = maximumDimensions[0] * maximumDimensions[1] * maximumDimensions[2];
    return {true, volume};
}

std::pair<bool, double> MagneticFilterArea::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double area = maximumDimensions[0] * maximumDimensions[2];
    return {true, area};
}

std::pair<bool, double> MagneticFilterHeight::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double height = maximumDimensions[1];
    return {true, height};
}

std::pair<bool, double> MagneticFilterTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);

    auto coreTemperatureModel = CoreTemperatureModel::factory(defaults.coreTemperatureModelDefault);

    auto coreTemperature = coreTemperatureModel->get_core_temperature(magnetic->get_core(), losses, defaults.ambientTemperature);
    double calculatedTemperature = coreTemperature.get_maximum_temperature();

    return {true, calculatedTemperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    MagneticFilterLosses fallback;
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES, fallback);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}

std::pair<bool, double> MagneticFilterVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    MagneticFilterLosses lossesFallback;
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES, lossesFallback);
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

} // namespace OpenMagnetics
