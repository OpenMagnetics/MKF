#include "Settings.h"
#include "WindingLosses.h"
#include "MAS.hpp"
#include "MagneticField.h"
#include "WindingOhmicLosses.h"
#include "WindingSkinEffectLosses.h"
#include "MagnetizingInductance.h"
#include "WindingProximityEffectLosses.h"

namespace OpenMagnetics {

WindingLossesPerElement combine_turn_losses_per_element(std::vector<WindingLossesPerElement> windingLossesPerTurn, std::vector<size_t> turnIndexesToCombine) {
    WindingLossesPerElement windingLossesThisElement;
    OhmicLosses ohmicLossesThisElement;
    WindingLossElement skinEffectLossesThisElement;
    WindingLossElement proximityEffectLossesThisElement;

    for (auto turnIndex : turnIndexesToCombine) {
        auto windingLossesThisTurn = windingLossesPerTurn[turnIndex];

        if (!windingLossesThisElement.get_ohmic_losses()) {
            ohmicLossesThisElement.set_losses(0);
            ohmicLossesThisElement.set_method_used(windingLossesThisTurn.get_ohmic_losses()->get_method_used());
            ohmicLossesThisElement.set_origin(windingLossesThisTurn.get_ohmic_losses()->get_origin());
            windingLossesThisElement.set_ohmic_losses(ohmicLossesThisElement);
        }
        ohmicLossesThisElement.get_mutable_losses() += windingLossesThisTurn.get_ohmic_losses()->get_losses();

        auto skinEffectHarmonicFrequencies = windingLossesThisTurn.get_skin_effect_losses()->get_harmonic_frequencies();
        auto skinEffectLossesPerHarmonic = windingLossesThisTurn.get_skin_effect_losses()->get_losses_per_harmonic();

        if (!windingLossesThisElement.get_skin_effect_losses()) {
            skinEffectLossesThisElement.set_method_used(windingLossesThisTurn.get_skin_effect_losses()->get_method_used());
            skinEffectLossesThisElement.set_origin(windingLossesThisTurn.get_skin_effect_losses()->get_origin());
            skinEffectLossesThisElement.set_harmonic_frequencies(skinEffectHarmonicFrequencies);
            skinEffectLossesThisElement.set_losses_per_harmonic(std::vector<double>(skinEffectLossesPerHarmonic.size(), 0));
            windingLossesThisElement.set_skin_effect_losses(skinEffectLossesThisElement);
        }
        for (size_t harmonicIndex = 0; harmonicIndex < skinEffectHarmonicFrequencies.size(); ++harmonicIndex) {
            skinEffectLossesThisElement.get_mutable_losses_per_harmonic()[harmonicIndex] += skinEffectLossesPerHarmonic[harmonicIndex];
        }

        auto proximityEffectHarmonicFrequencies = windingLossesThisTurn.get_proximity_effect_losses()->get_harmonic_frequencies();
        auto proximityEffectLossesPerHarmonic = windingLossesThisTurn.get_proximity_effect_losses()->get_losses_per_harmonic();

        if (!windingLossesThisElement.get_proximity_effect_losses()) {
            proximityEffectLossesThisElement.set_method_used(windingLossesThisTurn.get_proximity_effect_losses()->get_method_used());
            proximityEffectLossesThisElement.set_origin(windingLossesThisTurn.get_proximity_effect_losses()->get_origin());
            proximityEffectLossesThisElement.set_harmonic_frequencies(proximityEffectHarmonicFrequencies);
            proximityEffectLossesThisElement.set_losses_per_harmonic(std::vector<double>(proximityEffectLossesPerHarmonic.size(), 0));
            windingLossesThisElement.set_proximity_effect_losses(proximityEffectLossesThisElement);
        }
        for (size_t harmonicIndex = 0; harmonicIndex < proximityEffectHarmonicFrequencies.size(); ++harmonicIndex) {
            proximityEffectLossesThisElement.get_mutable_losses_per_harmonic()[harmonicIndex] += proximityEffectLossesPerHarmonic[harmonicIndex];
        }

    }
    windingLossesThisElement.set_ohmic_losses(ohmicLossesThisElement);
    windingLossesThisElement.set_skin_effect_losses(skinEffectLossesThisElement);
    windingLossesThisElement.set_proximity_effect_losses(proximityEffectLossesThisElement);

    return windingLossesThisElement;
}

WindingLossesOutput combine_turn_losses(WindingLossesOutput windingLossesOutput, CoilWrapper coil) {
    auto windingLossesPerTurn = windingLossesOutput.get_winding_losses_per_turn().value();
    auto layers = coil.get_layers_description_conduction();

    std::vector<WindingLossesPerElement> windingLossesPerLayer;
    for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex)
    {
        OhmicLosses ohmicLossesThisLayer;
        WindingLossElement skinEffectLossesThisLayer;
        WindingLossElement proximityEffectLossesThisLayer;

        auto layer = layers[layerIndex];
        auto turnsIndexes = coil.get_turns_indexes_by_layer(layer.get_name());
        WindingLossesPerElement windingLossesThisLayer = combine_turn_losses_per_element(windingLossesPerTurn, turnsIndexes);
        windingLossesPerLayer.push_back(windingLossesThisLayer);
    }
    windingLossesOutput.set_winding_losses_per_layer(windingLossesPerLayer);

    auto sections = coil.get_sections_description_conduction();
    std::vector<WindingLossesPerElement> windingLossesPerSection;
    for (size_t sectionIndex = 0; sectionIndex < sections.size(); ++sectionIndex)
    {
        OhmicLosses ohmicLossesThisSection;
        WindingLossElement skinEffectLossesThisSection;
        WindingLossElement proximityEffectLossesThisSection;

        auto section = sections[sectionIndex];
        auto turnsIndexes = coil.get_turns_indexes_by_section(section.get_name());
        WindingLossesPerElement windingLossesThisSection = combine_turn_losses_per_element(windingLossesPerTurn, turnsIndexes);
        windingLossesPerSection.push_back(windingLossesThisSection);
    }
    windingLossesOutput.set_winding_losses_per_section(windingLossesPerSection);

    std::vector<WindingLossesPerElement> windingLossesPerWinding;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        OhmicLosses ohmicLossesThisWinding;
        WindingLossElement skinEffectLossesThisWinding;
        WindingLossElement proximityEffectLossesThisWinding;

        auto winding = coil.get_functional_description()[windingIndex];
        auto turnsIndexes = coil.get_turns_indexes_by_winding(winding.get_name());
        WindingLossesPerElement windingLossesThisWinding = combine_turn_losses_per_element(windingLossesPerTurn, turnsIndexes);
        windingLossesPerWinding.push_back(windingLossesThisWinding);
    }
    windingLossesOutput.set_winding_losses_per_winding(windingLossesPerWinding);
    return windingLossesOutput;
}

WindingLossesOutput WindingLosses::calculate_losses(MagneticWrapper magnetic, OperatingPoint operatingPoint, double temperature) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
    windingLossesOutput = WindingSkinEffectLosses::calculate_skin_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, settings->get_harmonic_amplitude_threshold());

    MagneticField magneticField(_magneticFieldStrengthModel, OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN);

    int64_t totalNumberPhysicalTurns = 0;
    for (auto winding : magnetic.get_coil().get_functional_description()) {
        totalNumberPhysicalTurns += winding.get_number_turns() * winding.get_number_parallels();
    }

    auto previoustHarmonicAmplitudeThreshold = settings->get_harmonic_amplitude_threshold();
    if (settings->get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings->set_harmonic_amplitude_threshold(settings->get_harmonic_amplitude_threshold() * 2);
    }
    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);

    windingLossesOutput = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, windingWindowMagneticStrengthFieldOutput);

    if (settings->get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings->set_harmonic_amplitude_threshold(previoustHarmonicAmplitudeThreshold);
    }
    windingLossesOutput = combine_turn_losses(windingLossesOutput, magnetic.get_coil());

    return windingLossesOutput;
}

ResistanceMatrixAtFrequency WindingLosses::calculate_resistance_matrix(MagneticWrapper magnetic, double temperature, double frequency) {
    ResistanceMatrixAtFrequency resistanceMatrixAtFrequency;
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();

    std::vector<std::vector<DimensionWithTolerance>> resistanceMatrix;
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        resistanceMatrix.push_back(std::vector<DimensionWithTolerance>{});
        for (size_t secondWindingIndex = 0; secondWindingIndex < magnetic.get_coil().get_functional_description().size(); ++secondWindingIndex) {
            DimensionWithTolerance dimension;
            resistanceMatrix[windingIndex].push_back(dimension);
        }
    }

    auto settings = OpenMagnetics::Settings::GetInstance();
    auto previousSetting = settings->get_magnetic_field_include_fringing();
    settings->set_magnetic_field_include_fringing(false);

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrent = 1;

    // diagonal calculation. mask => {1, 0, 0, ...}, {0, 1, 0, ...}, ...
    for (size_t enabledWindingIndex = 0; enabledWindingIndex < magnetic.get_coil().get_functional_description().size(); ++enabledWindingIndex) {
        std::vector<double> currentMask;
        for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
            if (enabledWindingIndex == windingIndex) {
                currentMask.push_back(virtualCurrent * sqrt(2));
            }
            else {
                currentMask.push_back(0);
            }
        }
        auto operatingPoint = InputsWrapper::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
        double totalLossses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        double effectiveResistance = totalLossses / pow(virtualCurrent, 2);

        resistanceMatrix[enabledWindingIndex][enabledWindingIndex].set_nominal(effectiveResistance);
    }

    // rest of elements calculation. mask => {1, 1, 0, ...}, {1, 0, 1, ...}, ...
    for (size_t enabledWindingIndex = 0; enabledWindingIndex < magnetic.get_coil().get_functional_description().size(); ++enabledWindingIndex) {
        for (size_t secondEnabledWindingIndex = enabledWindingIndex + 1; secondEnabledWindingIndex < magnetic.get_coil().get_functional_description().size(); ++secondEnabledWindingIndex) {
            std::vector<double> currentMask;
            double virtualCurrent = 1;
            for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
                if (enabledWindingIndex == windingIndex || secondEnabledWindingIndex == windingIndex) {
                    currentMask.push_back(virtualCurrent * sqrt(2));
                }
                else {
                    currentMask.push_back(0);
                }
            }
            auto operatingPoint = InputsWrapper::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
            double totalLossses =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

            double mutualResistance = (totalLossses - (resolve_dimensional_values(resistanceMatrix[enabledWindingIndex][enabledWindingIndex]) * pow(virtualCurrent, 2) + resolve_dimensional_values(resistanceMatrix[secondEnabledWindingIndex][secondEnabledWindingIndex]) * pow(virtualCurrent, 2))) / (2 * virtualCurrent * virtualCurrent);
            resistanceMatrix[enabledWindingIndex][secondEnabledWindingIndex].set_nominal(mutualResistance);
            resistanceMatrix[secondEnabledWindingIndex][enabledWindingIndex].set_nominal(mutualResistance);
        }
    }

    resistanceMatrixAtFrequency.set_frequency(frequency);
    resistanceMatrixAtFrequency.set_matrix(resistanceMatrix);

    settings->set_magnetic_field_include_fringing(previousSetting);
    return resistanceMatrixAtFrequency;
}

double WindingLosses::calculate_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
{
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first;
}

double WindingLosses::calculate_effective_resistance_per_meter(const WireWrapper& wire, double effectiveFrequency, double temperature)
{
    return WindingOhmicLosses::calculate_effective_resistance_per_meter(wire, effectiveFrequency, temperature);
}


double WindingLosses::calculate_skin_effect_resistance_per_meter(WireWrapper wire, SignalDescriptor current, double temperature)
{
    if (!current.get_processed()->get_rms()) {
        throw std::runtime_error("Current processed is missing field RMS");
    }
    auto currentRms = current.get_processed()->get_rms().value();
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first / pow(currentRms, 2);
}


} // namespace OpenMagnetics
