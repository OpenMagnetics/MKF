#include "support/Settings.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/Inductance.h"
#include "MAS.hpp"
#include "physical_models/MagnetizingInductance.h"
#include "support/Exceptions.h"

namespace OpenMagnetics {

double WindingLosses::get_total_winding_losses(WindingLossesPerElement windingLossesThisElement) {
    double totalLossses = 0;
    if (windingLossesThisElement.get_ohmic_losses()) {
        totalLossses += windingLossesThisElement.get_ohmic_losses()->get_losses();
    }

    if (windingLossesThisElement.get_skin_effect_losses()) {
        auto skinEffectLossesPerHarmonic = windingLossesThisElement.get_skin_effect_losses()->get_losses_per_harmonic();
        for (size_t harmonicIndex = 0; harmonicIndex < skinEffectLossesPerHarmonic.size(); ++harmonicIndex) {
            totalLossses += skinEffectLossesPerHarmonic[harmonicIndex];
        }
    }

    if (windingLossesThisElement.get_proximity_effect_losses()) {
        auto proximityEffectLossesPerHarmonic = windingLossesThisElement.get_proximity_effect_losses()->get_losses_per_harmonic();
        for (size_t harmonicIndex = 0; harmonicIndex < proximityEffectLossesPerHarmonic.size(); ++harmonicIndex) {
            totalLossses += proximityEffectLossesPerHarmonic[harmonicIndex];
        }
    }
    return totalLossses;
}

double WindingLosses::get_total_winding_losses(std::vector<WindingLossesPerElement> windingLossesPerElement) {
    double totalLossses = 0;
    for (auto windingLossesThisElement : windingLossesPerElement) {
        totalLossses += WindingLosses::get_total_winding_losses(windingLossesThisElement);
    }
    return totalLossses;
}

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

WindingLossesOutput WindingLosses::combine_turn_losses(WindingLossesOutput windingLossesOutput, Coil coil) {
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
        windingLossesThisLayer.set_name(layer.get_name());
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
        windingLossesThisSection.set_name(section.get_name());
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
        windingLossesThisWinding.set_name(winding.get_name());
        windingLossesPerWinding.push_back(windingLossesThisWinding);
    }
    windingLossesOutput.set_winding_losses_per_winding(windingLossesPerWinding);
    return windingLossesOutput;
}

double WindingLosses::calculate_effective_resistance_of_winding(Magnetic magnetic, size_t windingIndex, double frequency, double temperature) {
    auto magnetizingInductanceModel = MagnetizingInductance();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    double virtualCurrentRms = 1;
    std::vector<double> currentMask = {virtualCurrentRms * sqrt(2)};
    for (auto turnsRatio : turnsRatios) {
        currentMask.push_back(virtualCurrentRms * sqrt(2) * turnsRatio);
    }

    auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
    auto windingLossesPerWinding =  WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses_per_winding().value();

    auto proximityLossesPerharmonic = windingLossesPerWinding[windingIndex].get_proximity_effect_losses()->get_losses_per_harmonic();
    double proximityLosses = std::accumulate(proximityLossesPerharmonic.begin(), proximityLossesPerharmonic.end(), 0.0);
    auto skinLossesPerharmonic = windingLossesPerWinding[windingIndex].get_skin_effect_losses()->get_losses_per_harmonic();
    double skinLosses = std::accumulate(skinLossesPerharmonic.begin(), skinLossesPerharmonic.end(), 0.0);
    double lossesThisWinding = windingLossesPerWinding[windingIndex].get_ohmic_losses()->get_losses() + proximityLosses + skinLosses;

    double effectiveResistance = lossesThisWinding / pow(operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_processed()->get_rms().value(), 2);

    return effectiveResistance;
}


WindingLossesOutput WindingLosses::calculate_losses(Magnetic magnetic, OperatingPoint operatingPoint, double temperature) {
    auto& settings = OpenMagnetics::Settings::GetInstance();

    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(magnetic.get_coil(), operatingPoint, temperature);
    windingLossesOutput = WindingSkinEffectLosses::calculate_skin_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, settings.get_harmonic_amplitude_threshold(), _models.skinEffectModel);

    auto isPlanar = true;
    auto isRectangular = true;
    auto isFoil = true;
    for (auto wire : magnetic.get_mutable_coil().get_wires()) {
        if (wire.get_type() != WireType::PLANAR) {
            isPlanar = false;
            break;
        }
    }
    for (auto wire : magnetic.get_mutable_coil().get_wires()) {
        if (wire.get_type() != WireType::RECTANGULAR) {
            isRectangular = false;
            break;
        }
    }
    for (auto wire : magnetic.get_mutable_coil().get_wires()) {
        if (wire.get_type() != WireType::FOIL) {
            isFoil = false;
            break;
        }
    }
    // Only auto-select model if not explicitly specified
    MagneticFieldStrengthModels modelToUse;
    if (_models.magneticFieldStrengthModel.has_value()) {
        modelToUse = _models.magneticFieldStrengthModel.value();
    }
    else if (isPlanar) {
        modelToUse = OpenMagnetics::MagneticFieldStrengthModels::WANG;
    }
    else if (isRectangular) {
        modelToUse = OpenMagnetics::MagneticFieldStrengthModels::BINNS_LAWRENSON;
    }
    else if (isFoil) {
        modelToUse = OpenMagnetics::MagneticFieldStrengthModels::WANG;
    }
    else {
        modelToUse = OpenMagnetics::MagneticFieldStrengthModels::LAMMERANER;
    }
    
    // Use fringing model from struct if specified, otherwise use default ROSHEN
    MagneticFieldStrengthFringingEffectModels fringingModelToUse = 
        _models.magneticFieldStrengthFringingEffectModel.value_or(OpenMagnetics::MagneticFieldStrengthFringingEffectModels::ROSHEN);
    MagneticField magneticField(modelToUse, fringingModelToUse);

    int64_t totalNumberPhysicalTurns = 0;
    for (auto winding : magnetic.get_coil().get_functional_description()) {
        totalNumberPhysicalTurns += winding.get_number_turns() * winding.get_number_parallels();
    }

    auto previoustHarmonicAmplitudeThreshold = settings.get_harmonic_amplitude_threshold();
    if (settings.get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings.set_harmonic_amplitude_threshold(settings.get_harmonic_amplitude_threshold() * 2);
    }
    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic);

    windingLossesOutput = WindingProximityEffectLosses::calculate_proximity_effect_losses(magnetic.get_coil(), temperature, windingLossesOutput, windingWindowMagneticStrengthFieldOutput, _models.proximityEffectModel);

    if (settings.get_harmonic_amplitude_threshold_quick_mode() && totalNumberPhysicalTurns > _quickModeForManyTurnsThreshold) {
        settings.set_harmonic_amplitude_threshold(previoustHarmonicAmplitudeThreshold);
    }
    windingLossesOutput = combine_turn_losses(windingLossesOutput, magnetic.get_coil());

    return windingLossesOutput;
}

ScalarMatrixAtFrequency WindingLosses::calculate_resistance_matrix(Magnetic magnetic, double temperature, double frequency) {
    ScalarMatrixAtFrequency scalarMatrixAtFrequency;
    scalarMatrixAtFrequency.set_frequency(frequency);
    
    auto& functionalDescription = magnetic.get_coil().get_functional_description();
    size_t numWindings = functionalDescription.size();
    auto turnsRatios = magnetic.get_mutable_coil().get_turns_ratios();

    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto previousSetting = settings.get_magnetic_field_include_fringing();
    settings.set_magnetic_field_include_fringing(false);

    auto magnetizingInductanceModel = MagnetizingInductance();
    auto magnetizingInductance = resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance());
    
    // Calculate inductance matrix to get sqrt(L_i/L_j) ratios for mutual resistance scaling
    // Per Spreen (1990), mutual resistance R_ij should be scaled using inductance ratio, not turns ratio
    Inductance inductanceModel;
    auto inductanceMatrix = inductanceModel.calculate_inductance_matrix(magnetic, frequency);
    auto& inductanceMagnitude = inductanceMatrix.get_magnitude();
    
    // Store self-inductances for later use in scaling
    std::vector<double> selfInductances(numWindings);
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName = functionalDescription[i].get_name();
        selfInductances[i] = inductanceMagnitude.at(windingName).at(windingName).get_nominal().value();
    }
    
    double virtualCurrent = 1;
    std::map<std::string, std::map<std::string, DimensionWithTolerance>> resistanceMagnitude;

    // Diagonal calculation: R_ii from losses with only winding i excited
    // Current mask => {I, 0, 0, ...}, {0, I, 0, ...}, ...
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName_i = functionalDescription[i].get_name();
        
        std::vector<double> currentMask(numWindings, 0.0);
        currentMask[i] = virtualCurrent * sqrt(2);
        
        auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
        double totalLosses = WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

        // R_ii = P / I^2
        double effectiveResistance = totalLosses / pow(virtualCurrent, 2);

        DimensionWithTolerance resistanceValue;
        resistanceValue.set_nominal(effectiveResistance);
        resistanceMagnitude[windingName_i][windingName_i] = resistanceValue;
    }

    // Off-diagonal calculation: R_ij from losses with windings i and j excited
    // Per Spreen (1990): P_total = R_11*I_1^2 + R_22*I_2^2 + 2*R_12*I_1*I_2
    // Using I_2 = sqrt(L_1/L_2) * I_1 for proper scaling based on inductance ratio
    for (size_t i = 0; i < numWindings; ++i) {
        std::string windingName_i = functionalDescription[i].get_name();
        
        for (size_t j = i + 1; j < numWindings; ++j) {
            std::string windingName_j = functionalDescription[j].get_name();
            
            // Use sqrt(L_i/L_j) as the current ratio, per Spreen (1990)
            // This properly accounts for the inductance-based coupling
            double inductanceRatio = std::sqrt(selfInductances[i] / selfInductances[j]);
            
            double I_i = virtualCurrent;
            double I_j = virtualCurrent * inductanceRatio;
            
            std::vector<double> currentMask(numWindings, 0.0);
            currentMask[i] = I_i * sqrt(2);
            currentMask[j] = I_j * sqrt(2);
            
            auto operatingPoint = Inputs::create_operating_point_with_sinusoidal_current_mask(frequency, magnetizingInductance, temperature, turnsRatios, currentMask);
            double totalLosses = WindingLosses().calculate_losses(magnetic, operatingPoint, temperature).get_winding_losses();

            // P_total = R_ii * I_i^2 + R_jj * I_j^2 + 2 * R_ij * I_i * I_j
            // Solving for R_ij:
            // R_ij = (P_total - R_ii * I_i^2 - R_jj * I_j^2) / (2 * I_i * I_j)
            double R_ii = resistanceMagnitude[windingName_i][windingName_i].get_nominal().value();
            double R_jj = resistanceMagnitude[windingName_j][windingName_j].get_nominal().value();
            
            double mutualResistance = (totalLosses - R_ii * pow(I_i, 2) - R_jj * pow(I_j, 2)) / (2 * I_i * I_j);
            
            DimensionWithTolerance resistanceValue;
            resistanceValue.set_nominal(mutualResistance);
            
            // Symmetric matrix: R_ij = R_ji
            resistanceMagnitude[windingName_i][windingName_j] = resistanceValue;
            resistanceMagnitude[windingName_j][windingName_i] = resistanceValue;
        }
    }

    scalarMatrixAtFrequency.set_magnitude(resistanceMagnitude);
    settings.set_magnetic_field_include_fringing(previousSetting);
    return scalarMatrixAtFrequency;
}

double WindingLosses::calculate_losses_per_meter(Wire wire, SignalDescriptor current, double temperature)
{
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first;
}

double WindingLosses::calculate_effective_resistance_per_meter(const Wire& wire, double effectiveFrequency, double temperature)
{
    return WindingOhmicLosses::calculate_effective_resistance_per_meter(wire, effectiveFrequency, temperature);
}


double WindingLosses::calculate_skin_effect_resistance_per_meter(Wire wire, SignalDescriptor current, double temperature)
{
    if (!current.get_processed()->get_rms()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Current processed is missing field RMS");
    }
    auto currentRms = current.get_processed()->get_rms().value();
    return WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first / pow(currentRms, 2);
}


} // namespace OpenMagnetics
