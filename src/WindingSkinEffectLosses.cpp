#include "WindingSkinEffectLosses.h"
#include "WindingOhmicLosses.h"
#include "Defaults.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include <functional>

namespace OpenMagnetics {

std::shared_ptr<WindingSkinEffectLossesModel>  WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels modelName){
    if (modelName == WindingSkinEffectLossesModels::WOJDA) {
        return std::make_shared<WindingSkinEffectLossesWojdaModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::ALBACH) {
        return std::make_shared<WindingSkinEffectLossesAlbachModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::PAYNE) {
        return std::make_shared<WindingSkinEffectLossesPayneModel>();
    }
    else
        throw std::runtime_error("Unknown wire skin effect losses mode, available options are: {DOWELL, WOJDA, ALBACH, PAYNE, NAN, VANDELAC_ZIOGAS, KAZIMIERCZUK, KUTKUT, FERREIRA, DIMITRAKAKIS, WANG, HOLGUIN, PERRY}");
}

std::shared_ptr<WindingSkinEffectLossesModel> WindingSkinEffectLosses::get_model(WireType wireType) {
    switch(wireType) {
        case WireType::ROUND: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        }
        case WireType::LITZ: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::WOJDA);
        }
        case WireType::RECTANGULAR: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::PAYNE);
        }
        case WireType::FOIL: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::PAYNE);
        }
        default:
            throw std::runtime_error("Unknown type of wire");
    }
}

double WindingSkinEffectLosses::calculate_skin_depth(WireMaterialDataOrNameUnion material, double frequency, double temperature) {
    WireMaterial wireMaterial;
    auto constants = Constants();
    if (std::holds_alternative<WireMaterial>(material)) {
        wireMaterial = std::get<WireMaterial>(material);
    }
    else if (std::holds_alternative<std::string>(material)) {
        wireMaterial = OpenMagnetics::find_wire_material_by_name(std::get<std::string>(material));
    }
    auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

    double skinDepth = sqrt(resistivity / (std::numbers::pi * frequency * constants.vacuumPermeability * wireMaterial.get_permeability()));
    return skinDepth;
};

double WindingSkinEffectLosses::calculate_skin_depth(WireWrapper wire, double frequency, double temperature) {
    return calculate_skin_depth(wire.resolve_material(), frequency, temperature);
};

std::pair<double, std::vector<std::pair<double, double>>> WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(WireWrapper wire, SignalDescriptor current, double temperature, double currentDivider) {
    auto defaults = Defaults();

    double dcResistancePerMeter = WindingOhmicLosses::calculate_dc_resistance_per_meter(wire, temperature);

    auto model = get_model(wire.get_type());

    if (!current.get_harmonics()) {
        throw std::runtime_error("Current is missing harmonics");
    }
    auto harmonics = current.get_harmonics().value();

    double totalSkinEffectLossesPerMeter = 0;
    double maximumHarmonicAmplitudeTimesRootFrequency = 0;
    for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
        if (harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]) > maximumHarmonicAmplitudeTimesRootFrequency) {
            maximumHarmonicAmplitudeTimesRootFrequency = harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]);
        }
    }

    std::vector<std::pair<double, double>> lossesPerHarmonic;

    for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex)
    {
        if ((harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex])) < maximumHarmonicAmplitudeTimesRootFrequency * defaults.windingLossesHarmonicAmplitudeThreshold) {
            continue;
        }

        auto harmonicRmsCurrent = harmonics.get_amplitudes()[harmonicIndex] / sqrt(2);  // Because a harmonic is always sinusoidal
        auto harmonicFrequency = harmonics.get_frequencies()[harmonicIndex];
        auto harmonicRmsCurrentInTurn = harmonicRmsCurrent * currentDivider;
        auto dcLossPerMeterThisHarmonic = pow(harmonicRmsCurrentInTurn, 2) * dcResistancePerMeter;

        auto turnLosses = model->calculate_turn_losses(wire, dcLossPerMeterThisHarmonic, harmonicFrequency, temperature);
        lossesPerHarmonic.push_back(std::pair<double, double>{turnLosses, harmonicFrequency});
        totalSkinEffectLossesPerMeter += turnLosses;
    }

    return {totalSkinEffectLossesPerMeter, lossesPerHarmonic};
}


WindingLossesOutput WindingSkinEffectLosses::calculate_skin_effect_losses(CoilWrapper coil, double temperature, WindingLossesOutput windingLossesOutput) {
    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding does not have turns description");
    }
    auto turns = coil.get_turns_description().value();
    auto currentDividerPerTurn = windingLossesOutput.get_current_divider_per_turn().value();
    auto operatingPoint = windingLossesOutput.get_current_per_winding().value();
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform() || 
        operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform()->get_data().size() == 0)
    {
        throw std::runtime_error("Input has no waveform. TODO: get waveform from processed data");
    }

    std::vector<std::shared_ptr<WindingSkinEffectLossesModel>> lossesModelPerWinding;

    auto windingLossesPerWinding = windingLossesOutput.get_winding_losses_per_winding().value();

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex)
    {
        auto model = get_model(coil.get_wire_type(windingIndex));
        lossesModelPerWinding.push_back(model);

        WindingLossElement skinEffectLosses;
        skinEffectLosses.set_method_used(model->method_name);
        skinEffectLosses.set_origin(ResultOrigin::SIMULATION);
        skinEffectLosses.get_mutable_harmonic_frequencies().push_back(0);
        skinEffectLosses.get_mutable_losses_per_harmonic().push_back(0);

        windingLossesPerWinding[windingIndex].set_skin_effect_losses(skinEffectLosses);
    }
    double totalSkinEffectLosses = 0;

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex)
    {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = coil.resolve_wire(windingIndex);
        double wireLength = turn.get_length();

        SignalDescriptor current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();

        auto lossesPerMeterPerHarmonic = calculate_skin_effect_losses_per_meter(wire, current, temperature, currentDividerPerTurn[turnIndex]).second;

        for (auto& lossesPerMeter : lossesPerMeterPerHarmonic) {
            auto skinEffectLosses = windingLossesPerWinding[windingIndex].get_skin_effect_losses().value();
            skinEffectLosses.get_mutable_harmonic_frequencies().push_back(lossesPerMeter.second);
            skinEffectLosses.get_mutable_losses_per_harmonic().push_back(lossesPerMeter.first * wireLength);
            totalSkinEffectLosses += lossesPerMeter.first * wireLength;
            windingLossesPerWinding[windingIndex].set_skin_effect_losses(skinEffectLosses);
        }

    }
    // throw std::runtime_error("totalSkinEffectLosses: " + Convert.ToString(totalSkinEffectLosses));
    windingLossesOutput.set_winding_losses_per_winding(windingLossesPerWinding);

    windingLossesOutput.set_method_used("AnalyticalModels");
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalSkinEffectLosses);
    return windingLossesOutput;
}

std::optional<double> WindingSkinEffectLossesModel::try_get_skin_factor(WireWrapper wire, double frequency, double temperature) {
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    std::size_t hash;
    if (!wire.get_name()) {
        hash = std::hash<double>{}(wire.get_number_conductors().value() * wire.get_maximum_outer_width() * wire.get_maximum_outer_height() );
    }
    else {
        hash = std::hash<std::string>{}(wire.get_name().value());
    }

    if (_skinFactorPerWirePerFrequencyPerTemperature.contains(hash)) {
        if (_skinFactorPerWirePerFrequencyPerTemperature[hash].contains(frequency)) {
            if (_skinFactorPerWirePerFrequencyPerTemperature[hash][frequency].contains(temperature)) {
                return _skinFactorPerWirePerFrequencyPerTemperature[hash][frequency][temperature];
            }
        }
    }

    return std::nullopt;
}


void WindingSkinEffectLossesModel::set_skin_factor(WireWrapper wire,  double frequency, double temperature, double skinFactor) {
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    std::size_t hash;
    if (!wire.get_name()) {
        hash = std::hash<double>{}(wire.get_number_conductors().value() * wire.get_maximum_outer_width() * wire.get_maximum_outer_height() );
    }
    else {
        hash = std::hash<std::string>{}(wire.get_name().value());
    }

    _skinFactorPerWirePerFrequencyPerTemperature[hash][frequency][temperature] = skinFactor;

}

double WindingSkinEffectLossesWojdaModel::calculate_penetration_ratio(WireWrapper wire, double frequency, double temperature) {
    Wire realWire;
    if (wire.get_type() == WireType::LITZ) {
        realWire = WireWrapper::resolve_strand(wire);
    }
    else {
        realWire = wire;
    }

    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(realWire, frequency, temperature);

    double penetrationRatio;
    switch(wire.get_type()) {
        case WireType::ROUND: {
            penetrationRatio = pow(std::numbers::pi / 4, 3 / 4) * resolve_dimensional_values(wire.get_conducting_diameter().value()) / skinDepth * sqrt(resolve_dimensional_values(wire.get_conducting_diameter().value()) / resolve_dimensional_values(wire.get_outer_diameter().value()));
            break;
        }
        case WireType::LITZ: {
            if (!wire.get_strand())
                throw std::runtime_error("Litz wire is missing strand information");

            auto strand = WireWrapper::resolve_strand(wire);
            penetrationRatio = pow(std::numbers::pi / 4, 3 / 4) * resolve_dimensional_values(strand.get_conducting_diameter().value()) / skinDepth * sqrt(resolve_dimensional_values(strand.get_conducting_diameter().value()) / resolve_dimensional_values(strand.get_outer_diameter().value()));
            break;
        }
        case WireType::RECTANGULAR: {
            penetrationRatio = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value())) / skinDepth;
            break;
        }
        case WireType::FOIL: {
            penetrationRatio = resolve_dimensional_values(wire.get_conducting_width().value()) / skinDepth;
            break;
        }
        default:
            throw std::runtime_error("Unknown type of wire");
    }

    return penetrationRatio;
}

double WindingSkinEffectLossesWojdaModel::calculate_skin_factor(WireWrapper wire, double frequency, double temperature) {
    auto penetrationRatio = calculate_penetration_ratio(wire, frequency, temperature);
    auto factor = penetrationRatio / 2 * (sinh(penetrationRatio) + sin(penetrationRatio)) / (cosh(penetrationRatio) - cos(penetrationRatio));
    return factor;
}

double WindingSkinEffectLossesWojdaModel::calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);

    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    }
    else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }

    auto turnLosses = dcLossTurn * (skinFactor - 1);
    return turnLosses;
}


double WindingSkinEffectLossesAlbachModel::calculate_skin_factor(WireWrapper wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double wireRadius;
    double wireOuterRadius;

    if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        wireRadius = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value())) / 2;
        wireOuterRadius = std::min(resolve_dimensional_values(wire.get_outer_width().value()), resolve_dimensional_values(wire.get_outer_height().value())) / 2;
    }
    else if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        wireRadius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        wireOuterRadius = resolve_dimensional_values(wire.get_outer_diameter().value()) / 2;
    }
    else {
        throw std::runtime_error("Unknown type of wire");
    }

    std::complex<double> alpha(1, 1);
    alpha *= wireRadius / skinDepth;
    double factor = 0.5 * (alpha * (modified_bessel_first_kind(0, alpha) / modified_bessel_first_kind(1, alpha) + wire.get_number_conductors().value() * (wire.get_number_conductors().value() - 1) * pow(wireRadius, 2) / pow(wireOuterRadius, 2) * modified_bessel_first_kind(1, alpha * wireRadius) / modified_bessel_first_kind(0, alpha * wireRadius))).real();

    return factor;
}


double WindingSkinEffectLossesAlbachModel::calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);

    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    }
    else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }

    auto turnLosses = dcLossTurn * (skinFactor - 1);
    return turnLosses;
}

double WindingSkinEffectLossesPayneModel::calculate_turn_losses(WireWrapper wire, double dcLossTurn, double frequency, double temperature) {
        double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
        double thinDimension;
        double thickDimension;

        if (resolve_dimensional_values(wire.get_conducting_height().value()) > resolve_dimensional_values(wire.get_conducting_width().value())) {
            thinDimension = resolve_dimensional_values(wire.get_conducting_width().value());
            thickDimension = resolve_dimensional_values(wire.get_conducting_height().value());
        }
        else {
            thinDimension = resolve_dimensional_values(wire.get_conducting_height().value());
            thickDimension = resolve_dimensional_values(wire.get_conducting_width().value());
        }
        double A = resolve_dimensional_values(wire.get_conducting_width().value()) * resolve_dimensional_values(wire.get_conducting_height().value()) * 1000000;  // en mm2, shame on you Payne...

        double p = pow(A, 0.5) / (1.26 * skinDepth * 1000);
        double Ff = 1.0 - exp(-0.026 * p);
        double Kc;
        try {
            Kc = 1.0 + Ff * (1.2 / exp(2.1 * thickDimension / thinDimension) + 1.2 / exp(2.1 * thinDimension / thickDimension));
        }
        catch (...) {
            Kc = 1.0;
        }
        double x = (2.0 * skinDepth / thickDimension * (1.0 + thickDimension / thinDimension) + 8.0 * pow(skinDepth / thickDimension, 3) / (thinDimension / thickDimension)) / (pow(thinDimension / thickDimension, 0.33) * exp(-3.5 * thickDimension / skinDepth) + 1.0);
        double acResistanceFactor =  (Kc / (1.0 - exp(-x))) - 1.0;
        auto turnLosses = dcLossTurn * acResistanceFactor;
        return turnLosses;
}

// double get_effective_current_density(Harmonics harmonics, WireWrapper wire, double temperature) {
//     for (size_t i=0; i<harmonics.get_amplitudes(); ++i) {
//         auto skinDepth = calculate_skin_depth(material, frequency, temperature);
//     }
// }

} // namespace OpenMagnetics

