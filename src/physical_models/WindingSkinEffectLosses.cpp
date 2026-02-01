#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "Defaults.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include <functional>
#include "support/Exceptions.h"

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
    else if (modelName == WindingSkinEffectLossesModels::FERREIRA) {
        return std::make_shared<WindingSkinEffectLossesFerreiraModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::LOTFI) {
        return std::make_shared<WindingSkinEffectLossesLotfiModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::KUTKUT) {
        return std::make_shared<WindingSkinEffectLossesKutkutModel>();
    }
    else
        throw ModelNotAvailableException("Unknown wire skin effect losses mode, available options are: {DOWELL, WOJDA, ALBACH, PAYNE, NAN, VANDELAC_ZIOGAS, KAZIMIERCZUK, KUTKUT, FERREIRA, DIMITRAKAKIS, WANG, HOLGUIN, PERRY}");
}

std::shared_ptr<WindingSkinEffectLossesModel> WindingSkinEffectLosses::get_model(WireType wireType, std::optional<WindingSkinEffectLossesModels> modelOverride) {
    // If an explicit model override is provided, use it
    if (modelOverride.has_value()) {
        return WindingSkinEffectLossesModel::factory(modelOverride.value());
    }
    
    // Otherwise, auto-select based on wire type
    switch(wireType) {
        case WireType::ROUND: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        }
        case WireType::LITZ: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        }
        case WireType::PLANAR: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        }
        case WireType::RECTANGULAR:
        {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::KUTKUT);
        }
        case WireType::FOIL: {
            return WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::KUTKUT);
        }
        default:
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
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

double WindingSkinEffectLosses::calculate_skin_depth(Wire wire, double frequency, double temperature) {
    return calculate_skin_depth(wire.resolve_material(), frequency, temperature);
};

std::pair<double, std::vector<std::pair<double, double>>> WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(Wire wire, SignalDescriptor current, double temperature, double currentDivider, std::optional<WindingSkinEffectLossesModels> modelOverride) {
    double dcResistancePerMeter = WindingOhmicLosses::calculate_dc_resistance_per_meter(wire, temperature);

    auto model = get_model(wire.get_type(), modelOverride);

    if (!current.get_harmonics()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing harmonics");
    }
    auto harmonics = current.get_harmonics().value();

    double totalSkinEffectLossesPerMeter = 0;

    std::vector<std::pair<double, double>> lossesPerHarmonic;

    for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex)
    {
        auto harmonicRmsCurrent = harmonics.get_amplitudes()[harmonicIndex] / sqrt(2);  // Because a harmonic is always sinusoidal
        auto harmonicFrequency = harmonics.get_frequencies()[harmonicIndex];
        auto harmonicRmsCurrentInTurn = harmonicRmsCurrent * currentDivider;
        auto dcLossPerMeterThisHarmonic = pow(harmonicRmsCurrentInTurn, 2) * dcResistancePerMeter;

        auto turnLosses = model->calculate_turn_losses(wire, dcLossPerMeterThisHarmonic, harmonicFrequency, temperature, harmonicRmsCurrentInTurn);
        lossesPerHarmonic.push_back(std::pair<double, double>{turnLosses, harmonicFrequency});
        totalSkinEffectLossesPerMeter += turnLosses;
    }

    return {totalSkinEffectLossesPerMeter, lossesPerHarmonic};
}


WindingLossesOutput WindingSkinEffectLosses::calculate_skin_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, double windingLossesHarmonicAmplitudeThreshold, std::optional<WindingSkinEffectLossesModels> modelOverride) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");
    }
    auto turns = coil.get_turns_description().value();
    auto currentDividerPerTurn = windingLossesOutput.get_current_divider_per_turn().value();
    auto operatingPoint = windingLossesOutput.get_current_per_winding().value();
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform() || 
        operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform()->get_data().size() == 0)
    {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Input has no waveform. TODO: get waveform from processed data");
    }
    operatingPoint = Inputs::prune_harmonics(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    auto windingLossesPerTurn = windingLossesOutput.get_winding_losses_per_turn().value();

    double totalSkinEffectLosses = 0;

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex)
    {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = coil.resolve_wire(windingIndex);
        double wireLength = turn.get_length();

        SignalDescriptor current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();

        auto lossesPerMeterPerHarmonic = calculate_skin_effect_losses_per_meter(wire, current, temperature, currentDividerPerTurn[turnIndex], modelOverride).second;

        WindingLossElement skinEffectLosses;
        auto model = get_model(coil.get_wire_type(windingIndex), modelOverride);
        skinEffectLosses.set_method_used(model->methodName);
        skinEffectLosses.set_origin(ResultOrigin::SIMULATION);
        skinEffectLosses.get_mutable_harmonic_frequencies().push_back(0);
        skinEffectLosses.get_mutable_losses_per_harmonic().push_back(0);

        for (auto& lossesPerMeter : lossesPerMeterPerHarmonic) {
            skinEffectLosses.get_mutable_harmonic_frequencies().push_back(lossesPerMeter.second);
            skinEffectLosses.get_mutable_losses_per_harmonic().push_back(lossesPerMeter.first * wireLength);
            totalSkinEffectLosses += lossesPerMeter.first * wireLength;
        }
        windingLossesPerTurn[turnIndex].set_skin_effect_losses(skinEffectLosses);

    }
    windingLossesOutput.set_winding_losses_per_turn(windingLossesPerTurn);

    windingLossesOutput.set_method_used("AnalyticalModels");
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalSkinEffectLosses);
    return windingLossesOutput;
}

std::optional<double> WindingSkinEffectLossesModel::try_get_skin_factor(Wire wire, double frequency, double temperature) {
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


void WindingSkinEffectLossesModel::set_skin_factor(Wire wire,  double frequency, double temperature, double skinFactor) {
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

/**
 * @brief Calculates the penetration ratio for winding resistance calculation.
 * 
 * Based on: "Analytical Optimization of Litz-Wire Windings Independent of Porosity Factor"
 * by Rafal P. Wojda and Marian K. Kazimierczuk
 * IEEE Transactions on Industry Applications, Vol. 54, No. 5, 2018
 * https://doi.org/10.1109/TIA.2018.2821647
 * 
 * The penetration ratio A (or Astr for litz wire) relates the conductor diameter
 * to the skin depth, accounting for the porosity factor.
 * 
 * From Eq. (4) in the paper:
 * 
 *   A = (π/4)^(3/4) * (d/δw) * √(d/p)
 * 
 * where:
 * - d = conducting diameter of wire (or strand for litz)
 * - δw = skin depth at operating frequency
 * - p = outer diameter (accounting for insulation)
 * 
 * The skin depth δw is given by Eq. (5):
 * 
 *   δw = √(ρw / (πμ₀f))
 * 
 * For litz wire, the effective number of layers is Nll = Nl * √k (Eq. 6),
 * where k = number of strands per bundle.
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Penetration ratio A (dimensionless)
 */
double WindingSkinEffectLossesWojdaModel::calculate_penetration_ratio(const Wire& wire, double frequency, double temperature) {
    Wire realWire;
    if (wire.get_type() == WireType::LITZ) {
        realWire = Wire(Wire::resolve_strand(wire));
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
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Litz wire is missing strand information");

            auto strand = Wire::resolve_strand(wire);
            penetrationRatio = pow(std::numbers::pi / 4, 3 / 4) * resolve_dimensional_values(strand.get_conducting_diameter()) / skinDepth * sqrt(resolve_dimensional_values(strand.get_conducting_diameter()) / resolve_dimensional_values(strand.get_outer_diameter().value()));
            break;
        }
        case WireType::PLANAR:
        case WireType::RECTANGULAR:
        {
            penetrationRatio = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value())) / skinDepth;
            break;
        }
        case WireType::FOIL: {
            penetrationRatio = resolve_dimensional_values(wire.get_conducting_width().value()) / skinDepth;
            break;
        }
        default:
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    return penetrationRatio;
}

/**
 * @brief Calculates the skin effect factor using Wojda's formulation of Dowell's equation.
 * 
 * Based on: "Analytical Optimization of Litz-Wire Windings Independent of Porosity Factor"
 * by Rafal P. Wojda and Marian K. Kazimierczuk
 * IEEE Transactions on Industry Applications, Vol. 54, No. 5, 2018
 * https://doi.org/10.1109/TIA.2018.2821647
 * 
 * The skin effect factor Fs is derived from Dowell's equation, expressing the
 * AC-to-DC resistance ratio increase due to non-uniform current distribution.
 * 
 * From Eq. (7) - the adapted Dowell's equation for a single conductor:
 * 
 *   Fs = A * [sinh(2A) + sin(2A)] / [cosh(2A) - cos(2A)]
 * 
 * where A is the penetration ratio from calculate_penetration_ratio().
 * 
 * For low frequencies (A << 1): Fs → 1 (uniform current distribution)
 * For high frequencies (A >> 1): Fs → A (current concentrated near surface)
 * 
 * Note: This is the skin effect portion only. The proximity effect is handled
 * separately by the proximity losses model.
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Skin factor Fs (ratio ≥ 1, returns 1 when Fs = 1 meaning no skin effect)
 */
double WindingSkinEffectLossesWojdaModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    auto penetrationRatio = calculate_penetration_ratio(wire, frequency, temperature);
    auto factor = penetrationRatio / 2 * (sinh(penetrationRatio) + sin(penetrationRatio)) / (cosh(penetrationRatio) - cos(penetrationRatio));
    return factor;
}

double WindingSkinEffectLossesWojdaModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
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


/**
 * @brief Calculates the skin effect factor for AC resistance increase.
 * 
 * Based on: "Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre 
 * parasitären Eigenschaften" by Manfred Albach
 * Springer Vieweg, 2017, ISBN 978-3-658-15081-5
 * Chapter 4: "Die Verluste in Luftspulen", Section 4.1 "Rms- und Skinverluste", pages 72-79
 * 
 * The skin effect causes current to flow predominantly near the conductor surface
 * at high frequencies, effectively reducing the conducting cross-section and
 * increasing AC resistance.
 * 
 * Key equations from Chapter 4:
 * 
 * - Eq. (4.10): Skin losses for sinusoidal current
 *   Ps = Pdc * Ks  where Ks is the skin factor
 * 
 * - Eq. (4.7): Skin factor using Bessel functions (for round wire)
 *   Ks = (ξ/2) * Re[I0(αrD) / I1(αrD)]
 *   where:
 *   - ξ = rD/δ (ratio of wire radius to skin depth)
 *   - α = (1+j)/δ (complex propagation constant)
 *   - I0, I1 = modified Bessel functions of first kind
 * 
 * - Skin depth δ from Eq. (1.39):
 *   δ = √(2ρ/(ωμ)) = √(ρ/(πfμ))
 *   At 100 kHz: δ_Cu ≈ 0.21 mm, δ_Al ≈ 0.27 mm
 * 
 * - For bundled conductors (litz wire), includes proximity between strands term:
 *   Factor += n(n-1) * (rD/rO)² * I1(αrD) / I0(αrD)
 *   where n = number of conductors, rO = outer radius
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Skin factor Ks (ratio of AC to DC resistance due to skin effect)
 */
double WindingSkinEffectLossesAlbachModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double wireRadius;
    double wireOuterRadius;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        wireRadius = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value())) / 2;
        wireOuterRadius = std::min(resolve_dimensional_values(wire.get_outer_width().value()), resolve_dimensional_values(wire.get_outer_height().value())) / 2;
    }
    else if (wire.get_type() == WireType::ROUND) {
        wireRadius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        wireOuterRadius = resolve_dimensional_values(wire.get_outer_diameter().value()) / 2;
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        wireOuterRadius = resolve_dimensional_values(strand.get_outer_diameter().value()) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    std::complex<double> alpha(1, 1);
    alpha *= wireRadius / skinDepth;
    double factor = 0.5 * (alpha * (modified_bessel_first_kind(0, alpha) / modified_bessel_first_kind(1, alpha) + wire.get_number_conductors().value() * (wire.get_number_conductors().value() - 1) * pow(wireRadius, 2) / pow(wireOuterRadius, 2) * modified_bessel_first_kind(1, alpha * wireRadius) / modified_bessel_first_kind(0, alpha * wireRadius))).real();

    return factor;
}


double WindingSkinEffectLossesAlbachModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
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

/**
 * @brief Calculates AC skin effect losses for rectangular conductors using Payne's empirical model.
 * 
 * Based on: "The AC Resistance of Rectangular Conductors" by Alan Payne
 * Application Note AP101, Issue 4, 2021
 * Available at: https://www.paynelectronics.co.uk/
 * 
 * Payne's model provides an empirical approximation for the AC resistance factor of
 * rectangular cross-section conductors, valid across a wide frequency range.
 * 
 * The model calculates:
 * 1. The parameter p = √A / (1.26·δ·1000), where A is cross-section area in mm²
 * 2. A frequency factor: Ff = 1 - exp(-0.026·p)
 * 3. The corner correction: Kc = 1 + Ff·(1.2/exp(2.1·a/b) + 1.2/exp(2.1·b/a))
 * 4. The AC/DC resistance ratio through an empirical exponential formula
 * 
 * where:
 * - A = conductor cross-sectional area (mm²)
 * - δ = skin depth (mm)
 * - a = thin dimension (height or width, whichever is smaller)
 * - b = thick dimension
 * 
 * Note: This model uses dimensions in mm (not SI) as per Payne's original formulation.
 * 
 * @param wire Wire object containing conductor geometry
 * @param dcLossTurn DC power loss for the turn [W]
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @param currentRms RMS current (unused in this model)
 * @return Turn losses including skin effect [W]
 */
double WindingSkinEffectLossesPayneModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
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
    // double acResistanceFactor =  (Kc / (1.0 - exp(-x))) - 1.0;
    double acResistanceFactor =  (Kc / (1.0 - exp(-x)));
    auto turnLosses = dcLossTurn * acResistanceFactor;
    return turnLosses;
}



/**
 * @brief Calculates skin factor using Ferreira's 1D approximation.
 * 
 * Based on: "Improved Analytical Modeling of Conductive Losses in Magnetic Components"
 * by J.A. Ferreira
 * IEEE Transactions on Power Electronics, Vol. 9, No. 1, January 1994
 * https://doi.org/10.1109/63.285503
 * 
 * Ferreira's model provides a closed-form skin factor using hyperbolic/trigonometric
 * functions, which is computationally simpler than Bessel function approaches.
 * 
 * For a conductor of thickness h (or diameter for round wire):
 * 
 *   Fs = (ξ/4) * [sinh(ξ) + sin(ξ)] / [cosh(ξ) - cos(ξ)]
 * 
 * where:
 * - ξ = h / δ (ratio of conductor dimension to skin depth)
 * - δ = √(ρ/(πfμ)) = skin depth
 * 
 * This formula assumes a 1D current distribution and is accurate when the
 * conductor width >> thickness. For conductors with aspect ratio close to 1,
 * 2D effects become significant.
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Skin factor Fs (dimensionless)
 */
double WindingSkinEffectLossesFerreiraModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double wireHeight;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        wireHeight = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        wireHeight = resolve_dimensional_values(wire.get_conducting_diameter().value());
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double xi = wireHeight / skinDepth;
    double factor = xi / 4 * (sinh(xi) + sin(xi)) / (cosh(xi) - cos(xi));

    return factor;
}


double WindingSkinEffectLossesFerreiraModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);

    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    }
    else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }

    auto turnLosses = dcLossTurn * skinFactor;
    return turnLosses;
}

/**
 * @brief Calculates AC skin effect losses using Lotfi's elliptic integral method.
 * 
 * Based on: "Winding Loss Measurement and Calculation for Magnetic Components"
 * by A.W. Lotfi, S.J. Davis, I. Farhadi
 * Colorado Power Electronics Center, University of Colorado, Boulder
 * 
 * This model uses complete elliptic integrals of the first kind to calculate
 * the AC resistance of conductors. For rectangular conductors, the physical
 * dimensions are mapped to an equivalent ellipse using conformal transformation.
 * 
 * The transformation for rectangular cross-sections:
 *   b' = max(height, width) / 2
 *   a' = min(height, width) / 2
 *   b = 2·b' / √π  (effective semi-major axis)
 *   a = a'·b/b'    (effective semi-minor axis)
 * 
 * The AC resistance is:
 *   Rac = ρ / (π²·δ·b) · K(c/b) · (1 - exp(-2a/δ))
 * 
 * where:
 * - K(k) = complete elliptic integral of the first kind
 * - c = √(b² - a²)
 * - δ = skin depth
 * - ρ = resistivity
 * 
 * @param wire Wire object containing conductor geometry
 * @param dcLossTurn DC power loss for the turn (unused, recalculated)
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @param currentRms RMS current for resistance calculation
 * @return Turn losses including skin effect [W]
 */
double WindingSkinEffectLossesLotfiModel::calculate_turn_losses(Wire wire, [[maybe_unused]] double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double b, a;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        double bPrima = std::max(resolve_dimensional_values(wire.get_conducting_height().value()), resolve_dimensional_values(wire.get_conducting_width().value())) / 2;
        double aPrima = std::min(resolve_dimensional_values(wire.get_conducting_height().value()), resolve_dimensional_values(wire.get_conducting_width().value())) / 2;
        b = 2 * bPrima / sqrt(std::numbers::pi);
        a = aPrima * b / bPrima;
    }
    else if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        b = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        a = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double c = sqrt(pow(b, 2) - pow(a, 2));
    auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    double acResistance = resistivity / (pow(std::numbers::pi, 2) * skinDepth * b) * comp_ellint_1(c / b) * (1 - exp(-2 * a / skinDepth));

    auto turnLosses = acResistance * pow(currentRms / sqrt(2), 2);
    return turnLosses;
}


/**
 * @brief Calculates skin effect losses using Kutkut's 2D edge effects correction.
 * 
 * Based on: "A Simple Technique to Evaluate Winding Losses Including Two-Dimensional
 * Edge Effects" by Nasser H. Kutkut
 * IEEE Transactions on Power Electronics, Vol. 13, No. 5, September 1998
 * https://doi.org/10.1109/63.712308
 * 
 * This model accounts for 2D edge effects in foil/rectangular windings that are not
 * captured by 1D analysis. The current density increases near conductor edges due to
 * perpendicular field components.
 * 
 * The AC resistance is modeled as a combination of low and high frequency asymptotes:
 * 
 * From Eq. (30) - combined asymptotic solution:
 * 
 *   Fr = [1 + (f/fl)^α + (f/fh)^β]^(1/γ)
 * 
 * where:
 * - fl = 3.22*ρ / (8*μ₀*b'*a') - low frequency corner (Eq. 31)
 * - fh = π²*ρ / (4*μ₀*a'²) * K(c/b)^(-2) - high frequency corner (Eq. 32)
 * - K = complete elliptic integral of first kind
 * - a', b' = semi-axes of equivalent ellipse (half of conductor dimensions)
 * - α = 2, β = 5.5, γ = 11 (empirical exponents)
 * 
 * At low frequency: Fr → 1 (DC behavior)
 * At high frequency: Fr → (f/fh)^(β/γ) ≈ f^0.5 (skin effect dominated)
 * 
 * The 2D correction factor (Eq. 35) can increase losses by 85% compared to 1D analysis.
 * 
 * @param wire Wire object (foil or rectangular conductor)
 * @param dcLossTurn DC losses for this turn [W]
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @param currentRms RMS current through the turn [A]
 * @return Additional skin effect losses [W]
 */
double WindingSkinEffectLossesKutkutModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]]double currentRms) {
    double bPrima, aPrima;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        bPrima = std::max(resolve_dimensional_values(wire.get_conducting_height().value()), resolve_dimensional_values(wire.get_conducting_width().value())) / 2;
        aPrima = std::min(resolve_dimensional_values(wire.get_conducting_height().value()), resolve_dimensional_values(wire.get_conducting_width().value())) / 2;
    }
    else if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        bPrima = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        aPrima = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }
    auto resistivityModel = OpenMagnetics::ResistivityModel::factory(OpenMagnetics::ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    double fl = 3.22 * resistivity / (8 * Constants().vacuumPermeability * bPrima * aPrima);
    double fh = pow(std::numbers::pi, 2) * resistivity / (4 * Constants().vacuumPermeability * pow(aPrima, 2)) * pow(comp_ellint_1(sqrt(1 - pow(aPrima, 2) / pow(bPrima, 2))), -2);
    double gamma = 11;
    double beta = 5.5;
    double alpha = 2;

    double acResistanceFactor = pow(1 + pow(frequency / fl, alpha) + pow(frequency / fh, beta), 1 / gamma);

    auto turnLosses = (acResistanceFactor - 1) * dcLossTurn;
    return turnLosses;
}

// double get_effective_current_density(Harmonics harmonics, Wire wire, double temperature) {
//     for (size_t i=0; i<harmonics.get_amplitudes(); ++i) {
//         auto skinDepth = calculate_skin_depth(material, frequency, temperature);
//     }
// }

} // namespace OpenMagnetics

