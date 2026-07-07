#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "constructive_models/Wire.h"
#include "Defaults.h"

#include <magic_enum.hpp>
#include <cmath>
#include <limits>
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

// Effective solid-conductor diameter of a round/litz wire. Litz wires carry NO
// wire-level conducting_diameter (it lives on the strand), so several skin/proximity
// models that shared a single `ROUND || LITZ` branch and read
// wire.get_conducting_diameter().value() threw std::bad_optional_access (opaque, no
// context) on every litz turn. Resolve the strand for litz; use the wire's own
// diameter for solid round; throw a specific error if it is genuinely absent.
inline double resolve_conducting_diameter(const Wire& wire) {
    if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        return resolve_dimensional_values(strand.get_conducting_diameter());
    }
    if (!wire.get_conducting_diameter()) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting diameter in round wire");
    }
    return resolve_dimensional_values(wire.get_conducting_diameter().value());
}

// PERF-003: Cached ResistivityModel to avoid repeated factory() calls
inline std::shared_ptr<ResistivityModel>& get_cached_resistivity_model() {
    static std::shared_ptr<ResistivityModel> m = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    return m;
}
// PERF-002: Composite hash to avoid collisions when wire has no name
inline std::size_t hash_combine_wire(int64_t n, double w, double h) {
    std::size_t seed = std::hash<int64_t>{}(n);
    seed ^= std::hash<double>{}(w) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= std::hash<double>{}(h) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

std::shared_ptr<WindingSkinEffectLossesModel> WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels modelName){
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
    else if (modelName == WindingSkinEffectLossesModels::DOWELL) {
        return std::make_shared<WindingSkinEffectLossesDowellModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::PERRY) {
        return std::make_shared<WindingSkinEffectLossesPerryModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::XI_NAN) {
        return std::make_shared<WindingSkinEffectLossesNanModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::KAZIMIERCZUK) {
        return std::make_shared<WindingSkinEffectLossesKazimierczukModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::DIMITRAKAKIS) {
        return std::make_shared<WindingSkinEffectLossesDimitrakakisModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::MUEHLETHALER) {
        return std::make_shared<WindingSkinEffectLossesMuehlethalerModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::WANG) {
        return std::make_shared<WindingSkinEffectLossesWangModel>();
    }
    else if (modelName == WindingSkinEffectLossesModels::HOLGUIN) {
        return std::make_shared<WindingSkinEffectLossesHolguinModel>();
    }
    else
        throw ModelNotAvailableException("Unknown wire skin effect losses mode, available options are: {DOWELL, WOJDA, ALBACH, PAYNE, NAN, KAZIMIERCZUK, KUTKUT, FERREIRA, LOTFI, DIMITRAKAKIS, MUEHLETHALER, WANG, HOLGUIN, PERRY}");
}


std::shared_ptr<WindingSkinEffectLossesModel> WindingSkinEffectLosses::get_model(WireType wireType, std::optional<WindingSkinEffectLossesModels> modelOverride) {
    // Resolve which model to use (unchanged selection logic): explicit override
    // first, then the user-selected model if manual selection is enabled, else
    // the wire-type default.
    WindingSkinEffectLossesModels resolvedModel;
    if (modelOverride.has_value()) {
        resolvedModel = modelOverride.value();
    }
    else if (Settings::GetInstance().get_coil_enable_user_winding_losses_models()) {
        resolvedModel = Settings::GetInstance().get_winding_skin_effect_losses_model();
    }
    else {
        switch(wireType) {
            case WireType::ROUND:
            case WireType::LITZ:
            case WireType::PLANAR:
                resolvedModel = WindingSkinEffectLossesModels::ALBACH;
                break;
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                resolvedModel = WindingSkinEffectLossesModels::KUTKUT;
                break;
            default:
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
        }
    }

    // Reuse one persistent instance per resolved model so its skin-factor cache
    // (_skinFactorPerWirePerFrequencyPerTemperature) survives across calls.
    // calculate_skin_factor (wire, frequency, temperature) is a pure, deterministic
    // function and the models hold no other mutable state, so a process-lifetime
    // per-model cache is correct and keeps each model's results separate.
    //
    // Previously get_model built a fresh model on every call, discarding the cache,
    // so the (expensive, Bessel-based) skin factor was recomputed for every wire ×
    // harmonic on every wire-advise pass. Profiled with callgrind: that path was
    // ~45% of CoreAdviser time. Single-threaded per the library's documented
    // threading contract (see Utils.h), consistent with the other model caches.
    // ABT #113: thread_local — the persistent models hold mutable per-wire
    // skin-factor caches, so they must not be shared between threads.
    static thread_local std::map<WindingSkinEffectLossesModels, std::shared_ptr<WindingSkinEffectLossesModel>> persistentModels;
    auto cached = persistentModels.find(resolvedModel);
    if (cached != persistentModels.end()) {
        return cached->second;
    }
    auto model = WindingSkinEffectLossesModel::factory(resolvedModel);
    persistentModels[resolvedModel] = model;
    return model;
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
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wireMaterial, temperature);

        // BUG-001 FIX: Guard against f<=0
    if (frequency <= 0) {
        return std::numeric_limits<double>::max();
    }
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

    size_t numAmplitudes = harmonics.get_amplitudes().size();
    size_t numFrequencies = harmonics.get_frequencies().size();
    if (numAmplitudes != numFrequencies) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Harmonics amplitudes and frequencies size mismatch: " + std::to_string(numAmplitudes) + " vs " + std::to_string(numFrequencies));
    }

    for (size_t harmonicIndex = 1; harmonicIndex < numAmplitudes; ++harmonicIndex)
    {
        auto harmonicRmsCurrent = harmonics.get_amplitudes()[harmonicIndex] / sqrt(2);  // Because a harmonic is always sinusoidal
        auto harmonicFrequency = harmonics.get_frequencies()[harmonicIndex];
        auto harmonicRmsCurrentInTurn = harmonicRmsCurrent * currentDivider;
        auto dcLossPerMeterThisHarmonic = pow(harmonicRmsCurrentInTurn, 2) * dcResistancePerMeter;

        if (std::isnan(harmonicFrequency) || harmonicFrequency <= 0) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Invalid frequency at harmonicIndex=" + std::to_string(harmonicIndex) + ": " + std::to_string(harmonicFrequency));
        }

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
    if (operatingPoint.get_excitations_per_winding().empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Operating point has no excitations for skin effect losses");
    }
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

        if (static_cast<size_t>(windingIndex) >= operatingPoint.get_excitations_per_winding().size() ||
            !operatingPoint.get_excitations_per_winding()[windingIndex].get_current()) {
            // Matches the ohmic-loss path's handling of windings without excitation
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing current excitation for winding " + std::to_string(windingIndex) + " in skin effect losses");
        }
        SignalDescriptor current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();

        auto lossesPerMeterPerHarmonic = calculate_skin_effect_losses_per_meter(wire, current, temperature, currentDividerPerTurn[turnIndex], modelOverride).second;

        LossElementPerHarmonic skinEffectLosses;
        auto model = get_model(coil.get_wire_type(windingIndex), modelOverride);
        skinEffectLosses.set_method_used(model->methodName);
        skinEffectLosses.set_origin(ResultOrigin::SIMULATION);
        skinEffectLosses.get_mutable_harmonic_frequencies().push_back(0);
        skinEffectLosses.get_mutable_losses_per_harmonic().push_back(0);

        for (auto& lossesPerMeter : lossesPerMeterPerHarmonic) {
            skinEffectLosses.get_mutable_harmonic_frequencies().push_back(lossesPerMeter.second);
            double turnLoss = lossesPerMeter.first * wireLength;
            skinEffectLosses.get_mutable_losses_per_harmonic().push_back(turnLoss);
            totalSkinEffectLosses += turnLoss;
        }
        windingLossesPerTurn[turnIndex].set_skin_effect_losses(skinEffectLosses);

    }
    windingLossesOutput.set_winding_losses_per_turn(windingLossesPerTurn);

    windingLossesOutput.set_method_used("AnalyticalModels");
    if (std::isnan(totalSkinEffectLosses)) {
        throw NaNResultException("NaN found in total skin effect losses");
    }
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalSkinEffectLosses);
    return windingLossesOutput;
}

std::optional<double> WindingSkinEffectLossesModel::try_get_skin_factor(Wire wire, double frequency, double temperature) {
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    std::size_t hash;
    if (!wire.get_name()) {
        hash = hash_combine_wire(wire.get_number_conductors().value(), wire.get_maximum_outer_width(), wire.get_maximum_outer_height());
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
        hash = hash_combine_wire(wire.get_number_conductors().value(), wire.get_maximum_outer_width(), wire.get_maximum_outer_height());
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
            penetrationRatio = pow(std::numbers::pi / 4.0, 3.0 / 4.0) * resolve_conducting_diameter(wire) / skinDepth * sqrt(resolve_conducting_diameter(wire) / resolve_dimensional_values(wire.get_outer_diameter().value()));
            break;
        }
        case WireType::LITZ: {
            if (!wire.get_strand())
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Litz wire is missing strand information");

            auto strand = Wire::resolve_strand(wire);
            penetrationRatio = pow(std::numbers::pi / 4.0, 3.0 / 4.0) * resolve_dimensional_values(strand.get_conducting_diameter()) / skinDepth * sqrt(resolve_dimensional_values(strand.get_conducting_diameter()) / resolve_dimensional_values(strand.get_outer_diameter().value()));
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
    // Cited Eq. (7): Fs = A*(sinh 2A + sin 2A)/(cosh 2A - cos 2A), HF limit
    // Fs -> A. The previous (A/2)*g(A) was the same function evaluated at half
    // the argument, giving a high-frequency factor of A/2 (2x underestimate).
    auto factor = penetrationRatio * (sinh(2 * penetrationRatio) + sin(2 * penetrationRatio)) / (cosh(2 * penetrationRatio) - cos(2 * penetrationRatio));
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
        if (!wire.get_conducting_width() || !wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting dimensions for planar/rectangular/foil wire");
        }
        wireRadius = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value())) / 2;
        if (!wire.get_outer_width() || !wire.get_outer_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing outer dimensions for planar/rectangular/foil wire");
        }
        wireOuterRadius = std::min(resolve_dimensional_values(wire.get_outer_width().value()), resolve_dimensional_values(wire.get_outer_height().value())) / 2;
    }
    else if (wire.get_type() == WireType::ROUND) {
        if (!wire.get_conducting_diameter()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting diameter for round wire");
        }
        wireRadius = resolve_conducting_diameter(wire) / 2;
        if (!wire.get_outer_diameter()) {
            wireOuterRadius = Wire::calculate_outer_diameter(wire, OpenMagnetics::DimensionalValues::NOMINAL) / 2;
        } else {
            wireOuterRadius = resolve_dimensional_values(wire.get_outer_diameter().value()) / 2;
        }
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        // The bundle term n(n-1)*(rD/rO)^2 needs rO = outer radius of the whole
        // BUNDLE (it encodes the strand distribution inside the bundle). Using
        // the strand's own outer radius makes rD/rO ~ 1 and inflates the term
        // by (rO_bundle/rO_strand)^2, independent of packing.
        if (!wire.get_outer_diameter()) {
            wireOuterRadius = Wire::calculate_outer_diameter(wire, OpenMagnetics::DimensionalValues::NOMINAL) / 2;
        } else {
            wireOuterRadius = resolve_dimensional_values(wire.get_outer_diameter().value()) / 2;
        }
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    std::complex<double> alpha(1, 1);
    alpha *= wireRadius / skinDepth;
    
    if (std::isnan(skinDepth) || std::isnan(wireRadius) || std::isnan(wireOuterRadius)) {
        throw NaNResultException("NaN in Albach model inputs: skinDepth=" + std::to_string(skinDepth) + 
                                 ", wireRadius=" + std::to_string(wireRadius) + 
                                 ", wireOuterRadius=" + std::to_string(wireOuterRadius));
    }

    if (wireOuterRadius <= 0) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_RESULT, 
            "wireOuterRadius must be positive in Albach model: " + std::to_string(wireOuterRadius));
    }

    // I1/I0 with a large-argument asymptotic branch (Utils helper). This replaces
    // the previous cap-at-100 + inf-substitution machinery, and fixes the bundle
    // term: its argument used to be alpha * wireRadius — alpha is ALREADY
    // (1+j)*r/delta (dimensionless), so the extra radius made the argument
    // meters-valued (~1e-7) and the n(n-1)(rD/rO)^2 strand-proximity term
    // vanished identically (Re[(1+j)^2] = 0). The documented bundle ratio
    // I1(alpha*rD)/I0(alpha*rD) uses the SAME normalized argument.
    std::complex<double> ratioI1I0 = modified_bessel_ratio_I1_I0(alpha);
    std::complex<double> besselRatio1 = 1.0 / ratioI1I0;
    std::complex<double> besselRatio2 = ratioI1I0;
    
    if (!wire.get_number_conductors()) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing number of conductors for wire");
    }
    double numConductors = wire.get_number_conductors().value();
    double geometricFactor = (numConductors * (numConductors - 1) * pow(wireRadius, 2)) / pow(wireOuterRadius, 2);
    
    double factor = 0.5 * (alpha * (besselRatio1 + geometricFactor * besselRatio2)).real();

    // Final safety check - cap to reasonable physical values
    // Skin factor should typically be between 1 (no skin effect) and ~10 (severe skin effect)
    const double MAX_SKIN_FACTOR = 100.0;
    if (std::isnan(factor) || std::isinf(factor)) {
        factor = MAX_SKIN_FACTOR;
    } else if (factor > MAX_SKIN_FACTOR) {
        factor = MAX_SKIN_FACTOR;
    } else if (factor < 1.0) {
        factor = 1.0;  // Minimum skin factor is 1 (no additional losses)
    }

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

    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        // Payne ("AC Resistance of Rectangular Conductors") is a rectangular model;
        // map round/litz to the equal-area square (side = d·√(π/4)), as the other skin
        // models do. Litz carries the diameter on the strand — resolve it (was an
        // unguarded get_conducting_height()/width().value() → bad_optional_access).
        double d;
        if (wire.get_type() == WireType::LITZ) {
            auto strand = wire.resolve_strand();
            d = resolve_dimensional_values(strand.get_conducting_diameter());
        }
        else {
            if (!wire.get_conducting_diameter()) {
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting diameter in round wire");
            }
            d = resolve_conducting_diameter(wire);
        }
        thinDimension = thickDimension = d * std::sqrt(std::numbers::pi / 4.0);
    }
    else {
        if (!wire.get_conducting_height() || !wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting width/height in wire");
        }
        double h = resolve_dimensional_values(wire.get_conducting_height().value());
        double w = resolve_dimensional_values(wire.get_conducting_width().value());
        if (h > w) { thinDimension = w; thickDimension = h; }
        else       { thinDimension = h; thickDimension = w; }
    }
    double A = thinDimension * thickDimension * 1000000;  // en mm2, shame on you Payne... (= equal-area square for round/litz)

    double p = pow(A, 0.5) / (1.26 * skinDepth * 1000);
    double Ff = 1.0 - exp(-0.026 * p);
    // Floating-point arithmetic never throws, so the previous try/catch around
    // this expression was dead code (and would have silently zeroed the AC
    // resistance correction had it ever fired).
    double Kc = 1.0 + Ff * (1.2 / exp(2.1 * thickDimension / thinDimension) + 1.2 / exp(2.1 * thinDimension / thickDimension));
    double x = (2.0 * skinDepth / thickDimension * (1.0 + thickDimension / thinDimension) + 8.0 * pow(skinDepth / thickDimension, 3) / (thinDimension / thickDimension)) / (pow(thinDimension / thickDimension, 0.33) * exp(-3.5 * thickDimension / skinDepth) + 1.0);
    // double acResistanceFactor =  (Kc / (1.0 - exp(-x))) - 1.0;
    double acResistanceFactor =  (Kc / (1.0 - exp(-x)));
    auto turnLosses = dcLossTurn * (acResistanceFactor - 1);
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
        wireHeight = resolve_conducting_diameter(wire);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double xi = wireHeight / skinDepth;
    // Standalone-foil Dowell skin factor: Fs = (xi/2)*(sinh xi + sin xi)/(cosh xi
    // - cos xi), which tends to 1 at DC (the previous xi/4 tended to 0.5,
    // fabricating loss at low frequency once combined with the dcLoss below)
    double factor = xi / 2 * (sinh(xi) + sin(xi)) / (cosh(xi) - cos(xi));

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

    // Return only the AC increment: the DC part is already accounted for by
    // WindingOhmicLosses on the full waveform RMS (matches every sibling model)
    auto turnLosses = dcLossTurn * (skinFactor - 1);
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
        b = resolve_conducting_diameter(wire) / 2;
        a = resolve_conducting_diameter(wire) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double c = sqrt(pow(b, 2) - pow(a, 2));
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    double acResistance = resistivity / (pow(std::numbers::pi, 2) * skinDepth * b) * comp_ellint_1(c / b) * (1 - exp(-2 * a / skinDepth));

    // Lotfi gives total Rac per meter. Convert to extra loss beyond DC:
    // P_skin = Rac * Irms^2 - Rdc * Irms^2 = (Rac - Rdc) * Irms^2
    // Since dcLossTurn = Rdc * Irms^2, and Rac * Irms^2 = acResistance * currentRms^2:
    auto turnLosses = acResistance * pow(currentRms, 2) - dcLossTurn;
    return (turnLosses > 0) ? turnLosses : 0.0;
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
        bPrima = resolve_conducting_diameter(wire) / 2;
        aPrima = resolve_conducting_diameter(wire) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
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
// ============================================================================
// DOWELL MODEL — Implementation
// Dowell 1966, Proc. IEE Vol. 113, No. 8
// FR_skin = M' = ζ·(sinh 2ζ + sin 2ζ)/(cosh 2ζ − cos 2ζ),  ζ = h/δ
// Round wire: h_eq = d·√(π/4)
// ============================================================================
double WindingSkinEffectLossesDowellModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double conductorHeight;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        conductorHeight = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                                   resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND) {
        // Dowell replaces round with equivalent square of same area: h = d·√(π/4)
        conductorHeight = resolve_conducting_diameter(wire) * sqrt(std::numbers::pi / 4.0);
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        conductorHeight = resolve_dimensional_values(strand.get_conducting_diameter()) * sqrt(std::numbers::pi / 4.0);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double zeta = conductorHeight / skinDepth;
    if (zeta < 1e-10) return 1.0;

    // Dowell Eq. (10), m=1: M' = ζ·(sinh 2ζ + sin 2ζ)/(cosh 2ζ − cos 2ζ)
    return zeta * (sinh(2 * zeta) + sin(2 * zeta)) / (cosh(2 * zeta) - cos(2 * zeta));
}

double WindingSkinEffectLossesDowellModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// PERRY MODEL — Implementation
// Perry 1979, IEEE Trans. PAS Vol. PAS-98, No. 1, pp. 116-123
// FR = Δ·(sinh 2Δ + sin 2Δ)/(cosh 2Δ − cos 2Δ),  Δ = T/δ  [Eq. 1 / 13a]
// Multi-layer proximity F₂'(Δ) term (Eq. 12) NOT included here.
// ============================================================================
double WindingSkinEffectLossesPerryModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double conductorThickness;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        conductorThickness = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                                       resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND) {
        conductorThickness = resolve_conducting_diameter(wire) * sqrt(std::numbers::pi / 4.0);
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        conductorThickness = resolve_dimensional_values(strand.get_conducting_diameter()) * sqrt(std::numbers::pi / 4.0);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double delta = conductorThickness / skinDepth;
    if (delta < 1e-10) return 1.0;

    // Perry Eq. (1) / (13a): FR = Δ·(sinh 2Δ + sin 2Δ)/(cosh 2Δ − cos 2Δ)
    return delta * (sinh(2 * delta) + sin(2 * delta)) / (cosh(2 * delta) - cos(2 * delta));
}

double WindingSkinEffectLossesPerryModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// DIMITRAKAKIS MODEL — Implementation
// Dimitrakakis et al. 2007, EPE Conference
// Skin (Eq. 5-6): FR_skin = Re[α·I₀(α)/I₁(α)]/2, α = (1+j)·r/δ
// Proximity r_prox (Eq. 6, ber₂/bei₂ terms) → proximity model.
// Rectangular fallback: Dowell 1D skin factor.
// ============================================================================
double WindingSkinEffectLossesDimitrakakisModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireRadius;
        if (wire.get_type() == WireType::LITZ) {
            auto strand = Wire::resolve_strand(wire);
            wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        } else {
            wireRadius = resolve_conducting_diameter(wire) / 2;
        }
        if (wireRadius / skinDepth < 1e-10) return 1.0;

        // Exact Bessel: FR = Re[α·I₀(α)/I₁(α)] / 2,  α = (1+j)·r/δ
        std::complex<double> alpha(1.0, 1.0);
        alpha *= wireRadius / skinDepth;
        // Ratio helper has a large-argument asymptotic branch; the raw truncated
        // series gave NEGATIVE skin factors for r/delta >~ 25
        return 0.5 * (alpha / modified_bessel_ratio_I1_I0(alpha)).real();
    }
    else {
        // Rectangular/foil/planar: Dowell 1D fallback
        double h = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                            resolve_dimensional_values(wire.get_conducting_height().value()));
        double zeta = h / skinDepth;
        if (zeta < 1e-10) return 1.0;
        return zeta * (sinh(2 * zeta) + sin(2 * zeta)) / (cosh(2 * zeta) - cos(2 * zeta));
    }
}

double WindingSkinEffectLossesDimitrakakisModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// MUEHLETHALER MODEL — Implementation
// Mühlethaler 2012, ETH Zurich DISS. ETH NO. 20217
//
// Round (Eq. 4.6-4.7, Appendix A.8.1):
//   FR = Re[(1+j)·(r/δ)·I₀((1+j)r/δ)/I₁((1+j)r/δ)] / 2
//
// Foil (Eq. 4.20, Appendix A.7.1):
//   FF = (Δ/4)·(sinh Δ + sin Δ)/(cosh Δ − cos Δ),  Δ = h/δ
//
// Proximity factors (Eq. 4.8-4.9) → proximity model.
// ============================================================================
double WindingSkinEffectLossesMuehlethalerModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireRadius;
        if (wire.get_type() == WireType::LITZ) {
            auto strand = Wire::resolve_strand(wire);
            wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        } else {
            wireRadius = resolve_conducting_diameter(wire) / 2;
        }
        if (wireRadius / skinDepth < 1e-10) return 1.0;

        // Bessel form of Eq. (4.7): FR = Re[α·I₀(α)/I₁(α)] / 2
        std::complex<double> alpha(1.0, 1.0);
        alpha *= wireRadius / skinDepth;
        // Ratio helper has a large-argument asymptotic branch; the raw truncated
        // series gave NEGATIVE skin factors for r/delta >~ 25
        return 0.5 * (alpha / modified_bessel_ratio_I1_I0(alpha)).real();
    }
    else {
        // Foil/rectangular/planar: standalone-foil Dowell skin factor
        // Fs = (Δ/2)·(sinh Δ + sin Δ)/(cosh Δ − cos Δ), which tends to 1 at DC.
        // The previous Δ/4 form tended to 0.5, so calculate_turn_losses below
        // (dcLoss·(Fs−1)) returned NEGATIVE loss at low frequency and half the
        // true skin loss at high frequency — the same convention bug fixed in the
        // Ferreira model (see the Δ/2 comment there); this foil branch was missed.
        double h = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                            resolve_dimensional_values(wire.get_conducting_height().value()));
        double delta = h / skinDepth;
        if (delta < 1e-10) return 1.0;
        return (delta / 2.0) * (sinh(delta) + sin(delta)) / (cosh(delta) - cos(delta));
    }
}

double WindingSkinEffectLossesMuehlethalerModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// NAN MODEL — Implementation
// Nan & Sullivan 2003, IEEE PESC, pp. 853-860
// FR = η·(sinh 2η + sin 2η)/(cosh 2η − cos 2η),  η = d_eq/(2δ)  [Eq. 2]
// Proximity correction (Eq. 7, porosity-corrected) → proximity model.
// ============================================================================
double WindingSkinEffectLossesNanModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double conductorDimension;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        conductorDimension = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                                       resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND) {
        conductorDimension = resolve_conducting_diameter(wire) * sqrt(std::numbers::pi / 4.0);
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        conductorDimension = resolve_dimensional_values(strand.get_conducting_diameter()) * sqrt(std::numbers::pi / 4.0);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    // Nan's η = d_eq / (2δ)  — half-thickness convention
    double eta = conductorDimension / (2.0 * skinDepth);
    if (eta < 1e-10) return 1.0;

    // Eq. (2): FR = η·(sinh 2η + sin 2η)/(cosh 2η − cos 2η)
    return eta * (sinh(2 * eta) + sin(2 * eta)) / (cosh(2 * eta) - cos(2 * eta));
}

double WindingSkinEffectLossesNanModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// KAZIMIERCZUK MODEL — Implementation
// Whitman & Kazimierczuk 2019, IEEE TPEL, DOI: 10.1109/TPEL.2019.2904582
// Base Dowell skin factor; cylindrical correction C_F requires winding radii.
// TODO: apply C_F ≈ 1 + h²/(12·R_in·R_out) when radii are available.
// ============================================================================
double WindingSkinEffectLossesKazimierczukModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double conductorHeight;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        conductorHeight = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                                   resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND) {
        conductorHeight = resolve_conducting_diameter(wire) * sqrt(std::numbers::pi / 4.0);
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = Wire::resolve_strand(wire);
        conductorHeight = resolve_dimensional_values(strand.get_conducting_diameter()) * sqrt(std::numbers::pi / 4.0);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    double zeta = conductorHeight / skinDepth;
    if (zeta < 1e-10) return 1.0;

    // Base Dowell skin factor (cylindrical correction C_F needs R_in, R_out)
    return zeta * (sinh(2 * zeta) + sin(2 * zeta)) / (cosh(2 * zeta) - cos(2 * zeta));
}

double WindingSkinEffectLossesKazimierczukModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}


// ============================================================================
// WANG MODEL — Implementation
// Wang et al. 2018, IEEE APEC, DOI: 10.1109/APEC.2018.8341586
//
// 2D model (Eq. 2-5), rectangular conductors:
//   f(z) = z·(sinh 2z+sin 2z)/(cosh 2z−cos 2z)
//   λ = 0.01·(c/h) + 0.66
//   Hy/Hx = c·[1/(λh) + 1/(c−λh)] / π
//   FR_2D = f(h/δ) + (Hy/Hx)²·(h/c)·f(c/δ)
//
// Round wire: Dowell 1D fallback.
// ============================================================================
double WindingSkinEffectLossesWangModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        // Fallback to Dowell for round/litz
        double d;
        if (wire.get_type() == WireType::LITZ) {
            auto strand = Wire::resolve_strand(wire);
            d = resolve_dimensional_values(strand.get_conducting_diameter()) * sqrt(std::numbers::pi / 4.0);
        } else {
            d = resolve_conducting_diameter(wire) * sqrt(std::numbers::pi / 4.0);
        }
        double zeta = d / skinDepth;
        if (zeta < 1e-10) return 0.0;
        double skinFactor = zeta * (sinh(2 * zeta) + sin(2 * zeta)) / (cosh(2 * zeta) - cos(2 * zeta));
        return dcLossTurn * (skinFactor - 1);
    }

    double c = resolve_dimensional_values(wire.get_conducting_width().value());
    double h = resolve_dimensional_values(wire.get_conducting_height().value());
    if (h > c) std::swap(c, h);

    double k = c / h;
    double lambda = 0.01 * k + 0.66;   // Eq. (5)

    // 1D skin function
    auto f = [](double z) -> double {
        if (z < 1e-10) return 1.0;
        return z * (sinh(2 * z) + sin(2 * z)) / (cosh(2 * z) - cos(2 * z));
    };

    double FR_x = f(h / skinDepth);   // 1-D skin R_ac/R_dc in the thin direction (→1 at DC)
    double FR_y = f(c / skinDepth);   // 1-D skin R_ac/R_dc in the wide (edge) direction

    double lambdaH = lambda * h;
    double Hy_over_Hx = 0.0;
    if (lambdaH > 0 && (c - lambdaH) > 0) {
        Hy_over_Hx = c * (1.0 / lambdaH + 1.0 / (c - lambdaH)) / std::numbers::pi;
    }

    // R_ac/R_dc = skin term + edge (2-D) term. The 2-D edge contribution is an
    // ADDED loss that must VANISH at DC (δ→∞): like every eddy/proximity effect it
    // grows from zero with frequency. The previous form used FR_y directly, but
    // f(c/δ)→1 at DC, so factor_2D → 1 + (Hy/Hx)²·(h/c) (~5 for a 10×0.2 mm foil) —
    // a fabricated ~4× DC loss at arbitrarily low frequency (ABT #115). Use the
    // excess skin factor (FR_y − 1), which → 0 at DC (restoring R_ac/R_dc → 1) and
    // ≈ FR_y at HF (unchanged high-frequency behaviour). NB: the absolute Eq.-2
    // magnitude (and the Hy/Hx normalisation) is not FEM-calibrated here — this
    // fixes the DC-limit pathology only.
    double factor_2D = FR_x + pow(Hy_over_Hx, 2) * (h / c) * (FR_y - 1.0);
    double turnLosses = dcLossTurn * (factor_2D - 1);
    return (turnLosses < 0) ? 0.0 : turnLosses;
}


// ============================================================================
// HOLGUIN MODEL — Implementation
// Holguin et al. 2014, IEEE APEC, pp. 2009-2014
// FR = Re[(1+j)·(r₀/δ)·I₀((1+j)r₀/δ)/I₁((1+j)r₀/δ)] / 2  [from Lammeraner]
// Gap-fringing proximity (Eq. 1-5) → proximity model.
// Rectangular/foil fallback: Dowell 1D skin factor.
// ============================================================================
double WindingSkinEffectLossesHolguinModel::calculate_skin_factor(const Wire& wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        double wireRadius;
        if (wire.get_type() == WireType::LITZ) {
            auto strand = Wire::resolve_strand(wire);
            wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        } else {
            wireRadius = resolve_conducting_diameter(wire) / 2;
        }
        if (wireRadius / skinDepth < 1e-10) return 1.0;

        // Exact Bessel solution: FR = Re[α·I₀(α)/I₁(α)] / 2,  α = (1+j)·r₀/δ
        std::complex<double> alpha(1.0, 1.0);
        alpha *= wireRadius / skinDepth;
        // Ratio helper has a large-argument asymptotic branch; the raw truncated
        // series gave NEGATIVE skin factors for r/delta >~ 25
        return 0.5 * (alpha / modified_bessel_ratio_I1_I0(alpha)).real();
    }
    else {
        // Rectangular/foil: Dowell 1D skin factor
        double h = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                            resolve_dimensional_values(wire.get_conducting_height().value()));
        double zeta = h / skinDepth;
        if (zeta < 1e-10) return 1.0;
        return zeta * (sinh(2 * zeta) + sin(2 * zeta)) / (cosh(2 * zeta) - cos(2 * zeta));
    }
}

double WindingSkinEffectLossesHolguinModel::calculate_turn_losses(Wire wire, double dcLossTurn, double frequency, double temperature, [[maybe_unused]] double currentRms) {
    double skinFactor;
    auto optionalSkinFactor = try_get_skin_factor(wire, frequency, temperature);
    if (optionalSkinFactor) {
        skinFactor = optionalSkinFactor.value();
    } else {
        skinFactor = calculate_skin_factor(wire, frequency, temperature);
        set_skin_factor(wire, frequency, temperature, skinFactor);
    }
    return dcLossTurn * (skinFactor - 1);
}

} // namespace OpenMagnetics

