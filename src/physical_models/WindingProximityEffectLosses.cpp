#include "physical_models/WindingSkinEffectLosses.h"
#include "WindingProximityEffectLosses.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Resistivity.h"
#include "Defaults.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"
#include "support/Logger.h"

namespace OpenMagnetics {


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

std::shared_ptr<WindingProximityEffectLossesModel>  WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels modelName){
    if (modelName == WindingProximityEffectLossesModels::ROSSMANITH) {
        return std::make_shared<WindingProximityEffectLossesRossmanithModel>();
    }
    else if (modelName == WindingProximityEffectLossesModels::WANG) {
        return std::make_shared<WindingProximityEffectLossesWangModel>();
    }
    else if (modelName == WindingProximityEffectLossesModels::FERREIRA) {
        return std::make_shared<WindingProximityEffectLossesFerreiraModel>();
    }
    else if (modelName == WindingProximityEffectLossesModels::ALBACH) {
        return std::make_shared<WindingProximityEffectLossesAlbachModel>();
    }
    else if (modelName == WindingProximityEffectLossesModels::LAMMERANER) {
        return std::make_shared<WindingProximityEffectLossesLammeranerModel>();
    }
    else if (modelName == WindingProximityEffectLossesModels::DOWELL) {
        return std::make_shared<WindingProximityEffectLossesDowellModel>();
    }
    else
        throw ModelNotAvailableException("Unknown wire proximity effect losses mode, available options are: {ROSSMANITH, WANG, FERREIRA, ALBACH, LAMMERANER, DOWELL}");
}

std::shared_ptr<WindingProximityEffectLossesModel> WindingProximityEffectLosses::get_model(WireType wireType, std::optional<WindingProximityEffectLossesModels> modelOverride) {
    // If an explicit model override is provided, use it
    if (modelOverride.has_value()) {
        return WindingProximityEffectLossesModel::factory(modelOverride.value());
    }

    // Check if user has enabled manual model selection
    auto& settings = Settings::GetInstance();
    if (settings.get_coil_enable_user_winding_losses_models()) {
        // Use the user's selected model from Settings
        auto userSelectedModel = settings.get_winding_proximity_effect_losses_model();
        return WindingProximityEffectLossesModel::factory(userSelectedModel);
    }

    // Otherwise, auto-select based on wire type (default behavior)
    switch(wireType) {
        case WireType::ROUND: {
            return WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels::FERREIRA);
        }
        case WireType::LITZ: {
            return WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels::FERREIRA);
        }
        case WireType::PLANAR: {
            return WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels::WANG);
        }
        case WireType::RECTANGULAR: {
            return WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels::WANG);
        }
        case WireType::FOIL: {
            return WindingProximityEffectLossesModel::factory(WindingProximityEffectLossesModels::WANG);
        }
        default:
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }
}

std::optional<double> WindingProximityEffectLossesModel::try_get_proximity_factor(Wire wire, double frequency, double temperature) {
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

    if (_proximityFactorPerWirePerFrequencyPerTemperature.contains(hash)) {
        if (_proximityFactorPerWirePerFrequencyPerTemperature[hash].contains(frequency)) {
            if (_proximityFactorPerWirePerFrequencyPerTemperature[hash][frequency].contains(temperature)) {
                return _proximityFactorPerWirePerFrequencyPerTemperature[hash][frequency][temperature];
            }
        }
    }

    return std::nullopt;
}

void WindingProximityEffectLossesModel::set_proximity_factor(Wire wire,  double frequency, double temperature, double proximityFactor) {
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

    _proximityFactorPerWirePerFrequencyPerTemperature[hash][frequency][temperature] = proximityFactor;

}

std::pair<double, std::vector<std::pair<double, double>>> WindingProximityEffectLosses::calculate_proximity_effect_losses_per_meter(Wire wire, double temperature, std::vector<ComplexField> fields, std::optional<WindingProximityEffectLossesModels> modelOverride) {
    auto model = get_model(wire.get_type(), modelOverride);
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }

    double totalProximityEffectLossesPerMeter = 0;
    std::vector<std::pair<double, double>> lossesPerHarmonic;

    for (auto& complexField : fields) {
        auto frequency = complexField.get_frequency();
        auto dataForThisTurn = complexField.get_data();

        auto turnLosses = model->calculate_turn_losses(wire, frequency, dataForThisTurn, temperature);

        if (std::isnan(turnLosses)) {
            throw NaNResultException("NaN found in proximity effect losses per meter");
        }
        lossesPerHarmonic.push_back(std::pair<double, double>{turnLosses, frequency});
        totalProximityEffectLossesPerMeter += turnLosses;
    }

    return {totalProximityEffectLossesPerMeter, lossesPerHarmonic};
}

WindingLossesOutput WindingProximityEffectLosses::calculate_proximity_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput, std::optional<WindingProximityEffectLossesModels> modelOverride) {
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");
    }
    auto turns = coil.get_turns_description().value();

    auto windingLossesPerTurn = windingLossesOutput.get_winding_losses_per_turn().value();

    double totalProximityEffectLosses = 0;

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex)
    {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = coil.resolve_wire(windingIndex);
        double wireLength = turn.get_length();

        std::vector<ComplexField> fields;
        for (auto& fieldPerHarmonic : windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()) {
            std::vector<ComplexFieldPoint> data;
            for (auto& fieldPoint : fieldPerHarmonic.get_data()) {
                if (!fieldPoint.get_turn_index()) {
                    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Missing turn index in field point");
                }
                if (fieldPoint.get_turn_index().value() == int(turnIndex)) {
                    data.push_back(fieldPoint);
                }
            }
            
            ComplexField complexField;
            complexField.set_data(data);
            complexField.set_frequency(fieldPerHarmonic.get_frequency());
            fields.push_back(complexField);
        }

        auto lossesPerHarmonicThisTurn = calculate_proximity_effect_losses_per_meter(wire, temperature, fields, modelOverride).second;


        WindingLossElement proximityEffectLossesThisTurn;
        auto model = get_model(coil.get_wire_type(windingIndex), modelOverride);
        proximityEffectLossesThisTurn.set_method_used(model->methodName);
        proximityEffectLossesThisTurn.set_origin(ResultOrigin::SIMULATION);
        proximityEffectLossesThisTurn.get_mutable_harmonic_frequencies().push_back(0);
        proximityEffectLossesThisTurn.get_mutable_losses_per_harmonic().push_back(0);

        for (auto& lossesThisHarmonic : lossesPerHarmonicThisTurn) {
            if (std::isnan(lossesThisHarmonic.first)) {
                throw NaNResultException("NaN found in proximity effect losses");
            }
            proximityEffectLossesThisTurn.get_mutable_harmonic_frequencies().push_back(lossesThisHarmonic.second);
            proximityEffectLossesThisTurn.get_mutable_losses_per_harmonic().push_back(lossesThisHarmonic.first * wireLength);

            totalProximityEffectLosses += lossesThisHarmonic.first * wireLength;

        }

        windingLossesPerTurn[turnIndex].set_proximity_effect_losses(proximityEffectLossesThisTurn);
    }
    windingLossesOutput.set_winding_losses_per_turn(windingLossesPerTurn);

    windingLossesOutput.set_method_used("AnalyticalModels");
    if (std::isnan(totalProximityEffectLosses)) {
        throw NaNResultException("NaN found in total proximity effect losses");
    }
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalProximityEffectLosses);
    return windingLossesOutput;
}

/**
 * @brief Calculates proximity factor using Rossmanith's model with Bessel functions.
 * 
 * Based on: H. Rossmanith et al., "Measurement and Characterization of High
 * Frequency Losses in Nonideal Litz Wires", IEEE Transactions on Power Electronics,
 * Vol. 26, No. 11, Nov. 2011
 * 
 * This model uses modified Bessel functions of the first kind to calculate the
 * proximity effect factor for round and litz wires. For rectangular conductors,
 * it uses the standard solution from: J. A. Ferreira, "Improved Analytical
 * Modeling of Conductive Losses in Magnetic Components", IEEE TPEL, Vol. 9, No. 1, 1994
 * 
 * For round wires, the proximity factor is:
 * 
 *   G_prox = 2π * α * I₁(α) / I₀(α)
 * 
 * where:
 * - α = (1+j) * r/δ, with r = wire radius, δ = skin depth
 * - I₀, I₁ = modified Bessel functions of the first kind
 * 
 * For rectangular/planar wires:
 *   G_prox = (w·h/δ) * [sinh(w/δ) - sin(w/δ)] / [cosh(w/δ) + cos(w/δ)]
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Proximity factor for loss calculation
 */
double WindingProximityEffectLossesRossmanithModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double factor;
    std::complex<double> alpha(1, 1);

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        double wireWidth = resolve_dimensional_values(wire.get_conducting_width().value());
        double wireHeight = resolve_dimensional_values(wire.get_conducting_height().value());

        factor = wireHeight * wireWidth / skinDepth * (sinh(wireWidth / skinDepth) - sin(wireWidth / skinDepth)) / (cosh(wireWidth / skinDepth) + cos(wireWidth / skinDepth));
    }
    else if (wire.get_type() == WireType::ROUND ) {
        double wireRadius;
        wireRadius = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
        alpha *= wireRadius / skinDepth;
        factor = 2 * std::numbers::pi * (alpha * modified_bessel_first_kind(1, alpha) / modified_bessel_first_kind(0, alpha)).real();
    }
    else if (wire.get_type() == WireType::LITZ) {
        double wireRadius;
        auto strand = wire.resolve_strand();
        wireRadius = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
        alpha *= wireRadius / skinDepth;
        factor = 2 * std::numbers::pi * (alpha * modified_bessel_first_kind(1, alpha) / modified_bessel_first_kind(0, alpha)).real();
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    return factor;
}

double WindingProximityEffectLossesRossmanithModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto optionalproximityFactor = try_get_proximity_factor(wire, frequency, temperature);

    if (optionalproximityFactor) {
        proximityFactor = optionalproximityFactor.value();
    }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }

    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    // BUG-002 FIX: Use RMS of |H|
    double He2_sum = 0;
    for (auto& datum : data) {
        He2_sum += pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
    }
    double He2_rms = He2_sum / data.size();
    double turnLosses = resistivity * He2_rms * proximityFactor;

    turnLosses *= wire.get_number_conductors().value();

    return turnLosses;
}

/**
 * @brief Calculates proximity effect losses using Wang's 2D field decomposition model.
 * 
 * Based on: Wang et al., "Improved Analytical Model for High-Frequency Winding Losses"
 * and the field decomposition approach for rectangular conductors.
 * 
 * This model decomposes the external magnetic field into components parallel and
 * perpendicular to the conductor faces, computing losses separately for each.
 * 
 * For a rectangular conductor with width c and height h:
 * 
 *   P_prox_x = c*h*ρ/δ * Hx_avg² * [sinh(h/δ) - sin(h/δ)] / [cosh(h/δ) + cos(h/δ)]
 *   P_prox_y = h*c*ρ/δ * Hy_avg² * [sinh(c/δ) - sin(c/δ)] / [cosh(c/δ) + cos(c/δ)]
 * 
 * where:
 * - Hx1, Hx2 = field at bottom and top surfaces (labeled 'bottom', 'top')
 * - Hy1, Hy2 = field at left and right surfaces (labeled 'left', 'right')
 * - Hx_avg = (Hx1 + Hx2) / 2, Hy_avg = (Hy1 + Hy2) / 2
 * 
 * This 2D approach captures the different penetration depths for fields parallel
 * to the short vs long dimension of the conductor.
 * 
 * For field components not aligned with the conductor faces, Ferreira's model
 * is used as a fallback.
 * 
 * @param wire Wire object (rectangular/foil conductor)
 * @param frequency Operating frequency [Hz]
 * @param data Field points at conductor surfaces (labeled top/bottom/left/right)
 * @param temperature Operating temperature [K]
 * @return Proximity effect losses for this turn [W]
 */
double WindingProximityEffectLossesWangModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double c = 0;
    double h = 0;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {

        if (!wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting height in wire");
        }
        c = resolve_dimensional_values(wire.get_conducting_width().value());
        h = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    else {
        throw NotImplementedException("Model not implemented for ROUND and LITZ");
    }

    double Hx1 = 0, Hx2 = 0, Hy1 = 0, Hy2 = 0;
    double nonPlanarHe = 0;
    for (auto& datum : data) {
        if (!datum.get_label()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing label in induced point");
        }
        else if (datum.get_label().value() == "top") {
            nonPlanarHe += datum.get_imaginary(); 
            Hx2 += datum.get_real();
        }
        else if (datum.get_label().value() == "bottom") {
            nonPlanarHe += datum.get_imaginary(); 
            Hx1 += datum.get_real();
        }
        else if (datum.get_label().value() == "right") {
            nonPlanarHe += datum.get_real(); 
            Hy2 += datum.get_imaginary();
        }
        else if (datum.get_label().value() == "left") {
            nonPlanarHe += datum.get_real(); 
            Hy1 += datum.get_imaginary();
        }
    }


    double turnLosses = 0;
    double hTerm = h / skinDepth;
    double cTerm = c / skinDepth;
    if (hTerm > 710) {
        // to avoid cosh overflow according to https://cpp-lang.net/docs/std/math/mathematical_functions/cosh/
        hTerm = 710;
    }
    if (cTerm > 710) {
        // to avoid cosh overflow according to https://cpp-lang.net/docs/std/math/mathematical_functions/cosh/
        cTerm = 710;
    }
    turnLosses += c * h * resistivity / skinDepth * pow((Hx2 + Hx1) / 2, 2) * (sinh(hTerm) - sin(hTerm)) / (cosh(hTerm) + cos(hTerm));
    turnLosses += h * c * resistivity / skinDepth * pow((Hy2 + Hy1) / 2, 2) * (sinh(cTerm) - sin(cTerm)) / (cosh(cTerm) + cos(cTerm));

    // BUG-003 FIX: Normalize nonPlanarHe by data.size()
    if (nonPlanarHe != 0 && !data.empty()) {
        nonPlanarHe /= data.size();
        double proximityFactor = WindingProximityEffectLossesFerreiraModel::calculate_proximity_factor(wire, frequency, temperature);
        turnLosses += proximityFactor * pow(nonPlanarHe, 2);
    }

    turnLosses *= wire.get_number_conductors().value();

    return turnLosses;
}

/**
 * @brief Calculates proximity effect factor using Ferreira's analytical method.
 * 
 * Based on: "Improved Analytical Modeling of Conductive Losses in Magnetic Components"
 * by J.A. Ferreira, IEEE Transactions on Power Electronics, Vol. 9, No. 1, January 1994
 * 
 * This implementation uses the orthogonality principle between skin and proximity effects
 * (Appendix A of the paper). The proximity losses per unit length are:
 * 
 *   P_prox = G · H_e²                                           (Eq. A7)
 * 
 * For round conductors, the factor G uses Kelvin functions (ber, bei):
 * 
 *   G = -2πγ · [ber₂(γ)·ber'(γ) + bei₂(γ)·bei'(γ)] / [ber²(γ) + bei²(γ)]   (Eq. A8)
 * 
 * where γ = d/(δ√2), d = wire diameter, δ = skin depth.
 * 
 * For rectangular conductors, the factor uses hyperbolic/trigonometric functions:
 * 
 *   G = (√π/2σ) · ξ · [sinh(ξ) - sin(ξ)] / [cosh(ξ) + cos(ξ)]   (Eq. A9)
 * 
 * where ξ = min(h,w)/δ is the normalized smaller dimension.
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Proximity factor G for loss calculation [Ω/m per (A/m)²]
 */
double WindingProximityEffectLossesFerreiraModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    double factor;
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (!wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting height in wire");
        }
        double w = resolve_dimensional_values(wire.get_conducting_width().value());
        double h = resolve_dimensional_values(wire.get_conducting_height().value());

        double xi = std::min(h, w) / skinDepth;

        factor = w * xi * resistivity * (sinh(xi) - sin(xi)) / (cosh(xi) + cos(xi));
        if (std::isnan(factor)) {
            throw NaNResultException("NaN found in Ferreira's proximity factor");
        }
    }
    else if (wire.get_type() == WireType::ROUND ||wire.get_type() == WireType::LITZ) {
        double wireDiameter;
        if (wire.get_type() == WireType::ROUND) {
            wireDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
        }
        else {
            auto strand = wire.resolve_strand();
            wireDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        }
        double gamma = wireDiameter / (skinDepth * sqrt(2));
        factor = - 2 * gamma * resistivity * (kelvin_function_real(2, gamma) * derivative_kelvin_function_real(0, gamma) + kelvin_function_imaginary(2, gamma) * derivative_kelvin_function_imaginary(0, gamma)) / (pow(kelvin_function_real(0, gamma), 2) + pow(kelvin_function_imaginary(0, gamma), 2));


    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    return factor;
}

double WindingProximityEffectLossesFerreiraModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto optionalproximityFactor = try_get_proximity_factor(wire, frequency, temperature);

    if (optionalproximityFactor) {
        proximityFactor = optionalproximityFactor.value();
    }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He = 0;

    for (auto& datum : data) {
        if (std::isnan(datum.get_real())) {
            throw NaNResultException("NaN found in Ferreira proximity losses calculation");
        }
        if (std::isnan(datum.get_imaginary())) {
            throw NaNResultException("NaN found in Ferreira proximity losses calculation");
        }
        He = std::max(He, hypot(datum.get_real(), datum.get_imaginary()));
    }

    double turnLosses = proximityFactor * pow(He, 2);
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    turnLosses *= wire.get_number_conductors().value();
        if (std::isnan(turnLosses)) {
            throw NaNResultException("NaN found in Ferreira proximity losses calculation: frequency=" + 
                std::to_string(frequency) + ", proximityFactor=" + std::to_string(proximityFactor) + 
                ", He=" + std::to_string(He));
        }

    return turnLosses;
}

/**
 * @brief Calculates proximity effect losses per turn using the Albach model.
 * 
 * Based on: "Induktivitäten in der Leistungselektronik: Spulen, Trafos und ihre 
 * parasitären Eigenschaften" by Manfred Albach
 * Springer Vieweg, 2017, ISBN 978-3-658-15081-5
 * Chapter 4: "Die Verluste in Luftspulen", Section 4.2 "Proximityverluste", pages 79-96
 * 
 * The Albach model calculates proximity losses based on the external magnetic field
 * causing eddy currents in nearby conductors. The model accounts for:
 * - Inhomogeneous external field distribution (Section 4.2.2)
 * - Proximity factor Ds,k for different field harmonics k (Eq. 4.48)
 * - Rückwirkung (reaction) of proximity currents on external field (Section 4.2.3)
 * 
 * Key equations from Chapter 4:
 * 
 * - Eq. (4.32): Proximity losses for homogeneous external field
 *   Pprox = (l/4) * |He|² * Ds
 *   where Ds is the proximity factor based on Bessel functions
 * 
 * - Eq. (4.48): Generalized proximity losses for inhomogeneous field
 *   Pprox = (l/4) * Σ[(|bk+ck|² + |ak-dk|²) * Ds,k]
 *   where ak, bk, ck, dk are Fourier coefficients of surface field
 * 
 * - Skin depth δ from Eq. (1.39):
 *   δ = √(2ρ/(ωμ))  where ρ = resistivity, ω = 2πf
 * 
 * - Complex propagation parameter α from Eq. (4.16):
 *   α = (1+j)/δ
 * 
 * For foil/rectangular conductors, the formula uses tanh for 1D approximation:
 *   Pprox = c * ρ * He² * Re(α*d * tanh(α*d/2))
 *   where c = conductor height, d = conductor width
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param data Vector of complex magnetic field points around the conductor
 * @param temperature Operating temperature [K]
 * @return Proximity effect losses for this turn [W]
 */
double WindingProximityEffectLossesAlbachModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double d = 0;
    double c = 0;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (!wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting height in wire");
        }
        d = resolve_dimensional_values(wire.get_conducting_width().value());
        c = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    else {
        d = resolve_dimensional_values(wire.get_conducting_diameter().value());
        c = resolve_dimensional_values(wire.get_conducting_diameter().value());
    }
    std::complex<double> alpha(1, 1);
    alpha /= skinDepth;

    // BUG-002 FIX: Use RMS of |H|
    double He2_sum = 0;
    for (auto& datum : data) {
        He2_sum += pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
    }
    double He2_rms = He2_sum / data.size();
    double turnLosses = c * resistivity * He2_rms * (alpha * d * tanh(alpha * d / 2.0)).real();

    turnLosses *= wire.get_number_conductors().value();

    if (std::isnan(turnLosses)) {
        throw NaNResultException("NaN found in Albach's model for proximity effect losses");
    }

    return turnLosses;
}

/**
 * @brief Calculates proximity factor using Lammeraner's low-frequency approximation.
 * 
 * Based on: "Eddy Currents" by J. Lammeraner and M. Stafl
 * Iliffe Books, London, 1966
 * 
 * Also referenced in: Kutkut 1998, "A Simple Technique to Evaluate Winding Losses
 * Including Two-Dimensional Edge Effects"
 * 
 * Lammeraner's model provides a simplified proximity factor valid for low frequencies
 * where the skin depth is larger than the conductor dimension (d/δ < 1).
 * 
 * The proximity factor scales with the fourth power of the ratio:
 * 
 *   F_prox = 2πρ * (r/2δ)⁴ / 4
 * 
 * where:
 * - r = conducting dimension (radius for round, min dimension for rectangular)
 * - δ = skin depth
 * - ρ = resistivity
 * 
 * The losses are: P_prox = (Hx² + Hy²) * F_prox
 * 
 * This model is computationally efficient but underestimates losses at high frequencies
 * where the skin effect becomes significant.
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Proximity factor for loss calculation
 */
double WindingProximityEffectLossesLammeranerModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double factor;
    double wireConductingDimension;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        wireConductingDimension = std::min(resolve_dimensional_values(wire.get_conducting_width().value()), resolve_dimensional_values(wire.get_conducting_height().value()));
    }
    else if (wire.get_type() == WireType::ROUND ) {
        wireConductingDimension = resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2;
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        wireConductingDimension = resolve_dimensional_values(strand.get_conducting_diameter()) / 2;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    // CRITICAL FIX: Lammeraner's formula is only valid for low frequencies where d/δ < 1
    // (Lammeraner & Stafl, 1966; Kutkut 1998). At high frequencies where skin depth is 
    // much smaller than conductor dimension, this formula overestimates losses by (d/δ)^4.
    // For rectangular wires at 1MHz with d=0.9mm and δ=66μm, d/δ=13.6, leading to
    // errors of ~34,000x (19,833,400% error observed in testing).
    double normalizedDimension = wireConductingDimension / skinDepth;
    
    if (normalizedDimension > 1.0) {
        // Model not valid at this frequency - d/δ must be < 1 for Lammeraner approximation
        throw InvalidInputException(ErrorCode::INVALID_INPUT, 
            "Lammeraner proximity model is only valid for low frequencies where conductor_dimension/skin_depth < 1. "
            "Current ratio is " + std::to_string(normalizedDimension) + 
            ". Use a different proximity model (e.g., WANG or FERREIRA) for high frequencies.");
    }

    factor = 2.0 * std::numbers::pi * resistivity * pow((wireConductingDimension / 2) / skinDepth, 4) / 4;

    return factor;
}

double WindingProximityEffectLossesLammeranerModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto optionalproximityFactor = try_get_proximity_factor(wire, frequency, temperature);

    if (optionalproximityFactor) {
        proximityFactor = optionalproximityFactor.value();
    }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }

    // BUG-004 FIX: Use mean of |H|^2 instead of averaging components separately
    double He2_sum = 0;
    for (auto& datum : data) {
        He2_sum += pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
    }
    double He2_mean = He2_sum / data.size();

    double turnLosses = He2_mean * proximityFactor;
    turnLosses *= wire.get_number_conductors().value();

    if (std::isnan(turnLosses)) {
        throw NaNResultException("NaN found in Lammeraner's model for proximity effect losses");
    }

    return turnLosses;
}


/**
 * @brief Calculates proximity factor using Dowell's analytical method for transformer windings.
 * 
 * Based on: "Effects of eddy currents in transformer windings" by P. L. Dowell
 * Proceedings of the IEE, Vol. 113, No. 8, August 1966
 * https://ieeexplore.ieee.org/document/5247417
 * 
 * Dowell's method models rectangular conductors in layered windings. The AC resistance
 * ratio FR (Eq. 10) is decomposed into skin effect (M') and proximity effect (D') parts:
 * 
 *   FR = M' + (m² - 1)/3 · D'                                    (Eq. 10)
 * 
 * where the proximity contribution uses the D function:
 * 
 *   M = αh · coth(αh)                                             (Eq. 32)
 *   D = 2αh · tanh(αh/2)                                          (Eq. 33)
 * 
 * with:
 * - α = (1+j)/δ, the complex propagation constant
 * - δ = skin depth = √(2ρ/(ωμ₀))
 * - h = conductor height (dimension in the direction of field penetration)
 * - M' and D' are the real parts of M and D, respectively
 * 
 * For round conductors, Dowell's approach replaces them with equivalent-area square
 * conductors (Section 3, Fig. 3 of the paper), yielding h = d·√(π/4).
 * 
 * The proximity factor for a single conductor of width a and height h exposed to
 * external field He gives losses per unit length:
 * 
 *   P_prox/l = (a · ρ / h) · He² · D'
 * 
 * @param wire Wire object containing conductor geometry
 * @param frequency Operating frequency [Hz]
 * @param temperature Operating temperature [K]
 * @return Proximity factor for loss calculation [Ω·m per (A/m)²]
 */
double WindingProximityEffectLossesDowellModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model(); // PERF-003: cached
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double factor;

    double h = 0; // conductor height (dimension in direction of field penetration)
    double a = 0; // conductor width (dimension perpendicular to field penetration)

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (!wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Missing conducting height in wire");
        }
        a = resolve_dimensional_values(wire.get_conducting_width().value());
        h = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    else if (wire.get_type() == WireType::ROUND) {
        // Dowell's approach: replace round conductor with equivalent square of same area
        // (Section 3, Fig. 3 of the paper)
        double d = resolve_dimensional_values(wire.get_conducting_diameter().value());
        a = d * sqrt(std::numbers::pi / 4.0);
        h = a;
    }
    else if (wire.get_type() == WireType::LITZ) {
        // Apply Dowell's round-to-square equivalence to each strand
        auto strand = wire.resolve_strand();
        double d = resolve_dimensional_values(strand.get_conducting_diameter());
        a = d * sqrt(std::numbers::pi / 4.0);
        h = a;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
    }

    // Dowell's complex propagation constant: α = (1+j)/δ  (derived from Eq. 27)
    std::complex<double> alpha(1.0, 1.0);
    alpha /= skinDepth;

    // Compute D = 2αh·tanh(αh/2) from Dowell's Eq. 33
    // D' = Re(D) is the proximity effect factor
    std::complex<double> ah = alpha * h;
    std::complex<double> D = 2.0 * ah * std::tanh(ah / 2.0);
    double D_prime = D.real();

    // Proximity factor: P_prox/l = factor · He²
    // Derived from Dowell's leakage impedance (Eq. 8) for a single conductor
    // of width a and height h exposed to external field He:
    //   factor = a · ρ · D' / h
    //
    // Dimensional analysis: [m] · [Ω·m] · [1] / [m] = [Ω·m]
    // Then: [Ω·m] · [A²/m²] = [W/m] ✓
    factor = a * resistivity * D_prime / h;

    if (std::isnan(factor)) {
        throw NaNResultException("NaN found in Dowell's proximity factor");
    }

    return factor;
}

double WindingProximityEffectLossesDowellModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto optionalproximityFactor = try_get_proximity_factor(wire, frequency, temperature);

    if (optionalproximityFactor) {
        proximityFactor = optionalproximityFactor.value();
    }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }

    // Use RMS of |H|² (consistent with Rossmanith, Albach, and Lammeraner models)
    double He2_sum = 0;
    for (auto& datum : data) {
        if (std::isnan(datum.get_real())) {
            throw NaNResultException("NaN found in Dowell proximity losses calculation");
        }
        if (std::isnan(datum.get_imaginary())) {
            throw NaNResultException("NaN found in Dowell proximity losses calculation");
        }
        He2_sum += pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
    }
    double He2_rms = He2_sum / data.size();

    double turnLosses = proximityFactor * He2_rms;

    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    turnLosses *= wire.get_number_conductors().value();

    if (std::isnan(turnLosses)) {
        throw NaNResultException("NaN found in Dowell's model for proximity effect losses: frequency=" +
            std::to_string(frequency) + ", proximityFactor=" + std::to_string(proximityFactor) +
            ", He2_rms=" + std::to_string(He2_rms));
    }

    return turnLosses;
}

// =============================================================================
// NAN & SULLIVAN 2003 — Modified Dowell for round conductors
// =============================================================================
double WindingProximityEffectLossesNanModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model();
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    double d = 0;
    if (wire.get_type() == WireType::ROUND) {
        d = resolve_dimensional_values(wire.get_conducting_diameter().value());
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        d = resolve_dimensional_values(strand.get_conducting_diameter());
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "Nan/Sullivan 2003 model only supports ROUND and LITZ wire");
    }

    double X = d / skinDepth;  // Eq. 16

    // Default fitted coefficients for square packing (v/d≈1, h/d≈1), Table I
    constexpr double k1 = 1.45, k2 = 0.33;
    constexpr double K  = 0.29, b  = 1.1, n = 3.0, w = 0.54;

    // Modified Dowell function G' (Eq. 15)
    double sqk2X = std::sqrt(k2) * X;
    double G_mod = k1 * std::sqrt(k2) * X
                 * (std::sinh(sqk2X) - std::sin(sqk2X))
                 / (std::cosh(sqk2X) + std::cos(sqk2X));

    // Smooth transition function g(X) (Eq. 17)
    double three_n = 3.0 * n;
    double g_smooth = K * X / std::pow(std::pow(X, -three_n) + std::pow(b, three_n), 1.0 / three_n);

    // Weighted combination (Eq. 18)
    double G = (1.0 - w) * G_mod + w * g_smooth;

    // Convert: P/l = G * H²/σ  →  factor [Ω·m] = G * ρ
    double factor = G * resistivity;

    if (std::isnan(factor))
        throw NaNResultException("NaN in Nan/Sullivan proximity factor");
    return factor;
}

double WindingProximityEffectLossesNanModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto opt = try_get_proximity_factor(wire, frequency, temperature);
    if (opt) { proximityFactor = opt.value(); }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He2_sum = 0;
    for (auto& d : data) {
        if (std::isnan(d.get_real()) || std::isnan(d.get_imaginary()))
            throw NaNResultException("NaN in field data for Nan proximity model");
        He2_sum += pow(d.get_real(), 2) + pow(d.get_imaginary(), 2);
    }
    double turnLosses = proximityFactor * He2_sum / data.size();
    if (!wire.get_number_conductors()) wire.set_number_conductors(1);
    turnLosses *= wire.get_number_conductors().value();
    if (std::isnan(turnLosses))
        throw NaNResultException("NaN in Nan model turn losses: f=" + std::to_string(frequency));
    return turnLosses;
}


// =============================================================================
// WOJDA & KAZIMIERCZUK 2012 — Field-averaging proximity model
// =============================================================================
double WindingProximityEffectLossesWojdaModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model();
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    auto constants = Constants();
    double mu0   = constants.vacuumPermeability;
    double omega  = 2.0 * std::numbers::pi * frequency;
    double factor;

    auto wt = wire.get_type();
    if (wt == WireType::FOIL || wt == WireType::PLANAR) {
        // Eq. 42: per-turn per-metre proximity factor
        // R_pe/l = ηh²·μ₀²·ω²·h / (12·ρ·b)  where ηh = h/p ≈ 1 for dense foil
        double h = resolve_dimensional_values(wire.get_conducting_height().value());
        double bw = resolve_dimensional_values(wire.get_conducting_width().value());
        // factor = μ₀²·ω²·h³ / (12·ρ·bw)
        factor = std::pow(mu0, 2) * std::pow(omega, 2) * std::pow(h, 3)
               / (12.0 * resistivity * bw);
    }
    else if (wt == WireType::RECTANGULAR) {
        // Same form as foil (Eq. 42/45):
        double h  = resolve_dimensional_values(wire.get_conducting_height().value());
        double bw = resolve_dimensional_values(wire.get_conducting_width().value());
        factor = std::pow(mu0, 2) * std::pow(omega, 2) * std::pow(h, 3)
               / (12.0 * resistivity * bw);
    }
    else if (wt == WireType::ROUND) {
        // Eq. 70: R_pe = ηb²·π²·μ₀²·ω²·Nl²·lT·d² / (576·ρ)
        // Per-conductor per-metre, single layer (Nl=1):
        //   factor = π²·μ₀²·ω²·d⁴ / (128·ρ)  (consistent with Sullivan SFD at low freq)
        double d = resolve_dimensional_values(wire.get_conducting_diameter().value());
        factor = std::pow(std::numbers::pi, 2) * std::pow(mu0, 2) * std::pow(omega, 2) * std::pow(d, 4)
               / (128.0 * resistivity);
    }
    else if (wt == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        double d = resolve_dimensional_values(strand.get_conducting_diameter());
        factor = std::pow(std::numbers::pi, 2) * std::pow(mu0, 2) * std::pow(omega, 2) * std::pow(d, 4)
               / (128.0 * resistivity);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown wire type for Wojda model");
    }

    if (std::isnan(factor))
        throw NaNResultException("NaN in Wojda proximity factor");
    return factor;
}

double WindingProximityEffectLossesWojdaModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto opt = try_get_proximity_factor(wire, frequency, temperature);
    if (opt) { proximityFactor = opt.value(); }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He2_sum = 0;
    for (auto& d : data) {
        if (std::isnan(d.get_real()) || std::isnan(d.get_imaginary()))
            throw NaNResultException("NaN in field data for Wojda proximity model");
        He2_sum += pow(d.get_real(), 2) + pow(d.get_imaginary(), 2);
    }
    double turnLosses = proximityFactor * He2_sum / data.size();
    if (!wire.get_number_conductors()) wire.set_number_conductors(1);
    turnLosses *= wire.get_number_conductors().value();
    if (std::isnan(turnLosses))
        throw NaNResultException("NaN in Wojda model turn losses: f=" + std::to_string(frequency));
    return turnLosses;
}


// =============================================================================
// SULLIVAN SFD 2001 — Squared-Field-Derivative method
// =============================================================================
double WindingProximityEffectLossesSullivanModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model();
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    auto constants = Constants();
    double mu0  = constants.vacuumPermeability;
    double omega = 2.0 * std::numbers::pi * frequency;

    double d = 0;
    if (wire.get_type() == WireType::ROUND) {
        d = resolve_dimensional_values(wire.get_conducting_diameter().value());
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        d = resolve_dimensional_values(strand.get_conducting_diameter());
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "Sullivan SFD model only supports ROUND and LITZ wire");
    }

    // Appendix A, Eq. (1): P_inst/l = (π·d⁴)/(128·ρ) · |dB/dt|²
    // For sinusoidal B=μ₀·H·sin(ωt): <|dB/dt|²> = μ₀²·ω²·H²_rms
    // → factor [Ω·m] = π·μ₀²·ω²·d⁴ / (128·ρ)
    double factor = std::numbers::pi * std::pow(mu0, 2) * std::pow(omega, 2) * std::pow(d, 4)
                  / (128.0 * resistivity);

    if (std::isnan(factor))
        throw NaNResultException("NaN in Sullivan SFD proximity factor");
    return factor;
}

double WindingProximityEffectLossesSullivanModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto opt = try_get_proximity_factor(wire, frequency, temperature);
    if (opt) { proximityFactor = opt.value(); }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He2_sum = 0;
    for (auto& d : data) {
        if (std::isnan(d.get_real()) || std::isnan(d.get_imaginary()))
            throw NaNResultException("NaN in field data for Sullivan SFD proximity model");
        He2_sum += pow(d.get_real(), 2) + pow(d.get_imaginary(), 2);
    }
    double turnLosses = proximityFactor * He2_sum / data.size();
    if (!wire.get_number_conductors()) wire.set_number_conductors(1);
    turnLosses *= wire.get_number_conductors().value();
    if (std::isnan(turnLosses))
        throw NaNResultException("NaN in Sullivan SFD model turn losses: f=" + std::to_string(frequency));
    return turnLosses;
}


// =============================================================================
// BARTOLI ET AL. 1996 — Litz-wire internal + external proximity
// =============================================================================
double WindingProximityEffectLossesBartoliModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model();
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    double d_s = 0;
    int    n_s = 1;

    if (wire.get_type() == WireType::ROUND) {
        d_s = resolve_dimensional_values(wire.get_conducting_diameter().value());
        n_s = 1;
    }
    else if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        d_s = resolve_dimensional_values(strand.get_conducting_diameter());
        n_s = wire.get_number_conductors().value_or(1);
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA,
            "Bartoli model only supports ROUND and LITZ wire");
    }

    // Kelvin-Bessel proximity variable: y_s = d_s·√2 / δ  (Eq. 4 adapted)
    double y_s = d_s * std::sqrt(2.0) / skinDepth;

    // Proximity Bessel factor K2(y):
    //   Low-freq  (y→0):    K2 ≈ y⁴/64
    //   High-freq (y→∞):    K2 ≈ y/(2√2)
    //   Smooth Padé bridge:
    double K2;
    if (y_s < 0.01) {
        K2 = std::pow(y_s, 4) / 64.0;
    }
    else {
        double low  = std::pow(y_s, 4) / 64.0;
        double high = y_s / (2.0 * std::sqrt(2.0));
        K2 = (low * high) / (low + high);   // harmonic mean → low for small y, high for large y
    }

    // External proximity factor per strand: factor_ext = 2π·ρ·K2 (Eq. 10 converted)
    double factor_ext = 2.0 * std::numbers::pi * resistivity * K2;

    double factor = factor_ext;

    // Internal proximity (strand-to-strand) for litz only (Eq. 13):
    //   Contribution ≈ k_s · K2 · factor_ext  with 50% reduction for twisted litz
    if (wire.get_type() == WireType::LITZ && n_s > 1) {
        // Packing factor k_s = n_s · d_s² / d_outer²
        double d_outer = 0;
        if (wire.get_conducting_diameter())
            d_outer = resolve_dimensional_values(wire.get_conducting_diameter().value());
        else if (wire.get_outer_diameter())
            d_outer = resolve_dimensional_values(wire.get_outer_diameter().value());

        if (d_outer > 0) {
            double k_s = n_s * std::pow(d_s, 2) / std::pow(d_outer, 2);
            // ×0.5 for twisted litz (Section II: 50% reduction of internal proximity)
            factor += 0.5 * k_s * factor_ext;
        }
    }

    if (std::isnan(factor))
        throw NaNResultException("NaN in Bartoli proximity factor");
    return factor;
}

double WindingProximityEffectLossesBartoliModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto opt = try_get_proximity_factor(wire, frequency, temperature);
    if (opt) { proximityFactor = opt.value(); }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He2_sum = 0;
    for (auto& d : data) {
        if (std::isnan(d.get_real()) || std::isnan(d.get_imaginary()))
            throw NaNResultException("NaN in field data for Bartoli proximity model");
        He2_sum += pow(d.get_real(), 2) + pow(d.get_imaginary(), 2);
    }
    double turnLosses = proximityFactor * He2_sum / data.size();
    if (!wire.get_number_conductors()) wire.set_number_conductors(1);
    turnLosses *= wire.get_number_conductors().value();
    if (std::isnan(turnLosses))
        throw NaNResultException("NaN in Bartoli model turn losses: f=" + std::to_string(frequency));
    return turnLosses;
}


// =============================================================================
// VANDELAC & ZIOGAS 1988 — F1/F2 proximity model (α=1 case)
// =============================================================================
double WindingProximityEffectLossesVandelacModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    auto& resistivityModel = get_cached_resistivity_model();
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    double h = 0, a = 0;
    auto wt = wire.get_type();

    if (wt == WireType::FOIL || wt == WireType::PLANAR || wt == WireType::RECTANGULAR) {
        a = resolve_dimensional_values(wire.get_conducting_width().value());
        h = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    else if (wt == WireType::ROUND) {
        // Eq. 12(b)/(c): round→square equivalence with same cross-section
        double d = resolve_dimensional_values(wire.get_conducting_diameter().value());
        a = d * std::sqrt(std::numbers::pi / 4.0);
        h = a;
    }
    else if (wt == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        double d = resolve_dimensional_values(strand.get_conducting_diameter());
        a = d * std::sqrt(std::numbers::pi / 4.0);
        h = a;
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown wire type for Vandelac model");
    }

    // p = h / (δ·√2)   (Eq. 27, using δ = δ_actual)
    double p = h / (skinDepth * std::sqrt(2.0));

    // F1(p) = [sinh(2p) + sin(2p)] / [cosh(2p) - cos(2p)]  (Eq. 29)
    // F2(p) = [sinh(p)·cos(p) + cosh(p)·sin(p)] / [cosh(2p) - cos(2p)]  (Eq. 30)
    double denom, F1, F2;
    if (p < 1e-4) {
        // Taylor series: F1 → 1, F2 → 1/2  as p→0
        F1 = 1.0;
        F2 = 0.5;
    }
    else {
        denom = std::cosh(2.0 * p) - std::cos(2.0 * p);
        F1 = (std::sinh(2.0 * p) + std::sin(2.0 * p)) / denom;
        F2 = (std::sinh(p) * std::cos(p) + std::cosh(p) * std::sin(p)) / denom;
    }

    // Proximity-only case (α=1, β=0), from Eq. 25 differentiated:
    //   Q_prox = H²_e·ρ / δ · [3·F1 - 4·F2] / 2
    // Per-unit-length for conductor width a:
    //   P_prox/l = Q_prox · a = H²_e · (a·ρ·(3·F1-4·F2)) / (2·δ)
    //
    // Note: at α=1 the total loss from Eq. 25 is:
    //   Q_total = H²/(2·σ·δ) · [(1+α²)·F1 - 2α·F2]
    //           = H²·ρ/δ · [(1+1)·F1 - 2·F2] / 2
    //           = H²·ρ/δ · [2·F1 - 2·F2] / 2
    //           = H²·ρ/δ · (F1 - F2)     [for field on BOTH surfaces]
    // But in our architecture H is the EXTERNAL field (one side only, proximity only):
    //   factor = a · ρ · (3·F1 - 4·F2) / (2·δ)
    double factor = a * resistivity * (3.0 * F1 - 4.0 * F2) / (2.0 * skinDepth);

    // Guard against negative factor at extreme p
    if (factor < 0) factor = 0;

    if (std::isnan(factor))
        throw NaNResultException("NaN in Vandelac proximity factor");
    return factor;
}

double WindingProximityEffectLossesVandelacModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    double proximityFactor;
    auto opt = try_get_proximity_factor(wire, frequency, temperature);
    if (opt) { proximityFactor = opt.value(); }
    else {
        proximityFactor = calculate_proximity_factor(wire, frequency, temperature);
        set_proximity_factor(wire, frequency, temperature, proximityFactor);
    }
    double He2_sum = 0;
    for (auto& d : data) {
        if (std::isnan(d.get_real()) || std::isnan(d.get_imaginary()))
            throw NaNResultException("NaN in field data for Vandelac proximity model");
        He2_sum += pow(d.get_real(), 2) + pow(d.get_imaginary(), 2);
    }
    double turnLosses = proximityFactor * He2_sum / data.size();
    if (!wire.get_number_conductors()) wire.set_number_conductors(1);
    turnLosses *= wire.get_number_conductors().value();
    if (std::isnan(turnLosses))
        throw NaNResultException("NaN in Vandelac model turn losses: f=" + std::to_string(frequency));
    return turnLosses;
}


} // namespace OpenMagnetics

