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
    else
        throw ModelNotAvailableException("Unknown wire proximity effect losses mode, available options are: {ROSSMANITH, WANG, FERREIRA, ALBACH, LAMMERANER}");
}

std::shared_ptr<WindingProximityEffectLossesModel> WindingProximityEffectLosses::get_model(WireType wireType, std::optional<WindingProximityEffectLossesModels> modelOverride) {
    if (modelOverride.has_value()) {
        return WindingProximityEffectLossesModel::factory(modelOverride.value());
    }
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


} // namespace OpenMagnetics

