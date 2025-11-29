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

namespace OpenMagnetics {

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
        throw std::runtime_error("Unknown wire proximity effect losses mode, available options are: {ROSSMANITH, WANG, FERREIRA, ALBACH, LAMMERANER}");
}

std::shared_ptr<WindingProximityEffectLossesModel> WindingProximityEffectLosses::get_model(WireType wireType) {
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
            throw std::runtime_error("Unknown type of wire");
    }
}

std::optional<double> WindingProximityEffectLossesModel::try_get_proximity_factor(Wire wire, double frequency, double temperature) {
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
        hash = std::hash<double>{}(wire.get_number_conductors().value() * wire.get_maximum_outer_width() * wire.get_maximum_outer_height() );
    }
    else {
        hash = std::hash<std::string>{}(wire.get_name().value());
    }

    _proximityFactorPerWirePerFrequencyPerTemperature[hash][frequency][temperature] = proximityFactor;

}

std::pair<double, std::vector<std::pair<double, double>>> WindingProximityEffectLosses::calculate_proximity_effect_losses_per_meter(Wire wire, double temperature, std::vector<ComplexField> fields) {
    auto model = get_model(wire.get_type());
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
            throw std::runtime_error("NaN found in proximity effect losses per meter");
        }
        lossesPerHarmonic.push_back(std::pair<double, double>{turnLosses, frequency});
        totalProximityEffectLossesPerMeter += turnLosses;
    }

    return {totalProximityEffectLossesPerMeter, lossesPerHarmonic};
}

WindingLossesOutput WindingProximityEffectLosses::calculate_proximity_effect_losses(Coil coil, double temperature, WindingLossesOutput windingLossesOutput, WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput) {
    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding does not have turns description");
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
                    throw std::runtime_error("Missing turn index in field point");
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

        auto lossesPerHarmonicThisTurn = calculate_proximity_effect_losses_per_meter(wire, temperature, fields).second;


        WindingLossElement proximityEffectLossesThisTurn;
        auto model = get_model(coil.get_wire_type(windingIndex));
        proximityEffectLossesThisTurn.set_method_used(model->methodName);
        proximityEffectLossesThisTurn.set_origin(ResultOrigin::SIMULATION);
        proximityEffectLossesThisTurn.get_mutable_harmonic_frequencies().push_back(0);
        proximityEffectLossesThisTurn.get_mutable_losses_per_harmonic().push_back(0);

        for (auto& lossesThisHarmonic : lossesPerHarmonicThisTurn) {
            if (std::isnan(lossesThisHarmonic.first)) {
                throw std::runtime_error("NaN found in proximity effect losses");
            }
            proximityEffectLossesThisTurn.get_mutable_harmonic_frequencies().push_back(lossesThisHarmonic.second);
            proximityEffectLossesThisTurn.get_mutable_losses_per_harmonic().push_back(lossesThisHarmonic.first * wireLength);

            totalProximityEffectLosses += lossesThisHarmonic.first * wireLength;

        }

        windingLossesPerTurn[turnIndex].set_proximity_effect_losses(proximityEffectLossesThisTurn);
    }
    windingLossesOutput.set_winding_losses_per_turn(windingLossesPerTurn);

    windingLossesOutput.set_method_used("AnalyticalModels");
    windingLossesOutput.set_winding_losses(windingLossesOutput.get_winding_losses() + totalProximityEffectLosses);
    return windingLossesOutput;
}

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
        throw std::runtime_error("Unknown type of wire");
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

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);

    double He = 0;
    for (auto& datum : data) {
        He += sqrt(pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2));
    }
    He /= data.size();
    double turnLosses = resistivity * pow(He, 2) * proximityFactor;

    turnLosses *= wire.get_number_conductors().value();

    return turnLosses;
}

double WindingProximityEffectLossesWangModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double c = 0;
    double h = 0;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {

        if (!wire.get_conducting_width()) {
            throw std::runtime_error("Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw std::runtime_error("Missing conducting height in wire");
        }
        c = resolve_dimensional_values(wire.get_conducting_width().value());
        h = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    else {
        throw std::runtime_error("Model not implemented for ROUND and LITZ");
    }

    double Hx1 = 0, Hx2 = 0, Hy1 = 0, Hy2 = 0;
    double nonPlanarHe = 0;
    for (auto& datum : data) {
        if (!datum.get_label()) {
            throw std::runtime_error("Missing label in induced point");
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


    double turnLosses = c * h * resistivity / skinDepth * pow((Hx2 + Hx1) / 2, 2) * (sinh(h / skinDepth) - sin(h / skinDepth)) / (cosh(h / skinDepth) + cos(h / skinDepth));
    turnLosses += h * c * resistivity / skinDepth * pow((Hy2 + Hy1) / 2, 2) * (sinh(c / skinDepth) - sin(c / skinDepth)) / (cosh(c / skinDepth) + cos(c / skinDepth));

    // For the H field not planned for in Wang's model, we use Ferreira's
    if (nonPlanarHe != 0) {
        double proximityFactor = WindingProximityEffectLossesFerreiraModel::calculate_proximity_factor(wire, frequency, temperature);
        turnLosses += proximityFactor * pow(nonPlanarHe, 2);
    }

    turnLosses *= wire.get_number_conductors().value();

    return turnLosses;
}


double WindingProximityEffectLossesFerreiraModel::calculate_proximity_factor(Wire wire, double frequency, double temperature) {
    double factor;
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (!wire.get_conducting_width()) {
            throw std::runtime_error("Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw std::runtime_error("Missing conducting height in wire");
        }
        double w = resolve_dimensional_values(wire.get_conducting_width().value());
        double h = resolve_dimensional_values(wire.get_conducting_height().value());

        double xi = std::min(h, w) / skinDepth;

        factor = w * xi * resistivity * (sinh(xi) - sin(xi)) / (cosh(xi) + cos(xi));
        if (std::isnan(factor)) {
            throw std::runtime_error("NaN found in Ferreira's proximity factor");
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
        throw std::runtime_error("Unknown type of wire");
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
            throw std::runtime_error("NaN found in Ferreira proximity losses calculation");
        }
        if (std::isnan(datum.get_imaginary())) {
            throw std::runtime_error("NaN found in Ferreira proximity losses calculation");
        }
        He = std::max(He, hypot(datum.get_real(), datum.get_imaginary()));
    }

    double turnLosses = proximityFactor * pow(He, 2);
    if (!wire.get_number_conductors()) {
        wire.set_number_conductors(1);
    }
    turnLosses *= wire.get_number_conductors().value();
        if (std::isnan(turnLosses)) {
            std::cout << "frequency:" << frequency << std::endl;
            std::cout << "proximityFactor:" << proximityFactor << std::endl;
            std::cout << "He:" << He << std::endl;
            throw std::runtime_error("NaN found in Ferreira proximity losses calculation");
        }

    return turnLosses;
}

double WindingProximityEffectLossesAlbachModel::calculate_turn_losses(Wire wire, double frequency, std::vector<ComplexFieldPoint> data, double temperature) {
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(wire.resolve_material(), temperature);
    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
    double d = 0;
    double c = 0;

    if (wire.get_type() == WireType::PLANAR || wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
        if (!wire.get_conducting_width()) {
            throw std::runtime_error("Missing conducting width in wire");
        }
        if (!wire.get_conducting_height()) {
            throw std::runtime_error("Missing conducting height in wire");
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

    double He = 0;

    for (auto& datum : data) {
        He += hypot(datum.get_real(), datum.get_imaginary());
    }
    He /= data.size();
    double turnLosses = c * resistivity * pow(He, 2) * (alpha * d * tanh(alpha * d / 2.0)).real();

    turnLosses *= wire.get_number_conductors().value();

    if (std::isnan(turnLosses)) {
        throw std::runtime_error("NaN found in Albach's model for proximity effect losses");
    }

    return turnLosses;
}

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
        throw std::runtime_error("Unknown type of wire");
    }

    auto resistivityModel = ResistivityModel::factory(ResistivityModels::WIRE_MATERIAL);
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

    double Hx = 0;
    double Hy = 0;
    for (auto& datum : data) {
        Hx += datum.get_real();
        Hy += datum.get_imaginary();
    }
    Hx /= data.size();
    Hy /= data.size();

    double turnLosses = (pow(Hx, 2) + pow(Hy, 2)) * proximityFactor;
    turnLosses *= wire.get_number_conductors().value();

    if (std::isnan(turnLosses)) {
        throw std::runtime_error("NaN found in Lammeraner's model for proximity effect losses");
    }

    return turnLosses;
}


} // namespace OpenMagnetics

