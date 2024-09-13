#include "MagneticField.h"
#include "WireWrapper.h"
#include "MagneticWrapper.h"
#include "CoilMesher.h"
#include "Models.h"
#include "Reluctance.h"
#include "MagneticSimulator.h"
#include "InitialPermeability.h"
#include "Settings.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {

SignalDescriptor MagneticField::calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                          double reluctance,
                                                          double numberTurns) {
    SignalDescriptor magneticFlux;
    Waveform magneticFluxWaveform;
    std::vector<double> magneticFluxData;
    auto compressedMagnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();

    if (InputsWrapper::is_waveform_sampled(compressedMagnetizingCurrentWaveform)) {
        compressedMagnetizingCurrentWaveform = InputsWrapper::compress_waveform(compressedMagnetizingCurrentWaveform);
    }

    for (auto& datum : compressedMagnetizingCurrentWaveform.get_data()) {
        magneticFluxData.push_back(datum * numberTurns / reluctance);
    }

    if (compressedMagnetizingCurrentWaveform.get_time()) {
        magneticFluxWaveform.set_time(compressedMagnetizingCurrentWaveform.get_time());
    }

    magneticFluxWaveform.set_data(magneticFluxData);
    magneticFlux.set_waveform(magneticFluxWaveform);
    if (magnetizingCurrent.get_harmonics()) {
        auto harmonics = magnetizingCurrent.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] *= numberTurns / reluctance;
        }
        magneticFlux.set_harmonics(harmonics);
    }
    return magneticFlux;
}
SignalDescriptor MagneticField::calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                  double area) {
    SignalDescriptor magneticFluxDensity;
    Waveform magneticFluxDensityWaveform;
    std::vector<double> magneticFluxDensityData;
    auto magneticFluxWaveform = magneticFlux.get_waveform().value();

    if (magneticFluxWaveform.get_time()) {
        magneticFluxDensityWaveform.set_time(magneticFluxWaveform.get_time());
    }

    for (auto& datum : magneticFluxWaveform.get_data()) {
        magneticFluxDensityData.push_back(datum / area);
    }

    magneticFluxDensityWaveform.set_data(magneticFluxDensityData);
    magneticFluxDensity.set_waveform(magneticFluxDensityWaveform);
    if (magneticFlux.get_harmonics()) {
        auto harmonics = magneticFlux.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= area;
        }
        magneticFluxDensity.set_harmonics(harmonics);
    }
    magneticFluxDensity.set_processed(
        InputsWrapper::calculate_basic_processed_data(magneticFluxDensityWaveform));

    return magneticFluxDensity;
}

SignalDescriptor MagneticField::calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability) {
    SignalDescriptor magneticFieldStrength;
    Waveform magneticFieldStrengthWaveform;
    std::vector<double> magneticFieldStrengthData;
    auto constants = Constants();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();

    if (magneticFluxDensityWaveform.get_time()) {
        magneticFieldStrengthWaveform.set_time(magneticFluxDensityWaveform.get_time());
    }

    for (auto& datum : magneticFluxDensityWaveform.get_data()) {
        magneticFieldStrengthData.push_back(datum / (initialPermeability * constants.vacuumPermeability));
    }

    magneticFieldStrengthWaveform.set_data(magneticFieldStrengthData);
    magneticFieldStrength.set_waveform(magneticFieldStrengthWaveform);
    if (magneticFluxDensity.get_harmonics()) {
        auto harmonics = magneticFluxDensity.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= (initialPermeability * constants.vacuumPermeability);
        }
        magneticFieldStrength.set_harmonics(harmonics);
    }
    magneticFieldStrength.set_processed(
        InputsWrapper::calculate_basic_processed_data(magneticFieldStrengthWaveform));

    return magneticFieldStrength;
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory(MagneticFieldStrengthModels modelName) {
    if (modelName == MagneticFieldStrengthModels::BINNS_LAWRENSON) {
        return std::make_shared<MagneticFieldStrengthBinnsLawrensonModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::LAMMERANER) {
        return std::make_shared<MagneticFieldStrengthLammeranerModel>();
    }
    else
        throw std::runtime_error("Unknown Magnetic Field Strength model, available options are: {BINNS_LAWRENSON, LAMMERANER}");
}

std::shared_ptr<MagneticFieldStrengthFringingEffectModel> MagneticField::factory(MagneticFieldStrengthFringingEffectModels modelName) {
    if (modelName == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        return std::make_shared<MagneticFieldStrengthAlbachModel>();
    }
    else if (modelName == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
        return std::make_shared<MagneticFieldStrengthRoshenModel>();
    }
    else
        throw std::runtime_error("Unknown Magnetic Field Strength Fringing Effect model, available options are: {ALBACH, ROSHEN}");
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory() {
    auto defaults = Defaults();
    return factory(defaults.magneticFieldStrengthModelDefault);
}

bool is_inside_inducing_turns(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, WireWrapper inducingWire) {
    double distanceX = fabs(inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0]);
    double distanceY = fabs(inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1]);
    if (inducingWire.get_type() == WireType::ROUND || inducingWire.get_type() == WireType::LITZ) {
        if (hypot(distanceX, distanceY) < inducingWire.get_maximum_outer_width() / 2) {
            return true;
        }
    }
    else {
        if (distanceX < inducingWire.get_maximum_outer_width() / 2 && distanceY < inducingWire.get_maximum_outer_height() / 2) {
            return true;
        }
    }

    return false;
}

bool is_inside_core(FieldPoint inducedFieldPoint, CoreWrapper core) {
    if (core.get_shape_family() != CoreShapeFamily::T) {
        return false;
    }
    double radius = sqrt(pow(inducedFieldPoint.get_point()[0], 2) + pow(inducedFieldPoint.get_point()[1], 2));
    double coreColumnWidth = core.get_columns()[0].get_width();
    auto processedDescription = core.get_processed_description().value();
    double coreWidth = processedDescription.get_width();

    if (radius * 1.05 > coreWidth / 2) {
        return false;
    }
    if (radius * 0.95 < (coreWidth / 2 - coreColumnWidth)) {
        return false;
    }
    return true;
}

double get_magnetic_field_strength_gap(OperatingPoint operatingPoint, MagneticWrapper magnetic, double frequency) {
    auto numberTurns = magnetic.get_mutable_coil().get_number_turns(0);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    OpenMagnetics::InitialPermeability initial_permeability;
    double initialPermeability = initial_permeability.get_initial_permeability(magnetic.get_mutable_core().resolve_material(), std::nullopt, std::nullopt, frequency);
    double reluctance = reluctanceModel->get_core_reluctance(magnetic.get_core(), initialPermeability).get_core_reluctance();
    if (!operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current()) {
        throw std::runtime_error("Magnetizing current is missing");
    }

    auto magnetizingCurrent = operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value();
    if (!magnetizingCurrent.get_waveform()) {
        throw std::runtime_error("Magnetizing current is missing waveform");
    }
    if (!magnetizingCurrent.get_waveform()->get_time()) {
        magnetizingCurrent = InputsWrapper::standarize_waveform(magnetizingCurrent, frequency);
    }
    auto magneticFlux = MagneticField::calculate_magnetic_flux(magnetizingCurrent, reluctance, numberTurns);
    auto magneticFluxDensity = MagneticField::calculate_magnetic_flux_density(magneticFlux, magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area());

    double magneticFieldStrengthGap = magneticFluxDensity.get_processed()->get_peak().value() / Constants().vacuumPermeability;
    return magneticFieldStrengthGap;
}

WindingWindowMagneticStrengthFieldOutput MagneticField::calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, std::optional<Field> externalInducedField, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto includeFringing = settings->get_magnetic_field_include_fringing();
    CoilMesher coilMesher; 
    std::vector<Field> inducingFields;

    std::vector<int8_t> currentDirectionPerWinding;
    if (!customCurrentDirectionPerWinding) {
        currentDirectionPerWinding.push_back(1);
        for (size_t windingIndex = 1; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
            currentDirectionPerWinding.push_back(-1);
        }
    }
    else {
        currentDirectionPerWinding = customCurrentDirectionPerWinding.value();
    }

    if (externalInducedField){
        auto aux = coilMesher.generate_mesh_inducing_coil(magnetic, operatingPoint, 0, currentDirectionPerWinding);
        // We only process the harmonic that comes from the external field
        for (auto field : aux) {
            if (field.get_frequency() == externalInducedField.value().get_frequency()) {
                inducingFields.push_back(field);
                break;
            }
        }
    }
    else {
        inducingFields = coilMesher.generate_mesh_inducing_coil(magnetic, operatingPoint, settings->get_harmonic_amplitude_threshold(), currentDirectionPerWinding);
    }

    if (!magnetic.get_coil().get_turns_description()) {
        throw std::runtime_error("Missing turns description in coil");
    }
    auto wirePerWinding = magnetic.get_mutable_coil().get_wires();
    auto turns = magnetic.get_coil().get_turns_description().value();

    std::vector<ComplexField> complexFieldPerHarmonic;
    for (auto& fieldPerHarmonic : inducingFields) {
        ComplexField field;
        field.set_frequency(fieldPerHarmonic.get_frequency());
        complexFieldPerHarmonic.push_back(field);
    }

    std::vector<Field> inducedFields;
    if (externalInducedField) {
        inducedFields.push_back(externalInducedField.value());
    }
    else {
        inducedFields = coilMesher.generate_mesh_induced_coil(magnetic, operatingPoint, settings->get_harmonic_amplitude_threshold());
    }

    if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        if (includeFringing) {
            if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {

                auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
                auto magnetizingCurrent = InputsWrapper::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                                       resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()));

                operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                // throw std::runtime_error("Operating point is missing magnetizing current");
            }
            if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()) {
                auto excitations = operatingPoint.get_excitations_per_winding();
                auto magnetizingCurrent = excitations[0].get_magnetizing_current().value();
                auto processed = InputsWrapper::calculate_basic_processed_data(magnetizingCurrent.get_waveform().value());
                magnetizingCurrent.set_processed(processed);
                excitations[0].set_magnetizing_current(magnetizingCurrent);
                operatingPoint.set_excitations_per_winding(excitations);
                // throw std::runtime_error("Operating point is missing magnetizing current processed data");

            }

            for (size_t harmonicIndex = 0; harmonicIndex < inducingFields.size(); ++harmonicIndex){

                if (inducingFields[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {

                    double frequency = inducingFields[harmonicIndex].get_frequency();
                    double magneticFieldStrengthGap = get_magnetic_field_strength_gap(operatingPoint, magnetic, frequency);
                    
                    for (auto& gap : magnetic.get_core().get_functional_description().get_gapping()) {
                        auto fieldPoint = _fringingEffectModel->get_equivalent_inducing_point_for_gap(gap, magneticFieldStrengthGap);
                        inducingFields[harmonicIndex].get_mutable_data().push_back(fieldPoint);
                    }
                }
            }
        }
    }

    for (size_t harmonicIndex = 0; harmonicIndex < inducingFields.size(); ++harmonicIndex){
        std::vector<ComplexFieldPoint> fieldPoints;

        if (inducedFields[harmonicIndex].get_data().size() == 0) {
            throw std::runtime_error("Empty complexField");
        }

        for (auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
            double totalInducedFieldX = 0;
            double totalInducedFieldY = 0;

            if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
                // For the main harmonic we calculate the fringing effect for each gap
                if (includeFringing && inducedFields[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {
                    if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
                        auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
                        auto magnetizingCurrent = InputsWrapper::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                                               resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()));

                        operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                        // throw std::runtime_error("Operating point is missing magnetizing current");
                    }
                    if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()) {
                        auto excitations = operatingPoint.get_excitations_per_winding();
                        auto magnetizingCurrent = excitations[0].get_magnetizing_current().value();
                        auto processed = InputsWrapper::calculate_basic_processed_data(magnetizingCurrent.get_waveform().value());
                        magnetizingCurrent.set_processed(processed);
                        excitations[0].set_magnetizing_current(magnetizingCurrent);
                        operatingPoint.set_excitations_per_winding(excitations);
                        // throw std::runtime_error("Operating point is missing magnetizing current processed data");
                    }
                    double frequency = inducingFields[harmonicIndex].get_frequency();
                    double magneticFieldStrengthGap = get_magnetic_field_strength_gap(operatingPoint, magnetic, frequency);
                    for (auto& gap : magnetic.get_core().get_functional_description().get_gapping()) {
                        if (gap.get_coordinates().value()[0] < 0) {
                            continue;
                        }
                        auto complexFieldPoint = _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, inducedFieldPoint);

                        totalInducedFieldX += complexFieldPoint.get_real();
                        totalInducedFieldY += complexFieldPoint.get_imaginary();
                    }
                }
            }

            for (auto& inducingFieldPoint : inducingFields[harmonicIndex].get_data()) {
                std::optional<WireWrapper> inducingWire;
                if (inducingFieldPoint.get_turn_index()) {
                    auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turns[inducingFieldPoint.get_turn_index().value()].get_winding());
                    inducingWire = wirePerWinding[windingIndex];
                    if (inducedFieldPoint.get_turn_index()) {
                        if (inducedFieldPoint.get_turn_index().value() == inducingFieldPoint.get_turn_index().value()) {
                            continue;
                        }
                    }
                    // else if (is_inside_inducing_turns(inducingFieldPoint, inducedFieldPoint, inducingWire.value())) {
                    //     continue;
                    // }
                    else if (is_inside_core(inducedFieldPoint, magnetic.get_core())) {
                        continue;
                    }
                }

                auto complexFieldPoint = _model->get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint, inducingWire);

                totalInducedFieldX += complexFieldPoint.get_real();
                totalInducedFieldY += complexFieldPoint.get_imaginary();
            }
            ComplexFieldPoint complexFieldPoint;
            complexFieldPoint.set_point(inducedFieldPoint.get_point());
            complexFieldPoint.set_real(totalInducedFieldX);
            complexFieldPoint.set_imaginary(totalInducedFieldY);
            if (inducedFieldPoint.get_turn_index()) {
                complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
            }
            if (inducedFieldPoint.get_label()) {
                complexFieldPoint.set_label(inducedFieldPoint.get_label().value());
            }
            fieldPoints.push_back(complexFieldPoint);
        }
        complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
    }

    WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput;
    windingWindowMagneticStrengthFieldOutput.set_field_per_frequency(complexFieldPerHarmonic);
    windingWindowMagneticStrengthFieldOutput.set_method_used(std::string{magic_enum::enum_name(_magneticFieldStrengthModel)});
    windingWindowMagneticStrengthFieldOutput.set_origin(ResultOrigin::SIMULATION);
    return windingWindowMagneticStrengthFieldOutput;
}

ComplexFieldPoint MagneticFieldStrengthBinnsLawrensonModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire) {
    double Hx;
    double Hy;
    if (!inducingWire || inducingWire.value().get_type() == WireType::ROUND || inducingWire.value().get_type() == WireType::LITZ) {

        double distanceX = inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0];
        double distanceY = inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1];
        if (hypot(distanceX, distanceY) < inducingWire.value().get_maximum_outer_width() / 2) {
            Hx = 0;
            Hy = 0;
        }
        else { 
            double divisor = 2 * std::numbers::pi * (pow(distanceY, 2) + pow(distanceX, 2));
            Hx = -inducingFieldPoint.get_value() * (distanceY) / divisor;
            Hy = inducingFieldPoint.get_value() * (distanceX) / divisor;
        }
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw std::runtime_error("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }
    else {
        auto wire = inducingWire.value();
        double a = resolve_dimensional_values(wire.get_conducting_width().value()) / 2;
        double b = resolve_dimensional_values(wire.get_conducting_height().value()) / 2;
        double x = inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0];
        double y = inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1];

        if (inducingFieldPoint.get_rotation()) {
            double modulo = hypot(x, y);
            double currentAngle = atan2(y, x);
            double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
            if (currentAngle < 0) {
                currentAngle += 2 * std::numbers::pi;
            }
            double totalAngle = currentAngle - turnAngle;

            x = modulo * cos(totalAngle);
            y = modulo * sin(totalAngle);
        }

        double r1 = hypot(y + b, x - a);
        double r2 = hypot(y + b, x + a);
        double r3 = hypot(y - b, x + a);
        double r4 = hypot(y - b, x - a);

        double tetha1 = atan((y + b) / (x - a));
        double tetha2 = atan((y + b) / (x + a));
        double tetha3 = atan((y - b) / (x + a));
        double tetha4 = atan((y - b) / (x - a));
        if (fabs(x) < a && fabs(y) < b) {
            Hx = 0;
            Hy = 0;
        }
        else if (std::isnan(tetha1) || std::isnan(tetha2) || std::isnan(tetha3) || std::isnan(tetha4)) {
            Hx = 0;
            Hy = 0;
        } 
        else {
            if (x == a) {
                if ((y + b) > 0) {
                    tetha1 = std::numbers::pi / 2;
                }
                else {
                    tetha1 = -std::numbers::pi / 2;
                }
                if ((y - b) > 0) {
                    tetha4 = std::numbers::pi / 2;
                }
                else {
                    tetha4 = -std::numbers::pi / 2;
                }
            }

            if (x > a && -b < y && y < b) {

            }
            else {
                if (x > a && y < -b) {
                    tetha1 += 2 * std::numbers::pi;
                }
                else if (x < a || y < -b) {
                    tetha1 += std::numbers::pi;
                }

                if (x > -a && y < -b) {
                    tetha2 += 2 * std::numbers::pi;
                }
                else if (x < -a || y < -b) {
                    tetha2 += std::numbers::pi;
                }

                if (x > -a && y < b) {
                    tetha3 += 2 * std::numbers::pi;
                }
                else if (x < -a || y < b) {
                    tetha3 += std::numbers::pi;
                }

                if (x > a && y < b) {
                    tetha4 += 2 * std::numbers::pi;
                }
                else if (x < a || y < b) {
                    tetha4 += std::numbers::pi;
                }
            }
        }

        double common_part = inducingFieldPoint.get_value() / (8.0 * std::numbers::pi * a * b);
        Hx = common_part * ((y + b) * (tetha1 - tetha2) - (y - b) * (tetha4 - tetha3) + (x + a) * log(r2 / r3) - (x - a) * log(r1 / r4));

        Hy = -common_part * ((x + a) * (tetha2 - tetha3) - (x - a) * (tetha1 - tetha4) + (y + b) * log(r2 / r1) - (y - b) * log(r3 / r4));
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw std::runtime_error("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }


    if (inducingFieldPoint.get_rotation()) {
        double modulo = hypot(Hx, Hy);
        double currentAngle = atan2(Hy, Hx);
        if (currentAngle < 0) {
            currentAngle += 2 * std::numbers::pi;
        }
        double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
        double totalAngle = currentAngle + turnAngle;
        Hx = modulo * cos(totalAngle);
        Hy = modulo * sin(totalAngle);
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw std::runtime_error("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;   
}

ComplexFieldPoint MagneticFieldStrengthLammeranerModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire) {
    double Hx;
    double Hy;
    if (!inducingWire || inducingWire.value().get_type() == WireType::ROUND || inducingWire.value().get_type() == WireType::LITZ) {
        double turnLength = 1;
        if (inducingFieldPoint.get_turn_length()) {
            turnLength = inducingFieldPoint.get_turn_length().value();
        }
        double distance = hypot(inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1], inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0]);
        double angle = atan2(inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0], inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1]);
        if (distance < inducingWire.value().get_maximum_outer_width() / 2) {
            Hx = 0;
            Hy = 0;
        }
        else {
            double ex = cos(angle - std::numbers::pi / 2);
            double ey = sin(angle - std::numbers::pi / 2);
            double magneticFiledStrengthModule = -inducingFieldPoint.get_value() / 2 / std::numbers::pi / distance * turnLength / hypot(turnLength, distance);
            Hx = magneticFiledStrengthModule * ex;
            Hy = magneticFiledStrengthModule * ey;
        }
    }
    else {
        throw std::runtime_error("Rectangular wires not implemented yet");
    }
    if (std::isnan(Hx) || std::isnan(Hy)) {
        throw std::runtime_error("NaN found in Lammeraner's model for magnetic field");
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;
}



FieldPoint MagneticFieldStrengthAlbachModel::get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap) {
    if (!gap.get_section_dimensions()) {
        throw std::runtime_error("Gap is missing section dimensions");
    }
    if (!gap.get_coordinates()) {
        throw std::runtime_error("Gap is missing coordinates");
    }
    double rc = gap.get_section_dimensions().value()[0] / 2;
    double xi = gap.get_length() / (2 * rc);
    double x = 1 - 1.05 * xi - 2.88 * pow(xi, 2) - 8.8 * pow(xi, 3);
    if (x < 0) {
        throw std::runtime_error("Something went wrong with Albach method with x");
    }
    double current = (magneticFieldStrengthGap * gap.get_length()) / (0.25 - 1.569 * xi + 4.34 * pow(xi, 2) - 7.042 * pow(xi, 3));
    if (current < 0) {
        throw std::runtime_error("Something went wrong with Albach method with current");
    }
    double eta = x * rc;

    if (eta > gap.get_section_dimensions().value()[0] / 2) {
        throw std::runtime_error("Something went wrong with Albach method with eta");
    }
    FieldPoint fieldPoint;
    if (gap.get_coordinates().value()[0] > 0) {
        fieldPoint.set_point({gap.get_coordinates().value()[0] - eta, gap.get_coordinates().value()[1]});
    }
    else {
        fieldPoint.set_point({gap.get_coordinates().value()[0] + eta, gap.get_coordinates().value()[1]});
    }
    fieldPoint.set_point({eta, gap.get_coordinates().value()[1]});
    fieldPoint.set_value(current);
    return fieldPoint;
}

ComplexFieldPoint MagneticFieldStrengthRoshenModel::get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) {
    double distanceFromCenterEdgeGapX;
    if (gap.get_coordinates().value()[0] == 0) {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] + gap.get_section_dimensions().value()[0] / 2);
    }
    else {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] - gap.get_section_dimensions().value()[0] / 2);
    }
    double distanceFromCenterEdgeGapY = inducedFieldPoint.get_point()[1] - gap.get_coordinates().value()[1];
    double halfGapLength = gap.get_length() / 2;

    double magneticIntensityXDividend = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY - halfGapLength, 2);
    double magneticIntensityXDivisor = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY + halfGapLength, 2);
    double Hx = -magneticFieldStrengthGap / 2 / std::numbers::pi * log(magneticIntensityXDividend / magneticIntensityXDivisor);

    double m;
    if (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) > pow(halfGapLength, 2)) {
        m = 0;
    }
    else {
        m = 1;
    }

    double x = distanceFromCenterEdgeGapX * halfGapLength / (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) - pow(halfGapLength, 2));
    double Hy = -magneticFieldStrengthGap / std::numbers::pi * (atan(x) + m * std::numbers::pi);

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}


} // namespace OpenMagnetics
