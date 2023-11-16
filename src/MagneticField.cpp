#include "MagneticField.h"
#include "WireWrapper.h"
#include "MagneticWrapper.h"
#include "CoilBreaker.h"

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
                                                          double numberTurns,
                                                          double frequency) {
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

    return magneticFlux;
}
SignalDescriptor MagneticField::calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                  double area,
                                                                  double frequency) {
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
    magneticFluxDensity.set_processed(
        InputsWrapper::calculate_basic_processed_data(magneticFluxDensityWaveform));

    return magneticFluxDensity;
}

SignalDescriptor MagneticField::calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability,
                                                                    double frequency) {
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
    magneticFieldStrength.set_processed(
        InputsWrapper::calculate_basic_processed_data(magneticFieldStrengthWaveform));

    return magneticFieldStrength;
}

bool is_inside_inducing_turns(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, WireWrapper inducingWire) {
    double distanceX = fabs(inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0]);
    double distanceY = fabs(inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1]);
    if (inducingWire.get_type() == WireType::ROUND || inducingWire.get_type() == WireType::LITZ) {
        if (hypot(distanceX, distanceY) < inducingWire.get_maximum_width() / 2) {
            return true;
        }
    }
    else {
        if (distanceX < inducingWire.get_maximum_width() / 2 && distanceY < inducingWire.get_maximum_height() / 2) {
            return true;
        }
    }

    return false;
}

WindingWindowMagneticStrengthFieldOutput MagneticField::calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, std::optional<std::vector<FieldPoint>> externalInducedFieldPoints) {

    CoilBreaker coilBreaker;
    coilBreaker.set_mirroring_dimension(_mirroringDimension);
    auto windingWindowCurrentFieldOutput = coilBreaker.breakdown_coil(magnetic, operatingPoint, 0.01);

    if (!magnetic.get_coil().get_turns_description()) {
        throw std::runtime_error("Missing turns description in coil");
    }
    auto wirePerWinding = magnetic.get_mutable_coil().get_wires();
    auto turns = magnetic.get_coil().get_turns_description().value();

    std::vector<ComplexField> complexFieldPerHarmonic;
    for (auto& fieldPerHarmonic : windingWindowCurrentFieldOutput.get_field_per_frequency()) {
        ComplexField field;
        field.set_frequency(fieldPerHarmonic.get_frequency());
        complexFieldPerHarmonic.push_back(field);
    }

    for (size_t harmonicIndex = 0; harmonicIndex < windingWindowCurrentFieldOutput.get_field_per_frequency().size(); ++harmonicIndex){
        std::vector<ComplexFieldPoint> fieldPoints;
        std::vector<FieldPoint> inducedFieldPoints;
        if (externalInducedFieldPoints) {
            inducedFieldPoints = externalInducedFieldPoints.value();
        }
        else {
            if (!magnetic.get_core().get_processed_description()) {
                throw std::runtime_error("Core is not processed");
            }
            WindingWindowElement windingWindow = magnetic.get_core().get_processed_description().value().get_winding_windows()[0]; // Hardcoded
            auto windingWindowOrigin = windingWindow.get_coordinates().value();
            double windingWindowWidth = windingWindow.get_width().value();
            double windingWindowHeight = windingWindow.get_height().value();
            for (auto& fieldPoint : windingWindowCurrentFieldOutput.get_field_per_frequency()[harmonicIndex].get_data()) {
                if (windingWindowOrigin[0] <= fieldPoint.get_point()[0] &&
                    fieldPoint.get_point()[0] <= (windingWindowOrigin[0] + windingWindowWidth) &&
                    (windingWindowOrigin[1] - windingWindowHeight / 2) <= fieldPoint.get_point()[1] &&
                    fieldPoint.get_point()[1] <= (windingWindowOrigin[1] + windingWindowHeight / 2)) {

                    inducedFieldPoints.push_back(fieldPoint);
                }
            }
        }

        std::cout << "inducedFieldPoints.size(): "  << inducedFieldPoints.size() << std::endl;

        for (auto& inducedFieldPoint : inducedFieldPoints) {
            double totalInducedFieldX = 0;
            double totalInducedFieldY = 0;

            // For the main harmonic we calculate the fringing effect for each gap
            if (_includeFringing && windingWindowCurrentFieldOutput.get_field_per_frequency()[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {
                if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
                    throw std::runtime_error("Operating point is missing magnetizing current");
                }
                if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()) {
                    throw std::runtime_error("Operating point is missing magnetizing current processed data");
                }
                auto numberTurns = magnetic.get_mutable_coil().get_number_turns(0);
                auto magnetizingCurrentPeak = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_peak().value();
                for (auto& gap : magnetic.get_core().get_functional_description().get_gapping()) {
                    // std::cout << gap.get_length() << std::endl;
                    if (gap.get_coordinates().value()[0] < 0) {
                        continue;
                    }
                    auto complexFieldPoint = _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, numberTurns, magnetizingCurrentPeak, inducedFieldPoint);

                    totalInducedFieldX += complexFieldPoint.get_real();
                    totalInducedFieldY += complexFieldPoint.get_imaginary();
                }
            }
// compobar si es espejado con las coordandas en los inducedFieldPoints
            for (auto& inducingFieldPoint : windingWindowCurrentFieldOutput.get_field_per_frequency()[harmonicIndex].get_data()) {
                std::optional<WireWrapper> inducingWire;
                if (inducingFieldPoint.get_turn_index()) {
                    auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turns[inducingFieldPoint.get_turn_index().value()].get_winding());
                    inducingWire = wirePerWinding[windingIndex];
                    if (inducedFieldPoint.get_turn_index()) {
                        if (inducedFieldPoint.get_turn_index().value() == inducingFieldPoint.get_turn_index().value()) {
                            continue;
                        }
                    }
                    else if (is_inside_inducing_turns(inducingFieldPoint, inducedFieldPoint, inducingWire.value())) {
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
            fieldPoints.push_back(complexFieldPoint);
        }
        complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
    }

    WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput;
    windingWindowMagneticStrengthFieldOutput.set_field_per_frequency(complexFieldPerHarmonic);
    return windingWindowMagneticStrengthFieldOutput;
}


std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory(MagneticFieldStrengthModels modelName) {
    if (modelName == MagneticFieldStrengthModels::BINNS_LAWRENSON) {
        return std::make_shared<MagneticFieldStrengthBinnsLawrensonModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::LAMMERANER) {
        return std::make_shared<MagneticFieldStrengthLammeranerModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::ROSHEN) {
        return std::make_shared<MagneticFieldStrengthRoshenModel>();
    }
    else
        throw std::runtime_error("Unknown Magnetic Field Strength mode, available options are: {BINNS_LAWRENSON, LAMMERANER, ROSHEN}");
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory() {
    auto defaults = Defaults();
    return factory(defaults.magneticFieldStrengthModelDefault);
}

ComplexFieldPoint MagneticFieldStrengthBinnsLawrensonModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire) {
    double Hx;
    double Hy;
    if (inducingWire.value().get_type() == WireType::ROUND || inducingWire.value().get_type() == WireType::LITZ) {
        double divisor = 2 * std::numbers::pi * (pow(inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1], 2) + pow(inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0], 2));
        Hx = inducingFieldPoint.get_value() * (inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0]) / divisor;
        Hy = -inducingFieldPoint.get_value() * (inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1]) / divisor;
    }
    else {
        auto wire = inducingWire.value();
        double a = resolve_dimensional_values(wire.get_conducting_width().value()) / 2;
        double b = resolve_dimensional_values(wire.get_conducting_height().value()) / 2;
        double x = inducedFieldPoint.get_point()[0];
        double y = inducedFieldPoint.get_point()[1];
        auto bottomLeftPoint = std::vector<double>{inducingFieldPoint.get_point()[0] - a, inducingFieldPoint.get_point()[1] - b};
        auto bottomRightPoint = std::vector<double>{inducingFieldPoint.get_point()[0] + a, inducingFieldPoint.get_point()[1] - b};
        auto topLeftPoint = std::vector<double>{inducingFieldPoint.get_point()[0] - a, inducingFieldPoint.get_point()[1] + b};
        auto topRightPoint = std::vector<double>{inducingFieldPoint.get_point()[0] + a, inducingFieldPoint.get_point()[1] + b};
        double tetha4 = atan((y - topRightPoint[1]) / (x - topRightPoint[0]));
        double tetha3 = atan((y - topLeftPoint[1]) / (x - topLeftPoint[0]));
        double tetha1 = atan((y - bottomRightPoint[1]) / (x - bottomRightPoint[0]));
        double tetha2 = atan((y - bottomLeftPoint[1]) / (x - bottomLeftPoint[0]));
        double r4 = hypot(y - topRightPoint[1], x - topRightPoint[0]);
        double r3 = hypot(y - topLeftPoint[1], x - topLeftPoint[0]);
        double r1 = hypot(y - bottomRightPoint[1], x - bottomRightPoint[0]);
        double r2 = hypot(y - bottomLeftPoint[1], x - bottomLeftPoint[0]);
        double common_part = inducingFieldPoint.get_value() / (8 * std::numbers::pi * a * b);
        Hx = -common_part * ((y + b) * (tetha1 - tetha2) - (y - b) * (tetha4 - tetha3) + (x + a) * log(r2 / r3) - (x - a) * log(r1 / r4));
        Hy = common_part * ((x + a) * (tetha2 - tetha3) - (x - a) * (tetha1 - tetha4) + (y + b) * log(r2 / r1) - (y - b) * log(r3 / r4));

        // throw std::runtime_error("Rectangular wires not implemented yet");
    }
    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index());
    complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length());
    return complexFieldPoint;   
}

ComplexFieldPoint MagneticFieldStrengthLammeranerModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<WireWrapper> inducingWire) {
    double Hx;
    double Hy;
    if (inducingWire.value().get_type() == WireType::ROUND || inducingWire.value().get_type() == WireType::LITZ) {
        if (!inducingFieldPoint.get_turn_length()) {
            throw std::runtime_error("Missing length in turn");
        }
        double distance = hypot(inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1], inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0]);
        double angle = atan2(inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0], inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1]);
        double ex = cos(angle - std::numbers::pi / 2);
        double ey = sin(angle - std::numbers::pi / 2);
        double magneticFiledStrengthModule = -inducingFieldPoint.get_value() / 2 / std::numbers::pi / distance * inducingFieldPoint.get_turn_length().value() / hypot(inducingFieldPoint.get_turn_length().value(), distance);
        Hx = magneticFiledStrengthModule * ex;
        Hy = magneticFiledStrengthModule * ey;
    }
    else {
        throw std::runtime_error("Rectangular wires not implemented yet");
    }
    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}

ComplexFieldPoint MagneticFieldStrengthRoshenModel::get_magnetic_field_strength_between_gap_and_point(CoreGap gap, int64_t numberTurns, double magnetizingCurrentPeak, FieldPoint inducedFieldPoint) {
    double distanceFromCenterEdgeGapX;
    if (gap.get_coordinates().value()[0] == 0) {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] + gap.get_section_dimensions().value()[0] / 2);
    }
    else {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] - gap.get_section_dimensions().value()[0] / 2);
    }
    double distanceFromCenterEdgeGapY = inducedFieldPoint.get_point()[1] - gap.get_coordinates().value()[1];
    double halfGapLength = gap.get_length() / 2;
    double magneticIntensityInGap = 0.9 * numberTurns * magnetizingCurrentPeak / gap.get_length();

    double magneticIntensityXDividend = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY - halfGapLength, 2);
    double magneticIntensityXDivisor = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY + halfGapLength, 2);
    double Hx = -magneticIntensityInGap / 2 / std::numbers::pi * log(magneticIntensityXDividend / magneticIntensityXDivisor);

    double m;
    if (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) > pow(halfGapLength, 2)) {
        m = 0;
    }
    else {
        m = 1;
    }

    double x = distanceFromCenterEdgeGapX * halfGapLength / (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) - pow(halfGapLength, 2));
    double Hy = -magneticIntensityInGap / std::numbers::pi * (atan(x) + m * std::numbers::pi);

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
