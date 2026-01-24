#include "physical_models/MagneticField.h"
#include "physical_models/StrayCapacitance.h"
#include "support/Painter.h"
#include "support/CoilMesher.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <chrono>
#include <thread>
#include <filesystem>
#include <set>
#include "support/Exceptions.h"


namespace OpenMagnetics {

ComplexField PainterInterface::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics()) {
            auto current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();
            if (!current.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from current");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(current.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            current.set_harmonics(harmonics);
            if (!current.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                current.set_processed(processed);
            }
            else {
                if (!current.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    current.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_current(current);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    bool mirroringDimension = settings.get_painter_mirroring_dimension();

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, true).first;

    MagneticField magneticField;
    settings.set_magnetic_field_include_fringing(includeFringing);
    settings.set_magnetic_field_mirroring_dimension(mirroringDimension);
    ComplexField field;
    {
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

    }
    auto turns = magnetic.get_coil().get_turns_description().value();

    if (turns[0].get_additional_coordinates()) {
        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
            if (turns[turnIndex].get_additional_coordinates()) {
                turns[turnIndex].set_coordinates(turns[turnIndex].get_additional_coordinates().value()[0]);
            }
        }
        magnetic.get_mutable_coil().set_turns_description(turns);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        auto additionalField = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        for (size_t pointIndex = 0; pointIndex < field.get_data().size(); ++pointIndex) {
            field.get_mutable_data()[pointIndex].set_real(field.get_mutable_data()[pointIndex].get_real() + additionalField.get_mutable_data()[pointIndex].get_real());
            field.get_mutable_data()[pointIndex].set_imaginary(field.get_mutable_data()[pointIndex].get_imaginary() + additionalField.get_mutable_data()[pointIndex].get_imaginary());
        }
    }


    return field;
}

Field PainterInterface::calculate_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    if (!operatingPoint.get_excitations_per_winding()[0].get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "voltage is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage()->get_harmonics()) {
            auto voltage = operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage().value();
            if (!voltage.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from voltage");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            voltage.set_harmonics(harmonics);
            if (!voltage.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                voltage.set_processed(processed);
            }
            else {
                if (!voltage.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    voltage.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_voltage(voltage);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    bool mirroringDimension = settings.get_painter_mirroring_dimension();

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    auto oldCoilMesherInsideTurnsFactor = settings.get_coil_mesher_inside_turns_factor();
    settings.set_coil_mesher_inside_turns_factor(1.2);
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, false, false).first;
    settings.set_coil_mesher_inside_turns_factor(oldCoilMesherInsideTurnsFactor);

    StrayCapacitance strayCapacitance;
    settings.set_magnetic_field_include_fringing(includeFringing);
    settings.set_magnetic_field_mirroring_dimension(mirroringDimension);

    auto [pixelXDimension, pixelYDimension] = Painter::get_pixel_dimensions(magnetic);

    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();


    auto capacitanceOutput = strayCapacitance.calculate_capacitance(coil);
    auto electricEnergyAmongTurns = capacitanceOutput.get_electric_energy_among_turns().value();

    std::set<std::pair<size_t, size_t>> turnsCombinations;
    for (size_t pointIndex = 0; pointIndex < inducedField.get_data().size(); ++pointIndex) {
        auto inducedFieldPoint = inducedField.get_data()[pointIndex];
        double fieldValue = 0;
        std::set<std::pair<std::string, std::string>> turnsCombinations;

        for (auto [firstTurnName, aux] : electricEnergyAmongTurns) {
            auto firstTurn = coil.get_turn_by_name(firstTurnName);
            for (auto [secondTurnName, energy] : aux) {
                auto key = std::make_pair(firstTurnName, secondTurnName);

                if (turnsCombinations.contains(key) || turnsCombinations.contains(std::make_pair(secondTurnName, firstTurnName))) {
                    continue;
                }
                turnsCombinations.insert(key);

                auto secondTurn = coil.get_turn_by_name(secondTurnName);
                double pixelArea = Painter::get_pixel_area_between_turns(firstTurn.get_coordinates(), firstTurn.get_dimensions().value(), firstTurn.get_cross_sectional_shape().value(), secondTurn.get_coordinates(), secondTurn.get_dimensions().value(), secondTurn.get_cross_sectional_shape().value(), inducedFieldPoint.get_point(), std::max(pixelXDimension, pixelYDimension));
                if (pixelArea > 0) {
                    double area = StrayCapacitance::calculate_area_between_two_turns(firstTurn, secondTurn);
                    double energyDensity = energy / area;
                    fieldValue += energyDensity * pixelArea;
                }
            }
        }
        inducedField.get_mutable_data()[pointIndex].set_value(fieldValue);
    }

    return inducedField;
}

std::shared_ptr<PainterInterface> Painter::factory(bool useAdvancedPainter, std::filesystem::path filepath, bool addProportionForColorBar, bool showTicks) {
    if (useAdvancedPainter) {
        return std::make_shared<AdvancedPainter>(filepath, addProportionForColorBar, showTicks);
    }
    else {
        return std::make_shared<BasicPainter>(filepath);
    }
}

void Painter::paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    _painter->paint_magnetic_field(operatingPoint, magnetic, harmonicIndex, inputField);
}

void Painter::paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<Field> inputField) {
    _painter->paint_electric_field(operatingPoint, magnetic, harmonicIndex, inputField);
}

void Painter::paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs, std::optional<OperatingPoint> operatingPoint, double temperature) {
    _painter->paint_wire_losses(magnetic, outputs, operatingPoint, temperature);
}

std::string Painter::export_svg() {
    return _painter->export_svg();
}

void Painter::export_png() {
    _painter->export_png();
}

void Painter::paint_core(Magnetic magnetic) {
    _painter->paint_core(magnetic);
}

void Painter::paint_bobbin(Magnetic magnetic) {
    _painter->paint_bobbin(magnetic);
}

void Painter::paint_coil_sections(Magnetic magnetic) {
    _painter->paint_coil_sections(magnetic);
}

void Painter::paint_coil_layers(Magnetic magnetic) {
    _painter->paint_coil_layers(magnetic);
}

void Painter::paint_coil_turns(Magnetic magnetic) {
    _painter->paint_coil_turns(magnetic);
}

void Painter::paint_wire(Wire wire) {
    _painter->paint_wire(wire);
}

void Painter::paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex) {
    _painter->paint_wire_with_current_density(wire, operatingPoint, windingIndex);
}

void Painter::paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature) {
    _painter->paint_wire_with_current_density(wire, current, frequency, temperature);
}

void Painter::paint_waveform(Waveform waveform) {
    paint_waveform(waveform.get_data(), waveform.get_time());
}

void Painter::paint_waveform(std::vector<double> data, std::optional<std::vector<double>> time) {
    _painter->paint_waveform(data, time);
}

void Painter::paint_curve(Curve2D curve2D, bool logScale) {
    _painter->paint_curve(curve2D, logScale);
}

void Painter::paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension) {
    _painter->paint_rectangle(xCoordinate, yCoordinate, xDimension, yDimension);
}

void Painter::paint_circle(double xCoordinate, double yCoordinate, double radius){
    _painter->paint_circle(xCoordinate, yCoordinate, radius);
}

double Painter::get_pixel_area_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    return dimension * dimension * get_pixel_proportion_between_turns(firstTurnCoordinates, firstTurnDimensions, firstTurncrossSectionalShape, secondTurnCoordinates, secondTurnDimensions, secondTurncrossSectionalShape, pixelCoordinates, dimension);
}

double Painter::get_pixel_proportion_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    // auto factor = Defaults().overlappingFactorSurroundingTurns;
    auto x1 = firstTurnCoordinates[0];
    auto y1 = firstTurnCoordinates[1];
    auto x2 = secondTurnCoordinates[0];
    auto y2 = secondTurnCoordinates[1];

    if (y2 == y1 && x2 == x1) {
        return 0;
    }


    double firstTurnMaximumDimension = 0;
    double secondTurnMaximumDimension = 0;

    if (firstTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        firstTurnMaximumDimension = hypot(firstTurnDimensions[0], firstTurnDimensions[1]);
    }
    else {
        firstTurnMaximumDimension = firstTurnDimensions[0];
    }
    if (secondTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        secondTurnMaximumDimension = hypot(secondTurnDimensions[0], secondTurnDimensions[1]);
    }
    else {
        secondTurnMaximumDimension = secondTurnDimensions[0];
    }

    double semiAverageDimensionOf12 = (firstTurnMaximumDimension + secondTurnMaximumDimension) / 4;

    auto x0 = pixelCoordinates[0];
    auto y0 = pixelCoordinates[1];

    auto distanceFrom0toLine12 = fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(pow(y2 - y1, 2) + pow(x2 - x1, 2));
    auto distanceFrom0toCenter1 = hypot(x1 - x0, y1 - y0);
    auto distanceFrom0toCenter2 = hypot(x2 - x0, y2 - y0);
    auto distanceFromCenter1toCenter2 = hypot(x2 - x1, y2 - y1);

    if (distanceFrom0toCenter1 > distanceFromCenter1toCenter2 || distanceFrom0toCenter2 > distanceFromCenter1toCenter2) {
        return 0;   
    }

    double proportion;
    if (distanceFrom0toLine12 - dimension / 2 > semiAverageDimensionOf12) {
        proportion = 0;
    }
    else if (distanceFrom0toLine12 + dimension / 2 < semiAverageDimensionOf12) {
        proportion = 1;
    }
    else {
        proportion = (semiAverageDimensionOf12 - (distanceFrom0toLine12 - dimension / 2)) / dimension;
    }

    proportion *= (semiAverageDimensionOf12 - distanceFrom0toLine12) / semiAverageDimensionOf12;

    return proportion;
}

std::pair<double, double> Painter::get_pixel_dimensions(Magnetic magnetic) {

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();
    double coreWindingWindowWidth = magnetic.get_mutable_core().get_winding_window().get_width().value();

    double pixelXDimension = coreWindingWindowWidth / numberPointsX;
    double pixelYDimension = coreColumnHeight / numberPointsY;

    return {pixelXDimension, pixelYDimension};
}


} // namespace OpenMagnetics
