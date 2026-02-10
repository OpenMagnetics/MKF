#include "physical_models/WindingLosses.h"
#include "support/Painter.h"
#include <cfloat>
#include <map>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "support/Exceptions.h"

namespace OpenMagnetics {

std::vector<SVG::Point> scale_points(std::vector<SVG::Point> points, double imageHeight, double scale = 1) {
    auto constants = Constants();
    std::vector<SVG::Point> scaledPoints;
    for (auto& point : points) {
        scaledPoints.push_back(SVG::Point(point.first * scale, (imageHeight / 2 - point.second) * scale));
    }
    return scaledPoints;
}

std::vector<double> BasicPainter::get_image_size(Magnetic magnetic) {
    auto core = magnetic.get_core();

    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    
    double showingCoreWidth;
    auto family = core.get_shape_family();
    switch (family) {
        case CoreShapeFamily::C:
        case CoreShapeFamily::U:
        case CoreShapeFamily::UR:
            _extraDimension = 1;
            showingCoreWidth = (processedDescription.get_width() - mainColumn.get_width() / 2);
            break;
        case CoreShapeFamily::T:
            _extraDimension = Coil::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());
            showingCoreWidth = processedDescription.get_width() * _extraDimension;
            break;
        default:
            _extraDimension = 1;
            showingCoreWidth = processedDescription.get_width() / 2;
            break;
    }

    double showingCoreHeight = processedDescription.get_height() * _extraDimension;
    return {showingCoreWidth, showingCoreHeight};
}


void BasicPainter::paint_round_wire(double xCoordinate, double yCoordinate, Wire wire, std::optional<std::string> label) {
    if (!wire.get_outer_diameter()) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing outerDiameter");
    }

    double outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
    double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
    double insulationThickness = (outerDiameter - conductingDiameter) / 2;
    auto coating = wire.resolve_coating();
    size_t numberLines = 0;
    double strokeWidth = 0;
    double lineRadiusIncrease = 0;
    double currentLineDiameter = conductingDiameter;
    auto coatingColor = settings.get_painter_color_insulation();
    if (coating) {
        CoatingInfo coatingInfo = process_coating(insulationThickness, coating.value());
        strokeWidth = coatingInfo.strokeWidth;
        numberLines = coatingInfo.numberLines;
        lineRadiusIncrease = coatingInfo.lineRadiusIncrease;
        coatingColor = coatingInfo.coatingColor;
    }
    coatingColor = std::regex_replace(std::string(coatingColor), std::regex("0x"), "#");

    SVG::Group* shapes = _root.add_child<SVG::Group>();

    double opacity = 1;
    // if (_fieldPainted) {
    //     opacity = 0.25;
    // }
    // Paint insulation
    {
        std::string cssClassName = generate_random_string();

        _root.style("." + cssClassName).set_attr("fill", coatingColor).set_attr("opacity", opacity);
        paint_circle(xCoordinate, yCoordinate, outerDiameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
    }

    // Paint copper
    {
        if (!wire.get_conducting_diameter()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing conducting diameter");
        }
        std::string colorClass;
        // if (_fieldPainted) {
        //     colorClass = "copper_translucent";
        // }
        // else {
            colorClass = "copper";
        // }
        paint_circle(xCoordinate, yCoordinate, conductingDiameter / 2, colorClass, shapes, 360, 0, {0, 0}, label);
    }

    // Paint layer separation lines
    {
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("opacity", opacity).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_lines()), std::regex("0x"), "#"));
        
        for (size_t i = 0; i < numberLines; ++i) {
            paint_circle(xCoordinate, yCoordinate, currentLineDiameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
            currentLineDiameter += lineRadiusIncrease;
        }
    }
}

void BasicPainter::paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire, std::optional<std::string> label) {
    if (!wire.get_outer_diameter()) {
        wire.set_nominal_value_outer_diameter(wire.calculate_outer_diameter());
    }

    bool simpleMode = settings.get_painter_simple_litz();
    auto coating = wire.resolve_coating();
    size_t numberConductors = wire.get_number_conductors().value();

    double outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
    auto strand = wire.resolve_strand();
    double strandOuterDiameter = resolve_dimensional_values(strand.get_outer_diameter().value());
    double conductingDiameter = 0;


    if (coating->get_type() == InsulationWireCoatingType::BARE) {
        conductingDiameter = outerDiameter;
        auto strandCoating = Wire::resolve_coating(strand);
        double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        auto conductingDiameterTheory = Wire::get_outer_diameter_bare_litz(strandConductingDiameter, numberConductors, strandCoating->get_grade().value());
        if (conductingDiameterTheory > conductingDiameter) {
            conductingDiameter = conductingDiameterTheory;
            outerDiameter = conductingDiameterTheory;
            wire.set_nominal_value_outer_diameter(outerDiameter);
            // set_image_size(wire);
        }

    }
    else if (coating->get_type() == InsulationWireCoatingType::SERVED) {
        if (!coating->get_number_layers()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Number layers missing in litz served");
        }
        auto strandCoating = Wire::resolve_coating(strand);
        double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        conductingDiameter = Wire::get_outer_diameter_bare_litz(strandConductingDiameter, numberConductors, strandCoating->get_grade().value());
        if (outerDiameter <= conductingDiameter) {
            double servedThickness = Wire::get_serving_thickness_from_standard(coating->get_number_layers().value(), outerDiameter);
            outerDiameter = conductingDiameter + 2 * servedThickness;
            wire.set_nominal_value_outer_diameter(outerDiameter);
            // set_image_size(wire);
        }
    }
    else if (coating->get_type() == InsulationWireCoatingType::INSULATED) {
        auto strandCoating = Wire::resolve_coating(strand);
        double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
        conductingDiameter = Wire::get_outer_diameter_bare_litz(strandConductingDiameter, numberConductors, strandCoating->get_grade().value());
        if (outerDiameter <= conductingDiameter) {
            double insulationThickness = coating->get_number_layers().value() * coating->get_thickness_layers().value();
            outerDiameter = conductingDiameter + 2 * insulationThickness;
            wire.set_nominal_value_outer_diameter(outerDiameter);
            // set_image_size(wire);
        }
    }
    else {
        throw NotImplementedException("Coating type not implemented for Litz yet");
    }

    double insulationThickness = (outerDiameter - conductingDiameter) / 2;
    if (insulationThickness < 0) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Insulation thickness cannot be negative");
    }

    size_t numberLines = 0;
    double strokeWidth = 0;
    double lineRadiusIncrease = 0;
    double currentLineDiameter = conductingDiameter;
    auto coatingColor = settings.get_painter_color_insulation();
    if (coating) {
        CoatingInfo coatingInfo = process_coating(insulationThickness, coating.value());
        strokeWidth = coatingInfo.strokeWidth;
        numberLines = coatingInfo.numberLines;
        lineRadiusIncrease = coatingInfo.lineRadiusIncrease;
        coatingColor = coatingInfo.coatingColor;
    }


    coatingColor = std::regex_replace(std::string(coatingColor), std::regex("0x"), "#");

    auto shapes = _root.add_child<SVG::Group>();

    // Paint insulation
    {
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", coatingColor);
        paint_circle(xCoordinate, yCoordinate, outerDiameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
    }
    // Paint layer separation lines
    {
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_lines()), std::regex("0x"), "#"));
        
        for (size_t i = 0; i < numberLines; ++i) {
            paint_circle(xCoordinate, yCoordinate, currentLineDiameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
            currentLineDiameter += lineRadiusIncrease;
        }
    }

    if (simpleMode) {
        paint_circle(xCoordinate, yCoordinate, conductingDiameter / 2, "copper", shapes, 360, 0, {0, 0}, label);
    }
    else {
        // Contour
        {
            paint_circle(xCoordinate, yCoordinate, conductingDiameter / 2, "white", shapes, 360, 0, {0, 0}, label);
        }

        auto coordinateFilePath = settings.get_painter_cci_coordinates_path();
        coordinateFilePath = coordinateFilePath.append("cci" + std::to_string(numberConductors) + ".txt");

        std::ifstream in(coordinateFilePath);
        std::vector<std::pair<double, double>> coordinates;
        bool advancedMode = settings.get_painter_advanced_litz();

        if (in) {
            std::string line;

            while (getline(in, line)) {
                std::stringstream sep(line);
                std::string field;

                std::vector<double> numbers;

                while (getline(sep, field, ' ')) {
                    try {
                        numbers.push_back(stod(field));
                    }
                    catch (...) {
                        continue;
                    }
                }

                coordinates.push_back({numbers[1], numbers[2]});
            }

            for (size_t i = 0; i < numberConductors; ++i) {
                double internalXCoordinate = conductingDiameter / 2 * coordinates[i].first;
                double internalYCoordinate = conductingDiameter / 2 * coordinates[i].second;

                if (advancedMode) {
                    BasicPainter::paint_round_wire(xCoordinate + internalXCoordinate, -(yCoordinate + internalYCoordinate), strand);
                }
                else {
                    paint_circle(xCoordinate + internalXCoordinate, yCoordinate - internalYCoordinate, strandOuterDiameter / 2, "copper", shapes, 360, 0, {0, 0}, label);
                }
            }
        }
        else {
            double currentRadius = 0;
            double currentAngle = 0;
            double angleCoveredThisLayer = 0;
            double strandAngle = 360;
            for (size_t i = 0; i < numberConductors; ++i) {
                double internalXCoordinate = currentRadius * cos(currentAngle / 180 * std::numbers::pi);
                double internalYCoordinate = currentRadius * sin(currentAngle / 180 * std::numbers::pi);

                if (advancedMode) {
                    BasicPainter::paint_round_wire(xCoordinate + internalXCoordinate, -(yCoordinate + internalYCoordinate), strand);
                }
                else {
                    paint_circle(xCoordinate + internalXCoordinate, yCoordinate - internalYCoordinate, strandOuterDiameter / 2, "copper", shapes, 360, 0, {0, 0}, label);
                }

                if (currentRadius > 0) {
                    strandAngle = wound_distance_to_angle(strandOuterDiameter, currentRadius);
                }

                if (angleCoveredThisLayer + strandAngle * 1.99 > 360) {
                    currentRadius += strandOuterDiameter;
                    if (currentRadius + strandOuterDiameter / 2 > conductingDiameter / 2) {
                        // We cut down some strands to avoid visual error, which should be done only at thousands of strands due to cci_coords files
                        break;
                    }
                    angleCoveredThisLayer = 0;
                }
                else {
                    currentAngle += strandAngle;
                    angleCoveredThisLayer += strandAngle;
                }
            }
        }

    }
}

void BasicPainter::paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension) {
    return paint_rectangle(xCoordinate, yCoordinate, xDimension, yDimension, "point");
}

void BasicPainter::paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string cssClassName, SVG::Group* group, double angle, std::vector<double> center, std::optional<std::string> label) {
    std::vector<SVG::Point> turnPoints = {};
    turnPoints.push_back(SVG::Point(xCoordinate - xDimension / 2, yCoordinate + yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate + xDimension / 2, yCoordinate + yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate + xDimension / 2, yCoordinate - yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate - xDimension / 2, yCoordinate - yDimension / 2));
    if (group == nullptr) {
        group = _root.add_child<SVG::Group>();
    }
    *group << SVG::Polygon(scale_points(turnPoints, 0, _scale));
    auto turnSvg = _root.get_children<SVG::Polygon>().back();
    turnSvg->set_attr("class", cssClassName);
    turnSvg->set_attr("transform", "rotate( " + std::to_string(-(angle)) + " " + std::to_string(center[0] * _scale) + " " + std::to_string(center[1] * _scale) + ") ");
    if (label) {
        turnSvg->add_child<SVG::Title>(label.value());
    }
}

void BasicPainter::paint_circle(double xCoordinate, double yCoordinate, double radius) {
    return paint_circle(xCoordinate, yCoordinate, radius, "point");
}

void BasicPainter::paint_circle(double xCoordinate, double yCoordinate, double radius, std::string cssClassName, SVG::Group* group, double fillAngle, double angle, std::vector<double> center, std::optional<std::string> label) {
    if (group == nullptr) {
        group = _root.add_child<SVG::Group>();
    }
    *group << SVG::Circle(xCoordinate * _scale, -yCoordinate * _scale, radius * _scale);
    auto turnSvg = _root.get_children<SVG::Circle>().back();
    turnSvg = _root.get_children<SVG::Circle>().back();
    turnSvg->set_attr("class", cssClassName);

    // auto group = _root.add_child<SVG::Group>();
    // *group << SVG::Circle(xCoordinate * _scale, -yCoordinate * _scale, radius * _scale);
    if (angle != 0) {
        turnSvg->set_attr("transform", "rotate( " + std::to_string(angle) + " " + std::to_string(center[0]) + " " + std::to_string(center[1]) + ")");
    }

    if (fillAngle < 360) {
        double circlePerimeter = std::numbers::pi * 2 * radius * _scale;
        double angleProportion = fillAngle / 360;
        std::string termination = angleProportion < 1? "butt" : "round";
        group->set_attr("stroke-linecap", termination);
        group->set_attr("stroke-dashoffset", "0");
        group->set_attr("stroke-dasharray", std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)));
    }
    if (label) {
        turnSvg->add_child<SVG::Title>(label.value());
    }
}

void BasicPainter::paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle, std::vector<double> center, std::optional<std::string> label) {
    double outerWidth = 0;
    double outerHeight = 0;
    if (wire.get_outer_width()) {
        outerWidth = resolve_dimensional_values(wire.get_outer_width().value());
    }
    else {
        if (!wire.get_conducting_width()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing both outerWidth and conductingWidth");
        }
        outerWidth = resolve_dimensional_values(wire.get_conducting_width().value());
    }
    if (wire.get_outer_height()) {
        outerHeight = resolve_dimensional_values(wire.get_outer_height().value());
    }
    else {
        if (!wire.get_conducting_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing both outerHeight and conductingHeight");
        }
        outerHeight = resolve_dimensional_values(wire.get_conducting_height().value());
    }
    double conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
    double conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
    double insulationThicknessInWidth = (outerWidth - conductingWidth) / 2;
    double insulationThicknessInHeight = (outerHeight - conductingHeight) / 2;
    auto coating = wire.resolve_coating();
    size_t numberLines = 0;
    double strokeWidth = 0;
    double lineWidthIncrease = 0;
    double lineHeightIncrease = 0;
    double currentLineWidth = conductingWidth;
    double currentLineHeight = conductingHeight;

    std::string coatingColor = settings.get_painter_color_insulation();
    if (coating) {
        InsulationWireCoatingType insulationWireCoatingType = coating->get_type().value();

        switch (insulationWireCoatingType) {
            case InsulationWireCoatingType::BARE:
                break;
            case InsulationWireCoatingType::ENAMELLED:
                if (!coating->get_grade()) {
                    throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Enamelled wire missing grade");
                }
                numberLines = coating->get_grade().value() + 1;
                lineWidthIncrease = insulationThicknessInWidth / coating->get_grade().value() * 2;
                lineHeightIncrease = insulationThicknessInHeight / coating->get_grade().value() * 2;
                coatingColor = settings.get_painter_color_enamel();
                break;
            default:
                json insulationWireCoatingTypeJson;
                to_json(insulationWireCoatingTypeJson, insulationWireCoatingType);
                throw NotImplementedException("Coating type plot not implemented yet: " + to_string(insulationWireCoatingTypeJson));
        }
        strokeWidth = std::min(lineWidthIncrease / 10 / numberLines, lineHeightIncrease / 10 / numberLines);
    }
    coatingColor = std::regex_replace(std::string(coatingColor), std::regex("0x"), "#");

    SVG::Group* shapes = _root.add_child<SVG::Group>();
    // Paint insulation
    {
        std::string cssClassName = generate_random_string();

        _root.style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", coatingColor);
        paint_rectangle(xCoordinate, yCoordinate, outerWidth, outerHeight, cssClassName, shapes, angle, center);
    }

    // Paint copper
    {                
        paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, "copper", shapes, angle, center);
    }

    // Paint layer separation lines

    {
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_lines()), std::regex("0x"), "#"));
        for (size_t i = 0; i < numberLines; ++i) {
            paint_rectangle(xCoordinate, yCoordinate, currentLineWidth, currentLineHeight, cssClassName, shapes, angle, center);
            currentLineWidth += lineWidthIncrease;
            currentLineHeight += lineHeightIncrease;
        }
    }
}

void BasicPainter::paint_two_piece_set_coil_sections(Magnetic magnetic) {
    auto constants = Constants();

    if (!magnetic.get_coil().get_sections_description()) {
        throw CoilNotProcessedException("Winding sections not created");
    }

    auto sections = magnetic.get_coil().get_sections_description().value();

    auto shapes = _root.add_child<SVG::Group>();
    for (size_t i = 0; i < sections.size(); ++i){
        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
            _root.style(".section_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleSections[i % constants.coilPainterColorsScaleSections.size()]);
            paint_rectangle(sections[i].get_coordinates()[0], sections[i].get_coordinates()[1], sections[i].get_dimensions()[0], sections[i].get_dimensions()[1], "section_" + std::to_string(i), shapes);
        }
        else {
            paint_rectangle(sections[i].get_coordinates()[0], sections[i].get_coordinates()[1], sections[i].get_dimensions()[0], sections[i].get_dimensions()[1], "insulation", shapes);
        }
    }
    paint_two_piece_set_margin(magnetic);
}

void BasicPainter::paint_toroidal_coil_sections(Magnetic magnetic) {

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    if (!magnetic.get_coil().get_sections_description()) {
        throw CoilNotProcessedException("Winding sections not created");
    }

    double initialRadius = processedDescription.get_width() / 2 - mainColumn.get_width();

    auto sections = magnetic.get_coil().get_sections_description().value();

    for (size_t i = 0; i < sections.size(); ++i){
        {
            double strokeWidth = sections[i].get_dimensions()[0];
            double circleDiameter = (initialRadius - sections[i].get_coordinates()[0]) * 2;
            double angleProportion = sections[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";

            std::string cssClassName = generate_random_string();
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_copper()), std::regex("0x"), "#"));
            }
            else {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_insulation()), std::regex("0x"), "#"));
            }

            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, sections[i].get_dimensions()[1], -(sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2), {0, 0});
        }
    }

    paint_toroidal_margin(magnetic);
}

void BasicPainter::paint_two_piece_set_coil_layers(Magnetic magnetic) {
    auto constants = Constants();
    Coil coil = magnetic.get_coil();
    if (!coil.get_layers_description()) {
        throw CoilNotProcessedException("Winding layers not created");
    }

    auto layers = coil.get_layers_description().value();

    auto shapes = _root.add_child<SVG::Group>();
    for (size_t i = 0; i < layers.size(); ++i){
        if (layers[i].get_type() == ElectricalType::CONDUCTION) {
            _root.style(".layer_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleLayers[i % constants.coilPainterColorsScaleLayers.size()]);
            paint_rectangle(layers[i].get_coordinates()[0], layers[i].get_coordinates()[1], layers[i].get_dimensions()[0], layers[i].get_dimensions()[1], "layer_" + std::to_string(i), shapes);
        }
        else {
            paint_rectangle(layers[i].get_coordinates()[0], layers[i].get_coordinates()[1], layers[i].get_dimensions()[0], layers[i].get_dimensions()[1], "insulation", shapes);
        }
    }
    paint_two_piece_set_margin(magnetic);
}

void BasicPainter::paint_toroidal_coil_layers(Magnetic magnetic) {
    Coil winding = magnetic.get_coil();
    Core core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core has not been processed");
    }

    if (!winding.get_layers_description()) {
        throw CoilNotProcessedException("Winding layers not created");
    }

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    auto layers = winding.get_layers_description().value();
    double initialRadius = processedDescription.get_width() / 2 - mainColumn.get_width();

    for (size_t i = 0; i < layers.size(); ++i){
        {
            double strokeWidth = layers[i].get_dimensions()[0];
            double circleDiameter = (initialRadius - layers[i].get_coordinates()[0]) * 2;
            double angleProportion = layers[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";

            std::string cssClassName = generate_random_string();
            if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_copper()), std::regex("0x"), "#"));
            }
            else {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_insulation()), std::regex("0x"), "#"));
            }
            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, layers[i].get_dimensions()[1], -(layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2), {0, 0});
        }
    }

    paint_toroidal_margin(magnetic);
}

void BasicPainter::paint_two_piece_set_coil_turns(Magnetic magnetic) {
    auto constants = Constants();
    Coil coil = magnetic.get_coil();
    auto wirePerWinding = coil.get_wires();

    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding turns not created");
    }

    auto turns = coil.get_turns_description().value();

    auto shapes = _root.add_child<SVG::Group>();

    WiringTechnology coilType = WiringTechnology::WOUND;

    if (coil.get_groups_description()) {
        coilType = coil.get_groups_description().value()[0].get_type(); // TODO: take into account more groups
    }

    auto layers = coil.get_layers_description().value();

    if (coilType == WiringTechnology::WOUND) {
        for (size_t i = 0; i < layers.size(); ++i){
            if (layers[i].get_type() == ElectricalType::INSULATION) {
                paint_rectangle(layers[i].get_coordinates()[0], layers[i].get_coordinates()[1], layers[i].get_dimensions()[0], layers[i].get_dimensions()[1], "insulation", shapes);
            }
        }
        paint_two_piece_set_margin(magnetic);
    }
    else if (coilType == WiringTechnology::PRINTED){
        std::string styleClass = "fr4";
        if (!_fieldPainted) {
            // styleClass = "fr4_translucent";
            auto group = coil.get_groups_description().value()[0]; // TODO: take into account more groups
            paint_rectangle(group.get_coordinates()[0], group.get_coordinates()[1], group.get_dimensions()[0], group.get_dimensions()[1], styleClass, shapes);
        }
    }

    for (size_t i = 0; i < turns.size(); ++i){

        auto windingIndex = coil.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND) {
            paint_round_wire(turns[i].get_coordinates()[0], turns[i].get_coordinates()[1], wirePerWinding[windingIndex], turns[i].get_name());
        }
        else if (wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            paint_litz_wire(turns[i].get_coordinates()[0], turns[i].get_coordinates()[1], wirePerWinding[windingIndex], turns[i].get_name());
        }
        else {
            {
                _root.style(".turn_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleTurns[turns[i].get_parallel() % constants.coilPainterColorsScaleTurns.size()]);
                double xCoordinate = turns[i].get_coordinates()[0];
                double yCoordinate = turns[i].get_coordinates()[1];
                double outerWidth = 0;
                double outerHeight = 0;
                if (wire.get_outer_width()) {
                    outerWidth = resolve_dimensional_values(wire.get_outer_width().value());
                }
                else {
                    if (!wire.get_conducting_width()) {
                        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing both outerWidth and conductingWidth");
                    }
                    outerWidth = resolve_dimensional_values(wire.get_conducting_width().value());
                }
                if (wire.get_outer_height()) {
                    outerHeight = resolve_dimensional_values(wire.get_outer_height().value());
                }
                else {
                    if (!wire.get_conducting_height()) {
                        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire is missing both outerHeight and conductingHeight");
                    }
                    outerHeight = resolve_dimensional_values(wire.get_conducting_height().value());
                }
                paint_rectangle(xCoordinate, yCoordinate, outerWidth, outerHeight, "turn_" + std::to_string(i), shapes, 0, {0, 0}, turns[i].get_name());
            }

            if (wire.get_conducting_width() && wire.get_conducting_height()) {
                double xCoordinate = turns[i].get_coordinates()[0];
                double yCoordinate = turns[i].get_coordinates()[1];
                double conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
                double conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
                paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, "copper", shapes, 0, {0, 0}, turns[i].get_name());
            }
        }
    }
}

void BasicPainter::paint_toroidal_coil_turns(Magnetic magnetic) {
    Coil winding = magnetic.get_coil();
    auto wirePerWinding = winding.get_wires();

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    double initialRadius = processedDescription.get_width() / 2 - mainColumn.get_width();

    if (!winding.get_turns_description()) {
        throw CoilNotProcessedException("Winding turns not created");
    }

    auto turns = winding.get_turns_description().value();

    for (size_t i = 0; i < turns.size(); ++i){

        if (!turns[i].get_coordinate_system()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Turn is missing coordinate system");
        }
        if (!turns[i].get_rotation()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Turn is missing rotation");
        }
        if (turns[i].get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Painter: Turn coordinates are not in cartesian");
        }

        auto windingIndex = winding.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        double xCoordinate = turns[i].get_coordinates()[0];
        double yCoordinate = turns[i].get_coordinates()[1];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND) {
            paint_round_wire(xCoordinate, yCoordinate, wire);
        }
        else if ( wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            paint_litz_wire(xCoordinate, yCoordinate, wire);
        }
        else {
            double turnAngle = turns[i].get_rotation().value();
            std::vector<double> turnCenter = {(xCoordinate), (-yCoordinate)};
            paint_rectangular_wire(xCoordinate, yCoordinate, wire, turnAngle, turnCenter);
        }

        if (turns[i].get_additional_coordinates()) {
            auto additionalCoordinates = turns[i].get_additional_coordinates().value();

            for (auto additionalCoordinate : additionalCoordinates){
                double xAdditionalCoordinate = additionalCoordinate[0];
                double yAdditionalCoordinate = additionalCoordinate[1];
                if (wirePerWinding[windingIndex].get_type() == WireType::ROUND) {
                    paint_round_wire(xAdditionalCoordinate, yAdditionalCoordinate, wire);
                }
                else if (wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                    paint_litz_wire(xAdditionalCoordinate, yAdditionalCoordinate, wire);
                }
                else {
                    double turnAngle = turns[i].get_rotation().value();
                    std::vector<double> turnCenter = {(xAdditionalCoordinate), (-yAdditionalCoordinate)};
                    paint_rectangular_wire(xAdditionalCoordinate, yAdditionalCoordinate, wire, turnAngle, turnCenter);
                }
            }
        }
    }

    auto layers = winding.get_layers_description().value();

    for (size_t i = 0; i < layers.size(); ++i){
        if (layers[i].get_type() == ElectricalType::INSULATION) {

            double strokeWidth = layers[i].get_dimensions()[0];
            double circleDiameter = (initialRadius - layers[i].get_coordinates()[0]) * 2;
            double angleProportion = layers[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";

            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_insulation()), std::regex("0x"), "#"));
            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, layers[i].get_dimensions()[1], -(layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2), {0, 0});

            if (layers[i].get_additional_coordinates()) {
                circleDiameter = (initialRadius - layers[i].get_additional_coordinates().value()[0][0]) * 2;
                paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, layers[i].get_dimensions()[1], -(layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2), {0, 0});
            }
        }
    }

    paint_toroidal_margin(magnetic);
    // _root.autoscale();
}

void BasicPainter::paint_two_piece_set_bobbin(Magnetic magnetic) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw CoilNotProcessedException("Bobbin has not been processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();

    std::vector<double> bobbinCoordinates = std::vector<double>({0, 0, 0});
    if (bobbinProcessedDescription.get_coordinates()) {
        bobbinCoordinates = bobbinProcessedDescription.get_coordinates().value();
    }

    auto shapes = _root.add_child<SVG::Group>();
    double bobbinOuterWidth = bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() + bobbinProcessedDescription.get_winding_windows()[0].get_width().value();
    double bobbinOuterHeight = bobbinProcessedDescription.get_wall_thickness();
    for (auto& windingWindow: bobbinProcessedDescription.get_winding_windows()) {
        bobbinOuterHeight += windingWindow.get_height().value();
        bobbinOuterHeight += bobbinProcessedDescription.get_wall_thickness();
    }

    std::vector<SVG::Point> bobbinPoints = {};
    bobbinPoints.push_back(SVG::Point(bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() - bobbinProcessedDescription.get_column_thickness(),
                                      bobbinCoordinates[1] + bobbinOuterHeight / 2));
    bobbinPoints.push_back(SVG::Point(bobbinOuterWidth,
                                      bobbinCoordinates[1] + bobbinOuterHeight / 2));
    bobbinPoints.push_back(SVG::Point(bobbinOuterWidth,
                                      bobbinCoordinates[1] + bobbinOuterHeight / 2 - bobbinProcessedDescription.get_wall_thickness()));
    bobbinPoints.push_back(SVG::Point(bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value(),
                                      bobbinCoordinates[1] + bobbinOuterHeight / 2 - bobbinProcessedDescription.get_wall_thickness()));
    bobbinPoints.push_back(SVG::Point(bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value(),
                                      bobbinCoordinates[1] - bobbinOuterHeight / 2 + bobbinProcessedDescription.get_wall_thickness()));
    bobbinPoints.push_back(SVG::Point(bobbinOuterWidth,
                                      bobbinCoordinates[1] - bobbinOuterHeight / 2 + bobbinProcessedDescription.get_wall_thickness()));
    bobbinPoints.push_back(SVG::Point(bobbinOuterWidth,
                                      bobbinCoordinates[1] - bobbinOuterHeight / 2));
    bobbinPoints.push_back(SVG::Point(bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() - bobbinProcessedDescription.get_column_thickness(),
                                      bobbinCoordinates[1] - bobbinOuterHeight / 2));

    *shapes << SVG::Polygon(scale_points(bobbinPoints, 0, _scale));

    auto sectionSvg = _root.get_children<SVG::Polygon>().back();
    // sectionSvg->set_attr("class", "bobbin");
    if (_fieldPainted) {
        sectionSvg->set_attr("class", "bobbin_translucent");
    }
    else {
        sectionSvg->set_attr("class", "bobbin");
    }
}

void BasicPainter::paint_two_piece_set_core(Core core) {
    std::vector<SVG::Point> topPiecePoints = {};
    std::vector<SVG::Point> bottomPiecePoints = {};
    std::vector<std::vector<SVG::Point>> gapChunks = {};
    CoreShape shape = core.resolve_shape();
    auto processedDescription = core.get_processed_description().value();
    auto rightColumn = core.find_closest_column_by_coordinates(std::vector<double>({processedDescription.get_width() / 2, 0, -processedDescription.get_depth() / 2}));
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto family = core.get_shape_family();
    double showingCoreWidth;
    double showingMainColumnWidth;
    switch (family) {
        case MAS::CoreShapeFamily::C:
        case MAS::CoreShapeFamily::U:
        case MAS::CoreShapeFamily::UR:
            showingMainColumnWidth = mainColumn.get_width() / 2;
            showingCoreWidth = processedDescription.get_width() - mainColumn.get_width() / 2;
            break;
        default:
            showingCoreWidth = processedDescription.get_width() / 2;
            showingMainColumnWidth = mainColumn.get_width() / 2;
    }

    double rightColumnWidth;
    if (rightColumn.get_minimum_width()) {
        rightColumnWidth = rightColumn.get_minimum_width().value();
    }
    else {
        rightColumnWidth = rightColumn.get_width();
    }

    auto gapsInMainColumn = core.find_gaps_by_column(mainColumn);
    std::sort(gapsInMainColumn.begin(), gapsInMainColumn.end(), [](const CoreGap& lhs, const CoreGap& rhs) { return lhs.get_coordinates().value()[1] > rhs.get_coordinates().value()[1];});

    auto gapsInRightColumn = core.find_gaps_by_column(rightColumn);
    std::sort(gapsInRightColumn.begin(), gapsInRightColumn.end(), [](const CoreGap& lhs, const CoreGap& rhs) { return lhs.get_coordinates().value()[1] > rhs.get_coordinates().value()[1];});

    double lowestHeightTopCoreMainColumn = 0;
    double lowestHeightTopCoreRightColumn = 0;
    double highestHeightBottomCoreMainColumn = 0;
    double highestHeightBottomCoreRightColumn = 0;
    double topCoreOffset = 0;
    double bottomCoreOffset = 0;
    if (gapsInMainColumn.size() == 0) {
        lowestHeightTopCoreMainColumn = 0;
        highestHeightBottomCoreMainColumn = 0;
    }
    else {
        if (gapsInRightColumn.front().get_type() != GapType::ADDITIVE) {
            lowestHeightTopCoreMainColumn = gapsInMainColumn.front().get_coordinates().value()[1] + gapsInMainColumn.front().get_length() / 2;
            highestHeightBottomCoreMainColumn = gapsInMainColumn.back().get_coordinates().value()[1] - gapsInMainColumn.back().get_length() / 2;
        }
        else {
            topCoreOffset = gapsInMainColumn.front().get_length() / 2;
            bottomCoreOffset = -gapsInMainColumn.front().get_length() / 2;
        }
    }
    if (gapsInRightColumn.size() == 0) {
        lowestHeightTopCoreRightColumn = 0;
        highestHeightBottomCoreRightColumn = 0;
    }
    else {
        if (gapsInRightColumn.front().get_type() != GapType::ADDITIVE) {
            lowestHeightTopCoreRightColumn = gapsInRightColumn.front().get_coordinates().value()[1] + gapsInRightColumn.front().get_length() / 2;
            highestHeightBottomCoreRightColumn = gapsInRightColumn.back().get_coordinates().value()[1] - gapsInRightColumn.back().get_length() / 2;
        }
        else {
        }
    }

    topPiecePoints.push_back(SVG::Point(0, topCoreOffset + processedDescription.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth, topCoreOffset + processedDescription.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth, topCoreOffset + lowestHeightTopCoreRightColumn));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, topCoreOffset + lowestHeightTopCoreRightColumn));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, topCoreOffset + rightColumn.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingMainColumnWidth, topCoreOffset + mainColumn.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingMainColumnWidth, topCoreOffset + lowestHeightTopCoreMainColumn));
    topPiecePoints.push_back(SVG::Point(0, topCoreOffset + lowestHeightTopCoreMainColumn));

    for (size_t i = 1; i < gapsInMainColumn.size(); ++i)
    {
        std::vector<SVG::Point> chunk;
        chunk.push_back(SVG::Point(0, gapsInMainColumn[i - 1].get_coordinates().value()[1] - gapsInMainColumn[i - 1].get_length() / 2));
        chunk.push_back(SVG::Point(showingMainColumnWidth, gapsInMainColumn[i - 1].get_coordinates().value()[1] - gapsInMainColumn[i - 1].get_length() / 2));
        chunk.push_back(SVG::Point(showingMainColumnWidth, gapsInMainColumn[i].get_coordinates().value()[1] + gapsInMainColumn[i].get_length() / 2));
        chunk.push_back(SVG::Point(0, gapsInMainColumn[i].get_coordinates().value()[1] + gapsInMainColumn[i].get_length() / 2));
        gapChunks.push_back(chunk);
    }
    for (size_t i = 1; i < gapsInRightColumn.size(); ++i)
    {
        std::vector<SVG::Point> chunk;
        chunk.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, gapsInRightColumn[i - 1].get_coordinates().value()[1] - gapsInRightColumn[i - 1].get_length() / 2));
        chunk.push_back(SVG::Point(showingCoreWidth, gapsInRightColumn[i - 1].get_coordinates().value()[1] - gapsInRightColumn[i - 1].get_length() / 2));
        chunk.push_back(SVG::Point(showingCoreWidth, gapsInRightColumn[i].get_coordinates().value()[1] + gapsInRightColumn[i].get_length() / 2));
        chunk.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, gapsInRightColumn[i].get_coordinates().value()[1] + gapsInRightColumn[i].get_length() / 2));
        gapChunks.push_back(chunk);
    }

    bottomPiecePoints.push_back(SVG::Point(0, bottomCoreOffset - processedDescription.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth, bottomCoreOffset - processedDescription.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth, bottomCoreOffset + highestHeightBottomCoreRightColumn));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, bottomCoreOffset + highestHeightBottomCoreRightColumn));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, bottomCoreOffset - rightColumn.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingMainColumnWidth, bottomCoreOffset - mainColumn.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingMainColumnWidth, bottomCoreOffset + highestHeightBottomCoreMainColumn));
    bottomPiecePoints.push_back(SVG::Point(0, bottomCoreOffset + highestHeightBottomCoreMainColumn));

    auto shapes = _root.add_child<SVG::Group>();
    *shapes << SVG::Polygon(scale_points(topPiecePoints, 0, _scale));
    auto topPiece = _root.get_children<SVG::Polygon>().back();
    topPiece->set_attr("class", "ferrite");
    *shapes << SVG::Polygon(scale_points(bottomPiecePoints, 0, _scale));
    auto bottomPiece = _root.get_children<SVG::Polygon>().back();
    bottomPiece->set_attr("class", "ferrite");
    for (auto& chunk : gapChunks) {
        *shapes << SVG::Polygon(scale_points(chunk, 0, _scale));
        auto chunkPiece = _root.get_children<SVG::Polygon>().back();
        chunkPiece->set_attr("class", "ferrite");
    }

    _root.autoscale();
}

void BasicPainter::paint_toroidal_core(Core core) {
    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});

    double strokeWidth = mainColumn.get_width();
    double circleDiameter = processedDescription.get_width() - strokeWidth;

    std::string cssClassName = generate_random_string();
    _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_ferrite()), std::regex("0x"), "#"));
    paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr);

    // _root.autoscale();
    _root.set_attr("width", _imageWidth * _scale).set_attr("height", _imageHeight * _scale);  // TODO remove
    _root.set_attr("viewBox", std::to_string(-_imageWidth / 2 * _scale) + " " + std::to_string(-_imageHeight / 2 * _scale) + " " + std::to_string(_imageWidth * _scale) + " " + std::to_string(_imageHeight * _scale));  // TODO remove
}

void BasicPainter::paint_two_piece_set_margin(Magnetic magnetic) {
    auto sections = magnetic.get_coil().get_sections_description().value();
    for (size_t i = 0; i < sections.size(); ++i){
        if (sections[i].get_margin()) {
            auto margins = Coil::resolve_margin(sections[i]);
            if (margins[0] > 0) {
                auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
                auto bobbinProcessedDescription = bobbin.get_processed_description().value();
                std::vector<double> bobbinCoordinates = std::vector<double>({0, 0, 0});
                if (bobbinProcessedDescription.get_coordinates()) {
                    bobbinCoordinates = bobbinProcessedDescription.get_coordinates().value();
                }
                auto windingWindowDimensions = bobbin.get_winding_window_dimensions();
                auto windingWindowCoordinates = bobbin.get_winding_window_coordinates();
                auto sectionsOrientation = bobbin.get_winding_window_sections_orientation();
                double xCoordinate;
                double yCoordinate;
                double marginWidth;
                double marginHeight;
                if (sectionsOrientation == WindingOrientation::OVERLAPPING) {
                    xCoordinate = sections[i].get_coordinates()[0];
                    yCoordinate = bobbinCoordinates[1] + windingWindowCoordinates[1] + windingWindowDimensions[1] / 2 - margins[0] / 2;
                    marginWidth = sections[i].get_dimensions()[0];
                    marginHeight = margins[0];
                }
                else {
                    xCoordinate = bobbinCoordinates[0] + windingWindowCoordinates[0] - windingWindowDimensions[0] / 2 + margins[0] / 2;
                    yCoordinate = sections[i].get_coordinates()[1];
                    marginWidth =  margins[0];
                    marginHeight = sections[i].get_dimensions()[1];
                }
                paint_rectangle(xCoordinate, yCoordinate, marginWidth, marginHeight, "margin");
            }
            if (margins[1] > 0) {
                auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
                auto bobbinProcessedDescription = bobbin.get_processed_description().value();
                std::vector<double> bobbinCoordinates = std::vector<double>({0, 0, 0});
                if (bobbinProcessedDescription.get_coordinates()) {
                    bobbinCoordinates = bobbinProcessedDescription.get_coordinates().value();
                }
                auto margins = Coil::resolve_margin(sections[i]);
                auto windingWindowDimensions = bobbin.get_winding_window_dimensions();
                auto windingWindowCoordinates = bobbin.get_winding_window_coordinates();
                auto sectionsOrientation = bobbin.get_winding_window_sections_orientation();
                std::vector<std::vector<double>> marginPoints = {};
                double xCoordinate;
                double yCoordinate;
                double marginWidth;
                double marginHeight;
                if (sectionsOrientation == WindingOrientation::OVERLAPPING) {
                    xCoordinate = sections[i].get_coordinates()[0];
                    yCoordinate = bobbinCoordinates[1] + windingWindowCoordinates[1] - windingWindowDimensions[1] / 2 + margins[1] / 2;
                    marginWidth = sections[i].get_dimensions()[0];
                    marginHeight = margins[1];
                }
                else {
                    xCoordinate = bobbinCoordinates[0] + windingWindowCoordinates[0] + windingWindowDimensions[0] / 2 - margins[1] / 2;
                    yCoordinate = sections[i].get_coordinates()[1];
                    marginWidth = margins[1];
                    marginHeight = sections[i].get_dimensions()[1];
                }
                paint_rectangle(xCoordinate, yCoordinate, marginWidth, marginHeight, "margin");
            }
        }
    }
}

void BasicPainter::paint_toroidal_margin(Magnetic magnetic) {
    bool drawSpacer = settings.get_painter_draw_spacer();
    auto sections = magnetic.get_coil().get_sections_description().value();

    if (sections.size() == 1) {
        return;
    }

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    if (!magnetic.get_coil().get_sections_description()) {
        throw CoilNotProcessedException("Winding sections not created");
    }

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();

    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    auto sectionOrientation = windingWindows[0].get_sections_orientation();
    double largestThickness = 0;

    double windingWindowRadialHeight = processedDescription.get_width() / 2 - mainColumn.get_width();
    for (size_t i = 0; i < sections.size(); ++i){
        if (sections[i].get_margin()) {
            auto margins = Coil::resolve_margin(sections[i]);

            if (sectionOrientation == WindingOrientation::CONTIGUOUS) {

                if (drawSpacer) {
                    size_t nextSectionIndex = 0;
                    if (i < sections.size() - 2) {
                        nextSectionIndex = i + 2;
                    }
                    double leftMargin = Coil::resolve_margin(sections[i])[1];
                    double rightMargin = Coil::resolve_margin(sections[nextSectionIndex])[0];
                    double rectangleThickness = leftMargin + rightMargin;
                    if (rectangleThickness == 0) {
                        continue;
                    }
                    double leftAngle = sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2;
                    double rightAngle = sections[nextSectionIndex].get_coordinates()[1] - sections[nextSectionIndex].get_dimensions()[1] / 2;
                    largestThickness = std::max(largestThickness, rectangleThickness);

                    if (i >= sections.size() - 2) {
                        rightAngle += 360;
                    }
                    double rectangleAngleInRadians = (leftAngle + rightAngle) / 2 / 180 * std::numbers::pi;
                    std::vector<double> centerRadialPoint = {windingWindowRadialHeight * cos(rectangleAngleInRadians), windingWindowRadialHeight * sin(rectangleAngleInRadians)};

                    std::vector<std::vector<double>> marginPoints = {};
                    double xCoordinate = 0;
                    double yCoordinate = windingWindowRadialHeight / 2;
                    double rectangleWidth = rectangleThickness;
                    double rectangleHeight = windingWindowRadialHeight;
                    paint_rectangle(xCoordinate, yCoordinate, rectangleWidth, rectangleHeight, "spacer", nullptr, -90 + rectangleAngleInRadians * 180 / std::numbers::pi, {0, 0});
                }
                else {
                    if (margins[0] > 0) {
                        double strokeWidth = sections[i].get_dimensions()[0];
                        double circleDiameter = (windingWindowRadialHeight - sections[i].get_coordinates()[0]) * 2;

                        double angle = wound_distance_to_angle(margins[0], circleDiameter / 2 - strokeWidth / 2);
                        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                            std::string cssClassName = generate_random_string();
                            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#"));
                            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2), {0, 0});
                        }
                    }

                    if (margins[1] > 0) {
                        double strokeWidth = sections[i].get_dimensions()[0];
                        double circleDiameter = (windingWindowRadialHeight - sections[i].get_coordinates()[0]) * 2;

                        double angle = wound_distance_to_angle(margins[1], circleDiameter / 2 - strokeWidth / 2);
                        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                            std::string cssClassName = generate_random_string();
                            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#"));
                            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2 + angle), {0, 0});
                        }
                    }
                }
            }
            else {
                if (margins[0] > 0) {
                    double strokeWidth = margins[0];
                    double circleDiameter = (windingWindowRadialHeight - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2) + margins[0] / 2) * 2;

                    double angle = wound_distance_to_angle(sections[i].get_dimensions()[1], circleDiameter / 2 + strokeWidth / 2);
                    if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                        std::string cssClassName = generate_random_string();
                        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#"));
                        paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2), {0, 0});
                    }
                }
                if (margins[1] > 0) {

                    double strokeWidth = margins[1];
                    double circleDiameter = (windingWindowRadialHeight - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2) - margins[1] / 2) * 2;

                    double angle = wound_distance_to_angle(sections[i].get_dimensions()[1], circleDiameter / 2 + strokeWidth / 2);
                    if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                        std::string cssClassName = generate_random_string();
                        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#"));
                        paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2 + angle), {0, 0});
                    }
                }
            }
        }
    }

    if (drawSpacer) {
        paint_circle(0, 0, largestThickness, "spacer", nullptr);
    }
}

void BasicPainter::set_image_size(Wire wire) {
    _extraDimension = 0.1;
    double margin = std::max(wire.get_maximum_outer_width() * _extraDimension, wire.get_maximum_outer_height() * _extraDimension);
    auto showingWireWidth = wire.get_maximum_outer_width();
    if (wire.get_maximum_outer_width() > wire.get_maximum_outer_height()) {
        showingWireWidth += margin * 2;
    }
    else {
        showingWireWidth += margin;
    }

    double showingWireHeight = wire.get_maximum_outer_height() + 2 * margin;

    _imageHeight = showingWireHeight;
    _imageWidth = showingWireWidth;
}

void BasicPainter::set_image_size(Magnetic magnetic) {
    auto aux = get_image_size(magnetic);

    _imageHeight = aux[0];
    _imageWidth = aux[1];
}

void BasicPainter::paint_wire(Wire wire) {
    set_image_size(wire);
    _scale = constants.coilPainterScale * 10;

    switch (wire.get_type()) {
        case WireType::ROUND:
            paint_round_wire(_imageWidth / 2, 0, wire);
            break;
        case WireType::LITZ:
            paint_litz_wire(_imageWidth / 2, 0, wire);
            break;
        case WireType::PLANAR:
        case WireType::FOIL:
        case WireType::RECTANGULAR:
            paint_rectangular_wire(_imageWidth / 2, 0, wire, 0, {0, 0});
            break;
        default:
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown error");
    }
    _root.autoscale();
    _root.set_attr("width", _imageWidth * _scale).set_attr("height", _imageHeight * _scale);  // TODO remove
}

void BasicPainter::paint_core(Magnetic magnetic) {
    Core core = magnetic.get_core();
    set_image_size(magnetic);
    CoreShape shape = core.resolve_shape();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return paint_toroidal_core(core);
            break;
        default:
            return paint_two_piece_set_core(core);
            break;
    }
}

void BasicPainter::paint_bobbin(Magnetic magnetic) {
    Core core = magnetic.get_core();
    _imageHeight = core.get_processed_description()->get_height();
    CoreShape shape = core.resolve_shape();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            break;
        default:
            return paint_two_piece_set_bobbin(magnetic);
            break;
    }
}

void BasicPainter::paint_coil_sections(Magnetic magnetic) {
    Core core = magnetic.get_core();
    CoreShape shape = core.resolve_shape();
    _imageHeight = core.get_processed_description()->get_height();
    auto windingWindows = core.get_winding_windows();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return paint_toroidal_coil_sections(magnetic);
            break;
        default:
            return paint_two_piece_set_coil_sections(magnetic);
            break;
    }
}

void BasicPainter::paint_coil_layers(Magnetic magnetic) {
    Core core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core has not been processed");
    }
    _imageHeight = core.get_processed_description()->get_height();

    CoreShape shape = core.resolve_shape();
    auto windingWindows = core.get_winding_windows();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return paint_toroidal_coil_layers(magnetic);
            break;
        default:
            return paint_two_piece_set_coil_layers(magnetic);
            break;
    }
}

void BasicPainter::paint_coil_turns(Magnetic magnetic) {
    Core core = magnetic.get_core();
    CoreShape shape = core.resolve_shape();
    auto windingWindows = core.get_winding_windows();
    _imageHeight = core.get_processed_description()->get_height();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return paint_toroidal_coil_turns(magnetic);
            break;
        default:
            return paint_two_piece_set_coil_turns(magnetic);
            break;
    }
}

std::string BasicPainter::get_color(double minimumValue, double maximumValue, std::string minimumColor, std::string maximumColor, double value) {
    // Clamp the value
    value = clamp(value, minimumValue, maximumValue);
    
    // Calculate interpolation factor (0.0 at minimumValue, 1.0 at maximumValue)
    double t = (value - minimumValue) / (maximumValue - minimumValue);
    
    // Linearly interpolate each channel
    uint32_t result = get_uint_color_from_ratio(t);
    
    // Convert back to hex string
    return uint_to_hex(result, "#");
}


void BasicPainter::paint_field_point(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string color, std::string label) {

    SVG::Group* shapes = _root.add_child<SVG::Group>();

    std::string cssClassName = generate_random_string();

    // _root.style("." + cssClassName).set_attr("opacity", 0.5).set_attr("fill", color).set_attr("stroke", color);
    _root.style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", color).set_attr("stroke", color);
    paint_rectangle(xCoordinate, yCoordinate, xDimension, yDimension, cssClassName, shapes, 0, {0, 0}, label);
}


void BasicPainter::paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    set_image_size(magnetic);
    double minimumModule = DBL_MAX;
    double maximumModule = 0;
    _fieldPainted = true;
    std::vector<double> modules;
    // settings.set_painter_number_points_x(4);
    // settings.set_painter_number_points_y(4);
    bool logarithmicScale = settings.get_painter_logarithmic_scale();

    ComplexField field;
    if (inputField) {
        field = inputField.value();
    }
    else {
        field = calculate_magnetic_field(operatingPoint, magnetic, harmonicIndex);
    }

    auto [pixelXDimension, pixelYDimension] = Painter::get_pixel_dimensions(magnetic);

    for (size_t i = 0; i < field.get_data().size(); ++i) {
        auto datum = field.get_data()[i];

        double value;
        if (logarithmicScale) {
            value = hypot(log10(fabs(datum.get_real())), log10(fabs(datum.get_imaginary())));
        }
        else {
            value = hypot(datum.get_real(), datum.get_imaginary());
        }
        modules.push_back(value);
    }
    std::sort(modules.begin(), modules.end());
    size_t index05 = static_cast<size_t>(0.02 * (modules.size() - 1));
    size_t index95 = static_cast<size_t>(0.98 * (modules.size() - 1));
    double percentile05Value = modules[index05];
    double percentile95Value = modules[index95];

    if (!settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = percentile95Value;
    }
    if (!settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = percentile05Value;
    }

    if (settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = settings.get_painter_maximum_value_colorbar().value();
    }
    if (settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = settings.get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }

    auto magneticFieldMinimumColor = settings.get_painter_color_magnetic_field_minimum();
    auto magneticFieldMaximumColor = settings.get_painter_color_magnetic_field_maximum();
    
    for (auto datum : field.get_data()) {
        double value;
        if (logarithmicScale) {
            value = hypot(log10(fabs(datum.get_real())), log10(fabs(datum.get_imaginary())));
        }
        else {
            value = hypot(datum.get_real(), datum.get_imaginary());
        }
        auto color = get_color(minimumModule, maximumModule, magneticFieldMinimumColor, magneticFieldMaximumColor, value);

        std::stringstream stream;
        stream << std::scientific << std::setprecision(1) << value;
        std::string label = stream.str() + " A/m";
        paint_field_point(datum.get_point()[0], datum.get_point()[1], pixelXDimension, pixelYDimension, color, label);
    }
}

void BasicPainter::paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<Field> inputField) {
    set_image_size(magnetic);
    double minimumModule = DBL_MAX;
    double maximumModule = 0;
    _fieldPainted = true;
    std::vector<double> modules;
    // settings.set_painter_number_points_x(4);
    // settings.set_painter_number_points_y(4);
    bool logarithmicScale = settings.get_painter_logarithmic_scale();

    Field field;
    if (inputField) {
        field = inputField.value();
    }
    else {
        field = calculate_electric_field(operatingPoint, magnetic, harmonicIndex);
    }

    auto [pixelXDimension, pixelYDimension] = Painter::get_pixel_dimensions(magnetic);

    for (size_t i = 0; i < field.get_data().size(); ++i) {
        auto datum = field.get_data()[i];

        double value;
        if (logarithmicScale) {
            value = log10(fabs(datum.get_value()));
        }
        else {
            value = datum.get_value();
        }
        modules.push_back(value);
    }
    std::sort(modules.begin(), modules.end());
    size_t index05 = static_cast<size_t>(0.02 * (modules.size() - 1));
    size_t index95 = static_cast<size_t>(0.98 * (modules.size() - 1));
    double percentile05Value = modules[index05];
    double percentile95Value = modules[index95];

    if (!settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = percentile95Value;
    }
    if (!settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = percentile05Value;
    }

    if (settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = settings.get_painter_maximum_value_colorbar().value();
    }
    if (settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = settings.get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }

    auto magneticFieldMinimumColor = settings.get_painter_color_magnetic_field_minimum();
    auto magneticFieldMaximumColor = settings.get_painter_color_magnetic_field_maximum();

    auto windingWindow = magnetic.get_mutable_core().get_winding_window();

    if (windingWindow.get_width()) {

        std::string cssClassName = generate_random_string();

        auto color = get_color(minimumModule, maximumModule, magneticFieldMinimumColor, magneticFieldMaximumColor, minimumModule);
        color = std::regex_replace(std::string(color), std::regex("0x"), "#");
        _root.style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", color);

        paint_rectangle(windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2, windingWindow.get_coordinates().value()[1], windingWindow.get_width().value(), windingWindow.get_height().value(), cssClassName);
    }
    else {
        throw NotImplementedException("Not implemented yet");
    }

    for (auto datum : field.get_data()) {
        double value;
        if (logarithmicScale) {
            value = log10(fabs(datum.get_value()));
        }
        else {
            value = datum.get_value();
        }
        auto color = get_color(minimumModule, maximumModule, magneticFieldMinimumColor, magneticFieldMaximumColor, value);

        std::stringstream stream;
        stream << std::scientific << std::setprecision(1) << value;
        std::string label = stream.str() + " V/m";
        paint_field_point(datum.get_point()[0], datum.get_point()[1], pixelXDimension, pixelYDimension, color, label);
    }
}

void BasicPainter::paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs, std::optional<OperatingPoint> operatingPoint, double temperature) {
    set_image_size(magnetic);
    auto constants = Constants();
    Coil coil = magnetic.get_coil();
    double minimumModule = DBL_MAX;
    double maximumModule = 0;
    std::vector<double> modules;
    std::vector<double> modulesToSort;
    bool logarithmicScale = settings.get_painter_logarithmic_scale();

    if (!outputs && !operatingPoint) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Missing both outputs and operatingPoint in paint_wire_losses");
    }
    if (!outputs && operatingPoint) {
        WindingLossesOutput windingLossesOutput = WindingLosses().calculate_losses(magnetic, operatingPoint.value(), temperature);
        outputs = Outputs();
        outputs->set_winding_losses(windingLossesOutput);
    }

    // Build a map from turn name to total losses for lookup
    std::map<std::string, double> lossesPerTurnByName;
    if (outputs && outputs->get_winding_losses() && outputs->get_winding_losses()->get_winding_losses_per_turn()) {
        auto windingLossesPerTurn = outputs->get_winding_losses()->get_winding_losses_per_turn().value();
        for (auto windingLossesthisTurn : windingLossesPerTurn) {
            double totalLoss = WindingLosses::get_total_winding_losses(windingLossesthisTurn);
            if (windingLossesthisTurn.get_name().has_value()) {
                lossesPerTurnByName[windingLossesthisTurn.get_name().value()] = totalLoss;
            }
            modules.push_back(totalLoss);
            modulesToSort.push_back(totalLoss);
        }
    }

    std::sort(modulesToSort.begin(), modulesToSort.end());
    size_t index05 = static_cast<size_t>(0.02 * (modulesToSort.size() - 1));
    size_t index95 = static_cast<size_t>(0.98 * (modulesToSort.size() - 1));
    double percentile05Value = modulesToSort[index05];
    double percentile95Value = modulesToSort[index95];

    if (!settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = percentile95Value;
    }
    if (!settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = percentile05Value;
    }

    if (settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumModule = settings.get_painter_maximum_value_colorbar().value();
    }
    if (settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumModule = settings.get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }

    // Swap min/max colors so low losses are blue (cold) and high losses are red (hot)
    auto windingLossesMaximumColor = settings.get_painter_color_magnetic_field_maximum();
    auto windingLossesMinimumColor = settings.get_painter_color_magnetic_field_minimum();

    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding turns not created");
    }

    auto turns = coil.get_turns_description().value();
    auto shapes = _root.add_child<SVG::Group>();

    for (size_t i = 0; i < turns.size(); ++i){
        std::string turnName = turns[i].get_name();
        
        // Look up loss value by turn name, fallback to index-based if not found
        double lossValue;
        auto it = lossesPerTurnByName.find(turnName);
        if (it != lossesPerTurnByName.end()) {
            lossValue = it->second;
        } else if (i < modules.size()) {
            lossValue = modules[i];
        } else {
            continue;  // Skip if no loss data available for this turn
        }

        double value;
        if (logarithmicScale) {
            value = log10(lossValue);
        }
        else {
            value = lossValue;
        }
        auto color = get_color(minimumModule, maximumModule, windingLossesMinimumColor, windingLossesMaximumColor, value);

        std::stringstream stream;
        stream << std::scientific << std::setprecision(2) << lossValue;
        std::string label = stream.str() + " W";

        if (turns[i].get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
            double xCoordinate = turns[i].get_coordinates()[0];
            double yCoordinate = turns[i].get_coordinates()[1];
            double diameter = turns[i].get_dimensions().value()[0];
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color);
            paint_circle(xCoordinate, yCoordinate, diameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
        }
        else {
            if (turns[i].get_dimensions().value()[0] && turns[i].get_dimensions().value()[1]) {
                double xCoordinate = turns[i].get_coordinates()[0];
                double yCoordinate = turns[i].get_coordinates()[1];
                double conductingWidth = turns[i].get_dimensions().value()[0];
                double conductingHeight = turns[i].get_dimensions().value()[1];
                std::string cssClassName = generate_random_string();
                _root.style("." + cssClassName).set_attr("fill", color);
                paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, cssClassName, shapes, 0, {0, 0}, label);
            }
        }
    }
}

void BasicPainter::paint_temperature_field(Magnetic magnetic, const std::map<std::string, double>& nodeTemperatures, bool showColorBar, ColorPalette palette, double ambientTemperature) {
    set_image_size(magnetic);
    
    if (nodeTemperatures.empty()) {
        return;
    }
    
    // Find temperature range for color mapping (skip NaN/inf values)
    double minimumTemperature = DBL_MAX;
    double maximumTemperature = -DBL_MAX;
    for (const auto& [name, temp] : nodeTemperatures) {
        if (!std::isfinite(temp)) continue;  // Skip invalid temperatures
        if (temp < minimumTemperature) minimumTemperature = temp;
        if (temp > maximumTemperature) maximumTemperature = temp;
    }
    
    // If no valid temperatures found, fall back to ambient
    if (minimumTemperature == DBL_MAX || maximumTemperature == -DBL_MAX) {
        minimumTemperature = ambientTemperature;
        maximumTemperature = ambientTemperature + 1.0;
    }
    
    // Use the provided ambient temperature as the color scale minimum
    // This ensures the color scale shows the full range from ambient to max temp
    minimumTemperature = ambientTemperature;
    
    // Apply colorbar settings if provided
    if (settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumTemperature = settings.get_painter_minimum_value_colorbar().value();
    }
    if (settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumTemperature = settings.get_painter_maximum_value_colorbar().value();
    }
    if (minimumTemperature == maximumTemperature || !std::isfinite(minimumTemperature) || !std::isfinite(maximumTemperature)) {
        minimumTemperature = ambientTemperature;
        maximumTemperature = ambientTemperature + 1.0;
    }
    
    // Helper lambda to get color from temperature using the selected palette
    auto getColorForTemperature = [&](double temp) -> std::string {
        if (!std::isfinite(temp)) {
            // Return gray for invalid temperatures
            return "#808080";
        }
        double ratio = (temp - minimumTemperature) / (maximumTemperature - minimumTemperature);
        ratio = clamp(ratio, 0.0, 1.0);
        uint32_t colorUint = get_uint_color_from_palette(ratio, palette);
        return uint_to_hex(colorUint, "#");
    };
    
    Coil coil = magnetic.get_coil();
    Core core = magnetic.get_mutable_core();
    
    auto shapes = _root.add_child<SVG::Group>();
    
    // Add ambient temperature background for entire SVG (centered at origin)
    auto ambientColor = uint_to_hex(get_uint_color_from_palette(0.0, palette), "#");
    std::string bgClassName = generate_random_string();
    _root.style("." + bgClassName).set_attr("fill", ambientColor).set_attr("opacity", "1.0");
    paint_rectangle(0, 0, 10000, 10000, bgClassName, shapes);
    
    // Helper function to find a node temperature by exact name match
    auto findNodeTemperatureExact = [&nodeTemperatures](const std::string& name) -> std::optional<double> {
        auto it = nodeTemperatures.find(name);
        if (it != nodeTemperatures.end()) {
            return it->second;
        }
        return std::nullopt;
    };
    
    // Helper function to find a node temperature by partial name match
    auto findNodeTemperature = [&nodeTemperatures](const std::string& prefix) -> std::optional<double> {
        for (const auto& [name, temp] : nodeTemperatures) {
            if (name.find(prefix) != std::string::npos) {
                return temp;
            }
        }
        return std::nullopt;
    };
    
    // Paint core regions with their temperatures
    // Look for core-related node names (e.g., "Core_Column_0", "Core_Top_Yoke", etc.)
    auto processedElements = core.get_processed_description().value();
    auto columns = processedElements.get_columns();
    auto family = core.get_shape_family();
    
    // Get the main column and calculate 2D view geometry (same as paint_two_piece_set_core)
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto rightColumn = core.find_closest_column_by_coordinates(
        std::vector<double>({processedElements.get_width() / 2, 0, -processedElements.get_depth() / 2}));
    
    double showingCoreWidth;
    double showingMainColumnWidth;
    switch (family) {
        case MAS::CoreShapeFamily::C:
        case MAS::CoreShapeFamily::U:
        case MAS::CoreShapeFamily::UR:
            showingMainColumnWidth = mainColumn.get_width() / 2;
            showingCoreWidth = processedElements.get_width() - mainColumn.get_width() / 2;
            break;
        default:
            showingCoreWidth = processedElements.get_width() / 2;
            showingMainColumnWidth = mainColumn.get_width() / 2;
    }
    
    double rightColumnWidth;
    if (rightColumn.get_minimum_width()) {
        rightColumnWidth = rightColumn.get_minimum_width().value();
    }
    else {
        rightColumnWidth = rightColumn.get_width();
    }
    
    double coreHeight = processedElements.get_height();
    double mainColumnHeight = mainColumn.get_height();
    
    // Paint columns with matching 2D view geometry
    // Skip for toroidal cores as they use segments instead
    if (family != MAS::CoreShapeFamily::T) {
        for (size_t i = 0; i < columns.size(); ++i) {
            auto column = columns[i];
            std::string nodeName = "Core_Column_" + std::to_string(i);
        
        auto tempOpt = findNodeTemperature(nodeName);
        if (!tempOpt) {
            tempOpt = findNodeTemperature("Core");
        }
        
        if (tempOpt) {
            double temp = tempOpt.value();
            auto color = getColorForTemperature(temp);
            
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << temp;
            std::string label = nodeName + ": " + stream.str() + " °C";
            
            // Determine column position in 2D view (matching paint_two_piece_set_core)
            double xCoord, colWidth, colHeight;
            if (column.get_type() == ColumnType::CENTRAL) {
                // Central column: x from 0 to mainColumnWidth/2
                xCoord = showingMainColumnWidth / 2;
                colWidth = showingMainColumnWidth;
                colHeight = column.get_height();
            } else {
                // Lateral column: x from showingCoreWidth-rightColumnWidth to showingCoreWidth
                xCoord = showingCoreWidth - rightColumnWidth / 2;
                colWidth = rightColumnWidth;
                colHeight = column.get_height();
            }
            
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 1.0);
            paint_rectangle(xCoord, 0, colWidth, colHeight, cssClassName, shapes, 0, {0, 0}, label);
        }
        }
    }
    
    // Paint top yoke
    // Skip for toroidal cores as they don't have yokes
    if (family != MAS::CoreShapeFamily::T) {
        auto topYokeTempOpt = findNodeTemperatureExact("Core_Top_Yoke");
        if (!topYokeTempOpt) {
            topYokeTempOpt = findNodeTemperature("Core");
        }
        if (topYokeTempOpt) {
        double temp = topYokeTempOpt.value();
        auto color = getColorForTemperature(temp);
        
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << temp;
        std::string label = "Core_Top_Yoke: " + stream.str() + " °C";
        
        // Top yoke: from y=mainColumnHeight/2 to y=coreHeight/2, full width
        double yokeHeight = (coreHeight - mainColumnHeight) / 2;
        double yokeY = mainColumnHeight / 2 + yokeHeight / 2;
        
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 1.0);
            paint_rectangle(showingCoreWidth / 2, yokeY, showingCoreWidth, yokeHeight, cssClassName, shapes, 0, {0, 0}, label);
        }
    }
    
    // Paint bottom yoke
    // Skip for toroidal cores as they don't have yokes
    if (family != MAS::CoreShapeFamily::T) {
    auto bottomYokeTempOpt = findNodeTemperatureExact("Core_Bottom_Yoke");
    if (!bottomYokeTempOpt) {
        bottomYokeTempOpt = findNodeTemperature("Core");
    }
    if (bottomYokeTempOpt) {
        double temp = bottomYokeTempOpt.value();
        auto color = getColorForTemperature(temp);
        
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << temp;
        std::string label = "Core_Bottom_Yoke: " + stream.str() + " °C";
        
        // Bottom yoke: from y=-coreHeight/2 to y=-mainColumnHeight/2, full width
        double yokeHeight = (coreHeight - mainColumnHeight) / 2;
        double yokeY = -(mainColumnHeight / 2 + yokeHeight / 2);
        
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 1.0);
            paint_rectangle(showingCoreWidth / 2, yokeY, showingCoreWidth, yokeHeight, cssClassName, shapes, 0, {0, 0}, label);
        }
    }
    
    // Paint toroidal core segments if toroidal shape
    if (family == MAS::CoreShapeFamily::T) {
        // Get toroidal core dimensions
        // For toroids, the processedElements.get_width() gives outer diameter
        // and mainColumn.get_width() gives the core cross-section thickness
        double outerDiameter = processedElements.get_width();
        double coreThickness = mainColumn.get_width();
        double innerDiameter = outerDiameter - 2 * coreThickness;
        double outerRadius = outerDiameter / 2;
        double innerRadius = innerDiameter / 2;
        
        // Paint winding window (inside toroid) with ambient temperature color
        auto ambientColor = getColorForTemperature(ambientTemperature);
        
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << ambientTemperature;
        std::string label = "Winding Window: " + stream.str() + " °C";
        
        // For toroidal winding window, paint inner circle
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("fill", ambientColor).set_attr("opacity", "1.0");
        paint_circle(0, 0, innerRadius, cssClassName, shapes, 360, 0, {0, 0}, label);
        
        // Find the maximum turn radius to determine outer ring extent
        double maxTurnRadius = innerRadius;
        if (coil.get_turns_description()) {
            auto turns = coil.get_turns_description().value();
            for (const auto& turn : turns) {
                double x = turn.get_coordinates()[0];
                double y = turn.get_coordinates()[1];
                double turnRadius = std::sqrt(x * x + y * y);
                if (turn.get_dimensions().value()[0]) {
                    double wireDiameter = turn.get_dimensions().value()[0];
                    turnRadius += wireDiameter / 2;
                }
                maxTurnRadius = std::max(maxTurnRadius, turnRadius);
                
                // Check additional coordinates
                if (turn.get_additional_coordinates()) {
                    auto addCoords = turn.get_additional_coordinates().value();
                    for (const auto& coord : addCoords) {
                        double xAdd = coord[0];
                        double yAdd = coord[1];
                        double addRadius = std::sqrt(xAdd * xAdd + yAdd * yAdd);
                        if (turn.get_dimensions().value()[0]) {
                            double wireDiameter = turn.get_dimensions().value()[0];
                            addRadius += wireDiameter / 2;
                        }
                        maxTurnRadius = std::max(maxTurnRadius, addRadius);
                    }
                }
            }
        }
        
        // Number of segments (should match Temperature.cpp)
        const int numSegments = 8;
        const double anglePerSegment = 2.0 * std::numbers::pi / numSegments;
        
        // Draw each segment as an arc
        for (int i = 0; i < numSegments; ++i) {
            std::string nodeName = "Core_Segment_" + std::to_string(i);
            auto tempOpt = findNodeTemperatureExact(nodeName);
            
            if (tempOpt) {
                double temp = tempOpt.value();
                auto color = getColorForTemperature(temp);
                
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << temp;
                std::string label = nodeName + ": " + stream.str() + " °C";
                
                // Calculate start and end angles (start from top, go clockwise)
                // In SVG, angles are measured clockwise from 3 o'clock
                // We want to start from top (12 o'clock = -90 degrees)
                double startAngle = -std::numbers::pi / 2.0 + i * anglePerSegment;
                double endAngle = startAngle + anglePerSegment;
                
                // Convert to degrees for SVG
                double startAngleDeg = startAngle * 180.0 / std::numbers::pi;
                double endAngleDeg = endAngle * 180.0 / std::numbers::pi;
                
                // Calculate points for the arc path
                // Start point on outer circle
                double x1 = outerRadius * std::cos(startAngle);
                double y1 = outerRadius * std::sin(startAngle);
                
                // End point on outer circle
                double x2 = outerRadius * std::cos(endAngle);
                double y2 = outerRadius * std::sin(endAngle);
                
                // End point on inner circle
                double x3 = innerRadius * std::cos(endAngle);
                double y3 = innerRadius * std::sin(endAngle);
                
                // Start point on inner circle
                double x4 = innerRadius * std::cos(startAngle);
                double y4 = innerRadius * std::sin(startAngle);
                
                // Create SVG path for the arc segment (annular sector)
                // M = move to, A = arc, L = line, Z = close path
                std::stringstream pathData;
                pathData << "M " << (x1 * _scale) << "," << (-y1 * _scale)  // Move to start on outer circle
                         << " A " << (outerRadius * _scale) << "," << (outerRadius * _scale)  // Arc on outer circle
                         << " 0 0 0 "  // x-axis-rotation, large-arc-flag, sweep-flag (counter-clockwise)
                         << (x2 * _scale) << "," << (-y2 * _scale)  // End point on outer circle
                         << " L " << (x3 * _scale) << "," << (-y3 * _scale)  // Line to inner circle
                         << " A " << (innerRadius * _scale) << "," << (innerRadius * _scale)  // Arc on inner circle
                         << " 0 0 1 "  // clockwise for proper annular sector
                         << (x4 * _scale) << "," << (-y4 * _scale)  // Back to start on inner circle
                         << " Z";  // Close path
                
                // Add the path to the SVG
                std::string cssClassName = generate_random_string();
                _root.style("." + cssClassName)
                    .set_attr("fill", color)
                    .set_attr("opacity", "1.0");
                auto* path = shapes->add_child<SVG::Path>();
                path->set_attr("d", pathData.str());
                path->set_attr("class", cssClassName);
                path->add_child<SVG::Title>(label);
            }
        }
    }
    
    // Paint bobbin if present - bobbin is a variant<Bobbin, string>
    // Skip bobbin for toroidal cores as they don't have a traditional bobbin structure
    auto bobbinVariant = magnetic.get_coil().get_bobbin();
    if (family != MAS::CoreShapeFamily::T && std::holds_alternative<Bobbin>(bobbinVariant)) {
        
        // Draw each bobbin wall separately with its actual dimensions
        // Look for bobbin nodes by their specific names
        for (const auto& [name, temp] : nodeTemperatures) {
            if (name.find("Bobbin_") == 0) {
                auto color = getColorForTemperature(temp);
                
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << temp;
                std::string label = name + ": " + stream.str() + " °C";
                
                std::string cssClassName = generate_random_string();
                _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", "1.0");
                
                // Get bobbin wall dimensions from the node temperature data
                // The node positions and dimensions are stored in the Temperature object
                // Use bobbin's own coordinate system (not winding window coordinates)
                auto windingWindow = core.get_winding_window();
                if (windingWindow.get_width()) {
                    double wwX = windingWindow.get_coordinates().value()[0];
                    double wwY = windingWindow.get_coordinates().value()[1];
                    double wwWidth = windingWindow.get_width().value();
                    double wwHeight = windingWindow.get_height().value();
                    
                    // Get bobbin processed description for actual dimensions
                    auto bobbinProcessed = std::get<Bobbin>(bobbinVariant).get_processed_description();
                    double wallThickness = 0.002;  // Default 2mm
                    double columnThickness = 0.002;  // Default 2mm
                    double columnWidth = 0.005;  // Default 5mm (bobbin depth)
                    if (bobbinProcessed) {
                        wallThickness = bobbinProcessed->get_wall_thickness();
                        columnThickness = bobbinProcessed->get_column_thickness();
                        if (bobbinProcessed->get_column_width()) {
                            columnWidth = bobbinProcessed->get_column_width().value();
                        }
                    }
                    
                    // Calculate bobbin column position
                    // The bobbin sits on top of the core column
                    // Column left edge is at: showingMainColumnWidth (edge of core column)
                    // Column right edge is at: showingMainColumnWidth + columnThickness
                    double columnLeftEdge = showingMainColumnWidth;
                    double columnRightEdge = showingMainColumnWidth + columnThickness;
                    
                    // Calculate bobbin yoke vertical positions
                    // Yokes sit on top/bottom of the core column (where core yokes end)
                    double coreColumnTopY = mainColumnHeight / 2;
                    double coreColumnBottomY = -mainColumnHeight / 2;
                    
                    // First pass: draw the bobbin column wall
                    if (name.find("CentralColumn") != std::string::npos) {
                        // Central column wall - vertical strip
                        // Positioned at the right edge of the core column
                        double xCoord = (columnLeftEdge + columnRightEdge) / 2;
                        double yCoord = wwY;
                        double width = columnThickness;
                        double height = wwHeight;
                        paint_rectangle(xCoord, yCoord, width, height, cssClassName, shapes, 0, {0, 0}, label);
                    }
                }
            }
        }
    }
    
    // Second pass: draw bobbin yokes on top of the column
    auto bobbinVariant2 = magnetic.get_coil().get_bobbin();
    if (family != MAS::CoreShapeFamily::T && std::holds_alternative<Bobbin>(bobbinVariant2)) {
        for (const auto& [name, temp] : nodeTemperatures) {
            if (name.find("Bobbin_") == 0 && (name.find("TopYoke") != std::string::npos || name.find("BottomYoke") != std::string::npos)) {
                auto color = getColorForTemperature(temp);
                
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << temp;
                std::string label = name + ": " + stream.str() + " °C";
                
                std::string cssClassName = generate_random_string();
                _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", "1.0");
                
                auto windingWindow = core.get_winding_window();
                if (windingWindow.get_width()) {
                    double wwX = windingWindow.get_coordinates().value()[0];
                    double wwY = windingWindow.get_coordinates().value()[1];
                    double wwWidth = windingWindow.get_width().value();
                    double wwHeight = windingWindow.get_height().value();
                    
                    auto bobbinProcessed = std::get<Bobbin>(bobbinVariant2).get_processed_description();
                    double wallThickness = 0.002;
                    double columnThickness = 0.002;
                    if (bobbinProcessed) {
                        wallThickness = bobbinProcessed->get_wall_thickness();
                        columnThickness = bobbinProcessed->get_column_thickness();
                    }
                    
                    double columnLeftEdge = showingMainColumnWidth;
                    double coreColumnTopY = mainColumnHeight / 2;
                    double coreColumnBottomY = -mainColumnHeight / 2;
                    
                    double xCoord, yCoord, width, height;
                    
                    if (name.find("TopYoke") != std::string::npos) {
                        // Top yoke - full winding window width, flush against core
                        xCoord = columnLeftEdge + wwWidth / 2;
                        yCoord = coreColumnTopY - wallThickness / 2;
                        width = wwWidth;
                        height = wallThickness;
                        paint_rectangle(xCoord, yCoord, width, height, cssClassName, shapes, 0, {0, 0}, label);
                    } else if (name.find("BottomYoke") != std::string::npos) {
                        // Bottom yoke - full winding window width, flush against core
                        xCoord = columnLeftEdge + wwWidth / 2;
                        yCoord = coreColumnBottomY + wallThickness / 2;
                        width = wwWidth;
                        height = wallThickness;
                        paint_rectangle(xCoord, yCoord, width, height, cssClassName, shapes, 0, {0, 0}, label);
                    }
                }
            }
        }
    }
    
    // Paint insulation layers with their temperatures
    // Look for insulation layer nodes (e.g., "L_0", "L_1" for concentric, "IL_0_0_i", "IL_0_0_o" for toroidal)
    for (const auto& [name, temp] : nodeTemperatures) {
        // Handle concentric core insulation layers (L_0, L_1, etc.)
        if (name.find("L_") == 0 && name.find("IL_") != 0) {
            // Extract layer index from node name (L_<index>)
            std::string idxStr = name.substr(2);  // Skip "L_"
            size_t layerIdx = std::stoull(idxStr);
            
            // Find the insulation layer by index
            if (coil.get_layers_description()) {
                auto layers = coil.get_layers_description().value();
                size_t insulationIdx = 0;
                for (const auto& layer : layers) {
                    if (layer.get_type() == MAS::ElectricalType::INSULATION) {
                        if (insulationIdx == layerIdx) {
                            auto color = getColorForTemperature(temp);
                            
                            std::stringstream stream;
                            stream << std::fixed << std::setprecision(1) << temp;
                            std::string label = name + ": " + stream.str() + " °C";
                            
                            // Get layer geometry
                            if (layer.get_coordinates().size() >= 2 && layer.get_dimensions().size() >= 2) {
                                double xCoord = layer.get_coordinates()[0];
                                double yCoord = layer.get_coordinates()[1];
                                double width = layer.get_dimensions()[0];
                                double height = layer.get_dimensions()[1];
                                
                                // Draw insulation layer filled with temperature color
                                std::string cssClassName = generate_random_string();
                                _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", "1.0");
                                paint_rectangle(xCoord, yCoord, width, height, cssClassName, shapes, 0, {0, 0}, label);
                            }
                            break;
                        }
                        insulationIdx++;
                    }
                }
            }
        }
        // Handle toroidal insulation layer nodes (IL_<layer>_<segment>_<i/o>)
        else if (name.find("IL_") == 0 && family == MAS::CoreShapeFamily::T) {
            auto color = getColorForTemperature(temp);
            
            // Parse node name: IL_<layer>_<segment>_<i/o>
            std::vector<std::string> parts;
            std::stringstream ss(name);
            std::string part;
            while (std::getline(ss, part, '_')) {
                parts.push_back(part);
            }
            
            if (parts.size() >= 4) {
                size_t layerIdx = std::stoull(parts[1]);
                size_t segmentIdx = std::stoull(parts[2]);
                std::string innerOuter = parts[3]; // "i" or "o"
                
                // Find the insulation layer
                if (coil.get_layers_description()) {
                    auto layers = coil.get_layers_description().value();
                    
                    // Get core dimensions to calculate initialRadius (core outer radius)
                    auto processedDescription = core.get_processed_description().value();
                    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
                    double initialRadius = processedDescription.get_width() / 2 - mainColumn.get_width();
                    
                    size_t insulationIdx = 0;
                    for (const auto& layer : layers) {
                        if (layer.get_type() == MAS::ElectricalType::INSULATION) {
                            if (insulationIdx == layerIdx) {
                                // Calculate radius based on inner (_i) or outer (_o) node
                                double layerRadialPos = layer.get_coordinates()[0]; // Distance from winding window edge
                                double radialThickness = layer.get_dimensions()[0];
                                if (radialThickness < 1e-6) {
                                    radialThickness = 0.00005; // Default 50 microns
                                }
                                
                                double radius;
                                if (innerOuter == "i") {
                                    // Inner node - at inner surface of insulation (closer to core)
                                    radius = initialRadius - layerRadialPos - radialThickness;
                                } else {
                                    // Outer node - at outer surface of insulation (farther from core)
                                    // Use additionalCoordinates if available (like in paint_toroidal_coil_turns)
                                    if (layer.get_additional_coordinates()) {
                                        double outerRadialPos = layer.get_additional_coordinates().value()[0][0];
                                        radius = initialRadius - outerRadialPos;
                                    } else {
                                        // Fallback: use inner position
                                        radius = initialRadius - layerRadialPos;
                                    }
                                }
                                
                                double circleDiameter = radius * 2;
                                double anglePerSegment = 360.0 / 8.0; // 8 segments
                                double startAngle = segmentIdx * anglePerSegment;
                                double endAngle = (segmentIdx + 1) * anglePerSegment;
                                
                                // Create style with all attributes at once
                                std::string cssClassName = generate_random_string();
                                _root.style("." + cssClassName)
                                    .set_attr("stroke", color)
                                    .set_attr("fill", "none")
                                    .set_attr("opacity", "1.0")
                                    .set_attr("stroke-width", std::to_string(radialThickness * _scale));
                                
                                paint_circle(0, 0, circleDiameter / 2, cssClassName, shapes, 
                                            endAngle - startAngle, -startAngle - (endAngle - startAngle) / 2, {0, 0});
                                break;
                            }
                            insulationIdx++;
                        }
                    }
                }
            }
        }
    }
    
    // Paint coil turns with their temperatures
    if (coil.get_turns_description()) {
        auto turns = coil.get_turns_description().value();
        
        for (size_t i = 0; i < turns.size(); ++i) {
            auto& turn = turns[i];
            std::string turnName = turn.get_name();
            
            // Try to find temperature for this specific turn
            std::optional<double> tempOpt;
            
            // First try exact match with turn name
            auto it = nodeTemperatures.find(turnName);
            if (it != nodeTemperatures.end()) {
                tempOpt = it->second;
            } else {
                // Try with _inner or _outer suffix for toroidal turns with split nodes
                // Note: Node names use capitalized "_Inner" and "_Outer"
                auto innerIt = nodeTemperatures.find(turnName + "_Inner");
                auto outerIt = nodeTemperatures.find(turnName + "_Outer");
                if (innerIt != nodeTemperatures.end() && outerIt != nodeTemperatures.end()) {
                    // Average the inner and outer temperatures
                    tempOpt = (innerIt->second + outerIt->second) / 2.0;
                } else if (innerIt != nodeTemperatures.end()) {
                    tempOpt = innerIt->second;
                } else if (outerIt != nodeTemperatures.end()) {
                    tempOpt = outerIt->second;
                } else {
                    // Fallback: try layer or coil-wide temperature
                    tempOpt = findNodeTemperature("Coil_Layer");
                    if (!tempOpt) {
                        tempOpt = findNodeTemperature("Coil");
                    }
                }
            }
            
            if (tempOpt) {
                double temp = tempOpt.value();
                auto color = getColorForTemperature(temp);
                
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << temp;
                std::string label = turnName + ": " + stream.str() + " °C";
                
                if (turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
                    double xCoordinate = turn.get_coordinates()[0];
                    double yCoordinate = turn.get_coordinates()[1];
                    double diameter = turn.get_dimensions().value()[0];
                    std::string cssClassName = generate_random_string();
                    _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", "1.0");
                    paint_circle(xCoordinate, yCoordinate, diameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
                } else {
                    if (turn.get_dimensions().value()[0] && turn.get_dimensions().value()[1]) {
                        double xCoordinate = turn.get_coordinates()[0];
                        double yCoordinate = turn.get_coordinates()[1];
                        double conductingWidth = turn.get_dimensions().value()[0];
                        double conductingHeight = turn.get_dimensions().value()[1];
                        
                        // Get rotation angle if available (for rectangular wires in toroidal cores)
                        double turnAngle = 0;
                        std::vector<double> turnCenter = {xCoordinate, -yCoordinate};
                        if (turn.get_rotation()) {
                            turnAngle = turn.get_rotation().value();
                        }
                        
                        std::string cssClassName = generate_random_string();
                        _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", "1.0");
                        paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, cssClassName, shapes, turnAngle, turnCenter, label);
                    }
                }
                
                // Paint additional coordinates if present (e.g., outer half of toroidal turns)
                if (turn.get_additional_coordinates()) {
                    auto additionalCoords = turn.get_additional_coordinates().value();
                    
                    for (size_t addIdx = 0; addIdx < additionalCoords.size(); ++addIdx) {
                        // For toroidal turns, additional coords have _Outer suffix in thermal nodes
                        std::string additionalNodeName = turnName + "_Outer";
                        auto addTempIt = nodeTemperatures.find(additionalNodeName);
                        
                        double additionalTemp = temp;  // Default to same as main turn
                        std::string additionalLabel = turnName + "_outer";
                        
                        if (addTempIt != nodeTemperatures.end()) {
                            additionalTemp = addTempIt->second;
                            std::stringstream addStream;
                            addStream << std::fixed << std::setprecision(1) << additionalTemp;
                            additionalLabel = turnName + "_outer: " + addStream.str() + " °C";
                        }
                        
                        auto additionalColor = getColorForTemperature(additionalTemp);
                        
                        if (turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
                            double xCoord = additionalCoords[addIdx][0];
                            double yCoord = additionalCoords[addIdx][1];
                            double diameter = turn.get_dimensions().value()[0];
                            std::string cssClassName = generate_random_string();
                            _root.style("." + cssClassName).set_attr("fill", additionalColor).set_attr("opacity", "1.0");
                            paint_circle(xCoord, yCoord, diameter / 2, cssClassName, shapes, 360, 0, {0, 0}, additionalLabel);
                        } else {
                            if (turn.get_dimensions().value()[0] && turn.get_dimensions().value()[1]) {
                                double xCoord = additionalCoords[addIdx][0];
                                double yCoord = additionalCoords[addIdx][1];
                                double conductingWidth = turn.get_dimensions().value()[0];
                                double conductingHeight = turn.get_dimensions().value()[1];
                                
                                // Get rotation angle if available (for rectangular wires in toroidal cores)
                                double turnAngle = 0;
                                std::vector<double> turnCenter = {xCoord, -yCoord};
                                if (turn.get_rotation()) {
                                    turnAngle = turn.get_rotation().value();
                                }
                                
                                std::string cssClassName = generate_random_string();
                                _root.style("." + cssClassName).set_attr("fill", additionalColor).set_attr("opacity", "1.0");
                                paint_rectangle(xCoord, yCoord, conductingWidth, conductingHeight, cssClassName, shapes, turnAngle, turnCenter, additionalLabel);
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Add temperature color bar on the right side
    if (showColorBar) {
        // Find the extent of all geometry
        double minX = 0, maxX = 0;
        double minY = 0, maxY = 0;
        
        // Check core geometry
        if (family == MAS::CoreShapeFamily::T) {
            // Toroidal: centered at origin, symmetric
            double halfWidth = processedElements.get_width() / 2;
            minX = -halfWidth;
            maxX = halfWidth;
            minY = -halfWidth;
            maxY = halfWidth;
        } else {
            // Concentric cores: start at x=0
            minX = 0;
            maxX = showingCoreWidth;
            // Add margin to top and bottom for concentric cores
            double verticalMargin = coreHeight * 0.1;  // 10% margin on top and bottom
            minY = -coreHeight / 2 - verticalMargin;
            maxY = coreHeight / 2 + verticalMargin;
        }
        
        // Check turn positions to find actual extent
        if (coil.get_turns_description()) {
            auto turns = coil.get_turns_description().value();
            for (const auto& turn : turns) {
                auto coords = turn.get_coordinates();
                if (coords.size() >= 2) {
                    double turnX = coords[0];
                    double turnY = coords[1];
                    double turnRadius = 0;
                    if (turn.get_dimensions().has_value() && !turn.get_dimensions().value().empty()) {
                        turnRadius = turn.get_dimensions().value()[0] / 2;
                    }
                    minX = std::min(minX, turnX - turnRadius);
                    maxX = std::max(maxX, turnX + turnRadius);
                    minY = std::min(minY, turnY - turnRadius);
                    maxY = std::max(maxY, turnY + turnRadius);
                    
                    // For toroidal turns, also check additional coordinates (outer half of turn)
                    if (family == MAS::CoreShapeFamily::T && turn.get_additional_coordinates()) {
                        auto addCoords = turn.get_additional_coordinates().value();
                        for (const auto& addCoord : addCoords) {
                            if (addCoord.size() >= 2) {
                                double addX = addCoord[0];
                                double addY = addCoord[1];
                                minX = std::min(minX, addX - turnRadius);
                                maxX = std::max(maxX, addX + turnRadius);
                                minY = std::min(minY, addY - turnRadius);
                                maxY = std::max(maxY, addY + turnRadius);
                            }
                        }
                    }
                }
            }
        }
        
        // Position color bar on the right with some margin
        double marginFactor = 0.15;  // 15% margin from rightmost geometry
        double barWidth = (maxX - 0) * 0.08;  // 8% of horizontal extent
        double barHeight = (maxY - minY) * 0.7;  // 70% of vertical extent
        double barX = maxX + (maxX - 0) * marginFactor;
        double barY = (maxY + minY) / 2;
        
        // Draw white background box for color bar area (very tight fit)
        double textHeight = barHeight * 0.08;
        double textWidthEstimate = textHeight * 3.5;  // ~3-4 chars for "100°C"
        double whiteBgPaddingX = barWidth * 0.03;  // Minimal horizontal padding
        double whiteBgPaddingY = textHeight * 0.2;  // Minimal vertical padding
        // Total width: bar + small gap + text + minimal padding on both sides
        double whiteBgWidth = barWidth + barWidth * 0.15 + textWidthEstimate + whiteBgPaddingX * 2;
        double whiteBgHeight = barHeight + whiteBgPaddingY * 2;
        // Position: left edge at barX - barWidth/2 (start of color bar), centered
        double whiteBgCenterX = (barX - barWidth / 2) + whiteBgWidth / 2;
        std::string whiteBgClass = generate_random_string();
        _root.style("." + whiteBgClass).set_attr("fill", "#FFFFFF").set_attr("opacity", "1.0");
        paint_rectangle(whiteBgCenterX, barY, whiteBgWidth, whiteBgHeight, whiteBgClass, shapes);
        
        size_t numSteps = 20;
        
        // Draw color bar rectangles
        for (size_t i = 0; i < numSteps; ++i) {
            double t = static_cast<double>(i) / numSteps;
            double temp = minimumTemperature + t * (maximumTemperature - minimumTemperature);
            auto color = getColorForTemperature(temp);
            
            double stepHeight = barHeight / numSteps;
            // Y goes from bottom (low temp at i=0) to top (high temp at i=numSteps-1)
            double stepY = barY - barHeight / 2 + (i + 0.5) * stepHeight;
            
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color);
            paint_rectangle(barX, stepY, barWidth, stepHeight, cssClassName, shapes);
        }
        
        // Add outline around color bar
        std::string outlineClass = generate_random_string();
        _root.style("." + outlineClass).set_attr("fill", "none").set_attr("stroke", "#000000").set_attr("stroke-width", "1");
        paint_rectangle(barX, barY, barWidth, barHeight, outlineClass, shapes);
        
        // Add temperature labels
        std::vector<std::pair<double, std::string>> labels;
        labels.push_back({1.0, std::to_string(static_cast<int>(std::round(maximumTemperature))) + "°C"});
        labels.push_back({0.5, std::to_string(static_cast<int>(std::round((minimumTemperature + maximumTemperature) / 2))) + "°C"});
        labels.push_back({0.0, std::to_string(static_cast<int>(std::round(minimumTemperature))) + "°C"});
        
        for (const auto& [position, label] : labels) {
            double labelY = barY - barHeight / 2 + position * barHeight;
            // Adjust top label position to stay within white box
            if (position == 1.0) {
                labelY -= textHeight * 0.8;  // Move top label down
            }
            // Text position needs to be in scaled coordinates for SVG::Text
            // Position text closer to the color bar
            double textX = (barX + barWidth * 1.7) * _scale;
            double textY = -labelY * _scale;  // Note: SVG Y is inverted
            auto* text = _root.add_child<SVG::Text>(textX, textY, label);
            text->set_attr("font-size", std::to_string(barHeight * _scale * 0.08));
            text->set_attr("fill", "#000000");
            text->set_attr("text-anchor", "start");
        }
        
        // Update viewBox to include the color bar
        double colorBarRight = barX - barWidth/2 + whiteBgWidth + whiteBgPaddingX;
        double colorBarTop = barY + barHeight / 2 + whiteBgPaddingY;
        double colorBarBottom = barY - barHeight / 2 - whiteBgPaddingY;
        
        // Find overall bounds
        double newMinX = minX;
        double newMaxX = std::max(maxX, colorBarRight);
        double newMinY = std::min(minY, colorBarBottom);
        double newMaxY = std::max(maxY, colorBarTop);
        
        double viewBoxWidth = newMaxX - newMinX;
        double viewBoxHeight = newMaxY - newMinY;
        
        // Update the image dimensions and viewBox
        _imageWidth = viewBoxWidth;
        _imageHeight = viewBoxHeight;
        _root.set_attr("width", _imageWidth * _scale).set_attr("height", _imageHeight * _scale);
        _root.set_attr("viewBox", std::to_string(newMinX * _scale) + " " + std::to_string(-newMaxY * _scale) + " " + std::to_string(viewBoxWidth * _scale) + " " + std::to_string(viewBoxHeight * _scale));
    }
}

std::string BasicPainter::export_svg() {
    if (!_filepath.empty()) {
        if (!std::filesystem::exists(_filepath)) {
            std::filesystem::create_directory(_filepath);
        }
        std::ofstream outfile(_filepath.replace_filename(_filename));
        outfile << std::string(_root);
    }
    return std::string(_root);
}

void BasicPainter::paint_waveform_svg(
    const Waveform& waveform,
    const std::string& name,
    const std::string& color,
    double xOffset,
    double yOffset,
    double plotWidth,
    double plotHeight) {
    
    auto data = waveform.get_data();
    if (data.empty()) return;
    
    std::optional<std::vector<double>> timeOpt = waveform.get_time();
    std::vector<double> time;
    if (timeOpt && !timeOpt->empty()) {
        time = timeOpt.value();
    } else {
        // Generate time axis from 0 to 1
        time.resize(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            time[i] = static_cast<double>(i) / static_cast<double>(data.size() - 1);
        }
    }
    
    // Find min/max for scaling
    double minVal = *std::min_element(data.begin(), data.end());
    double maxVal = *std::max_element(data.begin(), data.end());
    double minTime = *std::min_element(time.begin(), time.end());
    double maxTime = *std::max_element(time.begin(), time.end());
    
    // Add padding
    double valRange = maxVal - minVal;
    if (valRange < 1e-12) valRange = 1.0;
    double timeRange = maxTime - minTime;
    if (timeRange < 1e-12) timeRange = 1.0;
    
    minVal -= valRange * 0.1;
    maxVal += valRange * 0.1;
    valRange = maxVal - minVal;
    
    // Create path string
    std::stringstream pathData;
    for (size_t i = 0; i < data.size(); ++i) {
        double x = xOffset + ((time[i] - minTime) / timeRange) * plotWidth;
        double y = yOffset + plotHeight - ((data[i] - minVal) / valRange) * plotHeight;
        
        if (i == 0) {
            pathData << "M " << x << " " << y;
        } else {
            pathData << " L " << x << " " << y;
        }
    }
    
    // Add path element
    std::string cssClassName = generate_random_string();
    _root.style("." + cssClassName)
        .set_attr("fill", "none")
        .set_attr("stroke", color)
        .set_attr("stroke-width", "1.5");
    
    auto* path = _root.add_child<SVG::Path>();
    path->set_attr("d", pathData.str());
    path->set_attr("class", cssClassName);
    
    // Add label
    auto* text = _root.add_child<SVG::Text>(xOffset + 5, yOffset + 15, name);
    text->set_attr("font-size", "12");
    text->set_attr("fill", color);
    
    // Add axis lines
    std::string axisClass = generate_random_string();
    _root.style("." + axisClass)
        .set_attr("stroke", "#888888")
        .set_attr("stroke-width", "1")
        .set_attr("fill", "none");
    
    // SVG::Line constructor takes (x1, x2, y1, y2) - NOT (x1, y1, x2, y2)
    
    // X-axis (bottom) - horizontal line from left to right at bottom of plot
    auto* xAxis = _root.add_child<SVG::Line>(xOffset, xOffset + plotWidth, yOffset + plotHeight, yOffset + plotHeight);
    xAxis->set_attr("class", axisClass);
    
    // Y-axis (left) - vertical line from top to bottom at left of plot
    auto* yAxis = _root.add_child<SVG::Line>(xOffset, xOffset, yOffset, yOffset + plotHeight);
    yAxis->set_attr("class", axisClass);
    
    // X-axis top border
    auto* xAxisTop = _root.add_child<SVG::Line>(xOffset, xOffset + plotWidth, yOffset, yOffset);
    xAxisTop->set_attr("class", axisClass);
    
    // Y-axis right border  
    auto* yAxisRight = _root.add_child<SVG::Line>(xOffset + plotWidth, xOffset + plotWidth, yOffset, yOffset + plotHeight);
    yAxisRight->set_attr("class", axisClass);
    
    // Add value labels
    std::stringstream maxLabel, minLabel;
    maxLabel << std::scientific << std::setprecision(2) << maxVal;
    minLabel << std::scientific << std::setprecision(2) << minVal;
    
    auto* maxText = _root.add_child<SVG::Text>(xOffset - 60, yOffset + 12, maxLabel.str());
    maxText->set_attr("font-size", "10");
    maxText->set_attr("fill", "#666666");
    
    auto* minText = _root.add_child<SVG::Text>(xOffset - 60, yOffset + plotHeight, minLabel.str());
    minText->set_attr("font-size", "10");
    minText->set_attr("fill", "#666666");
}

std::string BasicPainter::paint_operating_point_waveforms(
    const OperatingPoint& operatingPoint,
    const std::string& title,
    double width,
    double height) {
    
    std::vector<Waveform> waveforms;
    std::vector<std::string> waveformNames;
    
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    
    // Collect all windings' voltage and current waveforms
    for (size_t i = 0; i < excitations.size(); ++i) {
        const auto& excitation = excitations[i];
        std::string windingName = excitation.get_name().value_or("Winding " + std::to_string(i));
        
        // Add voltage waveform if present
        if (excitation.get_voltage() && excitation.get_voltage()->get_waveform()) {
            waveforms.push_back(excitation.get_voltage()->get_waveform().value());
            waveformNames.push_back(windingName + " Voltage (V)");
        }
        
        // Add current waveform if present
        if (excitation.get_current() && excitation.get_current()->get_waveform()) {
            waveforms.push_back(excitation.get_current()->get_waveform().value());
            waveformNames.push_back(windingName + " Current (A)");
        }
    }
    
    if (waveforms.empty()) {
        return "";  // No waveforms to plot
    }
    
    // Reset root SVG
    _root = SVG::SVG();
    _root.set_attr("width", std::to_string(static_cast<int>(width)));
    _root.set_attr("height", std::to_string(static_cast<int>(height)));
    _root.set_attr("viewBox", "0 0 " + std::to_string(static_cast<int>(width)) + " " + std::to_string(static_cast<int>(height)));
    
    // Add white background
    auto* bg = _root.add_child<SVG::Rect>(0, 0, width, height);
    bg->set_attr("fill", "#ffffff");
    
    // Add title
    auto* titleText = _root.add_child<SVG::Text>(width / 2, 25, title);
    titleText->set_attr("font-size", "16");
    titleText->set_attr("font-weight", "bold");
    titleText->set_attr("text-anchor", "middle");
    titleText->set_attr("fill", "#333333");
    
    // Calculate layout
    size_t numPlots = waveforms.size();
    double margin = 80;
    double plotSpacing = 20;
    double availableHeight = height - margin - 40;  // Top and bottom margins
    double plotHeight = (availableHeight - (numPlots - 1) * plotSpacing) / numPlots;
    double plotWidth = width - 2 * margin;
    
    // Color palette
    std::vector<std::string> colors = {
        "#1f77b4",  // Blue
        "#ff7f0e",  // Orange
        "#2ca02c",  // Green
        "#d62728",  // Red
        "#9467bd",  // Purple
        "#8c564b",  // Brown
        "#e377c2",  // Pink
        "#7f7f7f",  // Gray
        "#bcbd22",  // Olive
        "#17becf"   // Cyan
    };
    
    // Draw each waveform
    for (size_t i = 0; i < numPlots; ++i) {
        double yOffset = 50 + i * (plotHeight + plotSpacing);
        std::string color = colors[i % colors.size()];
        
        paint_waveform_svg(waveforms[i], waveformNames[i], color, margin, yOffset, plotWidth, plotHeight);
    }
    
    return export_svg();
}

std::string BasicPainter::paint_thermal_circuit_schematic(
    const std::vector<ThermalNetworkNode>& nodes,
    const std::vector<ThermalResistanceElement>& resistances,
    double width,
    double height,
    bool showQuadrantLabels) {
    
    // Configurable colors for thermal visualization
    const std::string COLOR_AMBIENT = "#4CAF50";      // Green - exposed to ambient air
    const std::string COLOR_CONDUCTION = "#FF9800";   // Orange - thermal conduction
    const std::string COLOR_NO_CONNECTION = "#795548"; // Brown - no thermal path
    
    // Create a new SVG root for the schematic
    _root = SVG::SVG();
    _root.set_attr("width", std::to_string(static_cast<int>(width)));
    _root.set_attr("height", std::to_string(static_cast<int>(height)));
    _root.set_attr("viewBox", "0 0 " + std::to_string(static_cast<int>(width)) + " " + std::to_string(static_cast<int>(height)));
    
    // Add white background
    auto* bg = _root.add_child<SVG::Rect>(0, 0, width, height);
    bg->set_attr("fill", "#ffffff");
    
    if (nodes.empty()) {
        return export_svg();
    }
    
    // Detect if this is a concentric core (not toroidal) - needed for proper scaling
    bool isConcentricCore = false;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
            node.part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
            node.part == ThermalNodePartType::CORE_TOP_YOKE ||
            node.part == ThermalNodePartType::CORE_BOTTOM_YOKE) {
            isConcentricCore = true;
            break;
        }
    }
    
    // =========================================================================
    // Calculate minimum node separation and scale factor for schematic
    // This ensures nodes don't overlap, especially with quadrant visualization
    // =========================================================================
    double minNodeSeparation = std::numeric_limits<double>::max();
    
    // Find minimum distance between ANY two non-ambient nodes (not just connected ones)
    // This is critical for dense layouts like T36 with many turns
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i].isAmbient()) continue;
        // Skip insulation layer nodes for min separation calculation
        if (nodes[i].part == ThermalNodePartType::INSULATION_LAYER) continue;
        
        for (size_t j = i + 1; j < nodes.size(); j++) {
            if (nodes[j].isAmbient()) continue;
            // Skip insulation layer nodes for min separation calculation
            if (nodes[j].part == ThermalNodePartType::INSULATION_LAYER) continue;
            
            // Calculate 3D Euclidean distance between nodes
            double dx = nodes[i].physicalCoordinates.size() >= 1 && nodes[j].physicalCoordinates.size() >= 1 ? 
                        nodes[i].physicalCoordinates[0] - nodes[j].physicalCoordinates[0] : 0.0;
            double dy = nodes[i].physicalCoordinates.size() >= 2 && nodes[j].physicalCoordinates.size() >= 2 ? 
                        nodes[i].physicalCoordinates[1] - nodes[j].physicalCoordinates[1] : 0.0;
            double dz = nodes[i].physicalCoordinates.size() >= 3 && nodes[j].physicalCoordinates.size() >= 3 ? 
                        nodes[i].physicalCoordinates[2] - nodes[j].physicalCoordinates[2] : 0.0;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist > 1e-6) {  // Avoid zero distances
                minNodeSeparation = std::min(minNodeSeparation, dist);
            }
        }
    }
    
    // Define node visualization requirements (in SVG units)
    double nodeCircleRadius = 8.0;       // Radius of node circles (base)
    double nodeQuadrantRadius = nodeCircleRadius * 1.15;  // Extended radius with quadrants
    double minClearance = 12.0;          // Minimum clearance between node edges
    
    // Calculate required SVG space between node centers
    // We need: 2 node radii (with quadrants) + minimum clearance
    double minSvgSpaceNeeded = 2 * nodeQuadrantRadius + minClearance;
    
    // Calculate scale factor automatically from minimum node separation
    // scaleFactor = SVG_pixels / physical_meters
    // For the smallest separation, we need at least minSvgSpaceNeeded pixels
    double scaleFactor = 5000.0;  // Fallback if no valid distances found
    
    if (minNodeSeparation < std::numeric_limits<double>::max() && minNodeSeparation > 1e-6) {
        // Calculate exact scale needed for minimum separation
        scaleFactor = minSvgSpaceNeeded / minNodeSeparation;
        
        // Add 150% safety margin for comfortable spacing and clear resistor display
        scaleFactor *= 2.5;
        
        // Ensure minimum practical scale
        scaleFactor = std::max(1000.0, scaleFactor);
    }
    
    // Find bounding box of all nodes in physical coordinates
    // Skip nodes with invalid coordinates (outside reasonable range of -10 to 10 meters)
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    
    const double REASONABLE_COORD_LIMIT = 10.0;  // 10 meters is way larger than any magnetic component
    
    for (const auto& node : nodes) {
        if (!node.isAmbient() && node.physicalCoordinates.size() >= 2) {
            // Skip insulation layer nodes for bounding box calculation
            // They can be positioned far from turns and cause excessive scaling
            if (node.part == ThermalNodePartType::INSULATION_LAYER) {
                continue;
            }
            
            double x = node.physicalCoordinates[0];
            double y = node.physicalCoordinates[1];
            
            // Skip nodes with unreasonable coordinates (likely corrupted or uninitialized)
            if (std::abs(x) > REASONABLE_COORD_LIMIT || std::abs(y) > REASONABLE_COORD_LIMIT) {
                std::cerr << "Warning: Node '" << node.name << "' has invalid coordinates (" 
                          << x << ", " << y << "), skipping in schematic" << std::endl;
                continue;
            }
            
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
    }
    
    // Calculate required SVG canvas size based on physical dimensions
    double physicalWidth = (maxX - minX) * scaleFactor;
    double physicalHeight = (maxY - minY) * scaleFactor;
    
    // Add margins
    double margin = 80;
    double schematicWidth = physicalWidth + 2 * margin;
    double schematicHeight = physicalHeight + 2 * margin;  // Removed extra space for legend
    
    // Update SVG dimensions if using magnetic coordinates
    if (schematicWidth > 100 && schematicHeight > 100) {
        _root.set_attr("width", std::to_string(static_cast<int>(schematicWidth)));
        _root.set_attr("height", std::to_string(static_cast<int>(schematicHeight)));
        _root.set_attr("viewBox", "0 0 " + std::to_string(static_cast<int>(schematicWidth)) + " " + 
                                    std::to_string(static_cast<int>(schematicHeight)));
        
        // Update background
        bg->set_attr("width", std::to_string(schematicWidth));
        bg->set_attr("height", std::to_string(schematicHeight));
        
        // Use scaled magnetic coordinates for layout
        width = schematicWidth;
        height = schematicHeight;
    }
    
    // =========================================================================
    // End of scaling calculation
    // =========================================================================
    
    // Categorize nodes by part type
    std::vector<const ThermalNetworkNode*> coreNodes;
    std::map<int, std::vector<const ThermalNetworkNode*>> coilNodesByWinding;
    
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::AMBIENT) {
            // Ambient node - no visualization needed
            continue;
        }
        
        // Skip nodes with invalid coordinates
        if (node.physicalCoordinates.size() >= 2) {
            double x = node.physicalCoordinates[0];
            double y = node.physicalCoordinates[1];
            if (std::abs(x) > REASONABLE_COORD_LIMIT || std::abs(y) > REASONABLE_COORD_LIMIT) {
                continue;  // Skip invalid nodes
            }
        }
        
        if (node.part == ThermalNodePartType::TURN) {
            // Extract winding index from windingIndex field or name
            int windingIdx = 0;
            if (node.windingIndex.has_value()) {
                windingIdx = node.windingIndex.value();
            }
            coilNodesByWinding[windingIdx].push_back(&node);
        } else {
            // All other parts are core nodes
            coreNodes.push_back(&node);
        }
    }
    
    // =========================================================================
    // GEOMETRIC LAYOUT: Use actual physical coordinates with calculated scaling
    // =========================================================================
    
    // Map node coordinates to SVG positions
    auto mapNodeToSvg = [&](const ThermalNetworkNode& node) -> std::pair<double, double> {
        if (node.physicalCoordinates.size() < 2) {
            // Fallback for nodes without coordinates
            return {width / 2, height / 2};
        }
        
        // Check for unreasonable coordinates
        double x = node.physicalCoordinates[0];
        double y = node.physicalCoordinates[1];
        if (std::abs(x) > REASONABLE_COORD_LIMIT || std::abs(y) > REASONABLE_COORD_LIMIT) {
            // Return a position off-screen for invalid nodes
            return {-1000, -1000};
        }
        
        // Map physical coordinates to SVG space
        // Physical: (node.physicalCoordinates[0], node.physicalCoordinates[1]) in meters
        // SVG: (margin to width-margin, margin to height-margin)
        
        double svgX = margin + (x - minX) * scaleFactor;
        double svgY = margin + (y - minY) * scaleFactor;
        
        return {svgX, svgY};
    };
    
    // Layout parameters
    double defaultNodeRadius = 8;  // Default size for core/ambient nodes
    double groundY = height - 60;
    
    // === GEOMETRIC LAYOUT ===
    // Nodes positioned based on their physical coordinates
    // Resistors drawn as connections between node positions
    // Ambient ground at bottom
    
    // Ground bus and ground symbol removed for cleaner visualization
    
    // === DRAW THERMAL RESISTANCES AS CONNECTIONS ===
    auto* resistorGroup = _root.add_child<SVG::Group>();
    resistorGroup->set_attr("id", "resistances");
    
    // Helper to get limit coordinate for a specific quadrant face from a node
    auto getQuadrantLimitCoordinate = [&](const ThermalNetworkNode& node, ThermalNodeFace face) -> std::optional<std::array<double, 3>> {
        for (const auto& quadrant : node.quadrants) {
            if (quadrant.face == face) {
                return quadrant.limitCoordinates;
            }
        }
        return std::nullopt;
    };
    
    // Helper function to draw a resistor symbol between two points
    auto drawResistor = [&](double x1, double y1, double x2, double y2, double resistance) {
        // Calculate line direction
        double dx = x2 - x1;
        double dy = y2 - y1;
        double length = std::sqrt(dx*dx + dy*dy);
        
        if (length < 1e-6) return;  // Skip zero-length resistors
        
        // Normalize
        double ux = dx / length;
        double uy = dy / length;
        
        // Perpendicular unit vector
        double px = -uy;
        double py = ux;
        
        // Resistor zigzag parameters
        double zigzagLength = 16.0;
        double zigzagWidth = 4.0;
        int numZigs = 4;
        
        // Start and end clearance from nodes
        double clearance = defaultNodeRadius + 2.0;
        
        // Calculate where zigzag should be centered between the two nodes
        double availableLength = length - 2 * clearance;
        double zigzagStartOffset = clearance + (availableLength - zigzagLength) / 2;
        
        // Zigzag start and end positions
        double zigzagStartX = x1 + ux * zigzagStartOffset;
        double zigzagStartY = y1 + uy * zigzagStartOffset;
        double zigzagEndX = zigzagStartX + ux * zigzagLength;
        double zigzagEndY = zigzagStartY + uy * zigzagLength;
        
        // Draw line from node1 to zigzag start
        auto* wire1 = resistorGroup->add_child<SVG::Line>(x1, zigzagStartX, y1, zigzagStartY);
        wire1->set_attr("stroke", "#666666");
        wire1->set_attr("stroke-width", "1");
        
        // Draw resistor zigzag
        std::stringstream pathData;
        pathData << "M " << zigzagStartX << " " << zigzagStartY;
        
        for (int i = 1; i <= numZigs; ++i) {
            double t = static_cast<double>(i) / numZigs;
            double sign = (i % 2 == 1) ? 1.0 : -1.0;
            double midX = zigzagStartX + ux * zigzagLength * t + px * zigzagWidth * sign;
            double midY = zigzagStartY + uy * zigzagLength * t + py * zigzagWidth * sign;
            pathData << " L " << midX << " " << midY;
        }
        pathData << " L " << zigzagEndX << " " << zigzagEndY;
        
        auto* zigzag = resistorGroup->add_child<SVG::Path>();
        zigzag->set_attr("d", pathData.str());
        zigzag->set_attr("stroke", "#666666");
        zigzag->set_attr("stroke-width", "1.5");
        zigzag->set_attr("fill", "none");
        
        // Draw line from zigzag end to node2
        auto* wire2 = resistorGroup->add_child<SVG::Line>(zigzagEndX, x2, zigzagEndY, y2);
        wire2->set_attr("stroke", "#666666");
        wire2->set_attr("stroke-width", "1");
        
        // Label resistance value at midpoint
        double labelX = (x1 + x2) / 2 + px * 10;
        double labelY = (y1 + y2) / 2 + py * 10;
        
        std::stringstream rStr;
        if (resistance >= 1000.0) {
            rStr << std::fixed << std::setprecision(1) << resistance / 1000.0 << "k";
        } else if (resistance >= 1.0) {
            rStr << std::fixed << std::setprecision(1) << resistance;
        } else {
            rStr << std::fixed << std::setprecision(0) << resistance * 1000.0 << "m";
        }
        rStr << "K/W";
        
        auto* rLabel = resistorGroup->add_child<SVG::Text>(labelX, labelY, rStr.str());
        rLabel->set_attr("font-size", "6");
        rLabel->set_attr("fill", "#cc6600");
    };
    
    // Draw all resistances
    for (const auto& res : resistances) {
        const ThermalNetworkNode* node1 = nullptr;
        const ThermalNetworkNode* node2 = nullptr;
        if (res.nodeFromId < nodes.size()) node1 = &nodes[res.nodeFromId];
        if (res.nodeToId < nodes.size()) node2 = &nodes[res.nodeToId];
        
        if (!node1 || !node2) continue;
        if (node1->isAmbient() || node2->isAmbient()) continue;
        
        // Helper to get connection point for a node/quadrant
        auto getConnectionPoint = [&](const ThermalNetworkNode& node, ThermalNodeFace face) -> std::pair<double, double> {
            auto pos = mapNodeToSvg(node);
            
            // For nodes without specific faces, just return center position
            if (face == ThermalNodeFace::NONE) {
                return pos;
            }
            
            // Only TURN and INSULATION_LAYER nodes use limit coordinates for connection points
            bool useLimitCoords = (node.part == ThermalNodePartType::TURN || 
                                   node.part == ThermalNodePartType::INSULATION_LAYER);
            if (!useLimitCoords) {
                return pos;
            }
            
            // For turns and insulation layers, calculate offset based on face
            // For concentric: use horizontal/vertical offsets
            // For toroidal: use angled limit coordinates
            double offsetX = 0, offsetY = 0;
            
            if (isConcentricCore) {
                // Always try to use limit coordinates first (for both turns and insulation layers)
                auto limit = getQuadrantLimitCoordinate(node, face);
                if (limit && ((*limit)[0] != 0.0 || (*limit)[1] != 0.0)) {
                    double centerX = node.physicalCoordinates[0];
                    double centerY = node.physicalCoordinates[1];
                    offsetX = ((*limit)[0] - centerX) / 2.0 * scaleFactor;
                    offsetY = ((*limit)[1] - centerY) / 2.0 * scaleFactor;
                    return {pos.first + offsetX, pos.second + offsetY};
                }
                
                // Fallback: fixed horizontal/vertical offsets based on dimensions
                double nodeWidth = node.dimensions.width * scaleFactor / 2.0;
                double nodeHeight = node.dimensions.height * scaleFactor / 2.0;
                double radius = (nodeWidth + nodeHeight) / 4.0;  // Node radius in SVG coords
                
                switch (face) {
                    case ThermalNodeFace::RADIAL_INNER:   // Left (-X)
                        offsetX = -radius; offsetY = 0; break;
                    case ThermalNodeFace::RADIAL_OUTER:   // Right (+X)
                        offsetX = radius; offsetY = 0; break;
                    case ThermalNodeFace::TANGENTIAL_LEFT:  // Top (+Y)
                        offsetX = 0; offsetY = radius; break;
                    case ThermalNodeFace::TANGENTIAL_RIGHT: // Bottom (-Y)
                        offsetX = 0; offsetY = -radius; break;
                    default: break;
                }
                return {pos.first + offsetX, pos.second + offsetY};
            } else {
                // Toroidal: use angled limit coordinates from 3D geometry
                auto limit = getQuadrantLimitCoordinate(node, face);
                if (limit && ((*limit)[0] != 0.0 || (*limit)[1] != 0.0)) {
                    double centerX = node.physicalCoordinates[0];
                    double centerY = node.physicalCoordinates[1];
                    offsetX = ((*limit)[0] - centerX) / 2.0 * scaleFactor;
                    offsetY = ((*limit)[1] - centerY) / 2.0 * scaleFactor;
                    return {pos.first + offsetX, pos.second + offsetY};
                }
            }
            return pos;
        };
        
        // Get connection points
        auto [x1, y1] = getConnectionPoint(*node1, res.quadrantFrom);
        auto [x2, y2] = getConnectionPoint(*node2, res.quadrantTo);
        
        // Connection points calculated - proceed with drawing
        
        // Draw conduction as full resistor symbols
        if (res.type == HeatTransferType::CONDUCTION) {
            drawResistor(x1, y1, x2, y2, res.resistance);
        } else {
            // Use dashed line for convection/radiation
            auto* line = resistorGroup->add_child<SVG::Line>(x1, x2, y1, y2);
            if (res.type == HeatTransferType::NATURAL_CONVECTION || 
                res.type == HeatTransferType::FORCED_CONVECTION) {
                line->set_attr("stroke", "#3399ff");
                line->set_attr("stroke-dasharray", "3,2");
            } else {
                line->set_attr("stroke", "#ff6600");
                line->set_attr("stroke-dasharray", "2,1");
            }
            line->set_attr("stroke-width", "1");
            line->set_attr("opacity", "0.6");
        }
    }
    
    // === DRAW NODES ===
    auto* nodesGroup = _root.add_child<SVG::Group>();
    nodesGroup->set_attr("id", "nodes");
    
    // Draw core nodes (placeholder - will add quadrants after helper functions)
    std::vector<const ThermalNetworkNode*> coreNodesForQuadrants = coreNodes;
    
    // Helper to get surface coverage for a specific face from the node's quadrants
    auto getSurfaceCoverage = [&](const ThermalNetworkNode& node, ThermalNodeFace face) -> double {
        for (const auto& quadrant : node.quadrants) {
            if (quadrant.face == face) {
                return quadrant.surfaceCoverage;
            }
        }
        return 1.0;  // Default: fully exposed if quadrant not found
    };
    
    // Helper to check if a specific face has convection to ambient
    auto hasConvectionToAmbientOnFace = [&](size_t nodeIdx, ThermalNodeFace face, double nodeAngle) -> bool {
        const auto& node = nodes[nodeIdx];
        for (const auto& res : resistances) {
            if ((res.nodeFromId == nodeIdx || res.nodeToId == nodeIdx) &&
                (res.type == HeatTransferType::NATURAL_CONVECTION || 
                 res.type == HeatTransferType::FORCED_CONVECTION)) {
                size_t otherIdx = (res.nodeFromId == nodeIdx) ? res.nodeToId : res.nodeFromId;
                if (otherIdx < nodes.size() && nodes[otherIdx].isAmbient()) {
                    // For turn nodes: check which face connects to ambient
                    if (node.part == ThermalNodePartType::TURN) {
                        // Check which face actually has the convection connection
                        ThermalNodeFace connectedFace = (res.nodeFromId == nodeIdx) ? res.quadrantFrom : res.quadrantTo;
                        if (face == connectedFace) return true;  // Found matching face
                        // Otherwise continue checking other resistances
                    }
                    // For core nodes: check which faces connect to ambient AND have surface coverage > 0
                    else if (node.part == ThermalNodePartType::CORE_TOROIDAL_SEGMENT) {
                        // Toroidal core: only RADIAL_INNER and RADIAL_OUTER face ambient
                        // Check surface coverage - if fully covered (coverage = 0), no convection
                        if (face == ThermalNodeFace::RADIAL_INNER || face == ThermalNodeFace::RADIAL_OUTER) {
                            double coverage = getSurfaceCoverage(node, face);
                            return coverage > 0.01;  // Exposed if coverage > 1%
                        }
                        return false;  // Tangential faces don't have ambient convection
                    }
                    else if (node.part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
                             node.part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
                             node.part == ThermalNodePartType::CORE_TOP_YOKE ||
                             node.part == ThermalNodePartType::CORE_BOTTOM_YOKE) {
                        return true;  // Other core types: all exposed surfaces
                    }
                    // For insulation layer nodes: only radial faces connect to ambient
                    else if (node.part == ThermalNodePartType::INSULATION_LAYER) {
                        // Check which face actually has the convection connection
                        ThermalNodeFace connectedFace = (res.nodeFromId == nodeIdx) ? res.quadrantFrom : res.quadrantTo;
                        // Only RADIAL_INNER and RADIAL_OUTER faces should show convection
                        if (face == connectedFace && 
                            (face == ThermalNodeFace::RADIAL_INNER || face == ThermalNodeFace::RADIAL_OUTER)) {
                            return true;
                        }
                        // TANGENTIAL faces should not show convection (left unconnected)
                        return false;
                    }
                }
            }
        }
        return false;
    };
    
    // Helper to check if a specific face has conduction connection
    auto hasConductionConnectionOnFace = [&](size_t nodeIdx, ThermalNodeFace face, double nodeAngle) -> bool {
        for (const auto& res : resistances) {
            if (res.type == HeatTransferType::CONDUCTION) {
                // Check if this resistance connects to the given node on the given face
                if (res.nodeFromId == nodeIdx && res.quadrantFrom == face) {
                    return true;
                }
                if (res.nodeToId == nodeIdx && res.quadrantTo == face) {
                    return true;
                }
            }
        }
        return false;
    };
    
    // Helper to draw a pie slice (ring quadrant)
    auto drawPieSlice = [&](double cx, double cy, double radius, double startAngleDeg, double endAngleDeg, 
                            const std::string& color, double opacity) {
        double startAngle = startAngleDeg * M_PI / 180.0;
        double endAngle = endAngleDeg * M_PI / 180.0;
        
        double x1 = cx + radius * std::cos(startAngle);
        double y1 = cy + radius * std::sin(startAngle);
        double x2 = cx + radius * std::cos(endAngle);
        double y2 = cy + radius * std::sin(endAngle);
        
        std::stringstream pathData;
        pathData << "M " << cx << " " << cy;
        pathData << " L " << x1 << " " << y1;
        
        int largeArc = (endAngleDeg - startAngleDeg > 180) ? 1 : 0;
        pathData << " A " << radius << " " << radius << " 0 " << largeArc << " 1 " << x2 << " " << y2;
        pathData << " Z";
        
        auto* slice = nodesGroup->add_child<SVG::Path>();
        slice->set_attr("d", pathData.str());
        slice->set_attr("fill", color);
        slice->set_attr("opacity", std::to_string(opacity));
        slice->set_attr("stroke", "#222222");
        slice->set_attr("stroke-width", "0.5");
    };
    
    // Check if insulation layers exist in the model
    bool hasInsulationLayers = false;
    for (const auto& node : nodes) {
        if (node.part == ThermalNodePartType::INSULATION_LAYER) {
            hasInsulationLayers = true;
            break;
        }
    }
    
    // Draw core nodes with quadrants
    for (const auto* node : coreNodesForQuadrants) {
        // Find node index in full nodes array
        size_t nodeIdx = 0;
        for (size_t i = 0; i < nodes.size(); i++) {
            if (&nodes[i] == node) {
                nodeIdx = i;
                break;
            }
        }
        
        auto [x, y] = mapNodeToSvg(*node);
        
        // For concentric cores, use fixed orientation:
        // RADIAL_INNER: toward central column (-X direction) -> slice at 180°
        // RADIAL_OUTER: toward lateral column (+X direction) -> slice at 0°
        // TANGENTIAL_LEFT: toward top (+Y direction) -> slice at 90°
        // TANGENTIAL_RIGHT: toward bottom (-Y direction) -> slice at 270°
        bool isConcentric = (node->part == ThermalNodePartType::CORE_CENTRAL_COLUMN ||
                            node->part == ThermalNodePartType::CORE_LATERAL_COLUMN ||
                            node->part == ThermalNodePartType::CORE_TOP_YOKE ||
                            node->part == ThermalNodePartType::CORE_BOTTOM_YOKE ||
                            node->part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN ||
                            node->part == ThermalNodePartType::BOBBIN_TOP_YOKE ||
                            node->part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE);
        
        // Node circle (base) - white with connection colors on quadrants
        auto* circle = nodesGroup->add_child<SVG::Circle>(x, y, defaultNodeRadius);
        circle->set_attr("fill", "#ffffff");
        circle->set_attr("stroke", "#333333");
        circle->set_attr("stroke-width", "1.5");
        
        // Draw 4 slices on each core node
        std::array<ThermalNodeFace, 4> quadrants = {
            ThermalNodeFace::RADIAL_OUTER,
            ThermalNodeFace::RADIAL_INNER,
            ThermalNodeFace::TANGENTIAL_LEFT,
            ThermalNodeFace::TANGENTIAL_RIGHT
        };
        
        for (int q = 0; q < 4; q++) {
            double sliceStartAngle = 0;
            
            if (isConcentric) {
                // Fixed orientation for concentric cores
                switch (quadrants[q]) {
                    case ThermalNodeFace::RADIAL_OUTER:  // Toward lateral column (+X)
                        sliceStartAngle = -45;
                        break;
                    case ThermalNodeFace::RADIAL_INNER:  // Toward central column (-X)
                        sliceStartAngle = 135;
                        break;
                    case ThermalNodeFace::TANGENTIAL_LEFT:  // Toward top (+Y)
                        sliceStartAngle = 45;
                        break;
                    case ThermalNodeFace::TANGENTIAL_RIGHT:  // Toward bottom (-Y)
                        sliceStartAngle = 225;
                        break;
                    default:
                        sliceStartAngle = q * 90 - 45;
                        break;
                }
            } else {
                // Toroidal: Orient slices based on node's angular position
                double nodeAngle = std::atan2(node->physicalCoordinates.size() >= 2 ? node->physicalCoordinates[1] : 0,
                                             node->physicalCoordinates.size() >= 1 ? node->physicalCoordinates[0] : 1);
                switch (quadrants[q]) {
                    case ThermalNodeFace::RADIAL_OUTER:
                        sliceStartAngle = nodeAngle * 180.0 / M_PI - 45;
                        break;
                    case ThermalNodeFace::RADIAL_INNER:
                        sliceStartAngle = (nodeAngle + M_PI) * 180.0 / M_PI - 45;
                        break;
                    case ThermalNodeFace::TANGENTIAL_LEFT:
                        sliceStartAngle = (nodeAngle + M_PI/2) * 180.0 / M_PI - 45;
                        break;
                    case ThermalNodeFace::TANGENTIAL_RIGHT:
                        sliceStartAngle = (nodeAngle - M_PI/2) * 180.0 / M_PI - 45;
                        break;
                    default:
                        sliceStartAngle = q * 90 - 45;
                        break;
                }
            }
            double sliceEndAngle = sliceStartAngle + 90;
            
            // Determine color based on connections and surface coverage for this specific quadrant
            std::string sliceColor;
            
            // Get surface coverage for this quadrant (if available)
            double coverage = getSurfaceCoverage(*node, quadrants[q]);
            
            if (hasConductionConnectionOnFace(nodeIdx, quadrants[q], 0)) {
                sliceColor = COLOR_CONDUCTION;  // Orange: connected via conduction (priority)
            } else if (hasConvectionToAmbientOnFace(nodeIdx, quadrants[q], 0)) {
                sliceColor = COLOR_AMBIENT;  // Green: exposed to air
            } else if (!hasInsulationLayers && coverage > 0.5) {
                // Surface coverage > 50% means mostly exposed to air
                // Only use this when insulation layers are NOT present
                sliceColor = COLOR_AMBIENT;  // Green: mostly exposed
            } else {
                // No explicit connection and mostly covered
                sliceColor = COLOR_NO_CONNECTION;  // Brown: mostly covered
            }
            
            drawPieSlice(x, y, defaultNodeRadius * 1.15, sliceStartAngle, sliceEndAngle, sliceColor, 0.6);
        }
        
        // Temperature label
        std::stringstream tStr;
        tStr << std::fixed << std::setprecision(0) << node->temperature << "°";
        auto* tLabel = nodesGroup->add_child<SVG::Text>(x, y + 3, tStr.str());
        tLabel->set_attr("font-size", "7");
        tLabel->set_attr("text-anchor", "middle");
        tLabel->set_attr("fill", "#ffffff");
        tLabel->set_attr("font-weight", "bold");
        
        // Node name based on type
        std::string nodeLabel = "C";
        if (node->part == ThermalNodePartType::BOBBIN_CENTRAL_COLUMN ||
            node->part == ThermalNodePartType::BOBBIN_TOP_YOKE ||
            node->part == ThermalNodePartType::BOBBIN_BOTTOM_YOKE) {
            nodeLabel = "B";
        } else if (node->part == ThermalNodePartType::INSULATION_LAYER) {
            // Use full name for insulation layers (e.g., "IL_0_0_i")
            nodeLabel = node->name;
        }
        auto* nameLabel = nodesGroup->add_child<SVG::Text>(x, y + defaultNodeRadius + 10, nodeLabel);
        nameLabel->set_attr("font-size", "6");
        nameLabel->set_attr("text-anchor", "middle");
        nameLabel->set_attr("fill", "#666666");
    }
    
    // Draw insulation layer nodes with black dots at limit coordinates
    for (const auto* node : coreNodesForQuadrants) {
        if (node->part != ThermalNodePartType::INSULATION_LAYER) continue;
        
        auto [x, y] = mapNodeToSvg(*node);
        
        // Draw black dots at quadrant limit coordinates
        for (const auto& quadrant : node->quadrants) {
            if (quadrant.face == ThermalNodeFace::NONE) continue;
            
            // Calculate limit position in SVG coordinates
            double centerX = node->physicalCoordinates[0];
            double centerY = node->physicalCoordinates[1];
            double offsetX = quadrant.limitCoordinates[0] - centerX;
            double offsetY = quadrant.limitCoordinates[1] - centerY;
            
            // Scale to match visual node size (same as turn nodes)
            double limitX = margin + (centerX + offsetX / 2.0 - minX) * scaleFactor;
            double limitY = margin + (centerY + offsetY / 2.0 - minY) * scaleFactor;
            
            // Draw small black dot
            auto* limitDot = nodesGroup->add_child<SVG::Circle>(limitX, limitY, 2.0);
            limitDot->set_attr("fill", "#000000");
            limitDot->set_attr("stroke", "none");
        }
    }
    
    // Helper function to draw a triangle quadrant
    auto drawTriangleQuadrant = [&](double cx, double cy, double width, double height, 
                                     double angleRad, ThermalNodeFace face,
                                     const std::string& color, double opacity) {
        // For rectangular wire, draw pie-slice style quadrants
        // Each quadrant is a triangle from center to two adjacent corners
        double halfW = width / 2.0;
        double halfH = height / 2.0;
        
        // Local corner coordinates (before rotation):
        // top-left:    (-halfW, -halfH)
        // top-right:   (+halfW, -halfH)
        // bottom-right:(+halfW, +halfH)
        // bottom-left: (-halfW, +halfH)
        
        auto rotate = [&](double lx, double ly) -> std::pair<double, double> {
            double rx = cx + lx * std::cos(angleRad) - ly * std::sin(angleRad);
            double ry = cy + lx * std::sin(angleRad) + ly * std::cos(angleRad);
            return {rx, ry};
        };
        
        // Calculate all four corners in screen space
        auto tl = rotate(-halfW, -halfH);  // top-left
        auto tr = rotate(+halfW, -halfH);  // top-right
        auto br = rotate(+halfW, +halfH);  // bottom-right
        auto bl = rotate(-halfW, +halfH);  // bottom-left
        
        // For toroidal geometry, "inner" means toward the center of the toroid
        // and "outer" means away from the center. 
        // The node's angle (angleRad) points from the toroid center to the node.
        // So RADIAL_INNER is toward angleRad + PI, RADIAL_OUTER is toward angleRad.
        // We determine inner vs outer by comparing each corner's angle to the node's angle.
        
        auto cornerAngle = [&](const std::pair<double, double>& p) {
            return std::atan2(p.second - cy, p.first - cx);  // Angle from node center to corner
        };
        
        // Normalize angle to [0, 2*PI)
        auto normalize = [&](double a) {
            while (a < 0) a += 2 * M_PI;
            while (a >= 2 * M_PI) a -= 2 * M_PI;
            return a;
        };
        
        double nodeAngleNorm = normalize(angleRad);
        double innerDirection = normalize(nodeAngleNorm + M_PI);  // Toward center
        double outerDirection = normalize(nodeAngleNorm);          // Away from center
        
        // Calculate angular distance from inner/outer direction for each corner
        auto angularDist = [&](double a, double ref) {
            double d = std::abs(normalize(a - ref));
            if (d > M_PI) d = 2 * M_PI - d;
            return d;
        };
        
        struct Corner { std::pair<double, double> p; double angle; double innerDist; double outerDist; const char* name; };
        Corner corners[4] = {
            {tl, cornerAngle(tl), angularDist(cornerAngle(tl), innerDirection), angularDist(cornerAngle(tl), outerDirection), "tl"},
            {tr, cornerAngle(tr), angularDist(cornerAngle(tr), innerDirection), angularDist(cornerAngle(tr), outerDirection), "tr"},
            {br, cornerAngle(br), angularDist(cornerAngle(br), innerDirection), angularDist(cornerAngle(br), outerDirection), "br"},
            {bl, cornerAngle(bl), angularDist(cornerAngle(bl), innerDirection), angularDist(cornerAngle(bl), outerDirection), "bl"}
        };
        
        // Sort by how close each corner is to the inner direction (toward center)
        std::sort(corners, corners + 4, [](const Corner& a, const Corner& b) { return a.innerDist < b.innerDist; });
        
        // corners[0] and corners[1] are inner (closest to inner direction = toward center)
        // corners[2] and corners[3] are outer (farthest from inner direction = away from center)
        auto& inner1 = corners[0];
        auto& inner2 = corners[1];
        auto& outer1 = corners[2];
        auto& outer2 = corners[3];
        
        double x1, y1, x2, y2, x3, y3;
        
        switch (face) {
            case ThermalNodeFace::RADIAL_INNER:
                // RADIAL_INNER: the face pointing toward origin
                // Triangle: center -> inner corners
                x1 = cx; y1 = cy;
                x2 = inner1.p.first; y2 = inner1.p.second;
                x3 = inner2.p.first; y3 = inner2.p.second;
                break;
            case ThermalNodeFace::RADIAL_OUTER:
                // RADIAL_OUTER: the face pointing away from origin
                // Triangle: center -> outer corners
                x1 = cx; y1 = cy;
                x2 = outer1.p.first; y2 = outer1.p.second;
                x3 = outer2.p.first; y3 = outer2.p.second;
                break;
            case ThermalNodeFace::TANGENTIAL_LEFT:
            case ThermalNodeFace::TANGENTIAL_RIGHT:
                // For concentric cores (angleRad ≈ 0): use Y coordinates to determine top/bottom
                // For toroidal cores: use inner/outer corner sorting
                if (std::abs(angleRad) < 0.1 || std::abs(angleRad - M_PI) < 0.1 || 
                    std::abs(angleRad - 2*M_PI) < 0.1) {
                    // Concentric core: TANGENTIAL_LEFT = top, TANGENTIAL_RIGHT = bottom
                    // In SVG: top = smaller Y, bottom = larger Y
                    // tl and tr have smaller Y (top), bl and br have larger Y (bottom)
                    if (face == ThermalNodeFace::TANGENTIAL_LEFT) {
                        // Top face: use tl and tr
                        x1 = cx; y1 = cy;
                        x2 = tl.first; y2 = tl.second;
                        x3 = tr.first; y3 = tr.second;
                    } else {
                        // Bottom face: use bl and br
                        x1 = cx; y1 = cy;
                        x2 = bl.first; y2 = bl.second;
                        x3 = br.first; y3 = br.second;
                    }
                } else {
                    // Toroidal core: use inner/outer corner sorting
                    double d1 = std::sqrt(std::pow(inner1.p.first - outer1.p.first, 2) + 
                                          std::pow(inner1.p.second - outer1.p.second, 2));
                    double d2 = std::sqrt(std::pow(inner1.p.first - outer2.p.first, 2) + 
                                          std::pow(inner1.p.second - outer2.p.second, 2));
                    
                    auto& leftInner = inner1;
                    auto& rightInner = inner2;
                    auto& leftOuter = (d1 < d2) ? outer1 : outer2;
                    auto& rightOuter = (d1 < d2) ? outer2 : outer1;
                    
                    if (face == ThermalNodeFace::TANGENTIAL_LEFT) {
                        x1 = cx; y1 = cy;
                        x2 = leftInner.p.first; y2 = leftInner.p.second;
                        x3 = leftOuter.p.first; y3 = leftOuter.p.second;
                    } else {
                        x1 = cx; y1 = cy;
                        x2 = rightInner.p.first; y2 = rightInner.p.second;
                        x3 = rightOuter.p.first; y3 = rightOuter.p.second;
                    }
                }
                break;
            default:
                return;
        }
        
        // Create triangle path
        std::stringstream pathD;
        pathD << "M " << x1 << " " << y1 << " L " << x2 << " " << y2 << " L " << x3 << " " << y3 << " Z";
        
        auto* triangle = nodesGroup->add_child<SVG::Path>();
        triangle->set_attr("d", pathD.str());
        triangle->set_attr("fill", color);
        triangle->set_attr("opacity", std::to_string(opacity));
        triangle->set_attr("stroke", "#222222");
        triangle->set_attr("stroke-width", "0.5");
    };
    
    // Draw coil nodes
    for (auto& [wIdx, turnNodes] : coilNodesByWinding) {
        for (size_t nodeListIdx = 0; nodeListIdx < turnNodes.size(); nodeListIdx++) {
            const auto* node = turnNodes[nodeListIdx];
            
            // Find node index in full nodes array
            size_t nodeIdx = 0;
            for (size_t i = 0; i < nodes.size(); i++) {
                if (&nodes[i] == node) {
                    nodeIdx = i;
                    break;
                }
            }
            
            auto [x, y] = mapNodeToSvg(*node);
            
            // Calculate angular position for polar orientation
            // For concentric cores: turns are NOT angled (0 degrees), radial is purely horizontal
            double nodeAngle = 0.0;
            if (!isConcentricCore) {
                // Toroidal: calculate angle from position
                nodeAngle = std::atan2(node->physicalCoordinates.size() >= 2 ? node->physicalCoordinates[1] : 0,
                                      node->physicalCoordinates.size() >= 1 ? node->physicalCoordinates[0] : 1);
            }
            
            // Get node dimensions - scale all turns by scaleFactor/2, no clamping
            double nodeWidth = node->dimensions.width * scaleFactor / 2.0;
            double nodeHeight = node->dimensions.height * scaleFactor / 2.0;
            
            // Draw node shape based on cross-sectional shape
            if (node->crossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
                // Draw rotated rectangle for rectangular wire using a path
                // Calculate corner points
                double halfW = nodeWidth / 2.0;
                double halfH = nodeHeight / 2.0;
                
                // Corners before rotation
                double corners[4][2] = {
                    {x - halfW, y - halfH},  // top-left
                    {x + halfW, y - halfH},  // top-right
                    {x + halfW, y + halfH},  // bottom-right
                    {x - halfW, y + halfH}   // bottom-left
                };
                
                // Rotate corners by node angle
                std::string pathD = "M ";
                for (int i = 0; i < 4; i++) {
                    double dx = corners[i][0] - x;
                    double dy = corners[i][1] - y;
                    double rx = x + dx * std::cos(nodeAngle) - dy * std::sin(nodeAngle);
                    double ry = y + dx * std::sin(nodeAngle) + dy * std::cos(nodeAngle);
                    pathD += std::to_string(rx) + " " + std::to_string(ry);
                    if (i < 3) pathD += " L ";
                }
                pathD += " Z";
                
                auto* rectPath = nodesGroup->add_child<SVG::Path>();
                rectPath->set_attr("d", pathD);
                rectPath->set_attr("fill", "#ffffff");
                rectPath->set_attr("stroke", "#333333");
                rectPath->set_attr("stroke-width", "1.5");
            } else {
                // Draw circle for round wire - use average of width and height as diameter
                double diameter = (nodeWidth + nodeHeight) / 2.0;
                auto* circle = nodesGroup->add_child<SVG::Circle>(x, y, diameter / 2.0);
                circle->set_attr("fill", "#ffffff");
                circle->set_attr("stroke", "#333333");
                circle->set_attr("stroke-width", "1.5");
            }
            
            // Draw 4 triangle quadrants
            std::array<ThermalNodeFace, 4> quadrants = {
                ThermalNodeFace::RADIAL_OUTER,
                ThermalNodeFace::RADIAL_INNER,
                ThermalNodeFace::TANGENTIAL_LEFT,
                ThermalNodeFace::TANGENTIAL_RIGHT
            };
            
            for (int q = 0; q < 4; q++) {
                // Determine color based on connections for this specific quadrant
                std::string triColor;
                if (hasConductionConnectionOnFace(nodeIdx, quadrants[q], nodeAngle)) {
                    triColor = COLOR_CONDUCTION;  // Orange: connected via conduction (priority)
                } else if (hasConvectionToAmbientOnFace(nodeIdx, quadrants[q], nodeAngle)) {
                    triColor = COLOR_AMBIENT;  // Green: exposed to air
                } else {
                    // No explicit connection - not exposed to air
                    triColor = COLOR_NO_CONNECTION;  // Brown: not connected
                }
                
                // Draw quadrant shape based on cross-sectional shape
                if (node->crossSectionalShape == TurnCrossSectionalShape::ROUND) {
                    // For round wires: use pie slices (circular sectors)
                    double sliceStartAngle = 0;
                    if (isConcentricCore) {
                        // Concentric: fixed orientation (0 degrees, radial is horizontal)
                        switch (quadrants[q]) {
                            case ThermalNodeFace::RADIAL_OUTER:  // Right (+X)
                                sliceStartAngle = -45;
                                break;
                            case ThermalNodeFace::RADIAL_INNER:  // Left (-X)
                                sliceStartAngle = 135;
                                break;
                            case ThermalNodeFace::TANGENTIAL_LEFT:  // Up (+Y)
                                sliceStartAngle = 45;
                                break;
                            case ThermalNodeFace::TANGENTIAL_RIGHT:  // Down (-Y)
                                sliceStartAngle = 225;
                                break;
                            default:
                                sliceStartAngle = q * 90 - 45;
                                break;
                        }
                    } else {
                        // Toroidal: rotate based on position
                        switch (quadrants[q]) {
                            case ThermalNodeFace::RADIAL_OUTER:
                                sliceStartAngle = nodeAngle * 180.0 / M_PI - 45;
                                break;
                            case ThermalNodeFace::RADIAL_INNER:
                                sliceStartAngle = (nodeAngle + M_PI) * 180.0 / M_PI - 45;
                                break;
                            case ThermalNodeFace::TANGENTIAL_LEFT:
                                sliceStartAngle = (nodeAngle + M_PI/2) * 180.0 / M_PI - 45;
                                break;
                            case ThermalNodeFace::TANGENTIAL_RIGHT:
                                sliceStartAngle = (nodeAngle - M_PI/2) * 180.0 / M_PI - 45;
                                break;
                            default:
                                sliceStartAngle = q * 90 - 45;
                                break;
                        }
                    }
                    double sliceEndAngle = sliceStartAngle + 90;
                    double avgSize = (nodeWidth + nodeHeight) / 2.0;
                    drawPieSlice(x, y, avgSize * 0.6, sliceStartAngle, sliceEndAngle, triColor, 0.6);
                } else {
                    // For rectangular wires: use triangles
                    drawTriangleQuadrant(x, y, nodeWidth, nodeHeight, nodeAngle, quadrants[q], triColor, 0.6);
                }
            }
            
            // Draw limit coordinates as small black dots for turn quadrants only
            // Draw limit dots (black dots) for TURN and INSULATION_LAYER nodes
            if (node->part == ThermalNodePartType::TURN || 
                node->part == ThermalNodePartType::INSULATION_LAYER) {
                for (const auto& quadrant : node->quadrants) {
                    if (quadrant.face == ThermalNodeFace::NONE) continue;
                    
                    double limitX, limitY;
                    
                    if (isConcentricCore) {
                        // For concentric: use fixed horizontal/vertical offsets in SVG coordinates
                        // For rectangular wires (planar), use actual edge positions
                        // For round wires, use radius
                        double offsetX = 0, offsetY = 0;

                        if (node->crossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
                            // Rectangular wire: use actual dimensions (dots at edges)
                            switch (quadrant.face) {
                                case ThermalNodeFace::RADIAL_INNER:   // Left (-X)
                                    offsetX = -nodeWidth / 2.0; offsetY = 0; break;
                                case ThermalNodeFace::RADIAL_OUTER:   // Right (+X)
                                    offsetX = nodeWidth / 2.0; offsetY = 0; break;
                                case ThermalNodeFace::TANGENTIAL_LEFT: // Top (+Y in physics = -Y in SVG)
                                    offsetX = 0; offsetY = -nodeHeight / 2.0; break;
                                case ThermalNodeFace::TANGENTIAL_RIGHT: // Bottom (-Y in physics = +Y in SVG)
                                    offsetX = 0; offsetY = nodeHeight / 2.0; break;
                                default: break;
                            }
                        } else {
                            // Round wire: use averaged radius
                            double radius = (nodeWidth + nodeHeight) / 4.0;
                            switch (quadrant.face) {
                                case ThermalNodeFace::RADIAL_INNER:   // Left (-X)
                                    offsetX = -radius; offsetY = 0; break;
                                case ThermalNodeFace::RADIAL_OUTER:   // Right (+X)
                                    offsetX = radius; offsetY = 0; break;
                                case ThermalNodeFace::TANGENTIAL_LEFT: // Top (+Y in physics = -Y in SVG)
                                    offsetX = 0; offsetY = -radius; break;
                                case ThermalNodeFace::TANGENTIAL_RIGHT: // Bottom (-Y in physics = +Y in SVG)
                                    offsetX = 0; offsetY = radius; break;
                                default: break;
                            }
                        }
                        // Position directly in SVG coordinates (center is already scaled in x, y)
                        limitX = x + offsetX;
                        limitY = y + offsetY;
                    } else {
                        // For toroidal: use actual 3D limit coordinates (angled)
                        double centerX = node->physicalCoordinates[0];
                        double centerY = node->physicalCoordinates[1];
                        double offsetX = quadrant.limitCoordinates[0] - centerX;
                        double offsetY = quadrant.limitCoordinates[1] - centerY;
                        // Scale limit offset by 1/2 to match the visual node size
                        limitX = margin + (centerX + offsetX / 2.0 - minX) * scaleFactor;
                        limitY = margin + (centerY + offsetY / 2.0 - minY) * scaleFactor;
                    }
                    
                    // Draw small black dot
                    auto* limitDot = nodesGroup->add_child<SVG::Circle>(limitX, limitY, 2.0);
                    limitDot->set_attr("fill", "#000000");
                    limitDot->set_attr("stroke", "none");
                }
            }
            
            // Draw quadrant labels (initials) if enabled
            if (showQuadrantLabels && node->part == ThermalNodePartType::TURN) {
                auto getQuadrantInitial = [isConcentricCore](ThermalNodeFace face) -> std::string {
                    if (isConcentricCore) {
                        // Concentric cores: T, B, L, R mapping
                        switch (face) {
                            case ThermalNodeFace::RADIAL_INNER: return "L";      // LEFT (-X)
                            case ThermalNodeFace::RADIAL_OUTER: return "R";      // RIGHT (+X)
                            case ThermalNodeFace::TANGENTIAL_LEFT: return "T";   // TOP (+Y)
                            case ThermalNodeFace::TANGENTIAL_RIGHT: return "B";  // BOTTOM (-Y)
                            default: return "";
                        }
                    } else {
                        // Toroidal cores: RI, RO, TL, TR (note: TL and TR are swapped)
                        switch (face) {
                            case ThermalNodeFace::RADIAL_INNER: return "RI";
                            case ThermalNodeFace::RADIAL_OUTER: return "RO";
                            case ThermalNodeFace::TANGENTIAL_LEFT: return "TR";  // Swapped
                            case ThermalNodeFace::TANGENTIAL_RIGHT: return "TL"; // Swapped
                            default: return "";
                        }
                    }
                };
                
                for (const auto& quadrant : node->quadrants) {
                    if (quadrant.face == ThermalNodeFace::NONE) continue;
                    
                    std::string initial = getQuadrantInitial(quadrant.face);
                    if (initial.empty()) continue;
                    
                    double labelX, labelY;
                    
                    if (isConcentricCore) {
                        // For concentric: use fixed horizontal/vertical offsets in SVG coordinates
                        // Note: SVG Y increases downward, opposite of physics
                        double radius = (nodeWidth + nodeHeight) / 4.0;  // Node radius in SVG coordinates
                        double offsetX = 0, offsetY = 0;
                        switch (quadrant.face) {
                            case ThermalNodeFace::RADIAL_INNER:   // Left (-X)
                                offsetX = -radius; offsetY = 0; break;
                            case ThermalNodeFace::RADIAL_OUTER:   // Right (+X)
                                offsetX = radius; offsetY = 0; break;
                            case ThermalNodeFace::TANGENTIAL_LEFT: // Top (+Y in physics = -Y in SVG)
                                offsetX = 0; offsetY = -radius; break;
                            case ThermalNodeFace::TANGENTIAL_RIGHT: // Bottom (-Y in physics = +Y in SVG)
                                offsetX = 0; offsetY = radius; break;
                            default: break;
                        }
                        // Position label at 50% of the way from center to limit (inside the quadrant)
                        labelX = x + offsetX * 0.5;
                        labelY = y + offsetY * 0.5;
                    } else {
                        // For toroidal: use actual 3D limit coordinates (angled)
                        double centerX = node->physicalCoordinates[0];
                        double centerY = node->physicalCoordinates[1];
                        double offsetX = quadrant.limitCoordinates[0] - centerX;
                        double offsetY = quadrant.limitCoordinates[1] - centerY;
                        // Position label at 25% of the way from center to limit (inside the quadrant)
                        labelX = margin + (centerX + offsetX * 0.25 - minX) * scaleFactor;
                        labelY = margin + (centerY + offsetY * 0.25 - minY) * scaleFactor;
                    }
                    
                    // Draw quadrant initial
                    auto* qLabel = nodesGroup->add_child<SVG::Text>(labelX, labelY + 2, initial);
                    qLabel->set_attr("font-size", "4");
                    qLabel->set_attr("text-anchor", "middle");
                    qLabel->set_attr("fill", "#ffffff");
                    qLabel->set_attr("font-weight", "bold");
                }
            }
            
            // Temperature label
            std::stringstream tStr;
            tStr << std::fixed << std::setprecision(0) << node->temperature << "°";
            auto* tLabel = nodesGroup->add_child<SVG::Text>(x, y + 3, tStr.str());
            tLabel->set_attr("font-size", "7");
            tLabel->set_attr("text-anchor", "middle");
            tLabel->set_attr("fill", "#ffffff");
            tLabel->set_attr("font-weight", "bold");
            
            // Turn number below with inner/outer suffix for toroidal cores
            // Format: W_i_T_j_i or W_i_T_j_o (winding index, turn index, inner/outer)
            size_t windingNum = 0;
            size_t turnNum = 0;
            if (node->windingIndex.has_value()) {
                windingNum = node->windingIndex.value();
            }
            if (node->turnIndex.has_value()) {
                turnNum = node->turnIndex.value();
            }
            
            // Check if this is an inner or outer turn node (toroidal)
            std::string turnLabel = "W" + std::to_string(windingNum) + "_T" + std::to_string(turnNum);
            std::string nameLower = node->name;
            std::transform(nameLower.begin(), nameLower.end(), nameLower.begin(), ::tolower);
            if (nameLower.find("_inner") != std::string::npos) {
                turnLabel += "_i";
            } else if (nameLower.find("_outer") != std::string::npos) {
                turnLabel += "_o";
            }
            
            double labelOffset = std::max(nodeHeight, nodeWidth) / 2.0 + 12;
            auto* numLabel = nodesGroup->add_child<SVG::Text>(x, y + labelOffset, turnLabel);
            numLabel->set_attr("font-size", "6");
            numLabel->set_attr("text-anchor", "middle");
            numLabel->set_attr("fill", "#666666");
        }
    }
    
    // Legend removed for cleaner visualization
    
    return export_svg();
}

} // namespace OpenMagnetics
