#include "physical_models/WindingLosses.h"
#include "support/Painter.h"
#include <cfloat>
#include <map>
#include <algorithm>
#include <sstream>
#include <iomanip>
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

void BasicPainter::paint_temperature_field(Magnetic magnetic, const std::map<std::string, double>& nodeTemperatures, bool showColorBar) {
    set_image_size(magnetic);
    
    if (nodeTemperatures.empty()) {
        return;
    }
    
    // Find temperature range for color mapping
    double minimumTemperature = DBL_MAX;
    double maximumTemperature = -DBL_MAX;
    for (const auto& [name, temp] : nodeTemperatures) {
        if (temp < minimumTemperature) minimumTemperature = temp;
        if (temp > maximumTemperature) maximumTemperature = temp;
    }
    
    // Apply colorbar settings if provided
    if (settings.get_painter_minimum_value_colorbar().has_value()) {
        minimumTemperature = settings.get_painter_minimum_value_colorbar().value();
    }
    if (settings.get_painter_maximum_value_colorbar().has_value()) {
        maximumTemperature = settings.get_painter_maximum_value_colorbar().value();
    }
    if (minimumTemperature == maximumTemperature) {
        minimumTemperature = maximumTemperature - 1;
    }
    
    // Use blue (cold) to red (hot) color scale
    auto coldColor = settings.get_painter_color_magnetic_field_minimum();  // Typically blue
    auto hotColor = settings.get_painter_color_magnetic_field_maximum();   // Typically red
    
    Coil coil = magnetic.get_coil();
    Core core = magnetic.get_mutable_core();
    
    auto shapes = _root.add_child<SVG::Group>();
    
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
    for (size_t i = 0; i < columns.size(); ++i) {
        auto column = columns[i];
        std::string nodeName = "Core_Column_" + std::to_string(i);
        
        auto tempOpt = findNodeTemperature(nodeName);
        if (!tempOpt) {
            tempOpt = findNodeTemperature("Core");
        }
        
        if (tempOpt) {
            double temp = tempOpt.value();
            auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
            
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
            _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 0.8);
            paint_rectangle(xCoord, 0, colWidth, colHeight, cssClassName, shapes, 0, {0, 0}, label);
        }
    }
    
    // Paint top yoke
    auto topYokeTempOpt = findNodeTemperatureExact("Core_Top_Yoke");
    if (!topYokeTempOpt) {
        topYokeTempOpt = findNodeTemperature("Core");
    }
    if (topYokeTempOpt) {
        double temp = topYokeTempOpt.value();
        auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
        
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << temp;
        std::string label = "Core_Top_Yoke: " + stream.str() + " °C";
        
        // Top yoke: from y=mainColumnHeight/2 to y=coreHeight/2, full width
        double yokeHeight = (coreHeight - mainColumnHeight) / 2;
        double yokeY = mainColumnHeight / 2 + yokeHeight / 2;
        
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 0.8);
        paint_rectangle(showingCoreWidth / 2, yokeY, showingCoreWidth, yokeHeight, cssClassName, shapes, 0, {0, 0}, label);
    }
    
    // Paint bottom yoke
    auto bottomYokeTempOpt = findNodeTemperatureExact("Core_Bottom_Yoke");
    if (!bottomYokeTempOpt) {
        bottomYokeTempOpt = findNodeTemperature("Core");
    }
    if (bottomYokeTempOpt) {
        double temp = bottomYokeTempOpt.value();
        auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
        
        std::stringstream stream;
        stream << std::fixed << std::setprecision(1) << temp;
        std::string label = "Core_Bottom_Yoke: " + stream.str() + " °C";
        
        // Bottom yoke: from y=-coreHeight/2 to y=-mainColumnHeight/2, full width
        double yokeHeight = (coreHeight - mainColumnHeight) / 2;
        double yokeY = -(mainColumnHeight / 2 + yokeHeight / 2);
        
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("fill", color).set_attr("opacity", 0.8);
        paint_rectangle(showingCoreWidth / 2, yokeY, showingCoreWidth, yokeHeight, cssClassName, shapes, 0, {0, 0}, label);
    }
    
    // Paint bobbin if present - bobbin is a variant<Bobbin, string>
    auto bobbinVariant = magnetic.get_coil().get_bobbin();
    if (std::holds_alternative<Bobbin>(bobbinVariant)) {
        
        // Try to find bobbin temperature - node names are "Bobbin_Inner" or "Bobbin_Outer"
        auto bobbinTempOpt = findNodeTemperature("Bobbin");
        if (bobbinTempOpt) {
            double temp = bobbinTempOpt.value();
            auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
            
            std::stringstream stream;
            stream << std::fixed << std::setprecision(1) << temp;
            std::string label = "bobbin: " + stream.str() + " °C";
            
            auto windingWindow = core.get_winding_window();
            if (windingWindow.get_width()) {
                double xCoord = windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2;
                double yCoord = windingWindow.get_coordinates().value()[1];
                double width = windingWindow.get_width().value();
                double height = windingWindow.get_height().value();
                
                // Draw bobbin outline as a rectangle
                std::string cssClassName = generate_random_string();
                _root.style("." + cssClassName).set_attr("fill", "none").set_attr("stroke", color).set_attr("stroke-width", "3");
                paint_rectangle(xCoord, yCoord, width, height, cssClassName, shapes, 0, {0, 0}, label);
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
            
            // First try exact match
            auto it = nodeTemperatures.find(turnName);
            if (it != nodeTemperatures.end()) {
                tempOpt = it->second;
            } else {
                // Try to find by layer or section name - names are "Coil_Layer_0", "Coil", etc.
                tempOpt = findNodeTemperature("Coil_Layer");
                if (!tempOpt) {
                    tempOpt = findNodeTemperature("Coil");
                }
            }
            
            if (tempOpt) {
                double temp = tempOpt.value();
                auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
                
                std::stringstream stream;
                stream << std::fixed << std::setprecision(1) << temp;
                std::string label = turnName + ": " + stream.str() + " °C";
                
                if (turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
                    double xCoordinate = turn.get_coordinates()[0];
                    double yCoordinate = turn.get_coordinates()[1];
                    double diameter = turn.get_dimensions().value()[0];
                    std::string cssClassName = generate_random_string();
                    _root.style("." + cssClassName).set_attr("fill", color);
                    paint_circle(xCoordinate, yCoordinate, diameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
                } else {
                    if (turn.get_dimensions().value()[0] && turn.get_dimensions().value()[1]) {
                        double xCoordinate = turn.get_coordinates()[0];
                        double yCoordinate = turn.get_coordinates()[1];
                        double conductingWidth = turn.get_dimensions().value()[0];
                        double conductingHeight = turn.get_dimensions().value()[1];
                        std::string cssClassName = generate_random_string();
                        _root.style("." + cssClassName).set_attr("fill", color);
                        paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, cssClassName, shapes, 0, {0, 0}, label);
                    }
                }
            }
        }
    }
    
    // Add temperature legend/colorbar inside the winding window area (optional)
    if (showColorBar) {
        // Get winding window coordinates for positioning
        auto windingWindow = core.get_winding_window();
        double windowWidth = windingWindow.get_width() ? windingWindow.get_width().value() : _imageWidth * 0.3;
        double windowHeight = windingWindow.get_height() ? windingWindow.get_height().value() : _imageHeight * 0.6;
        auto windowCoords = windingWindow.get_coordinates();
        double windowX = windowCoords ? windowCoords.value()[0] : 0;
        double windowY = windowCoords ? windowCoords.value()[1] : 0;
        
        // Position legend in the right side of the winding window
        double legendWidth = windowWidth * 0.08;
        double legendHeight = windowHeight * 0.6;
        double legendX = windowX + windowWidth * 0.85;
        size_t numSteps = 10;
        
        for (size_t i = 0; i <= numSteps; ++i) {
            double t = static_cast<double>(i) / numSteps;
            double temp = minimumTemperature + t * (maximumTemperature - minimumTemperature);
            auto color = get_color(minimumTemperature, maximumTemperature, coldColor, hotColor, temp);
            
            double stepHeight = legendHeight / numSteps;
            // Y goes from bottom (low temp) to top (high temp)
            double stepY = windowY - legendHeight / 2 + (numSteps - i - 0.5) * stepHeight;
            
            std::string cssClassName = generate_random_string();
            _root.style("." + cssClassName).set_attr("fill", color);
            paint_rectangle(legendX, stepY, legendWidth, stepHeight, cssClassName, shapes);
            
            // Add temperature labels at key positions (min, mid, max)
            if (i == 0 || i == numSteps / 2 || i == numSteps) {
                std::stringstream stream;
                stream << std::fixed << std::setprecision(0) << temp;
                // Text position needs to be in scaled coordinates for SVG::Text
                double textX = (legendX + legendWidth * 0.7) * _scale;
                double textY = -stepY * _scale;  // Note: SVG Y is inverted
                auto* text = _root.add_child<SVG::Text>(textX, textY, stream.str() + "C");
                text->set_attr("font-size", std::to_string(windowHeight * _scale * 0.03));
                text->set_attr("fill", "#000000");
            }
        }
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
    const std::vector<ThermalNode>& nodes,
    const std::vector<ThermalResistanceElement>& resistances,
    double width,
    double height) {
    
    // Create a new SVG root for the schematic
    _root = SVG::SVG();
    _root.set_attr("width", std::to_string(static_cast<int>(width)));
    _root.set_attr("height", std::to_string(static_cast<int>(height)));
    _root.set_attr("viewBox", "0 0 " + std::to_string(static_cast<int>(width)) + " " + std::to_string(static_cast<int>(height)));
    
    // Add white background
    auto* bg = _root.add_child<SVG::Rect>(0, 0, width, height);
    bg->set_attr("fill", "#ffffff");
    
    // Add title
    auto* titleText = _root.add_child<SVG::Text>(width / 2, 25, "Thermal Equivalent Circuit");
    titleText->set_attr("font-size", "18");
    titleText->set_attr("font-weight", "bold");
    titleText->set_attr("text-anchor", "middle");
    titleText->set_attr("fill", "#333333");
    
    if (nodes.empty()) {
        return export_svg();
    }
    
    // Layout parameters
    double margin = 60;
    double nodeRadius = 25;
    double powerSourceRadius = 12;
    
    // Categorize nodes by type for layout
    std::vector<size_t> coreNodeIds;
    std::vector<size_t> coilNodeIds;
    std::vector<size_t> bobbinNodeIds;
    size_t ambientNodeId = 0;
    
    for (const auto& node : nodes) {
        switch (node.type) {
            case ThermalNodeType::CORE_CENTRAL_COLUMN:
            case ThermalNodeType::CORE_LATERAL_COLUMN:
            case ThermalNodeType::CORE_TOP_YOKE:
            case ThermalNodeType::CORE_BOTTOM_YOKE:
                coreNodeIds.push_back(node.id);
                break;
            case ThermalNodeType::COIL_SECTION:
            case ThermalNodeType::COIL_LAYER:
            case ThermalNodeType::COIL_TURN:
                coilNodeIds.push_back(node.id);
                break;
            case ThermalNodeType::BOBBIN_INNER:
            case ThermalNodeType::BOBBIN_OUTER:
                bobbinNodeIds.push_back(node.id);
                break;
            case ThermalNodeType::AMBIENT:
                ambientNodeId = node.id;
                break;
        }
    }
    
    // Calculate positions for each node
    std::map<size_t, std::pair<double, double>> nodePositions;
    
    // Limit displayed coil nodes to avoid overcrowding
    size_t maxCoilNodesToShow = 12;
    std::vector<size_t> displayedCoilNodeIds;
    if (coilNodeIds.size() <= maxCoilNodesToShow) {
        displayedCoilNodeIds = coilNodeIds;
    } else {
        // Sample evenly
        for (size_t i = 0; i < maxCoilNodesToShow; ++i) {
            size_t idx = i * coilNodeIds.size() / maxCoilNodesToShow;
            displayedCoilNodeIds.push_back(coilNodeIds[idx]);
        }
    }
    
    // Layout: 
    // - Core nodes on the left (vertical stack)
    // - Bobbin nodes in the middle-left
    // - Coil nodes in the center-right (grid layout)
    // - Ambient at the bottom center
    
    double coreX = margin + 80;
    double bobbinX = width * 0.35;
    double coilStartX = width * 0.45;
    double coilEndX = width - margin - 80;
    double topY = margin + 80;
    double bottomY = height - margin - 60;
    double ambientY = height - margin - 30;
    
    // Position core nodes vertically on the left
    if (!coreNodeIds.empty()) {
        double coreSpacing = (bottomY - topY - 100) / std::max(1.0, static_cast<double>(coreNodeIds.size() - 1));
        for (size_t i = 0; i < coreNodeIds.size(); ++i) {
            nodePositions[coreNodeIds[i]] = {coreX, topY + i * coreSpacing};
        }
    }
    
    // Position bobbin nodes
    if (!bobbinNodeIds.empty()) {
        double bobbinSpacing = (bottomY - topY - 100) / std::max(1.0, static_cast<double>(bobbinNodeIds.size() - 1));
        for (size_t i = 0; i < bobbinNodeIds.size(); ++i) {
            nodePositions[bobbinNodeIds[i]] = {bobbinX, topY + 50 + i * bobbinSpacing};
        }
    }
    
    // Position coil nodes in a grid layout
    if (!displayedCoilNodeIds.empty()) {
        size_t coilCols = std::min(static_cast<size_t>(4), displayedCoilNodeIds.size());
        size_t coilRows = (displayedCoilNodeIds.size() + coilCols - 1) / coilCols;
        
        double colSpacing = (coilEndX - coilStartX) / std::max(1.0, static_cast<double>(coilCols));
        double rowSpacing = (bottomY - topY - 100) / std::max(1.0, static_cast<double>(coilRows));
        
        for (size_t i = 0; i < displayedCoilNodeIds.size(); ++i) {
            size_t col = i % coilCols;
            size_t row = i / coilCols;
            nodePositions[displayedCoilNodeIds[i]] = {
                coilStartX + col * colSpacing + colSpacing / 2,
                topY + row * rowSpacing
            };
        }
    }
    
    // Position ambient node at bottom center
    nodePositions[ambientNodeId] = {width / 2, ambientY};
    
    // Find temperature range for color mapping
    double minTemp = std::numeric_limits<double>::max();
    double maxTemp = std::numeric_limits<double>::lowest();
    for (const auto& node : nodes) {
        if (!node.isAmbient()) {
            minTemp = std::min(minTemp, node.temperature);
            maxTemp = std::max(maxTemp, node.temperature);
        }
    }
    if (minTemp == maxTemp) {
        minTemp = maxTemp - 10;
    }
    
    // Helper to check if a node should be displayed
    auto isNodeDisplayed = [&nodePositions](size_t id) -> bool {
        return nodePositions.find(id) != nodePositions.end();
    };
    
    // Draw resistances (connections between nodes)
    auto* resistanceGroup = _root.add_child<SVG::Group>();
    resistanceGroup->set_attr("id", "resistances");
    
    for (const auto& res : resistances) {
        if (!isNodeDisplayed(res.nodeFromId) || !isNodeDisplayed(res.nodeToId)) {
            continue;  // Skip resistances to non-displayed nodes
        }
        
        auto [x1, y1] = nodePositions[res.nodeFromId];
        auto [x2, y2] = nodePositions[res.nodeToId];
        
        // Calculate direction and length
        double dx = x2 - x1;
        double dy = y2 - y1;
        double len = std::sqrt(dx * dx + dy * dy);
        
        if (len < 1e-6) continue;
        
        // Normalize direction
        double nx = dx / len;
        double ny = dy / len;
        
        // Start and end points (offset from node center)
        double startX = x1 + nx * nodeRadius;
        double startY = y1 + ny * nodeRadius;
        double endX = x2 - nx * nodeRadius;
        double endY = y2 - ny * nodeRadius;
        
        // Recalculate actual resistor length
        double actualLen = std::sqrt((endX - startX) * (endX - startX) + (endY - startY) * (endY - startY));
        
        // Draw resistor symbol using zigzag pattern
        // First, draw the connecting lines
        double resistorStart = actualLen * 0.25;
        double resistorEnd = actualLen * 0.75;
        
        // Line from start to resistor start
        auto* line1 = resistanceGroup->add_child<SVG::Line>(
            startX, startY,
            startX + nx * resistorStart, startY + ny * resistorStart);
        line1->set_attr("stroke", "#333333");
        line1->set_attr("stroke-width", "2");
        
        // Line from resistor end to node
        auto* line2 = resistanceGroup->add_child<SVG::Line>(
            startX + nx * resistorEnd, startY + ny * resistorEnd,
            endX, endY);
        line2->set_attr("stroke", "#333333");
        line2->set_attr("stroke-width", "2");
        
        // Draw zigzag resistor
        // Perpendicular direction for zigzag
        double px = -ny;
        double py = nx;
        
        double zigzagLen = resistorEnd - resistorStart;
        int numZigs = 5;
        double zigWidth = zigzagLen / (numZigs * 2);
        double zigHeight = 8;
        
        std::string pathD = "M " + std::to_string(startX + nx * resistorStart) + " " + 
                           std::to_string(startY + ny * resistorStart);
        
        for (int z = 0; z < numZigs; ++z) {
            double baseOffset = resistorStart + (z * 2 + 0.5) * zigWidth;
            double nextOffset = resistorStart + (z * 2 + 1.5) * zigWidth;
            
            // Zig up
            pathD += " L " + std::to_string(startX + nx * baseOffset + px * zigHeight) + " " +
                     std::to_string(startY + ny * baseOffset + py * zigHeight);
            // Zag down
            pathD += " L " + std::to_string(startX + nx * nextOffset - px * zigHeight) + " " +
                     std::to_string(startY + ny * nextOffset - py * zigHeight);
        }
        
        // Final point
        pathD += " L " + std::to_string(startX + nx * resistorEnd) + " " +
                 std::to_string(startY + ny * resistorEnd);
        
        auto* zigzag = resistanceGroup->add_child<SVG::Path>();
        zigzag->set_attr("d", pathD);
        zigzag->set_attr("stroke", "#333333");
        zigzag->set_attr("stroke-width", "2");
        zigzag->set_attr("fill", "none");
        
        // Add resistance value label
        double labelX = startX + nx * (resistorStart + zigzagLen / 2) + px * 20;
        double labelY = startY + ny * (resistorStart + zigzagLen / 2) + py * 20;
        
        std::stringstream resStr;
        resStr << std::fixed << std::setprecision(1) << res.resistance << " K/W";
        
        auto* resLabel = resistanceGroup->add_child<SVG::Text>(labelX, labelY, resStr.str());
        resLabel->set_attr("font-size", "9");
        resLabel->set_attr("text-anchor", "middle");
        resLabel->set_attr("fill", "#666666");
    }
    
    // Draw nodes
    auto* nodeGroup = _root.add_child<SVG::Group>();
    nodeGroup->set_attr("id", "nodes");
    
    for (const auto& node : nodes) {
        if (!isNodeDisplayed(node.id)) continue;
        
        auto [x, y] = nodePositions[node.id];
        
        if (node.isAmbient()) {
            // Draw ambient as ground symbol
            auto* groundGroup = nodeGroup->add_child<SVG::Group>();
            
            // Vertical line
            auto* gndLine = groundGroup->add_child<SVG::Line>(x, y - 15, x, y);
            gndLine->set_attr("stroke", "#333333");
            gndLine->set_attr("stroke-width", "2");
            
            // Three horizontal lines (decreasing width)
            for (int i = 0; i < 3; ++i) {
                double lineWidth = 20 - i * 5;
                auto* hLine = groundGroup->add_child<SVG::Line>(
                    x - lineWidth / 2, y + i * 5,
                    x + lineWidth / 2, y + i * 5);
                hLine->set_attr("stroke", "#333333");
                hLine->set_attr("stroke-width", "2");
            }
            
            // Label
            std::stringstream tempStr;
            tempStr << std::fixed << std::setprecision(1) << node.temperature << "°C";
            auto* tempLabel = nodeGroup->add_child<SVG::Text>(x, y + 25, "Ambient: " + tempStr.str());
            tempLabel->set_attr("font-size", "11");
            tempLabel->set_attr("text-anchor", "middle");
            tempLabel->set_attr("fill", "#333333");
        } else {
            // Draw node as circle with temperature-based color
            double tempRatio = (node.temperature - minTemp) / (maxTemp - minTemp);
            tempRatio = std::max(0.0, std::min(1.0, tempRatio));
            
            // Blue to red gradient
            int r = static_cast<int>(tempRatio * 255);
            int b = static_cast<int>((1 - tempRatio) * 255);
            int g = static_cast<int>((1 - std::abs(tempRatio - 0.5) * 2) * 128);
            
            std::stringstream colorStr;
            colorStr << "rgb(" << r << "," << g << "," << b << ")";
            
            auto* circle = nodeGroup->add_child<SVG::Circle>(x, y, nodeRadius);
            circle->set_attr("fill", colorStr.str());
            circle->set_attr("stroke", "#333333");
            circle->set_attr("stroke-width", "2");
            
            // Temperature label inside circle
            std::stringstream tempStr;
            tempStr << std::fixed << std::setprecision(0) << node.temperature << "°C";
            auto* tempLabel = nodeGroup->add_child<SVG::Text>(x, y + 4, tempStr.str());
            tempLabel->set_attr("font-size", "10");
            tempLabel->set_attr("text-anchor", "middle");
            tempLabel->set_attr("fill", "#ffffff");
            tempLabel->set_attr("font-weight", "bold");
            
            // Node name label below
            std::string shortName = node.name;
            if (shortName.length() > 15) {
                shortName = shortName.substr(0, 12) + "...";
            }
            auto* nameLabel = nodeGroup->add_child<SVG::Text>(x, y + nodeRadius + 12, shortName);
            nameLabel->set_attr("font-size", "9");
            nameLabel->set_attr("text-anchor", "middle");
            nameLabel->set_attr("fill", "#333333");
            
            // Draw power source if node has power dissipation
            if (node.powerDissipation > 0.001) {  // > 1 mW
                double psX = x + nodeRadius + 5;
                double psY = y - nodeRadius + 5;
                
                // Draw small circle with P symbol
                auto* psCircle = nodeGroup->add_child<SVG::Circle>(psX, psY, powerSourceRadius);
                psCircle->set_attr("fill", "#ffcc00");
                psCircle->set_attr("stroke", "#cc9900");
                psCircle->set_attr("stroke-width", "1.5");
                
                // P symbol
                auto* pLabel = nodeGroup->add_child<SVG::Text>(psX, psY + 4, "P");
                pLabel->set_attr("font-size", "10");
                pLabel->set_attr("text-anchor", "middle");
                pLabel->set_attr("fill", "#333333");
                pLabel->set_attr("font-weight", "bold");
                
                // Power value
                std::stringstream powerStr;
                if (node.powerDissipation >= 1.0) {
                    powerStr << std::fixed << std::setprecision(2) << node.powerDissipation << "W";
                } else {
                    powerStr << std::fixed << std::setprecision(0) << node.powerDissipation * 1000 << "mW";
                }
                auto* powerLabel = nodeGroup->add_child<SVG::Text>(psX, psY + powerSourceRadius + 10, powerStr.str());
                powerLabel->set_attr("font-size", "8");
                powerLabel->set_attr("text-anchor", "middle");
                powerLabel->set_attr("fill", "#666666");
            }
        }
    }
    
    // Add legend
    double legendX = margin;
    double legendY = height - 50;
    
    auto* legendGroup = _root.add_child<SVG::Group>();
    legendGroup->set_attr("id", "legend");
    
    // Temperature scale bar
    double scaleWidth = 150;
    double scaleHeight = 15;
    
    // Draw gradient bar
    for (int i = 0; i < 20; ++i) {
        double ratio = static_cast<double>(i) / 19.0;
        int r = static_cast<int>(ratio * 255);
        int b = static_cast<int>((1 - ratio) * 255);
        int g = static_cast<int>((1 - std::abs(ratio - 0.5) * 2) * 128);
        
        std::stringstream colorStr;
        colorStr << "rgb(" << r << "," << g << "," << b << ")";
        
        auto* rect = legendGroup->add_child<SVG::Rect>(
            legendX + i * scaleWidth / 20, legendY,
            scaleWidth / 20 + 1, scaleHeight);
        rect->set_attr("fill", colorStr.str());
    }
    
    // Scale labels
    std::stringstream minTempStr, maxTempStr;
    minTempStr << std::fixed << std::setprecision(0) << minTemp << "°C";
    maxTempStr << std::fixed << std::setprecision(0) << maxTemp << "°C";
    
    auto* minLabel = legendGroup->add_child<SVG::Text>(legendX, legendY + scaleHeight + 12, minTempStr.str());
    minLabel->set_attr("font-size", "10");
    minLabel->set_attr("text-anchor", "start");
    minLabel->set_attr("fill", "#333333");
    
    auto* maxLabel = legendGroup->add_child<SVG::Text>(legendX + scaleWidth, legendY + scaleHeight + 12, maxTempStr.str());
    maxLabel->set_attr("font-size", "10");
    maxLabel->set_attr("text-anchor", "end");
    maxLabel->set_attr("fill", "#333333");
    
    // Legend title
    auto* legendTitle = legendGroup->add_child<SVG::Text>(legendX + scaleWidth / 2, legendY - 5, "Temperature");
    legendTitle->set_attr("font-size", "10");
    legendTitle->set_attr("text-anchor", "middle");
    legendTitle->set_attr("fill", "#333333");
    
    // Node count info (if truncated)
    if (coilNodeIds.size() > maxCoilNodesToShow) {
        std::stringstream infoStr;
        infoStr << "Showing " << displayedCoilNodeIds.size() << " of " << coilNodeIds.size() << " coil nodes";
        auto* infoLabel = legendGroup->add_child<SVG::Text>(width - margin, legendY + scaleHeight + 12, infoStr.str());
        infoLabel->set_attr("font-size", "9");
        infoLabel->set_attr("text-anchor", "end");
        infoLabel->set_attr("fill", "#999999");
    }
    
    return export_svg();
}

} // namespace OpenMagnetics
