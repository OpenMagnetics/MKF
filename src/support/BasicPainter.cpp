#include "support/Painter.h"
#include <cfloat>

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
        throw std::runtime_error("Wire is missing outerDiameter");
    }

    double outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
    double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
    double insulationThickness = (outerDiameter - conductingDiameter) / 2;
    auto coating = wire.resolve_coating();
    size_t numberLines = 0;
    double strokeWidth = 0;
    double lineRadiusIncrease = 0;
    double currentLineDiameter = conductingDiameter;
    auto coatingColor = settings->get_painter_color_insulation();
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
        // opacity = 0.25;
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
            throw std::runtime_error("Wire is missing conducting diameter");
        }
        std::string colorClass;
        if (_fieldPainted) {
            colorClass = "copper_translucent";
        }
        else {
            colorClass = "copper";
        }
        paint_circle(xCoordinate, yCoordinate, conductingDiameter / 2, colorClass, shapes, 360, 0, {0, 0}, label);
    }

    // Paint layer separation lines
    {
        std::string cssClassName = generate_random_string();
        _root.style("." + cssClassName).set_attr("opacity", opacity).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
        
        for (size_t i = 0; i < numberLines; ++i) {
            paint_circle(xCoordinate, yCoordinate, currentLineDiameter / 2, cssClassName, shapes, 360, 0, {0, 0}, label);
            currentLineDiameter += lineRadiusIncrease;
        }
    }
}

void BasicPainter::paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire, std::optional<std::string> label) {
    if (!wire.get_outer_diameter()) {
        throw std::runtime_error("Wire is missing outerDiameter");
    }
    bool simpleMode = settings->get_painter_simple_litz();
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
            throw std::runtime_error("Number layers missing in litz served");
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
        throw std::runtime_error("Coating type not implemented for Litz yet");
    }

    double insulationThickness = (outerDiameter - conductingDiameter) / 2;
    if (insulationThickness < 0) {
        throw std::runtime_error("Insulation thickness cannot be negative");
    }

    size_t numberLines = 0;
    double strokeWidth = 0;
    double lineRadiusIncrease = 0;
    double currentLineDiameter = conductingDiameter;
    auto coatingColor = settings->get_painter_color_insulation();
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
        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
        
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

        auto coordinateFilePath = settings->get_painter_cci_coordinates_path();
        coordinateFilePath = coordinateFilePath.append("cci" + std::to_string(numberConductors) + ".txt");

        std::ifstream in(coordinateFilePath);
        std::vector<std::pair<double, double>> coordinates;
        bool advancedMode = settings->get_painter_advanced_litz();

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
            throw std::runtime_error("Wire is missing both outerWidth and conductingWidth");
        }
        outerWidth = resolve_dimensional_values(wire.get_conducting_width().value());
    }
    if (wire.get_outer_height()) {
        outerHeight = resolve_dimensional_values(wire.get_outer_height().value());
    }
    else {
        if (!wire.get_conducting_height()) {
            throw std::runtime_error("Wire is missing both outerHeight and conductingHeight");
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

    std::string coatingColor = settings->get_painter_color_insulation();
    if (coating) {
        InsulationWireCoatingType insulationWireCoatingType = coating->get_type().value();

        switch (insulationWireCoatingType) {
            case InsulationWireCoatingType::BARE:
                break;
            case InsulationWireCoatingType::ENAMELLED:
                if (!coating->get_grade()) {
                    throw std::runtime_error("Enamelled wire missing grade");
                }
                numberLines = coating->get_grade().value() + 1;
                lineWidthIncrease = insulationThicknessInWidth / coating->get_grade().value() * 2;
                lineHeightIncrease = insulationThicknessInHeight / coating->get_grade().value() * 2;
                coatingColor = settings->get_painter_color_enamel();
                break;
            default:
                throw std::runtime_error("Coating type plot not implemented yet");
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
        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
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
        throw std::runtime_error("Winding sections not created");
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
        throw std::runtime_error("Winding sections not created");
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
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_copper()), std::regex("0x"), "#"));
            }
            else {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_insulation()), std::regex("0x"), "#"));
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
        throw std::runtime_error("Winding layers not created");
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
        throw std::runtime_error("Core has not being processed");
    }

    if (!winding.get_layers_description()) {
        throw std::runtime_error("Winding layers not created");
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
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_copper()), std::regex("0x"), "#"));
            }
            else {
                _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_insulation()), std::regex("0x"), "#"));
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
        throw std::runtime_error("Winding turns not created");
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
        auto group = coil.get_groups_description().value()[0]; // TODO: take into account more groups
        paint_rectangle(group.get_coordinates()[0], group.get_coordinates()[1], group.get_dimensions()[0], group.get_dimensions()[1], "fr4", shapes);
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
                        throw std::runtime_error("Wire is missing both outerWidth and conductingWidth");
                    }
                    outerWidth = resolve_dimensional_values(wire.get_conducting_width().value());
                }
                if (wire.get_outer_height()) {
                    outerHeight = resolve_dimensional_values(wire.get_outer_height().value());
                }
                else {
                    if (!wire.get_conducting_height()) {
                        throw std::runtime_error("Wire is missing both outerHeight and conductingHeight");
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
        throw std::runtime_error("Winding turns not created");
    }

    auto turns = winding.get_turns_description().value();

    for (size_t i = 0; i < turns.size(); ++i){

        if (!turns[i].get_coordinate_system()) {
            throw std::runtime_error("Turn is missing coordinate system");
        }
        if (!turns[i].get_rotation()) {
            throw std::runtime_error("Turn is missing rotation");
        }
        if (turns[i].get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            throw std::runtime_error("Painter: Turn coordinates are not in cartesian");
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
            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_insulation()), std::regex("0x"), "#"));
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
        throw std::runtime_error("Bobbin has not being processed");
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
    _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_ferrite()), std::regex("0x"), "#"));
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
    bool drawSpacer = settings->get_painter_draw_spacer();
    auto sections = magnetic.get_coil().get_sections_description().value();

    if (sections.size() == 1) {
        return;
    }

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    if (!magnetic.get_coil().get_sections_description()) {
        throw std::runtime_error("Winding sections not created");
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
                            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_margin()), std::regex("0x"), "#"));
                            paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2), {0, 0});
                        }
                    }

                    if (margins[1] > 0) {
                        double strokeWidth = sections[i].get_dimensions()[0];
                        double circleDiameter = (windingWindowRadialHeight - sections[i].get_coordinates()[0]) * 2;

                        double angle = wound_distance_to_angle(margins[1], circleDiameter / 2 - strokeWidth / 2);
                        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                            std::string cssClassName = generate_random_string();
                            _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_margin()), std::regex("0x"), "#"));
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
                        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_margin()), std::regex("0x"), "#"));
                        paint_circle(0, 0, circleDiameter / 2, cssClassName, nullptr, angle, -(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2), {0, 0});
                    }
                }
                if (margins[1] > 0) {

                    double strokeWidth = margins[1];
                    double circleDiameter = (windingWindowRadialHeight - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2) - margins[1] / 2) * 2;

                    double angle = wound_distance_to_angle(sections[i].get_dimensions()[1], circleDiameter / 2 + strokeWidth / 2);
                    if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                        std::string cssClassName = generate_random_string();
                        _root.style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_margin()), std::regex("0x"), "#"));
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
            throw std::runtime_error("Unknown error");
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
        throw std::runtime_error("Core has not being processed");
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
    auto minColor = hex_to_uint(minimumColor);
    auto maxColor = hex_to_uint(maximumColor);
    
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
    // settings->set_painter_number_points_x(4);
    // settings->set_painter_number_points_y(4);
    auto mode = settings->get_painter_mode();
    bool logarithmicScale = settings->get_painter_logarithmic_scale();


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

    if (!settings->get_painter_maximum_value_colorbar()) {
        maximumModule = percentile95Value;
    }
    if (!settings->get_painter_minimum_value_colorbar()) {
        minimumModule = percentile05Value;
    }

    if (settings->get_painter_maximum_value_colorbar()) {
        maximumModule = settings->get_painter_maximum_value_colorbar().value();
    }
    if (settings->get_painter_minimum_value_colorbar()) {
        minimumModule = settings->get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }

    auto magneticFieldMinimumColor = settings->get_painter_color_magnetic_field_minimum();
    auto magneticFieldMaximumColor = settings->get_painter_color_magnetic_field_maximum();
    
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
        stream << std::fixed << std::setprecision(1) << value;
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
    // settings->set_painter_number_points_x(4);
    // settings->set_painter_number_points_y(4);
    auto mode = settings->get_painter_mode();
    bool logarithmicScale = settings->get_painter_logarithmic_scale();


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

    if (!settings->get_painter_maximum_value_colorbar()) {
        maximumModule = percentile95Value;
    }
    if (!settings->get_painter_minimum_value_colorbar()) {
        minimumModule = percentile05Value;
    }

    if (settings->get_painter_maximum_value_colorbar()) {
        maximumModule = settings->get_painter_maximum_value_colorbar().value();
    }
    if (settings->get_painter_minimum_value_colorbar()) {
        minimumModule = settings->get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }

    auto magneticFieldMinimumColor = settings->get_painter_color_magnetic_field_minimum();
    auto magneticFieldMaximumColor = settings->get_painter_color_magnetic_field_maximum();

    auto windingWindow = magnetic.get_mutable_core().get_winding_window();
    json mierda;
    to_json(mierda, windingWindow);

    if (windingWindow.get_width()) {

        std::string cssClassName = generate_random_string();

        auto color = get_color(minimumModule, maximumModule, magneticFieldMinimumColor, magneticFieldMaximumColor, minimumModule);
        color = std::regex_replace(std::string(color), std::regex("0x"), "#");
        _root.style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", color);
        std::cout << "bool(windingWindow.get_coordinates())" << bool(windingWindow.get_coordinates()) << std::endl;
        std::cout << "bool(windingWindow.get_width())" << bool(windingWindow.get_width()) << std::endl;
        std::cout << "bool(windingWindow.get_height())" << bool(windingWindow.get_height()) << std::endl;
        paint_rectangle(windingWindow.get_coordinates().value()[0] + windingWindow.get_width().value() / 2, windingWindow.get_coordinates().value()[1], windingWindow.get_width().value(), windingWindow.get_height().value(), cssClassName);
    }
    else {
        throw std::runtime_error("Not implemented yet");
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
        stream << std::fixed << std::setprecision(1) << value;
        std::string label = stream.str() + " V/m";
        paint_field_point(datum.get_point()[0], datum.get_point()[1], pixelXDimension, pixelYDimension, color, label);
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

} // namespace OpenMagnetics