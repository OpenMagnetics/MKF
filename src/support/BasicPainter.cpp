#include "support/Painter.h"

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


void BasicPainter::paint_round_wire(double xCoordinate, double yCoordinate, Wire wire) {
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

    auto shapes = _root->add_child<SVG::Group>();

    std::string cssClassName = generate_random_string();

    // Paint insulation
    {
        *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, outerDiameter / 2 * _scale);

        auto turnSvg = _root->get_children<SVG::Circle>().back();
        _root->style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", coatingColor);
        turnSvg->set_attr("class", cssClassName);

    }

    // Paint copper
    {
        if (!wire.get_conducting_diameter()) {
            throw std::runtime_error("Wire is missing conducting diameter");
        }
        *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, conductingDiameter / 2 * _scale);
        auto turnSvg = _root->get_children<SVG::Circle>().back();
        turnSvg = _root->get_children<SVG::Circle>().back();
        turnSvg->set_attr("class", "copper");
    }

    // Paint layer separation lines
    {
        std::string cssClassName = generate_random_string();
        _root->style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
        
        for (size_t i = 0; i < numberLines; ++i) {

            *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, currentLineDiameter / 2 * _scale);
            auto turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg->set_attr("class", cssClassName);
            currentLineDiameter += lineRadiusIncrease;
        }
    }
}

void BasicPainter::paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire) {
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

    auto shapes = _root->add_child<SVG::Group>();
    std::string cssClassName = generate_random_string();

    // Paint insulation
    {
        *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, outerDiameter / 2 * _scale);

        auto turnSvg = _root->get_children<SVG::Circle>().back();
        _root->style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", coatingColor);
        turnSvg->set_attr("class", cssClassName);
    }
    // Paint layer separation lines
    {
        std::string cssClassName = generate_random_string();
        _root->style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
        
        for (size_t i = 0; i < numberLines; ++i) {

            *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, currentLineDiameter / 2 * _scale);
            auto turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg->set_attr("class", cssClassName);
            currentLineDiameter += lineRadiusIncrease;
        }
    }

    if (simpleMode) {
        *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, conductingDiameter / 2 * _scale);
        auto turnSvg = _root->get_children<SVG::Circle>().back();
        turnSvg = _root->get_children<SVG::Circle>().back();
        turnSvg->set_attr("class", "copper");
    }
    else {
        // Contour
        {
            *shapes << SVG::Circle(xCoordinate * _scale, (_imageHeight / 2 - yCoordinate) * _scale, conductingDiameter / 2 * _scale);
            auto turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg = _root->get_children<SVG::Circle>().back();
            turnSvg->set_attr("class", "white");
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
                    BasicPainter::paint_round_wire(xCoordinate + internalXCoordinate, yCoordinate + internalYCoordinate, strand);
                }
                else {
                    *shapes << SVG::Circle((xCoordinate + internalXCoordinate) * _scale, (_imageHeight / 2 - yCoordinate - internalYCoordinate) * _scale, strandOuterDiameter / 2 * _scale);
                    auto turnSvg = _root->get_children<SVG::Circle>().back();
                    turnSvg = _root->get_children<SVG::Circle>().back();
                    turnSvg->set_attr("class", "copper");
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
                    BasicPainter::paint_round_wire(xCoordinate + internalXCoordinate, yCoordinate + internalYCoordinate, strand);
                }
                else {
                    *shapes << SVG::Circle((xCoordinate + internalXCoordinate) * _scale, (_imageHeight / 2 - yCoordinate - internalYCoordinate) * _scale, strandOuterDiameter / 2 * _scale);
                    auto turnSvg = _root->get_children<SVG::Circle>().back();
                    turnSvg = _root->get_children<SVG::Circle>().back();
                    turnSvg->set_attr("class", "copper");
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

void BasicPainter::paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string cssClassName) {
    std::vector<SVG::Point> turnPoints = {};
    turnPoints.push_back(SVG::Point(xCoordinate - xDimension / 2, yCoordinate + yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate + xDimension / 2, yCoordinate + yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate + xDimension / 2, yCoordinate - yDimension / 2));
    turnPoints.push_back(SVG::Point(xCoordinate - xDimension / 2, yCoordinate - yDimension / 2));

    auto shapes = _root->add_child<SVG::Group>();
    *shapes << SVG::Polygon(scale_points(turnPoints, _imageHeight, _scale));
    auto turnSvg = _root->get_children<SVG::Polygon>().back();
    turnSvg->set_attr("class", cssClassName);
}

void BasicPainter::paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle, std::vector<double> center) {
 
    if (!wire.get_outer_width()) {
        throw std::runtime_error("Wire is missing outerWidth");
    }
    if (!wire.get_outer_height()) {
        throw std::runtime_error("Wire is missing outerHeight");
    }

    double outerWidth = resolve_dimensional_values(wire.get_outer_width().value());
    double outerHeight = resolve_dimensional_values(wire.get_outer_height().value());
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
    auto shapes = _root->add_child<SVG::Group>();

    // Paint insulation

    {
        std::string cssClassName = generate_random_string();
        _root->style("." + cssClassName).set_attr("opacity", _opacity).set_attr("fill", coatingColor);
        paint_rectangle(xCoordinate, yCoordinate, outerWidth, outerHeight, cssClassName);
    }

    // Paint copper
    {                
        paint_rectangle(xCoordinate, yCoordinate, conductingWidth, conductingHeight, "copper");
    }

    // Paint layer separation lines

    {
        std::string cssClassName = generate_random_string();
        _root->style("." + cssClassName).set_attr("stroke-width", strokeWidth * _scale).set_attr("fill", "none").set_attr("stroke", std::regex_replace(std::string(settings->get_painter_color_lines()), std::regex("0x"), "#"));
        for (size_t i = 0; i < numberLines; ++i) {
            paint_rectangle(xCoordinate, yCoordinate, currentLineWidth, currentLineHeight, cssClassName);
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

    auto shapes = _root->add_child<SVG::Group>();
    for (size_t i = 0; i < sections.size(); ++i){

        std::vector<SVG::Point> sectionPoints = {};
        sectionPoints.push_back(SVG::Point(sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
        sectionPoints.push_back(SVG::Point(sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2));
        sectionPoints.push_back(SVG::Point(sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));
        sectionPoints.push_back(SVG::Point(sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2));

        *shapes << SVG::Polygon(scale_points(sectionPoints, _imageHeight, _scale));

        auto sectionSvg = _root->get_children<SVG::Polygon>().back();
        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
            _root->style(".section_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleSections[i % constants.coilPainterColorsScaleSections.size()]);
        }
        else {
            _root->style(".section_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
        }
        sectionSvg->set_attr("class", "section_" + std::to_string(i));
    }
    paint_two_piece_set_margin(magnetic);
}

void BasicPainter::paint_two_piece_set_coil_layers(Magnetic magnetic) {
    auto constants = Constants();
    Coil coil = magnetic.get_coil();
    if (!coil.get_layers_description()) {
        throw std::runtime_error("Winding layers not created");
    }

    auto layers = coil.get_layers_description().value();

    auto shapes = _root->add_child<SVG::Group>();
    for (size_t i = 0; i < layers.size(); ++i){



        std::vector<SVG::Point> layerPoints = {};
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));

        *shapes << SVG::Polygon(scale_points(layerPoints, _imageHeight, _scale));

        auto layerSvg = _root->get_children<SVG::Polygon>().back();
        if (layers[i].get_type() == ElectricalType::CONDUCTION) {
            _root->style(".layer_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleLayers[i % constants.coilPainterColorsScaleLayers.size()]);
        }
        else {
            _root->style(".layer_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
        }
        layerSvg->set_attr("class", "layer_" + std::to_string(i));
    }
    paint_two_piece_set_margin(magnetic);
}

void BasicPainter::paint_two_piece_set_coil_turns(Magnetic magnetic) {
    auto constants = Constants();
    Coil coil = magnetic.get_coil();
    auto wirePerWinding = coil.get_wires();

    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding turns not created");
    }

    auto turns = coil.get_turns_description().value();

    auto shapes = _root->add_child<SVG::Group>();
    for (size_t i = 0; i < turns.size(); ++i){

        auto windingIndex = coil.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND) {
            paint_round_wire(turns[i].get_coordinates()[0], turns[i].get_coordinates()[1], wirePerWinding[windingIndex]);
        }
        else if (wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            paint_litz_wire(turns[i].get_coordinates()[0], turns[i].get_coordinates()[1], wirePerWinding[windingIndex]);
        }
        else {
            {
                std::vector<SVG::Point> turnPoints = {};
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2));

                *shapes << SVG::Polygon(scale_points(turnPoints, _imageHeight, _scale));

                auto turnSvg = _root->get_children<SVG::Polygon>().back();
                _root->style(".turn_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleTurns[turns[i].get_parallel() % constants.coilPainterColorsScaleTurns.size()]);
                turnSvg->set_attr("class", "turn_" + std::to_string(i));
            }

            if (wire.get_conducting_width() && wire.get_conducting_height()) {
                std::vector<SVG::Point> turnPoints = {};
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2));

                *shapes << SVG::Polygon(scale_points(turnPoints, _imageHeight, _scale));

                auto turnSvg = _root->get_children<SVG::Polygon>().back();
                turnSvg->set_attr("class", "copper");
            }
            
        }
    }

    auto layers = coil.get_layers_description().value();

    for (size_t i = 0; i < layers.size(); ++i){
        if (layers[i].get_type() == ElectricalType::INSULATION) {

            std::vector<SVG::Point> layerPoints = {};
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));

            *shapes << SVG::Polygon(scale_points(layerPoints, _imageHeight, _scale));

            auto layerSvg = _root->get_children<SVG::Polygon>().back();

            _root->style(".layer_" + std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
            layerSvg->set_attr("class", "layer_" + std::to_string(i));
        }
    }
    paint_two_piece_set_margin(magnetic);
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

    auto shapes = _root->add_child<SVG::Group>();
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

    *shapes << SVG::Polygon(scale_points(bobbinPoints, _imageHeight, _scale));

    auto sectionSvg = _root->get_children<SVG::Polygon>().back();
    sectionSvg->set_attr("class", "bobbin");
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

    double lowestHeightTopCoreMainColumn;
    double lowestHeightTopCoreRightColumn;
    double highestHeightBottomCoreMainColumn;
    double highestHeightBottomCoreRightColumn;
    if (gapsInMainColumn.size() == 0) {
        lowestHeightTopCoreMainColumn = 0;
        highestHeightBottomCoreMainColumn = 0;
    }
    else {
        lowestHeightTopCoreMainColumn = gapsInMainColumn.front().get_coordinates().value()[1] + gapsInMainColumn.front().get_length() / 2;
        highestHeightBottomCoreMainColumn = gapsInMainColumn.back().get_coordinates().value()[1] - gapsInMainColumn.back().get_length() / 2;
    }
    if (gapsInRightColumn.size() == 0) {
        lowestHeightTopCoreRightColumn = 0;
        highestHeightBottomCoreRightColumn = 0;
    }
    else {
        lowestHeightTopCoreRightColumn = gapsInRightColumn.front().get_coordinates().value()[1] + gapsInRightColumn.front().get_length() / 2;
        highestHeightBottomCoreRightColumn = gapsInRightColumn.back().get_coordinates().value()[1] - gapsInRightColumn.back().get_length() / 2;
    }

    topPiecePoints.push_back(SVG::Point(0, processedDescription.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth, processedDescription.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth, lowestHeightTopCoreRightColumn));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, lowestHeightTopCoreRightColumn));
    topPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, rightColumn.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingMainColumnWidth, mainColumn.get_height() / 2));
    topPiecePoints.push_back(SVG::Point(showingMainColumnWidth, lowestHeightTopCoreMainColumn));
    topPiecePoints.push_back(SVG::Point(0, lowestHeightTopCoreMainColumn));

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

    bottomPiecePoints.push_back(SVG::Point(0, -processedDescription.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth, -processedDescription.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth, highestHeightBottomCoreRightColumn));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, highestHeightBottomCoreRightColumn));
    bottomPiecePoints.push_back(SVG::Point(showingCoreWidth - rightColumnWidth, -rightColumn.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingMainColumnWidth, -mainColumn.get_height() / 2));
    bottomPiecePoints.push_back(SVG::Point(showingMainColumnWidth, highestHeightBottomCoreMainColumn));
    bottomPiecePoints.push_back(SVG::Point(0, highestHeightBottomCoreMainColumn));

    auto shapes = _root->add_child<SVG::Group>();
    *shapes << SVG::Polygon(scale_points(topPiecePoints, _imageHeight, _scale));
    auto topPiece = _root->get_children<SVG::Polygon>().back();
    topPiece->set_attr("class", "ferrite");
    *shapes << SVG::Polygon(scale_points(bottomPiecePoints, _imageHeight, _scale));
    auto bottomPiece = _root->get_children<SVG::Polygon>().back();
    bottomPiece->set_attr("class", "ferrite");
    for (auto& chunk : gapChunks) {
        *shapes << SVG::Polygon(scale_points(chunk, _imageHeight, _scale));
        auto chunkPiece = _root->get_children<SVG::Polygon>().back();
        chunkPiece->set_attr("class", "ferrite");
    }

    _root->autoscale();
    _root->set_attr("width", _imageWidth / 2 * _scale).set_attr("height", _imageHeight * _scale);  // TODO remove
}

void BasicPainter::paint_two_piece_set_margin(Magnetic magnetic) {
    auto sections = magnetic.get_coil().get_sections_description().value();
    for (size_t i = 0; i < sections.size(); ++i){
        if (sections[i].get_margin()) {
            auto margins = sections[i].get_margin().value();
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
                auto margins = sections[i].get_margin().value();
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
    _root->autoscale();
    _root->set_attr("width", _imageWidth * _scale).set_attr("height", _imageHeight * _scale);  // TODO remove
}

void BasicPainter::paint_core(Magnetic magnetic) {
    Core core = magnetic.get_core();
    _imageHeight = core.get_processed_description()->get_height();
    _imageWidth = core.get_processed_description()->get_width();
    _scale = constants.coilPainterScale;
    CoreShape shape = core.resolve_shape();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
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
            throw std::runtime_error("T shapes does not have bobbins");
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
            throw std::runtime_error("T shapes not implemented yet");
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
            throw std::runtime_error("T shapes not implemented yet");
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
    _scale = constants.coilPainterScale;
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
            break;
        default:
            return paint_two_piece_set_coil_turns(magnetic);
            break;
    }
}

std::string BasicPainter::export_svg() {
    if (!_filepath.empty()) {
        if (!std::filesystem::exists(_filepath)) {
            std::filesystem::create_directory(_filepath);
        }
        std::ofstream outfile(_filepath.replace_filename(_filename));
        outfile << std::string(*_root);
    }
    return std::string(*_root);
}

} // namespace OpenMagnetics