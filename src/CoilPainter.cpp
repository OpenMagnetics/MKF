#include "CoilWrapper.h"
#include "CoilPainter.h"
#include "MAS.hpp"
#include "Utils.h"
#include "json.hpp"

namespace OpenMagnetics {

std::vector<SVG::Point> scale_points(std::vector<SVG::Point> points, double imageHeight) {
    auto constants = Constants();
    std::vector<SVG::Point> scaledPoints;
    for (auto& point : points) {
        scaledPoints.push_back(SVG::Point(point.first * constants.coilPainterScale, (imageHeight / 2 - point.second) * constants.coilPainterScale));
    }
    return scaledPoints;
}

SVG::SVG* CoilPainter::paint_two_piece_set_winding_sections(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    double imageHeight = core.get_processed_description()->get_height();

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

        *shapes << SVG::Polygon(scale_points(sectionPoints, imageHeight));

        auto sectionSvg = _root->get_children<SVG::Polygon>().back();
        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
            _root->style(".section_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleSections[i % constants.coilPainterColorsScaleSections.size()]);
        }
        else {
            _root->style(".section_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
        }
        sectionSvg->set_attr("class", "section_" +  std::to_string(i));
    }

    if (!std::filesystem::exists(_filepath)) {
        std::filesystem::create_directory(_filepath);
    }
    std::ofstream outfile(_filepath.replace_filename(_filename));
    outfile << std::string(*_root);
    return _root;
}

SVG::SVG* CoilPainter::paint_two_piece_set_winding_layers(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core has not being processed");
    }
    double imageHeight = core.get_processed_description()->get_height();

    if (!winding.get_layers_description()) {
        throw std::runtime_error("Winding layers not created");
    }

    auto layers = winding.get_layers_description().value();

    auto shapes = _root->add_child<SVG::Group>();
    for (size_t i = 0; i < layers.size(); ++i){



        std::vector<SVG::Point> layerPoints = {};
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));
        layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));

        *shapes << SVG::Polygon(scale_points(layerPoints, imageHeight));

        auto layerSvg = _root->get_children<SVG::Polygon>().back();
        if (layers[i].get_type() == ElectricalType::CONDUCTION) {
            _root->style(".layer_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleLayers[i % constants.coilPainterColorsScaleLayers.size()]);
        }
        else {
            _root->style(".layer_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
        }
        layerSvg->set_attr("class", "layer_" +  std::to_string(i));
    }

    if (!std::filesystem::exists(_filepath)) {
        std::filesystem::create_directory(_filepath);
    }
    std::ofstream outfile(_filepath.replace_filename(_filename));
    outfile << std::string(*_root);
    return _root;
}

SVG::SVG* CoilPainter::paint_two_piece_set_winding_turns(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    double imageHeight = core.get_processed_description()->get_height();
    auto wirePerWinding = winding.get_wires();

    if (!winding.get_turns_description()) {
        throw std::runtime_error("Winding turns not created");
    }

    auto turns = winding.get_turns_description().value();

    auto shapes = _root->add_child<SVG::Group>();
    for (size_t i = 0; i < turns.size(); ++i){

        auto windingIndex = winding.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND) {
            *shapes << SVG::Circle(turns[i].get_coordinates()[0] * constants.coilPainterScale,
                                   (imageHeight / 2 - turns[i].get_coordinates()[1]) * constants.coilPainterScale,
                                   resolve_dimensional_values(wire.get_outer_diameter().value()) / 2 * constants.coilPainterScale);

            auto turnSvg = _root->get_children<SVG::Circle>().back();
            _root->style(".turn_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleTurns[windingIndex % constants.coilPainterColorsScaleTurns.size()]);
            // _root->style(".turn_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleTurns[turns[i].get_parallel() % constants.coilPainterColorsScaleTurns.size()]);
            turnSvg->set_attr("class", "turn_" +  std::to_string(i));

            if (wire.get_conducting_diameter()) {
                *shapes << SVG::Circle(turns[i].get_coordinates()[0] * constants.coilPainterScale,
                                       (imageHeight / 2 - turns[i].get_coordinates()[1]) * constants.coilPainterScale,
                                       resolve_dimensional_values(wire.get_conducting_diameter().value()) / 2 * constants.coilPainterScale);

                turnSvg = _root->get_children<SVG::Circle>().back();
                turnSvg->set_attr("class", "copper");
            }
        }
        else {
            {
                std::vector<SVG::Point> turnPoints = {};
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2));

                *shapes << SVG::Polygon(scale_points(turnPoints, imageHeight));

                auto turnSvg = _root->get_children<SVG::Polygon>().back();
                _root->style(".turn_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", constants.coilPainterColorsScaleTurns[turns[i].get_parallel() % constants.coilPainterColorsScaleTurns.size()]);
                turnSvg->set_attr("class", "turn_" +  std::to_string(i));
            }

            if (wire.get_conducting_width() && wire.get_conducting_height()) {
                std::vector<SVG::Point> turnPoints = {};
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2));
                turnPoints.push_back(SVG::Point(turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2));

                *shapes << SVG::Polygon(scale_points(turnPoints, imageHeight));

                auto turnSvg = _root->get_children<SVG::Polygon>().back();
                turnSvg->set_attr("class", "copper");
            }
            
        }
    }

    auto layers = winding.get_layers_description().value();

    for (size_t i = 0; i < layers.size(); ++i){
        if (layers[i].get_type() == ElectricalType::INSULATION) {

            std::vector<SVG::Point> layerPoints = {};
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));
            layerPoints.push_back(SVG::Point(layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2));

            *shapes << SVG::Polygon(scale_points(layerPoints, imageHeight));

            auto layerSvg = _root->get_children<SVG::Polygon>().back();

            _root->style(".layer_" +  std::to_string(i)).set_attr("opacity", _opacity).set_attr("fill", "#E37E00");
            layerSvg->set_attr("class", "layer_" +  std::to_string(i));
        }
    }

    if (!std::filesystem::exists(_filepath)) {
        std::filesystem::create_directory(_filepath);
    }
    std::ofstream outfile(_filepath.replace_filename(_filename));
    outfile << std::string(*_root);
    return _root;
}

SVG::SVG* CoilPainter::paint_two_piece_set_bobbin(MagneticWrapper magnetic) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin has not being processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    CoreWrapper core = magnetic.get_core();
    double imageHeight = core.get_processed_description()->get_height();

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

    *shapes << SVG::Polygon(scale_points(bobbinPoints, imageHeight));

    auto sectionSvg = _root->get_children<SVG::Polygon>().back();
    sectionSvg->set_attr("class", "bobbin");

    if (!std::filesystem::exists(_filepath)) {
        std::filesystem::create_directory(_filepath);
    }
    std::ofstream outfile(_filepath.replace_filename(_filename));
    outfile << std::string(*_root);
    return _root;
}

SVG::SVG* CoilPainter::paint_two_piece_set_core(CoreWrapper core) {
    std::vector<SVG::Point> topPiecePoints = {};
    std::vector<SVG::Point> bottomPiecePoints = {};
    std::vector<std::vector<SVG::Point>> gapChunks = {};
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto processedDescription = core.get_processed_description().value();
    auto rightColumn = core.find_closest_column_by_coordinates(std::vector<double>({processedDescription.get_width() / 2, 0, -processedDescription.get_depth() / 2}));
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto family = core.get_shape_family();
    double showingCoreWidth;
    double showingMainColumnWidth;
    switch (family) {
        case OpenMagnetics::CoreShapeFamily::U:
        case OpenMagnetics::CoreShapeFamily::UR:
            showingMainColumnWidth = mainColumn.get_width() / 2;
            showingCoreWidth = processedDescription.get_width() - mainColumn.get_width() / 2;
            break;
        default:
            showingCoreWidth = processedDescription.get_width() / 2;
            showingMainColumnWidth = mainColumn.get_width() / 2;
    }

    double imageWidth = processedDescription.get_width();
    double imageHeight = processedDescription.get_height();
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
    *shapes << SVG::Polygon(scale_points(topPiecePoints, imageHeight));
    auto topPiece = _root->get_children<SVG::Polygon>().back();
    topPiece->set_attr("class", "ferrite");
    *shapes << SVG::Polygon(scale_points(bottomPiecePoints, imageHeight));
    auto bottomPiece = _root->get_children<SVG::Polygon>().back();
    bottomPiece->set_attr("class", "ferrite");
    for (auto& chunk : gapChunks) {
        *shapes << SVG::Polygon(scale_points(chunk, imageHeight));
        auto chunkPiece = _root->get_children<SVG::Polygon>().back();
        chunkPiece->set_attr("class", "ferrite");
    }

    _root->autoscale();
    _root->set_attr("width", imageWidth / 2 * 100000).set_attr("height", imageHeight * 100000);  // TODO remove

    if (!std::filesystem::exists(_filepath)) {
        std::filesystem::create_directory(_filepath);
    }
    std::ofstream outfile(_filepath.replace_filename(_filename));
    outfile << std::string(*_root);

    return _root;
}

SVG::SVG* CoilPainter::paint_core(MagneticWrapper magnetic) {
    CoreWrapper core = magnetic.get_core();
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
            break;
        default:
            return paint_two_piece_set_core(core);
            break;
    }
}

SVG::SVG* CoilPainter::paint_bobbin(MagneticWrapper magnetic) {
    CoreWrapper core = magnetic.get_core();
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes does not have bobbins");
            break;
        default:
            return paint_two_piece_set_bobbin(magnetic);
            break;
    }
}

SVG::SVG* CoilPainter::paint_winding_sections(MagneticWrapper magnetic) {
    CoreWrapper core = magnetic.get_core();
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto windingWindows = core.get_winding_windows();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
            break;
        default:
            return paint_two_piece_set_winding_sections(magnetic);
            break;
    }
}

SVG::SVG* CoilPainter::paint_winding_layers(MagneticWrapper magnetic) {
    CoreWrapper core = magnetic.get_core();
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto windingWindows = core.get_winding_windows();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
            break;
        default:
            return paint_two_piece_set_winding_layers(magnetic);
            break;
    }
}

SVG::SVG* CoilPainter::paint_winding_turns(MagneticWrapper magnetic) {
    CoreWrapper core = magnetic.get_core();
    CoreShape shape = std::get<CoreShape>(core.get_functional_description().get_shape());
    auto windingWindows = core.get_winding_windows();
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            throw std::runtime_error("T shapes not implemented yet");
            break;
        default:
            return paint_two_piece_set_winding_turns(magnetic);
            break;
    }
}

} // namespace OpenMagnetics
