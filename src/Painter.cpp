#include "MagneticField.h"
#include "Painter.h"
#include "CoilMesher.h"
#include "MAS.hpp"
#include "Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <../tests/TestingUtils.h>

namespace OpenMagnetics {

ComplexField Painter::calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex) {
    double bobbinWidthStart = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_coordinates().value()[0] - magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value() / 2;
    double bobbinWidth = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value();
    double coreColumnWidth = magnetic.get_mutable_core().get_columns()[0].get_width();
    double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    std::vector<double> bobbinPointsX = matplot::linspace(coreColumnWidth / 2, bobbinWidthStart + bobbinWidth, _numberPointsX);
    std::vector<double> bobbinPointsY = matplot::linspace(-coreColumnHeight / 2, coreColumnHeight / 2, _numberPointsY);
    std::vector<OpenMagnetics::FieldPoint> points;
    for (size_t j = 0; j < bobbinPointsY.size(); ++j) {
        for (size_t i = 0; i < bobbinPointsX.size(); ++i) {
            OpenMagnetics::FieldPoint fieldPoint;
            fieldPoint.set_point(std::vector<double>{bobbinPointsX[i], bobbinPointsY[j]});
            points.push_back(fieldPoint);
        }
    }

    Field inducedField;
    inducedField.set_data(points);
    inducedField.set_frequency(frequency);

    OpenMagnetics::MagneticField magneticField;
    magneticField.set_fringing_effect(_includeFringing);
    magneticField.set_mirroring_dimension(_mirroringDimension);
    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);

    auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
    return field;
}

void Painter::paint_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    matplot::gcf()->quiet_mode(true);
    matplot::cla();
    double bobbinWidthStart = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_coordinates().value()[0] - magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value() / 2;
    double bobbinWidth = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value();
    double bobbinHeight = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_height().value();
    double minimumModule = DBL_MAX;
    double maximumModule = 0;

    ComplexField field;
    if (inputField) {
        field = inputField.value();
    }
    else {
        field = calculate_magnetic_field(operatingPoint, magnetic, harmonicIndex);
    }

    if (_mode == PainterModes::CONTOUR) {

        std::vector<std::vector<double>> X(_numberPointsY, std::vector<double>(_numberPointsX, 0));
        std::vector<std::vector<double>> Y(_numberPointsY, std::vector<double>(_numberPointsX, 0));
        std::vector<std::vector<double>> M(_numberPointsY, std::vector<double>(_numberPointsX, 0));

        for (size_t j = 0; j < _numberPointsY; ++j) { 
            for (size_t i = 0; i < _numberPointsX; ++i) { 
                X[j][i] = field.get_data()[_numberPointsX * j + i].get_point()[0];
                Y[j][i] = field.get_data()[_numberPointsX * j + i].get_point()[1];
                if (_logarithmicScale) {
                    M[j][i] = hypot(log(fabs(field.get_data()[_numberPointsX * j + i].get_real())), log(fabs(field.get_data()[_numberPointsX * j + i].get_imaginary())));
                }
                else {
                    M[j][i] = hypot(field.get_data()[_numberPointsX * j + i].get_real(), field.get_data()[_numberPointsX * j + i].get_imaginary());
                }
                minimumModule = std::min(minimumModule, M[j][i]);
                maximumModule = std::max(maximumModule, M[j][i]);
            }
        }
        matplot::gcf()->quiet_mode(true);
        matplot::contourf(X, Y, M);
    }
    else if (_mode == PainterModes::QUIVER) {

        std::vector<double> X(_numberPointsY * _numberPointsX, 0);
        std::vector<double> Y(_numberPointsY * _numberPointsX, 0);
        std::vector<double> U(_numberPointsY * _numberPointsX, 0);
        std::vector<double> V(_numberPointsY * _numberPointsX, 0);
        std::vector<double> M(_numberPointsY * _numberPointsX, 0);

        for (size_t i = 0; i < _numberPointsY * _numberPointsX; ++i) {
            X[i] = field.get_data()[i].get_point()[0];
            Y[i] = field.get_data()[i].get_point()[1];
            if (_logarithmicScale) {
                U[i] = log(fabs(field.get_data()[i].get_real()));
                if (field.get_data()[i].get_real() < 0) {
                    U[i] = -U[i];
                }
                V[i] = log(fabs(field.get_data()[i].get_imaginary()));
                if (field.get_data()[i].get_imaginary() < 0) {
                    V[i] = -V[i];
                }
                M[i] = hypot(V[i], U[i]);
            }
            else {
                U[i] = field.get_data()[i].get_real();
                V[i] = field.get_data()[i].get_imaginary();
                M[i] = hypot(field.get_data()[i].get_real(), field.get_data()[i].get_imaginary());
            }
            minimumModule = std::min(minimumModule, M[i]);
            maximumModule = std::max(maximumModule, M[i]);
        }

        matplot::gcf()->quiet_mode(true);
        auto q = matplot::quiver(X, Y, U, V, M, 0.0001)->normalize(true).line_width(1.5);
        matplot::colorbar();
    }
    else {
        throw std::runtime_error("Unknown field painter mode");
    }

    if (_maximumValueColorbar) {
        maximumModule = _maximumValueColorbar.value();
    }
    if (_minimumValueColorbar) {
        minimumModule = _minimumValueColorbar.value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }
    matplot::colormap(matplot::palette::jet());
    matplot::colorbar().label("Magnetic Field Strength (A/m)");
    matplot::gca()->cblim(std::array<double, 2>{minimumModule, maximumModule});
    matplot::gcf()->size(bobbinWidth * _scale / 0.7, bobbinHeight * _scale);
    matplot::xlim({bobbinWidthStart, bobbinWidthStart + bobbinWidth});
    matplot::ylim({-bobbinHeight / 2, bobbinHeight / 2});
    matplot::gca()->position({0.0f, 0.0f, 0.7f, 1.0f});
    matplot::gca()->cb_inside(false);
    matplot::xticks({});
    matplot::yticks({});
}

void Painter::export_svg() {
    auto outFile = _filepath;
    matplot::save(outFile);
}
void Painter::export_png() {
    auto outFile = _filepath;
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile),std::regex(".svg"), ".png"));
    matplot::save(outFile);
}


void Painter::paint_core(MagneticWrapper magnetic) {
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

void Painter::paint_bobbin(MagneticWrapper magnetic) {
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

void Painter::paint_coil_sections(MagneticWrapper magnetic) {
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

void Painter::paint_coil_layers(MagneticWrapper magnetic) {
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

void Painter::paint_coil_turns(MagneticWrapper magnetic) {
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

void Painter::paint_two_piece_set_core(CoreWrapper core) {
    std::vector<std::vector<double>> topPiecePoints = {};
    std::vector<std::vector<double>> bottomPiecePoints = {};
    std::vector<std::vector<std::vector<double>>> gapChunks = {};
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

    double coreWidth = processedDescription.get_width();
    double coreHeight = processedDescription.get_height();


    matplot::gcf()->size(coreWidth / 2 * _scale, coreHeight * _scale);
    matplot::xlim({0, coreWidth / 2});
    matplot::ylim({-coreHeight / 2, coreHeight / 2});
    matplot::gca()->cb_inside(true);
    matplot::gca()->cb_position({0.05f, 0.05f, 0.05f, 0.9f});
    matplot::gca()->position({0.0f, 0.0f, 1.0f, 1.0f});

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



    topPiecePoints.push_back(std::vector<double>({0, processedDescription.get_height() / 2}));
    topPiecePoints.push_back(std::vector<double>({showingCoreWidth, processedDescription.get_height() / 2}));
    topPiecePoints.push_back(std::vector<double>({showingCoreWidth, lowestHeightTopCoreRightColumn}));
    topPiecePoints.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, lowestHeightTopCoreRightColumn}));
    topPiecePoints.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, rightColumn.get_height() / 2}));
    topPiecePoints.push_back(std::vector<double>({showingMainColumnWidth, mainColumn.get_height() / 2}));
    topPiecePoints.push_back(std::vector<double>({showingMainColumnWidth, lowestHeightTopCoreMainColumn}));
    topPiecePoints.push_back(std::vector<double>({0, lowestHeightTopCoreMainColumn}));

    for (size_t i = 1; i < gapsInMainColumn.size(); ++i)
    {
        std::vector<std::vector<double>> chunk;
        chunk.push_back(std::vector<double>({0, gapsInMainColumn[i - 1].get_coordinates().value()[1] - gapsInMainColumn[i - 1].get_length() / 2}));
        chunk.push_back(std::vector<double>({showingMainColumnWidth, gapsInMainColumn[i - 1].get_coordinates().value()[1] - gapsInMainColumn[i - 1].get_length() / 2}));
        chunk.push_back(std::vector<double>({showingMainColumnWidth, gapsInMainColumn[i].get_coordinates().value()[1] + gapsInMainColumn[i].get_length() / 2}));
        chunk.push_back(std::vector<double>({0, gapsInMainColumn[i].get_coordinates().value()[1] + gapsInMainColumn[i].get_length() / 2}));
        gapChunks.push_back(chunk);
    }
    for (size_t i = 1; i < gapsInRightColumn.size(); ++i)
    {
        std::vector<std::vector<double>> chunk;
        chunk.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, gapsInRightColumn[i - 1].get_coordinates().value()[1] - gapsInRightColumn[i - 1].get_length() / 2}));
        chunk.push_back(std::vector<double>({showingCoreWidth, gapsInRightColumn[i - 1].get_coordinates().value()[1] - gapsInRightColumn[i - 1].get_length() / 2}));
        chunk.push_back(std::vector<double>({showingCoreWidth, gapsInRightColumn[i].get_coordinates().value()[1] + gapsInRightColumn[i].get_length() / 2}));
        chunk.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, gapsInRightColumn[i].get_coordinates().value()[1] + gapsInRightColumn[i].get_length() / 2}));
        gapChunks.push_back(chunk);
    }

    bottomPiecePoints.push_back(std::vector<double>({0, -processedDescription.get_height() / 2}));
    bottomPiecePoints.push_back(std::vector<double>({showingCoreWidth, -processedDescription.get_height() / 2}));
    bottomPiecePoints.push_back(std::vector<double>({showingCoreWidth, highestHeightBottomCoreRightColumn}));
    bottomPiecePoints.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, highestHeightBottomCoreRightColumn}));
    bottomPiecePoints.push_back(std::vector<double>({showingCoreWidth - rightColumnWidth, -rightColumn.get_height() / 2}));
    bottomPiecePoints.push_back(std::vector<double>({showingMainColumnWidth, -mainColumn.get_height() / 2}));
    bottomPiecePoints.push_back(std::vector<double>({showingMainColumnWidth, highestHeightBottomCoreMainColumn}));
    bottomPiecePoints.push_back(std::vector<double>({0, highestHeightBottomCoreMainColumn}));

    {
        std::vector<double> x, y;
        for (auto& point : topPiecePoints) {
            x.push_back(point[0]);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorFerrite));
    }

    {
        std::vector<double> x, y;
        for (auto& point : bottomPiecePoints) {
            x.push_back(point[0]);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorFerrite));
    }

    for (auto& chunk : gapChunks){
        std::vector<double> x, y;
        for (auto& point : chunk) {
            x.push_back(point[0]);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorFerrite));
    }
}

void Painter::paint_two_piece_set_bobbin(MagneticWrapper magnetic) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin has not being processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    CoreWrapper core = magnetic.get_core();

    std::vector<double> bobbinCoordinates = std::vector<double>({0, 0, 0});
    if (bobbinProcessedDescription.get_coordinates()) {
        bobbinCoordinates = bobbinProcessedDescription.get_coordinates().value();
    }

    double bobbinOuterWidth = bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() + bobbinProcessedDescription.get_winding_windows()[0].get_width().value();
    double bobbinOuterHeight = bobbinProcessedDescription.get_wall_thickness();
    for (auto& windingWindow: bobbinProcessedDescription.get_winding_windows()) {
        bobbinOuterHeight += windingWindow.get_height().value();
        bobbinOuterHeight += bobbinProcessedDescription.get_wall_thickness();
    }

    std::vector<std::vector<double>> bobbinPoints = {};
    bobbinPoints.push_back(std::vector<double>({bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() - bobbinProcessedDescription.get_column_thickness(),
                                          bobbinCoordinates[1] + bobbinOuterHeight / 2}));
    bobbinPoints.push_back(std::vector<double>({bobbinOuterWidth,
                                          bobbinCoordinates[1] + bobbinOuterHeight / 2}));
    bobbinPoints.push_back(std::vector<double>({bobbinOuterWidth,
                                          bobbinCoordinates[1] + bobbinOuterHeight / 2 - bobbinProcessedDescription.get_wall_thickness()}));
    bobbinPoints.push_back(std::vector<double>({bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value(),
                                          bobbinCoordinates[1] + bobbinOuterHeight / 2 - bobbinProcessedDescription.get_wall_thickness()}));
    bobbinPoints.push_back(std::vector<double>({bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value(),
                                          bobbinCoordinates[1] - bobbinOuterHeight / 2 + bobbinProcessedDescription.get_wall_thickness()}));
    bobbinPoints.push_back(std::vector<double>({bobbinOuterWidth,
                                          bobbinCoordinates[1] - bobbinOuterHeight / 2 + bobbinProcessedDescription.get_wall_thickness()}));
    bobbinPoints.push_back(std::vector<double>({bobbinOuterWidth,
                                          bobbinCoordinates[1] - bobbinOuterHeight / 2}));
    bobbinPoints.push_back(std::vector<double>({bobbinCoordinates[0] + bobbinProcessedDescription.get_column_width().value() - bobbinProcessedDescription.get_column_thickness(),
                                          bobbinCoordinates[1] - bobbinOuterHeight / 2}));

    std::vector<double> x, y;
    for (auto& point : bobbinPoints) {
        x.push_back(point[0]);
        y.push_back(point[1]);
    }

    matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorBobbin));
}

void Painter::paint_two_piece_set_winding_sections(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();

    if (!magnetic.get_coil().get_sections_description()) {
        throw std::runtime_error("Winding sections not created");
    }

    auto sections = magnetic.get_coil().get_sections_description().value();

    for (size_t i = 0; i < sections.size(); ++i){

        std::vector<std::vector<double>> sectionPoints = {};
        sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
        sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
        sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
        sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));

        std::vector<double> x, y;
        for (auto& point : sectionPoints) {
            x.push_back(point[0]);
            y.push_back(point[1]);
        }
        if (sections[i].get_type() == ElectricalType::CONDUCTION) {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorCopper));
        }
        else {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorInsulation));
        }
    }
}

void Painter::paint_two_piece_set_winding_layers(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core has not being processed");
    }

    if (!winding.get_layers_description()) {
        throw std::runtime_error("Winding layers not created");
    }

    auto layers = winding.get_layers_description().value();

    for (size_t i = 0; i < layers.size(); ++i){



        std::vector<std::vector<double>> layerPoints = {};
        layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2}));
        layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2}));
        layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2}));
        layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2}));

        std::vector<double> x, y;
        for (auto& point : layerPoints) {
            x.push_back(point[0]);
            y.push_back(point[1]);
        }
        if (layers[i].get_type() == ElectricalType::CONDUCTION) {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorCopper));
        }
        else {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorInsulation));
        }
    }
}

void Painter::paint_two_piece_set_winding_turns(MagneticWrapper magnetic) {
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    auto wirePerWinding = winding.get_wires();

    if (!winding.get_turns_description()) {
        throw std::runtime_error("Winding turns not created");
    }

    auto turns = winding.get_turns_description().value();

    for (size_t i = 0; i < turns.size(); ++i){

        auto windingIndex = winding.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            double outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
            matplot::ellipse(turns[i].get_coordinates()[0] - outerDiameter / 2, turns[i].get_coordinates()[1] - outerDiameter / 2, outerDiameter, outerDiameter)->fill(true).color(matplot::to_array(_colorInsulation));

            if (wire.get_conducting_diameter()) {
                double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
                matplot::ellipse(turns[i].get_coordinates()[0] - conductingDiameter / 2, turns[i].get_coordinates()[1] - conductingDiameter / 2, conductingDiameter, conductingDiameter)->fill(true).color(matplot::to_array(_colorCopper));
            }
        }
        else {
            {
                std::vector<std::vector<double>> turnPoints = {};
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_outer_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_outer_height().value()) / 2}));

                std::vector<double> x, y;
                for (auto& point : turnPoints) {
                    x.push_back(point[0]);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorInsulation));
            }

            if (wire.get_conducting_width() && wire.get_conducting_height()) {
                std::vector<std::vector<double>> turnPoints = {};
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));

                std::vector<double> x, y;
                for (auto& point : turnPoints) {
                    x.push_back(point[0]);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorCopper));
            }
            
        }
    }

    auto layers = winding.get_layers_description().value();

    for (size_t i = 0; i < layers.size(); ++i){
        if (layers[i].get_type() == ElectricalType::INSULATION) {

            std::vector<std::vector<double>> layerPoints = {};
            layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2}));
            layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] + layers[i].get_dimensions()[1] / 2}));
            layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] + layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2}));
            layerPoints.push_back(std::vector<double>({layers[i].get_coordinates()[0] - layers[i].get_dimensions()[0] / 2, layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2}));

            std::vector<double> x, y;
            for (auto& point : layerPoints) {
                x.push_back(point[0]);
                y.push_back(point[1]);
            }
            if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorCopper));
            }
            else {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(_colorInsulation));
            }
        }
    }

}

} // namespace OpenMagnetics
