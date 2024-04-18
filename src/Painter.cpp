#include "MagneticField.h"
#include "Painter.h"
#include "CoilMesher.h"
#include "MAS.hpp"
#include "Utils.h"
#include "Settings.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <chrono>
#include <thread>

namespace OpenMagnetics {

std::string uint_to_hex(uint32_t i) {
  std::stringstream stream;
  stream << "0x" << std::setfill ('0') << std::setw(8) << std::hex << i;
  return stream.str();
}

std::string key_to_rgb_color(uint32_t i) {
  std::stringstream stream;
  stream << "rgb(" << std::setfill (' ') << std::setw(3) << ((i >> 16) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 8) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 0) & 0xFF) << ")";
  return stream.str();
}

std::string replace_key(std::string key, std::string line, std::string replacement) {
    size_t pos = 0;
    std::string token;
    std::string newLine;
    while ((pos = line.find(key)) != std::string::npos) {
        token = line.substr(0, pos);
        newLine += token;
        newLine += replacement;
        line.erase(0, pos + key.length());
    }
    newLine += line;
    return newLine;
}

ComplexField Painter::calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings->get_painter_include_fringing();
    bool mirroringDimension = settings->get_painter_mirroring_dimension();

    _extraDimension = CoilWrapper::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());

    size_t numberPointsX = settings->get_painter_number_points_x();
    size_t numberPointsY = settings->get_painter_number_points_y();
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY).first;

    OpenMagnetics::MagneticField magneticField;
    settings->set_magnetic_field_include_fringing(includeFringing);
    settings->set_magnetic_field_mirroring_dimension(mirroringDimension);
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

bool is_inside_inducing_turns(std::vector<double> point, CoilWrapper coil) {
    auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();

    for (auto turn : turns) {
        double turnRadius = 0;

        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        if (wires[windingIndex].get_type() == WireType::ROUND || wires[windingIndex].get_type() == WireType::LITZ) {
            turnRadius = wires[windingIndex].get_maximum_outer_width() / 2;
        }
        else {
            // Not implemented here for rectangulars
            return false;
        }

        double distanceX = fabs(point[0] - turn.get_coordinates()[0]);
        double distanceY = fabs(point[1] - turn.get_coordinates()[1]);
        if (hypot(distanceX, distanceY) < turnRadius) {
            return true;
        }

        if (turn.get_additional_coordinates()) {
            auto additionalCoordinates = turn.get_additional_coordinates().value();
            for (auto additionalCoordinate : additionalCoordinates) {
                double distanceX = fabs(point[0] - additionalCoordinate[0]);
                double distanceY = fabs(point[1] - additionalCoordinate[1]);
                if (hypot(distanceX, distanceY) < turnRadius) {
                    return true;
                }

            }
        }
    }

    return false;
}

void Painter::paint_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    matplot::gcf()->quiet_mode(true);
    matplot::cla();
    get_image_size(magnetic);

    double minimumModule = DBL_MAX;
    double maximumModule = 0;

    size_t numberPointsX = settings->get_painter_number_points_x();
    size_t numberPointsY = settings->get_painter_number_points_y();
    auto mode = settings->get_painter_mode();
    bool logarithmicScale = settings->get_painter_logarithmic_scale();

    ComplexField field;
    if (inputField) {
        field = inputField.value();
    }
    else {
        field = calculate_magnetic_field(operatingPoint, magnetic, harmonicIndex);
    }

    if (mode == PainterModes::CONTOUR) {

        std::vector<std::vector<double>> X(numberPointsY, std::vector<double>(numberPointsX, 0));
        std::vector<std::vector<double>> Y(numberPointsY, std::vector<double>(numberPointsX, 0));
        std::vector<std::vector<double>> M(numberPointsY, std::vector<double>(numberPointsX, 0));

        for (size_t j = 0; j < numberPointsY; ++j) { 
            for (size_t i = 0; i < numberPointsX; ++i) {
                X[j][i] = _offsetForColorBar + field.get_data()[numberPointsX * j + i].get_point()[0];
                Y[j][i] = field.get_data()[numberPointsX * j + i].get_point()[1];

                if (logarithmicScale) {
                    M[j][i] = hypot(log(fabs(field.get_data()[numberPointsX * j + i].get_real())), log(fabs(field.get_data()[numberPointsX * j + i].get_imaginary())));
                }
                else {
                    M[j][i] = hypot(field.get_data()[numberPointsX * j + i].get_real(), field.get_data()[numberPointsX * j + i].get_imaginary());
                }
                minimumModule = std::min(minimumModule, M[j][i]);
                maximumModule = std::max(maximumModule, M[j][i]);

            }
        }
        matplot::gcf()->quiet_mode(true);
        auto c = matplot::contourf(X, Y, M);
        c->font_size(99);
    }
    else if (mode == PainterModes::QUIVER || mode == PainterModes::SCATTER) {

        std::vector<double> X(field.get_data().size(), 0);
        std::vector<double> Y(field.get_data().size(), 0);
        std::vector<double> U(field.get_data().size(), 0);
        std::vector<double> V(field.get_data().size(), 0);
        std::vector<double> M(field.get_data().size(), 0); 

        for (size_t i = 0; i < field.get_data().size(); ++i) {
            if (is_inside_inducing_turns(field.get_data()[i].get_point(), magnetic.get_coil())) {
                X[i] = std::numeric_limits<double>::quiet_NaN();;
                Y[i] = std::numeric_limits<double>::quiet_NaN();;
                U[i] = std::numeric_limits<double>::quiet_NaN();;
                V[i] = std::numeric_limits<double>::quiet_NaN();;
                M[i] = std::numeric_limits<double>::quiet_NaN();;
                continue;
            }

            X[i] = _offsetForColorBar + field.get_data()[i].get_point()[0];
            Y[i] = field.get_data()[i].get_point()[1];

            if (logarithmicScale) {
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

        for (int i = field.get_data().size() - 1; i >= 0; --i) {
            if (std::isnan(X[i])) {
                X.erase(X.begin() + i);
                Y.erase(Y.begin() + i);
                U.erase(U.begin() + i);
                V.erase(V.begin() + i);
                M.erase(M.begin() + i);
            }
        }

        matplot::gcf()->quiet_mode(true);
        if (mode == PainterModes::QUIVER) {
            auto q = matplot::quiver(X, Y, U, V, M, 0.0001)->normalize(true).line_width(1.5);
        }
        else {
            auto c = matplot::scatter(X, Y, 30, M);
            c->marker_face(true);
        }
        matplot::colorbar();
    }
    else {
        throw std::runtime_error("Unknown field painter mode");
    }

    if (settings->get_painter_maximum_value_colorbar()) {
        maximumModule = settings->get_painter_maximum_value_colorbar().value();
    }
    if (settings->get_painter_maximum_value_colorbar()) {
        minimumModule = settings->get_painter_minimum_value_colorbar().value();
    }
    if (minimumModule == maximumModule) {
        minimumModule = maximumModule - 1;
    }
    matplot::colormap(matplot::palette::jet());

    int maximumDecimals = ceil(log10(maximumModule));

    std::vector<double> tickValues = {minimumModule};
    std::ostringstream oss;
    oss << std::setprecision(0) << std::fixed << roundFloat(minimumModule, 0);
    std::vector<std::string> tickLabels = {oss.str() + " A/m"};

    double lastValue = 0;
    int precision = 0;
    while (tickLabels.size() < 3) {
        tickLabels = {oss.str() + " A/m"};
        for (size_t value = minimumModule; value < maximumModule; value+=pow(10, maximumDecimals - 1)) {
                if (value == 0) {
                    continue;
                }
                lastValue = value;
                tickValues.push_back(roundFloat(value, -maximumDecimals + 1));
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(precision) << std::scientific << (roundFloat(value, -maximumDecimals + 1));
                auto label = replace_key("+0", oss.str(), "");
                if (label == "0e0") {
                    label = "0";
                }
                tickLabels.push_back(label + " A/m");
        }
        maximumDecimals -= 1;
        precision += 1;
        if (maximumDecimals <=0) {
            break;
        }
    }
    auto cb = matplot::colorbar().tick_values({tickValues}).ticklabels({tickLabels});
    matplot::gca()->cblim(std::array<double, 2>{minimumModule, maximumModule * 0.99});
    matplot::xticks({});
    matplot::yticks({});
}

void Painter::export_svg() {
    auto outFile = std::string {_filepath.string()};
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\\\"), std::string("/"))).string();

    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex(" "), "_")).string();
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\,"), "_")).string();
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\."), "_")).string();
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\:"), "_")).string();
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\/"), "_")).string();

    auto tempFile = outFile + "_temp.svg";
    matplot::save(tempFile);
    // TODO add more complex polling, check WebBackend
    std::this_thread::sleep_for(std::chrono::milliseconds(300));

    std::ofstream finalFile(outFile);

    std::string svgFile;
    std::string line;

    std::ifstream ifTempFile(tempFile);

    std::vector<std::string> lines;

    while (getline (ifTempFile, line)) {
        lines.push_back(line);
    }

    for (size_t lineIndex = 1; lineIndex < lines.size(); ++lineIndex) {
        std::string currentLine = lines[lineIndex];
        std::string newLine = lines[lineIndex];
        for (auto [key, changedLine] : _postProcessingChanges){
            if (currentLine.contains(key)) {
                if (lines[lineIndex - 1].contains(R"(stroke-linecap="butt")")) {
                    lines[lineIndex - 1] = replace_key(R"(stroke-linecap="butt")", lines[lineIndex - 1], _postProcessingChanges[key]);
                }

                lines[lineIndex] = replace_key(key, currentLine, _postProcessingColors[key]);

                if (lines[lineIndex].contains(R"(polygon)")) {
                    lines[lineIndex - 2] = replace_key(R"(<g )", lines[lineIndex - 2], _postProcessingChanges[key]);
                }

            }
        }

        if (lines[lineIndex].contains(R"(font-size="10.00")")) {
            lines[lineIndex] = replace_key(R"(font-size="10.00")", lines[lineIndex], R"(font-size=")" + std::to_string(_fontSize) + "\"");
        }

        if (lines[lineIndex].contains(R"(font-size="11.00")")) {
            lines[lineIndex] = replace_key(R"(font-size="11.00")", lines[lineIndex], R"(font-size=")" + std::to_string(_fontSize) + "\"");
        }


    }


    for (size_t lineIndex = 0; lineIndex < lines.size(); ++lineIndex) {
        finalFile << lines[lineIndex] << std::endl;
    }
    finalFile.close();
}

void Painter::export_png() {
    auto outFile = std::string {_filepath.string()};
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex("\\\\"), std::string("/"))).string();
    outFile = std::filesystem::path(std::regex_replace(std::string(outFile), std::regex(".svg"), std::string(".png")))
                  .string();
    matplot::save(outFile);
}

void Painter::paint_core(MagneticWrapper magnetic) {
    CoreShape shape = std::get<CoreShape>(magnetic.get_core().get_functional_description().get_shape());
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return paint_toroidal_core(magnetic);
            break;
        default:
            return paint_two_piece_set_core(magnetic);
            break;
    }
}

void Painter::paint_bobbin(MagneticWrapper magnetic) {
    CoreShape shape = std::get<CoreShape>(magnetic.get_core().get_functional_description().get_shape());
    switch(shape.get_family()) {
        case CoreShapeFamily::T:
            return;
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
            return paint_toroidal_winding_sections(magnetic);
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
            return paint_toroidal_winding_layers(magnetic);
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
            return paint_toroidal_winding_turns(magnetic);
            break;
        default:
            return paint_two_piece_set_winding_turns(magnetic);
            break;
    }
}

std::vector<double> Painter::get_image_size(MagneticWrapper magnetic) {
    auto core = magnetic.get_core();


    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    
    double showingCoreWidth;
    auto family = core.get_shape_family();
    switch (family) {
        case OpenMagnetics::CoreShapeFamily::U:
        case OpenMagnetics::CoreShapeFamily::UR:
            _extraDimension = 1;
            showingCoreWidth = (processedDescription.get_width() - mainColumn.get_width() / 2);
            if (_addProportionForColorBar) {
                double proportionForColorBar = 0.25;
                showingCoreWidth *= (1 + proportionForColorBar);
            }
            _offsetForColorBar = (showingCoreWidth - (processedDescription.get_width() - mainColumn.get_width() / 2)) / 2;
            break;
        case OpenMagnetics::CoreShapeFamily::T:
            _extraDimension = CoilWrapper::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());
            showingCoreWidth = processedDescription.get_width() * _extraDimension;
            if (_addProportionForColorBar) {
                double proportionForColorBar = 0.25;
                showingCoreWidth *= (1 + proportionForColorBar);
            }
            _offsetForColorBar = (showingCoreWidth - processedDescription.get_width() * _extraDimension) / 2;
            break;
        default:
            _extraDimension = 1;
            showingCoreWidth = processedDescription.get_width() / 2;
            if (_addProportionForColorBar) {
                double proportionForColorBar = 0.4;
                showingCoreWidth *= (1 + proportionForColorBar);
            }
            _offsetForColorBar = (showingCoreWidth - processedDescription.get_width() / 2);
            break;
    }


    double showingCoreHeight = processedDescription.get_height() * _extraDimension;

    _fontSize = std::max(1.0, 10 * showingCoreHeight / 0.01);
    _fontSize = std::min(100.0, _fontSize);

    return {showingCoreWidth, showingCoreHeight};
}

void Painter::paint_toroidal_core(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    auto aux = get_image_size(magnetic);
    double imageWidth = aux[0];
    double imageHeight = aux[1];

    matplot::gcf()->size(imageWidth * _scale, imageHeight * _scale);
    matplot::xlim({-imageWidth / 2, imageWidth / 2});
    matplot::ylim({-imageHeight / 2, imageHeight / 2});
    matplot::gca()->cb_inside(true);
    matplot::gca()->cb_position({0.01f, 0.05f, 0.05f, 0.9f});
    matplot::gca()->position({0.0f, 0.0f, 1.0f, 1.0f});


    auto core = magnetic.get_core();
    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});

    double strokeWidth = mainColumn.get_width();
    double circleDiameter = processedDescription.get_width() - strokeWidth;

    auto currentMapIndex = uint_to_hex(_currentMapIndex);
    auto key = key_to_rgb_color(_currentMapIndex);
    _currentMapIndex++;

    matplot::ellipse(_offsetForColorBar -circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
    _postProcessingChanges[key] = R"(stroke-linecap="round")";
    _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_ferrite(), 0, 16));
}

void Painter::paint_two_piece_set_core(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto core = magnetic.get_core();
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

    double coreHeight = processedDescription.get_height();

    auto aux = get_image_size(magnetic);
    double imageWidth = aux[0];
    double imageHeight = aux[1];

    matplot::gcf()->size(imageWidth * _scale, imageHeight * _scale);
    matplot::xlim({0, imageWidth});
    matplot::ylim({-imageHeight / 2, imageHeight / 2});
    matplot::gca()->cb_inside(true);
    matplot::gca()->cb_position({0.01f, 0.05f, 0.05f, 0.9f});
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
            x.push_back(point[0] + _offsetForColorBar);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_ferrite()));
    }

    {
        std::vector<double> x, y;
        for (auto& point : bottomPiecePoints) {
            x.push_back(point[0] + _offsetForColorBar);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_ferrite()));
    }

    for (auto& chunk : gapChunks){
        std::vector<double> x, y;
        for (auto& point : chunk) {
            x.push_back(point[0] + _offsetForColorBar);
            y.push_back(point[1]);
        }

        matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_ferrite()));
    }
}

void Painter::paint_two_piece_set_bobbin(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin has not being processed");
    }
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();

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
        x.push_back(point[0] + _offsetForColorBar);
        y.push_back(point[1]);
    }

    matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_bobbin()));
}

void Painter::paint_two_piece_set_margin(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
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
                std::vector<std::vector<double>> marginPoints = {};
                if (sectionsOrientation == WindingOrientation::OVERLAPPING) {
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] + windingWindowDimensions[1] / 2}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] + windingWindowDimensions[1] / 2}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] + windingWindowDimensions[1] / 2 - margins[0]}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] + windingWindowDimensions[1] / 2 - margins[0]}));
                }
                else if (sectionsOrientation == WindingOrientation::CONTIGUOUS) {
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] - windingWindowDimensions[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] - windingWindowDimensions[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] - windingWindowDimensions[0] / 2 + margins[0], sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] - windingWindowDimensions[0] / 2 + margins[0], sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
                }

                std::vector<double> x, y;
                for (auto& point : marginPoints) {
                    x.push_back(point[0] + _offsetForColorBar);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_margin()));
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
                if (sectionsOrientation == WindingOrientation::OVERLAPPING) {
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] - windingWindowDimensions[1] / 2}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] - windingWindowDimensions[1] / 2}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] - windingWindowDimensions[1] / 2 + margins[1]}));
                    marginPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, bobbinCoordinates[1] + windingWindowCoordinates[1] - windingWindowDimensions[1] / 2 + margins[1]}));
                }
                else if (sectionsOrientation == WindingOrientation::CONTIGUOUS) {
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] + windingWindowDimensions[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] + windingWindowDimensions[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] + windingWindowDimensions[0] / 2 - margins[1], sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
                    marginPoints.push_back(std::vector<double>({bobbinCoordinates[0] + windingWindowCoordinates[0] + windingWindowDimensions[0] / 2 - margins[1], sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
                }

                std::vector<double> x, y;
                for (auto& point : marginPoints) {
                    x.push_back(point[0] + _offsetForColorBar);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_margin()));
            }
        }
    }
}

void Painter::paint_toroidal_margin(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    bool drawSpacer = settings->get_painter_painter_draw_spacer();
    auto sections = magnetic.get_coil().get_sections_description().value();

    if (sections.size() == 1) {
        return;
    }

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    double coreWidth = processedDescription.get_width();
    double coreHeight = processedDescription.get_height();
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
            auto margins = sections[i].get_margin().value();

            if (sectionOrientation == WindingOrientation::CONTIGUOUS) {

                if (drawSpacer) {
                    size_t nextSectionIndex = 0;
                    if (i < sections.size() - 2) {
                        nextSectionIndex = i + 2;
                    }
                    double leftMargin = sections[i].get_margin().value()[1];
                    double rightMargin = sections[nextSectionIndex].get_margin().value()[0];
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
                    marginPoints.push_back({centerRadialPoint[0] - rectangleThickness / 2 * sin(rectangleAngleInRadians), centerRadialPoint[1] + rectangleThickness / 2 * cos(rectangleAngleInRadians)});
                    marginPoints.push_back({centerRadialPoint[0] + rectangleThickness / 2 * sin(rectangleAngleInRadians), centerRadialPoint[1] - rectangleThickness / 2 * cos(rectangleAngleInRadians)});
                    marginPoints.push_back({0 + rectangleThickness / 2 * sin(rectangleAngleInRadians), 0 - rectangleThickness / 2 * cos(rectangleAngleInRadians)});
                    marginPoints.push_back({0 - rectangleThickness / 2 * sin(rectangleAngleInRadians), 0 + rectangleThickness / 2 * cos(rectangleAngleInRadians)});

                    std::vector<double> x, y;
                    for (auto& point : marginPoints) {
                        x.push_back(point[0] + _offsetForColorBar);
                        y.push_back(point[1]);
                    }
                    matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_spacer()));

                }
                else {
                    if (margins[0] > 0) {
                        auto currentMapIndex = uint_to_hex(_currentMapIndex);
                        auto key = key_to_rgb_color(_currentMapIndex);
                        _currentMapIndex++;
                        double strokeWidth = sections[i].get_dimensions()[0];
                        double circleDiameter = (windingWindowRadialHeight - sections[i].get_coordinates()[0]) * 2;
                        double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

                        double angle = wound_distance_to_angle(margins[0], circleDiameter / 2 - strokeWidth / 2);
                        double angleProportion = angle / 360;
                        std::string termination = angleProportion < 1? "butt" : "round";

                        if (sections[i].get_type() == ElectricalType::CONDUCTION) {

                            matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
                            _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2 - angle)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                                            R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
                            _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_margin(), 0, 16));
                 
                        }
                    }

                    if (margins[1] > 0) {
                        double strokeWidth = sections[i].get_dimensions()[0];
                        double circleDiameter = (windingWindowRadialHeight - sections[i].get_coordinates()[0]) * 2;
                        double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

                        auto currentMapIndex = uint_to_hex(_currentMapIndex);
                        auto key = key_to_rgb_color(_currentMapIndex);
                        _currentMapIndex++;

                        double angle = wound_distance_to_angle(margins[1], circleDiameter / 2 - strokeWidth / 2);
                        double angleProportion = angle / 360;
                        std::string termination = angleProportion < 1? "butt" : "round";
                        if (sections[i].get_type() == ElectricalType::CONDUCTION) {

                            matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
                            _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                                            R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
                            _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_margin(), 0, 16));
                 
                        }
                    }
                }
            }
            else {
                if (margins[0] > 0) {
                    double strokeWidth = margins[0];
                    double circleDiameter = (windingWindowRadialHeight - (sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2) + margins[0] / 2) * 2;
                    double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

                    auto currentMapIndex = uint_to_hex(_currentMapIndex);
                    auto key = key_to_rgb_color(_currentMapIndex);
                    _currentMapIndex++;

                    double angle = wound_distance_to_angle(sections[i].get_dimensions()[1], circleDiameter / 2 + strokeWidth / 2);
                    double angleProportion = angle / 360;
                    std::string termination = angleProportion < 1? "butt" : "round";
                    if (sections[i].get_type() == ElectricalType::CONDUCTION) {

                        matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
                        _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2 - angle)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                                        R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
                        _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_margin(), 0, 16));
                    }
                }
                if (margins[1] > 0) {

                    double strokeWidth = margins[1];
                    double circleDiameter = (windingWindowRadialHeight - (sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2) - margins[1] / 2) * 2;
                    double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

                    auto currentMapIndex = uint_to_hex(_currentMapIndex);
                    auto key = key_to_rgb_color(_currentMapIndex);
                    _currentMapIndex++;

                    double angle = wound_distance_to_angle(sections[i].get_dimensions()[1], circleDiameter / 2 + strokeWidth / 2);
                    double angleProportion = angle / 360;
                    std::string termination = angleProportion < 1? "butt" : "round";
                    if (sections[i].get_type() == ElectricalType::CONDUCTION) {

                        matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
                        _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2 - angle)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                                        R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
                        _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_margin(), 0, 16));
                    }
                }
            }
        }
    }

    if (drawSpacer) {
        matplot::ellipse(_offsetForColorBar -largestThickness * 2 / 2, -largestThickness * 2 / 2, largestThickness * 2, largestThickness * 2)->fill(true).color(matplot::to_array(settings->get_painter_color_spacer()));
    }
}

void Painter::paint_two_piece_set_winding_sections(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();

    if (!magnetic.get_coil().get_sections_description()) {
        throw std::runtime_error("Winding sections not created");
    }

    auto sections = magnetic.get_coil().get_sections_description().value();

    for (size_t i = 0; i < sections.size(); ++i){

        {
            std::vector<std::vector<double>> sectionPoints = {};
            sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
            sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] + sections[i].get_dimensions()[1] / 2}));
            sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] + sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));
            sectionPoints.push_back(std::vector<double>({sections[i].get_coordinates()[0] - sections[i].get_dimensions()[0] / 2, sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2}));

            std::vector<double> x, y;
            for (auto& point : sectionPoints) {
                x.push_back(point[0] + _offsetForColorBar);
                y.push_back(point[1]);
            }
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_copper()));
            }
            else {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_insulation()));
            }
        }
    }

    paint_two_piece_set_margin(magnetic);
}

void Painter::paint_toroidal_winding_sections(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    double coreWidth = processedDescription.get_width();
    double coreHeight = processedDescription.get_height();
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
            double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

            auto currentMapIndex = uint_to_hex(_currentMapIndex);
            auto key = key_to_rgb_color(_currentMapIndex);
            _currentMapIndex++;

            double angleProportion = sections[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";
            matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
            _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(sections[i].get_coordinates()[1] - sections[i].get_dimensions()[1] / 2)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                                R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
            if (sections[i].get_type() == ElectricalType::CONDUCTION) {
                _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_copper(), 0, 16));
            }
            else {
                _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_insulation(), 0, 16));
            }
        }
    }

    paint_toroidal_margin(magnetic);
}

void Painter::paint_two_piece_set_winding_layers(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
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
            x.push_back(point[0] + _offsetForColorBar);
            y.push_back(point[1]);
        }
        if (layers[i].get_type() == ElectricalType::CONDUCTION) {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_copper()));
        }
        else {
            matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_insulation()));
        }
    }

    paint_two_piece_set_margin(magnetic);
}

void Painter::paint_toroidal_winding_layers(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    CoreWrapper core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core has not being processed");
    }

    if (!winding.get_layers_description()) {
        throw std::runtime_error("Winding layers not created");
    }

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    double coreWidth = processedDescription.get_width();
    double coreHeight = processedDescription.get_height();
    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});

    auto layers = winding.get_layers_description().value();
    double initialRadius = processedDescription.get_width() / 2 - mainColumn.get_width();

    for (size_t i = 0; i < layers.size(); ++i){
        {
            double strokeWidth = layers[i].get_dimensions()[0];
            double circleDiameter = (initialRadius - layers[i].get_coordinates()[0]) * 2;
            double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

            auto currentMapIndex = uint_to_hex(_currentMapIndex);
            auto key = key_to_rgb_color(_currentMapIndex);
            _currentMapIndex++;

            double angleProportion = layers[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";

            matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
            _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2)) + " " + std::to_string(coreWidth / 2 * _scale * _extraDimension) + " " + std::to_string(coreHeight / 2 * _scale * _extraDimension) + ")\" " + 
                                            R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
            if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_copper(), 0, 16));
            }
            else {
                _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_insulation(), 0, 16));
            }
        }
    }

    paint_toroidal_margin(magnetic);
}

void Painter::paint_round_wire(double xCoordinate, double yCoordinate, WireWrapper wire) {
    if (!wire.get_outer_diameter()) {
        throw std::runtime_error("Wire is missing outerDiameter");
    }
    auto settings = OpenMagnetics::Settings::GetInstance();

    double outerDiameter = resolve_dimensional_values(wire.get_outer_diameter().value());
    matplot::ellipse(_offsetForColorBar + xCoordinate - outerDiameter / 2, yCoordinate - outerDiameter / 2, outerDiameter, outerDiameter)->fill(true).color(matplot::to_array(settings->get_painter_color_insulation()));

    if (wire.get_conducting_diameter()) {
        double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
        matplot::ellipse(_offsetForColorBar + xCoordinate - conductingDiameter / 2, yCoordinate - conductingDiameter / 2, conductingDiameter, conductingDiameter)->fill(true).color(matplot::to_array(settings->get_painter_color_copper()));
    }
}

void Painter::paint_rectangular_wire(double xCoordinate, double yCoordinate, WireWrapper wire, double angle, std::vector<double> center) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    std::vector<std::vector<double>> turnPoints = {};
    turnPoints.push_back(std::vector<double>({xCoordinate - resolve_dimensional_values(wire.get_outer_width().value()) / 2, yCoordinate + resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
    turnPoints.push_back(std::vector<double>({xCoordinate + resolve_dimensional_values(wire.get_outer_width().value()) / 2, yCoordinate + resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
    turnPoints.push_back(std::vector<double>({xCoordinate + resolve_dimensional_values(wire.get_outer_width().value()) / 2, yCoordinate - resolve_dimensional_values(wire.get_outer_height().value()) / 2}));
    turnPoints.push_back(std::vector<double>({xCoordinate - resolve_dimensional_values(wire.get_outer_width().value()) / 2, yCoordinate - resolve_dimensional_values(wire.get_outer_height().value()) / 2}));

    auto currentMapIndex = uint_to_hex(_currentMapIndex);
    auto key = key_to_rgb_color(_currentMapIndex);
    _currentMapIndex++;

    std::vector<double> x, y;
    for (auto& point : turnPoints) {
        x.push_back(point[0] + _offsetForColorBar);
        y.push_back(point[1]);
    }
    matplot::fill(x, y)->fill(true).color(matplot::to_array(currentMapIndex));

    _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_insulation(), 0, 16));
    _postProcessingChanges[key] = R"(<g transform="rotate( )" + std::to_string(-(angle)) + " " + std::to_string(center[0] * _scale) + " " + std::to_string(center[1] * _scale) + ")\" ";

    if (wire.get_conducting_width() && wire.get_conducting_height()) {
        std::vector<std::vector<double>> turnPoints = {};
        turnPoints.push_back(std::vector<double>({xCoordinate - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, yCoordinate + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
        turnPoints.push_back(std::vector<double>({xCoordinate + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, yCoordinate + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
        turnPoints.push_back(std::vector<double>({xCoordinate + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, yCoordinate - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
        turnPoints.push_back(std::vector<double>({xCoordinate - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, yCoordinate - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));

        auto currentMapIndex = uint_to_hex(_currentMapIndex);
        auto key = key_to_rgb_color(_currentMapIndex);
        _currentMapIndex++;

        std::vector<double> x, y;
        for (auto& point : turnPoints) {
            x.push_back(point[0] + _offsetForColorBar);
            y.push_back(point[1]);
        }
        matplot::fill(x, y)->fill(true).color(matplot::to_array(currentMapIndex));

        _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_copper(), 0, 16));
        _postProcessingChanges[key] = R"(<g transform="rotate( )" + std::to_string(-(angle)) + " " + std::to_string(center[0] * _scale) + " " + std::to_string(center[1] * _scale) + ")\" ";

    }
}

void Painter::paint_two_piece_set_winding_turns(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    auto wirePerWinding = winding.get_wires();

    if (!winding.get_turns_description()) {
        throw std::runtime_error("Winding turns not created");
    }

    auto turns = winding.get_turns_description().value();

    for (size_t i = 0; i < turns.size(); ++i){

        auto windingIndex = winding.get_winding_index_by_name(turns[i].get_winding());
        auto wire = wirePerWinding[windingIndex];
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            paint_round_wire(turns[i].get_coordinates()[0], turns[i].get_coordinates()[1], wire);
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
                    x.push_back(point[0] + _offsetForColorBar);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_insulation()));
            }

            if (wire.get_conducting_width() && wire.get_conducting_height()) {
                std::vector<std::vector<double>> turnPoints = {};
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] + resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] + resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));
                turnPoints.push_back(std::vector<double>({turns[i].get_coordinates()[0] - resolve_dimensional_values(wire.get_conducting_width().value()) / 2, turns[i].get_coordinates()[1] - resolve_dimensional_values(wire.get_conducting_height().value()) / 2}));

                std::vector<double> x, y;
                for (auto& point : turnPoints) {
                    x.push_back(point[0] + _offsetForColorBar);
                    y.push_back(point[1]);
                }
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_copper()));
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
                x.push_back(point[0] + _offsetForColorBar);
                y.push_back(point[1]);
            }
            if (layers[i].get_type() == ElectricalType::CONDUCTION) {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_copper()));
            }
            else {
                matplot::fill(x, y)->fill(true).color(matplot::to_array(settings->get_painter_color_insulation()));
            }
        }
    }

    paint_two_piece_set_margin(magnetic);
}

void Painter::paint_toroidal_winding_turns(MagneticWrapper magnetic) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto constants = Constants();
    CoilWrapper winding = magnetic.get_coil();
    auto wirePerWinding = winding.get_wires();

    auto processedDescription = magnetic.get_core().get_processed_description().value();

    auto mainColumn = magnetic.get_mutable_core().find_closest_column_by_coordinates({0, 0, 0});
    auto aux = get_image_size(magnetic);
    double imageWidth = aux[0];
    double imageHeight = aux[1];

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
        if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
            paint_round_wire(xCoordinate, yCoordinate, wire);
        }
        else {
            double turnAngle = turns[i].get_rotation().value();
            std::vector<double> turnCenter = {(_offsetForColorBar + imageWidth / 2 + xCoordinate), (imageHeight / 2 - yCoordinate)};
            paint_rectangular_wire(xCoordinate, yCoordinate, wire, turnAngle, turnCenter);
        }

        if (turns[i].get_additional_coordinates()) {
            auto additionalCoordinates = turns[i].get_additional_coordinates().value();

            for (auto additionalCoordinate : additionalCoordinates){
                double xAdditionalCoordinate = additionalCoordinate[0];
                double yAdditionalCoordinate = additionalCoordinate[1];
                if (wirePerWinding[windingIndex].get_type() == WireType::ROUND || wirePerWinding[windingIndex].get_type() == WireType::LITZ) {
                    paint_round_wire(xAdditionalCoordinate, yAdditionalCoordinate, wire);
                }
                else {
                    double turnAngle = turns[i].get_rotation().value();
                    std::vector<double> turnCenter = {(_offsetForColorBar + imageWidth / 2 + xAdditionalCoordinate), (imageHeight / 2 - yAdditionalCoordinate)};
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
            double circlePerimeter = std::numbers::pi * circleDiameter * _scale;

            auto currentMapIndex = uint_to_hex(_currentMapIndex);
            auto key = key_to_rgb_color(_currentMapIndex);
            _currentMapIndex++;

            double angleProportion = layers[i].get_dimensions()[1] / 360;
            std::string termination = angleProportion < 1? "butt" : "round";

            matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
            _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2)) + " " + std::to_string((_offsetForColorBar + imageWidth / 2) * _scale) + " " + std::to_string(imageHeight / 2 * _scale) + ")\" " + 
                                            R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
            _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_insulation(), 0, 16));
            
            if (layers[i].get_additional_coordinates()) {
                circleDiameter = (initialRadius - layers[i].get_additional_coordinates().value()[0][0]) * 2;
                matplot::ellipse(_offsetForColorBar - circleDiameter / 2, -circleDiameter / 2, circleDiameter, circleDiameter)->line_width(strokeWidth * _scale).color(matplot::to_array(currentMapIndex));
                _postProcessingChanges[key] = R"( transform="rotate( )" + std::to_string(-(layers[i].get_coordinates()[1] - layers[i].get_dimensions()[1] / 2)) + " " + std::to_string((_offsetForColorBar + imageWidth / 2) * _scale) + " " + std::to_string(imageHeight / 2 * _scale) + ")\" " + 
                                                R"(stroke-linecap=")" + termination + R"(" stroke-dashoffset="0" stroke-dasharray=")" + std::to_string(circlePerimeter * angleProportion) + " " + std::to_string(circlePerimeter * (1 - angleProportion)) + "\"";
                _postProcessingColors[key] = key_to_rgb_color(stoi(settings->get_painter_color_insulation(), 0, 16));
            }


        }
    }

    paint_toroidal_margin(magnetic);
}

void Painter::paint_waveform(Waveform waveform) {
    std::vector<double> x, y;
    if (waveform.get_time()) {
        x = waveform.get_time().value(); 
    }
    else {
        x = matplot::linspace(0, 2 * std::numbers::pi, waveform.get_data().size());
    }
    y = waveform.get_data();

    matplot::plot(x, y);
}

} // namespace OpenMagnetics
