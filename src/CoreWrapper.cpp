#include "CoreWrapper.h"
#include "WireWrapper.h"
#include "Resistivity.h"
#include "Reluctance.h"
#include "InitialPermeability.h"
#include "spline.h"

#include "Defaults.h"
#include "Constants.h"
#include "json.hpp"
#include "Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

using json = nlohmann::json;

namespace OpenMagnetics {


CoreWrapper::CoreWrapper(const json& j, bool includeMaterialData, bool includeProcessedDescription, bool includeGeometricalDescription) {
    _includeMaterialData = includeMaterialData;
    from_json(j, *this);
    
    if (includeProcessedDescription) {
        process_data();
        process_gap();
    }

    if (!get_geometrical_description() && includeGeometricalDescription) {
        auto geometricalDescription = create_geometrical_description();

        set_geometrical_description(geometricalDescription);
    }
}

CoreWrapper::CoreWrapper(const MagneticCore core) {
    set_functional_description(core.get_functional_description());

    if (core.get_geometrical_description()) {
        set_geometrical_description(core.get_geometrical_description());
    }
    if (core.get_processed_description()) {
        set_processed_description(core.get_processed_description());
    }
    if (core.get_distributors_info()) {
        set_distributors_info(core.get_distributors_info());
    }
    if (core.get_manufacturer_info()) {
        set_manufacturer_info(core.get_manufacturer_info());
    }
    if (core.get_name()) {
        set_name(core.get_name());
    }
}

double CoreWrapper::get_depth() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_depth();
}
double CoreWrapper::get_height() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_height();
}
double CoreWrapper::get_width() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_width();
}

double roundFloat(double value, int64_t decimals) {
    return round(value * pow(10, decimals)) / pow(10, decimals);
}

double interp(std::vector<std::pair<double, double>> data, double temperature) {
    double result;

    if (data.size() == 0) {
        throw std::runtime_error("Data cannot be empty");
    }
    else if (data.size() == 1) {
        result = data[0].second;
    }
    else if (data.size() == 2) {
        result = data[0].second -
        (data[0].first - temperature) *
            (data[0].second - data[1].second) /
            (data[0].first - data[1].first);
    }
    else {
        std::sort(data.begin(), data.end(), [](const std::pair<double, double>& b1, const std::pair<double, double>& b2) {
                return b1.first < b2.first;
        });
        size_t n = data.size();
        std::vector<double> x, y;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || data[i].first != x.back()) {
                x.push_back(data[i].first);
                y.push_back(data[i].second);
            }
        }
        tk::spline interp(x, y, tk::spline::cspline_hermite, true);
        result = interp(temperature);
    }
    return result;
}

template<int decimals> double roundFloat(double value) {
    if (value < 0)
        return floor(value * pow(10, decimals)) / pow(10, decimals);
    else
        return ceil(value * pow(10, decimals)) / pow(10, decimals);
}

CoreShape flatten_dimensions(CoreShape shape) {
    CoreShape flattenedShape(shape);
    std::map<std::string, Dimension> dimensions = shape.get_dimensions().value();
    std::map<std::string, Dimension> flattenedDimensions;
    for (auto& dimension : dimensions) {
        double value = resolve_dimensional_values(dimension.second);
        flattenedDimensions[dimension.first] = value;
    }
    flattenedShape.set_dimensions(flattenedDimensions);
    return flattenedShape;
}

std::map<std::string, double> flatten_dimensions(std::map<std::string, Dimension> dimensions) {
    std::map<std::string, double> flattenedDimensions;
    for (auto& dimension : dimensions) {
        double value = resolve_dimensional_values(dimension.second);
        flattenedDimensions[dimension.first] = value;
    }
    return flattenedDimensions;
}

void CorePiece::process() {
    
    process_winding_window();
    process_columns();
    process_extra_data();

    auto [c1, c2, minimumArea] = get_shape_constants();

    if (c1 <= 0 || c2 <= 0 || minimumArea <= 0) {
        throw std::runtime_error("Shape constants cannot be negative or 0");
    }

    EffectiveParameters pieceEffectiveParameters;
    pieceEffectiveParameters.set_effective_length(pow(c1, 2) / c2);
    pieceEffectiveParameters.set_effective_area(c1 / c2);
    pieceEffectiveParameters.set_effective_volume(pow(c1, 3) / pow(c2, 2));
    pieceEffectiveParameters.set_minimum_area(minimumArea);
    set_partial_effective_parameters(pieceEffectiveParameters);
}

class CorePieceE : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["C"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 -
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = dimensions["B"] - dimensions["D"];
        double q = dimensions["C"];
        double s = dimensions["F"] / 2;
        double p = (dimensions["A"] - dimensions["E"]) / 2;

        lengths.push_back(dimensions["D"]);
        lengths.push_back((dimensions["E"] - dimensions["F"]) / 2);
        lengths.push_back(dimensions["D"]);
        lengths.push_back(std::numbers::pi / 8 * (p + h));
        lengths.push_back(std::numbers::pi / 8 * (s + h));

        areas.push_back(2 * q * p);
        areas.push_back(2 * q * h);
        areas.push_back(2 * s * q);
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceEtd : public CorePieceE {
  public:
    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double tetha;
        double aperture;
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            tetha = asin(dimensions["C"] / dimensions["E"]);
            aperture = dimensions["E"] / 2 * cos(tetha);
        }
        else {
            if (dimensions["G"] > 0) {
                aperture = dimensions["G"] / 2;
                tetha = acos(aperture / (dimensions["E"] / 2));
            }
            else {
                tetha = asin(dimensions["C"] / dimensions["E"]);
                aperture = dimensions["E"] / 2 * cos(tetha);
            }
        }
        double segmentArea = pow(dimensions["E"] / 2, 2) / 2 * (2 * tetha - sin(2 * tetha));
        double area =
            dimensions["C"] * (dimensions["A"] / 2 - aperture) - segmentArea;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_minimum_width(roundFloat<6>(dimensions["A"] / 2 - dimensions["E"] / 2));
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_width(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_depth()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = dimensions["B"] - dimensions["D"];
        double q = dimensions["C"];
        double s = dimensions["F"] / 2;
        double s1 = 0.5959 * s;
        double p = get_lateral_leg_area() / dimensions["C"];

        lengths.push_back(dimensions["D"]);
        lengths.push_back((dimensions["E"] - dimensions["F"]) / 2);
        lengths.push_back(dimensions["D"]);
        lengths.push_back(std::numbers::pi / 8 * (p + h));
        lengths.push_back(std::numbers::pi / 8 * (2 * s1 + h));

        areas.push_back(2 * q * p);
        areas.push_back(2 * q * h);
        areas.push_back(std::numbers::pi * pow(s, 2));
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {


            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceEl : public CorePieceE {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::OBLONG);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F2"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2) )+
                          (dimensions["F2"] - dimensions["F"]) *
                              dimensions["F"]);
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 -
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double a = dimensions["A"];
        double b = dimensions["B"];
        double c = dimensions["C"];
        double d = dimensions["D"];
        double e = dimensions["E"];
        double f = dimensions["F"];
        double f2 = dimensions["F2"];
        double r = 0;
        if ((dimensions.find("R") != dimensions.end()) && (dimensions["R"] != 0)) {
            r = dimensions["R"];
        }

        double a21 = (b - d) * c;
        double a23 = (f2 - f + std::numbers::pi * f / 2) * (b - d);
        double a3 = 1. / 2 * (1. / 4 * std::numbers::pi * pow(f, 2) + (f2 - f) * f);

        lengths.push_back(d);
        lengths.push_back(e / 2 - f / 2);
        lengths.push_back(d);
        lengths.push_back(std::numbers::pi / 8 * (a / 2 - e / 2 + b - d));
        lengths.push_back(std::numbers::pi / 8 * (a3 / f2 + b - d));

        areas.push_back(1. / 2 * (a - e) * c - 4 * (pow(r, 2) - 1. / 4 * std::numbers::pi * pow(r, 2)));
        areas.push_back(1. / 2 * (c + f2 - f + std::numbers::pi * f / 2) * (b - d));
        areas.push_back(a3);
        areas.push_back((areas[0] + a21) / 2);
        areas.push_back((a23 + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i] / 2;
            c2 += lengths[i] / (2 * pow(areas[i], 2)) / 2;
        }
        auto minimumArea = 2 * (*min_element(areas.begin(), areas.end()));

        return {c1, c2, minimumArea};
    }
};

class CorePieceEfd : public CorePieceEl {
  public:
    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F2"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 -
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"] + std::max(0., dimensions["K"]));
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double a = dimensions["A"];
        double b = dimensions["B"];
        double c = dimensions["C"];
        double d = dimensions["D"];
        double e = dimensions["E"];
        double f = dimensions["F"];
        double f2 = dimensions["F2"];
        double k = dimensions["K"];
        double q = dimensions["q"];

        lengths.push_back(d);
        lengths.push_back((e - f) / 2);
        lengths.push_back(d);
        lengths.push_back(std::numbers::pi / 8 * ((a - e) / 2 + b - d));
        lengths.push_back(std::numbers::pi / 4 * (f / 4 + sqrt(pow((c - f2 - 2 * k) / 2, 2) + pow((b - d) / 2, 2))));

        areas.push_back(c * (a - e) / 2);
        areas.push_back(c * (b - d));
        areas.push_back((f * f2 - 2 * pow(q, 2)) / 2);
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i] / 2;
            c2 += lengths[i] / (2 * pow(areas[i], 2)) / 2;
        }
        auto minimumArea = 2 * (*min_element(areas.begin(), areas.end()));

        return {c1, c2, minimumArea};
    }
};

class CorePieceEr : public CorePieceEtd {};

class CorePiecePlanarEr : public CorePieceEtd {};

class CorePiecePlanarE : public CorePieceE {};

class CorePiecePlanarEl : public CorePieceEl {};

class CorePieceEc : public CorePieceEtd {
  public:
    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double tetha = asin(dimensions["C"] / dimensions["E"]);
        double aperture = dimensions["E"] / 2 * cos(tetha);
        double segmentArea = pow(dimensions["E"] / 2, 2) / 2 * (2 * tetha - sin(2 * tetha));
        double clipHoleArea = std::numbers::pi * pow(dimensions["s"], 2) / 2;
        double area = dimensions["C"] * (dimensions["A"] / 2 - aperture) -
                      segmentArea - clipHoleArea;
        return area;
    }
};

class CorePieceEq : public CorePieceEtd {};

class CorePieceEp : public CorePieceE {
  public:
    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());

        double baseArea;
        double windingArea;
        double apertureArea;
        double k;
        if ((dimensions.find("K") == dimensions.end()) || (dimensions["K"] == 0)) {
            k = dimensions["F"] / 2;
        }
        else {
            k = dimensions["K"];
        }
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            baseArea = dimensions["A"] * dimensions["C"];
            windingArea = k * dimensions["E"] +
                          1 / 2 * std::numbers::pi * pow(dimensions["E"] / 2, 2);
            apertureArea = 0;
        }
        else {
            double aperture = dimensions["G"] / 2;
            double tetha = asin(aperture / (dimensions["E"] / 2));
            double segmentArea = (pow(dimensions["E"] / 2, 2) / 2 * (2 * tetha - sin(2 * tetha))) / 2;
            double apertureMaximumDepth =
                dimensions["C"] - k - dimensions["E"] / 2 * cos(tetha);
            apertureArea = aperture * apertureMaximumDepth - segmentArea;
            baseArea = dimensions["A"] / 2 * dimensions["C"];
            windingArea = k * dimensions["E"] / 2 +
                          1 / 4 * std::numbers::pi * pow(dimensions["E"] / 2, 2);
        }
        double area = baseArea - windingArea - apertureArea;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            lateralColumn.set_depth(roundFloat<6>(dimensions["C"] - dimensions["E"] / 2 ) -
                              dimensions["K"]);
            lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
            lateralColumn.set_minimum_width(roundFloat<6>(dimensions["A"] / 2 - dimensions["E"] / 2));
            lateralColumn.set_width(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_depth()));
            lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
            lateralColumn.set_coordinates({0, 0, roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_depth() / 2)});
            windingWindows.push_back(lateralColumn);
        }
        else {
            lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
            lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
            lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
            lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
            lateralColumn.set_coordinates({
                roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
            lateralColumn.set_coordinates({
                roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
        }
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double h1 = 2 * dimensions["B"];
        double h2 = 2 * dimensions["D"];
        double d1 = dimensions["E"];
        double d2 = dimensions["F"];
        double a = dimensions["A"];
        double b = dimensions["C"];
        double k;
        if ((dimensions.find("K") == dimensions.end()) || (dimensions["K"] == 0)) {
            k = dimensions["F"] / 2;
        }
        else {
            k = dimensions["K"];
        }
        double pi = std::numbers::pi;
        double a1 = a * b - pi * pow(d1, 2) / 8 - d1 * k;
        double a3 = pi * pow(d2 / 2, 2) + (k - d2 / 2) * d2;
        double alpha = atan(dimensions["E"] / 2 / k);
        double gamma = sqrt(((pi - alpha) * pow(d1, 2) + 2 * a1) / (4 * (pi - alpha)));
        double l4 = pi / 2 * (gamma - d1 / 2 + (h1 - h2) / 4);
        double a4 = 1. / 2 * (a * b - pi / 8 * pow(d1, 2) - d1 * d2 / 2 + (pi - alpha) * d1 * (h1 / 2 - h2 / 2));
        double l5 = pi / 2 * (0.29289 * (d2 / 2 + k) / 2 + (h1 - h2) / 4);
        double a5 = pi / 2 * (pow((d2 / 2 + k), 2) / 4 + (d2 / 2 + k) / 2 * (h1 - h2));

        areas.push_back(a1);
        areas.push_back(a3);
        areas.push_back(a4);
        areas.push_back(a5);

        lengths_areas.push_back(h2 / a1);
        lengths_areas.push_back(2 / (pi - alpha) / (h1 - h2) * log(d1 / (d2 / 2 + k)));
        lengths_areas.push_back(h2 / a3);
        lengths_areas.push_back(l4 / a4);
        lengths_areas.push_back(l5 / a5);

        lengths_areas_2.push_back(h2 / pow(a1, 2));
        lengths_areas_2.push_back(4 * (d1 - (d2 / 2 + k)) / pow(pi - alpha, 2) / pow(h1 - h2, 2) / d1 / (d2 / 2 + k));
        lengths_areas_2.push_back(h2 / pow(a3, 2));
        lengths_areas_2.push_back(l4 / pow(a4, 2));
        lengths_areas_2.push_back(l5 / pow(a5, 2));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths_areas.size(); ++i) {
            c1 += lengths_areas[i] / 2;
            c2 += lengths_areas_2[i] / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceLp : public CorePieceEp {
  public:
    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }
};

class CorePieceEpx : public CorePieceEp {
  public:
    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::OBLONG);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]) / 2 + roundFloat<6>(dimensions["K"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2) )+
                          (dimensions["K"] - dimensions["F"] / 2) *
                              dimensions["F"]);
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            lateralColumn.set_depth(roundFloat<6>(dimensions["C"] - dimensions["E"] / 2 )-
                              dimensions["K"]);
            lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
            lateralColumn.set_minimum_width(roundFloat<6>(dimensions["A"] / 2 - dimensions["E"] / 2));
            lateralColumn.set_width(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_depth()));
            lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
            lateralColumn.set_coordinates({
                0, 0,
                roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_depth() / 2 -
                              (dimensions["K"] - dimensions["F"] / 2) / 2)});
            windingWindows.push_back(lateralColumn);
        }
        else {
            lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
            lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
            lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
            lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
            lateralColumn.set_coordinates({
                roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
            lateralColumn.set_coordinates({
                roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
        }
        set_columns(windingWindows);
    }
};

class CorePieceRm : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["E"]);
    }

    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());

        double d2 = dimensions["E"];
        double a = dimensions["J"];
        double e = dimensions["G"];
        double p = sqrt(2) * dimensions["J"] - dimensions["A"];
        double pi = std::numbers::pi;
        double alpha = pi / 2;
        double beta = alpha - asin(e / d2);

        double a1 = 1. / 2 * pow(a, 2) * (1 + tan(beta - pi / 4)) - beta / 2 * pow(d2, 2) - 1. / 2 * pow(p, 2);
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double d2 = dimensions["E"];
        double d3 = dimensions["F"];
        double d4 = dimensions["H"];
        double a = dimensions["J"];
        double c = dimensions["C"];
        double e = dimensions["G"];
        double h = dimensions["B"] - dimensions["D"];
        double p = sqrt(2) * dimensions["J"] - dimensions["A"];
        double b = 0;
        double pi = std::numbers::pi;
        double alpha = pi / 2;
        double gamma = pi / 2;
        double beta = alpha - asin(e / d2);
        double lmin = (dimensions["E"] - dimensions["F"]) / 2;
        double lmax;
        double a7;
        double a8 = alpha / 8 * (pow(d2, 2) - pow(d3, 2));
        if (familySubtype == "1") {
            lmax = sqrt(1. / 4 * (pow(d2, 2) + pow(d3, 2)) - 1. / 2 * d2 * d3 * cos(alpha - beta));
            a7 = 1. / 4 *
                 (beta / 2 * pow(d2, 2) + 1. / 2 * pow(e, 2) * tan(beta) - 1. / 2 * pow(e, 2) * tan(alpha - gamma / 2) -
                  pi / 4 * pow(d3, 2));
        }
        else if (familySubtype == "2") {
            lmax = sqrt(1. / 4 * (pow(d2, 2) + pow(d3, 2)) - 1. / 2 * d2 * d3 * cos(alpha - beta)) -
                   b / (2 * sin(gamma / 2));
            a7 = 1. / 4 *
                 (beta / 2 * pow(d2, 2) - pi / 4 * pow(d3, 2) +
                  1. / 2 * (pow(b, 2) - pow(e, 2)) * tan(alpha - gamma / 2) + 1. / 2 * pow(e, 2) * tan(beta));
        }
        else if (familySubtype == "3") {
            lmax = e / 2 + 1. / 2 * (1 - sin(gamma / 2)) * (d2 - c);
            a7 = 1. / 4 * (beta / 2 * pow(d2, 2) - pi / 4 * pow(d3, 2) + 1. / 2 * pow(c, 2) * tan(alpha - beta));
        }
        else if (familySubtype == "4") {
            lmax = sqrt(1. / 4 * (pow(d2, 2) + pow(d3, 2)) - 1. / 2 * d2 * d3 * cos(alpha - beta));
            a7 = 1. / 4 *
                 (beta / 2 * pow(d2, 2) + 1. / 2 * d2 * d3 * sin(alpha - beta) +
                  1. / 2 * pow(c - d3, 2) * tan(gamma / 2) - pi / 4 * pow(d3, 2));
        }
        else {
            lmax = 0;
            a7 = 0;
        }

        double f = (lmin + lmax) / (2 * lmin);
        double D = a7 / a8;

        double l1 = 2 * dimensions["D"];
        double a1 = 1. / 2 * pow(a, 2) * (1 + tan(beta - pi / 4)) - beta / 2 * pow(d2, 2) - 1. / 2 * pow(p, 2);

        double l3 = 2 * dimensions["D"];
        double a3 = pi / 4 * (pow(d3, 2) - pow(d4, 2));

        double l4 = pi / 4 * (h + a / 2 - d2 / 2);
        double a4 = 1. / 2 * (a1 + 2 * beta * d2 * h);
        double l5 = pi / 4 * (d3 + h - sqrt(1. / 2 * (pow(d3, 2) + pow(d4, 2))));
        double a5 = 1. / 2 * (pi / 4 * (pow(d3, 2) - pow(d4, 2)) + 2 * alpha * d3 * h);

        areas.push_back(a1);
        areas.push_back(a3);
        areas.push_back(a4);
        areas.push_back(a5);

        lengths_areas.push_back(l1 / a1);
        lengths_areas.push_back(log(d2 / d3) * f / (D * pi * h));
        lengths_areas.push_back(l3 / a3);
        lengths_areas.push_back(l4 / a4);
        lengths_areas.push_back(l5 / a5);

        lengths_areas_2.push_back(l1 / pow(a1, 2));
        lengths_areas_2.push_back((1 / d3 - 1 / d2) * f / pow(D * pi * h, 2));
        lengths_areas_2.push_back(l3 / pow(a3, 2));
        lengths_areas_2.push_back(l4 / pow(a4, 2));
        lengths_areas_2.push_back(l5 / pow(a5, 2));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths_areas.size(); ++i) {
            c1 += lengths_areas[i] / 2;
            c2 += lengths_areas_2[i] / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePiecePq : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());

        double A = dimensions["A"];
        double C = dimensions["C"];
        double E = dimensions["E"];
        double G = dimensions["G"];

        double beta = acos(G / E);
        double I = E * sin(beta);

        double a1 = C * (A - G) - beta * pow(E, 2) / 2 + 1. / 2 * G * I;
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_depth(dimensions["C"]);
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_minimum_width(roundFloat<6>(dimensions["A"] / 2 - dimensions["E"] / 2));
        lateralColumn.set_width(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_depth()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double A = dimensions["A"];
        double B = dimensions["B"];
        double C = dimensions["C"];
        double D = dimensions["D"];
        double E = dimensions["E"];
        double F = dimensions["F"];
        double G = dimensions["G"];
        double J;
        double L;
        if ((dimensions.find("J") == dimensions.end()) || (dimensions["J"] == 0)) {
            J = dimensions["F"] / 2; // Totally made up base on drawings
            L = F + (C - F) / 3; // Totally made up base on drawings
        }
        else {
            J = dimensions["J"];
            L = dimensions["L"];
        }

        double pi = std::numbers::pi;
        double beta = acos(G / E);
        double alpha = atan(L / J);
        double I = E * sin(beta);
        double a7 = 1. / 8 * (beta * pow(E, 2) - alpha * pow(F, 2) + G * L - J * I);
        double a8 = pi / 16 * (pow(E, 2) - pow(F, 2));
        double a9 = 2 * alpha * F * (B - D);
        double a10 = 2 * beta * E * (B - D);
        double lmin = (dimensions["E"] - dimensions["F"]) / 2;
        double lmax = sqrt(pow(E, 2) + pow(F, 2) - 2 * E * F * cos(alpha - beta)) / 2;
        double f = (lmin + lmax) / (2 * lmin);
        double K = a7 / a8;

        double l1 = 2 * D;
        double a1 = C * (A - G) - beta * pow(E, 2) / 2 + 1. / 2 * G * I;
        double a2 = pi * K * E * F * (B - D) / (E - F) * log(E / F);
        double l2 = f * E * F / (E - F) * pow(log(E / F), 2);

        double l3 = 2 * D;
        double a3 = pi / 4 * (pow(F, 2));

        double l4 = pi / 4 * ((B - D) + A / 2 - E / 2);
        double a4 = 1. / 2 * (a1 + a10);
        double l5 = pi / 4 * ((B - D) + (1 - 1. / sqrt(2)) * F);
        double a5 = 1. / 2 * (a3 + a9);

        areas.push_back(a1);
        areas.push_back(a3);
        areas.push_back(a2);
        areas.push_back(a4);
        areas.push_back(a5);

        lengths_areas.push_back(l1 / a1);
        lengths_areas.push_back(l2 / a2);
        lengths_areas.push_back(l3 / a3);
        lengths_areas.push_back(l4 / a4);
        lengths_areas.push_back(l5 / a5);

        lengths_areas_2.push_back(l1 / pow(a1, 2));
        lengths_areas_2.push_back(l2 / pow(a2, 2));
        lengths_areas_2.push_back(l3 / pow(a3, 2));
        lengths_areas_2.push_back(l4 / pow(a4, 2));
        lengths_areas_2.push_back(l5 / pow(a5, 2));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths_areas.size(); ++i) {
            c1 += lengths_areas[i] / 2;
            c2 += lengths_areas_2[i] / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePiecePm : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["E"]);
    }

    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());

        double d1 = dimensions["A"];
        double d2 = dimensions["E"];
        double f = dimensions["G"];
        double b = dimensions["b"];
        double t = dimensions["t"];
        double pi = std::numbers::pi;

        double alpha = pi / 2;
        double beta = alpha - asin(f / d2);

        double a1 = beta / 2 * (pow(d1, 2) - pow(d2, 2)) - 2 * b * t;
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();

        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;
        double pi = std::numbers::pi;

        if (dimensions.find("alpha") == dimensions.end()) {
            if (familySubtype == "1") {
                dimensions["alpha"] = 120.;
            }
            else {
                dimensions["alpha"] = 90.;
            }
        }

        double d1 = dimensions["A"];
        double h1 = 2 * dimensions["B"];
        double h2 = 2 * dimensions["D"];
        double d2 = dimensions["E"];
        double d3 = dimensions["F"];
        double f = dimensions["G"];
        double d4 = dimensions["H"];
        double gamma = dimensions["alpha"] / 180 * pi;
        double b = dimensions["b"];
        double t = dimensions["t"];

        double alpha = pi / 2;
        double beta = alpha - asin(f / d2);
        double lmin = (dimensions["E"] - dimensions["F"]) / 2;
        double lmax = sqrt(1. / 4 * (pow(d2, 2) + pow(d3, 2)) - 1. / 2 * d2 * d3 * cos(alpha - beta));
        double g = (lmin + lmax) / (2 * lmin);
        double a7 = beta / 8 * pow(d2, 2) + 1. / 8 * pow(f, 2) * tan(beta) -
                    1. / 8 * pow(f, 2) * tan(alpha - gamma / 2) - pi / 16 * pow(d3, 2);
        double a8 = alpha / 8 * (pow(d2, 2) - pow(d3, 2));
        double D = a7 / a8;

        double a1 = beta / 2 * (pow(d1, 2) - pow(d2, 2)) - 2 * b * t;
        double l1 = h2;

        double l3 = h2;
        double a3 = pi / 4 * (pow(d3, 2) - pow(d4, 2));

        double l4 = pi / 8 * (h1 - h2 + d1 - d2);
        double a4 = 1. / 2 * (a1 + 2 * beta * d2 * (h1 - h2));
        double l5 = pi / 4 * (d3 + h1 - h2 - sqrt(1. / 2 * (pow(d3, 2) + pow(d4, 2))));
        double a5 = pi / 8 * (pow(d3, 2) - pow(d4, 2)) + alpha * d3 * (h1 - h2);

        areas.push_back(a1);
        areas.push_back(a3);
        areas.push_back(a4);
        areas.push_back(a5);

        lengths_areas.push_back(l1 / a1);
        lengths_areas.push_back(log(d2 / d3) * g / (D * pi * (h1 - h2) / 2));
        lengths_areas.push_back(l3 / a3);
        lengths_areas.push_back(l4 / a4);
        lengths_areas.push_back(l5 / a5);

        lengths_areas_2.push_back(l1 / pow(a1, 2));
        lengths_areas_2.push_back((1 / d3 - 1 / d2) * g / pow(D * pi * (h1 - h2) / 2, 2));
        lengths_areas_2.push_back(l3 / pow(a3, 2));
        lengths_areas_2.push_back(l4 / pow(a4, 2));
        lengths_areas_2.push_back(l5 / pow(a5, 2));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths_areas.size(); ++i) {
            c1 += lengths_areas[i] / 2;
            c2 += lengths_areas_2[i] / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceP : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["A"]);
    }

    double get_lateral_leg_area() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();
        double pi = std::numbers::pi;
        double d1 = dimensions["A"];
        double d2 = dimensions["E"];
        double b = dimensions["G"];
        double tetha = asin(2 * b / (d1 + d2));
        double n;
        if (familySubtype == "1" || familySubtype == "2") {
            n = 2;
        }
        else {
            n = 0;
        }

        double a1 = 1. / 4 * (pi - n * tetha) * (pow(d1, 2) - pow(d2, 2));
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        mainColumn.set_width(roundFloat<6>(dimensions["F"]));
        mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat<6>(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat<6>(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat<6>(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat<6>(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;
        double pi = std::numbers::pi;

        double r4 = dimensions["A"] / 2;
        double r3 = dimensions["E"] / 2;
        double r2 = dimensions["F"] / 2;
        double r1 = dimensions["H"] / 2;
        double h = dimensions["B"] - dimensions["D"];
        double h2 = 2 * dimensions["D"];
        double b = dimensions["G"];

        double s1 = r2 - sqrt((pow(r1, 2) + pow(r2, 2)) / 2);
        double s2 = sqrt((pow(r3, 2) + pow(r4, 2)) / 2) - r3;
        double n;
        if (familySubtype == "1" || familySubtype == "2") {
            n = 2;
        }
        else {
            n = 0;
        }

        double k1 = n * b * (r4 - r3);
        double k2 = 1 / (1 - n * b / (2 * pi * r3));
        double k3 = 1 - n * b / (pi * (r3 + r4));

        double a1 = pi * (r4 - r3) * (r4 + r3) - k1;
        double l1 = h2;

        double a3 = pi * (r2 - r1) * (r2 + r1);
        double l3 = h2;

        double l4 = pi / 4 * (2 * s2 + h);
        double a4 = pi / 2 * (pow(r4, 2) - pow(r3, 2) + 2 * r3 * h) * k3;
        double l5 = pi / 4 * (2 * s1 + h);
        double a5 = pi / 2 * (pow(r2, 2) - pow(r1, 2) + 2 * r2 * h);

        areas.push_back(a1);
        areas.push_back(a3);
        areas.push_back(a4);
        areas.push_back(a5);

        lengths_areas.push_back(l1 / a1);
        lengths_areas.push_back(1. / (pi * h) * log(r3 / r2) * k2);
        lengths_areas.push_back(l3 / a3);
        lengths_areas.push_back(l4 / a4);
        lengths_areas.push_back(l5 / a5);

        lengths_areas_2.push_back(l1 / pow(a1, 2));
        lengths_areas_2.push_back(1 / (2 * pow(pi * h, 2)) * (r3 - r2) / (r3 * r2) * k2);
        lengths_areas_2.push_back(l3 / pow(a3, 2));
        lengths_areas_2.push_back(l4 / pow(a4, 2));
        lengths_areas_2.push_back(l5 / pow(a5, 2));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths_areas.size(); ++i) {
            c1 += lengths_areas[i] / 2;
            c2 += lengths_areas_2[i] / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceU : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(dimensions["F"]) == 0)) {
                windingWindowWidth = dimensions["A"] - dimensions["C"] -
                                     dimensions["H"];
            }
            else {
                windingWindowWidth = dimensions["A"] - dimensions["F"] -
                                     dimensions["H"];
            }
        }
        else {
            windingWindowWidth = dimensions["E"];
        }

        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width(windingWindowWidth);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({(dimensions["A"] - windingWindowWidth) / 2 + windingWindowWidth / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(dimensions["H"]) == 0)) {
            mainColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        }
        else {
            mainColumn.set_width(roundFloat<6>(dimensions["H"]));
        }
        mainColumn.set_depth(roundFloat<6>(dimensions["C"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        lateralColumn.set_width(mainColumn.get_width());
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat<6>((dimensions["A"] + dimensions["E"]) / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = dimensions["B"] - dimensions["D"];
        double q = dimensions["C"];
        double s;
        double p;
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(dimensions["H"]) == 0)) {
            s = (dimensions["A"] - dimensions["E"]) / 2;
            p = (dimensions["A"] - dimensions["E"]) / 2;
        }
        else {
            s = dimensions["H"];
            p = dimensions["A"] - dimensions["E"] -
                dimensions["H"];
        }

        lengths.push_back(2 * dimensions["D"]);
        lengths.push_back(2 * dimensions["E"]);
        lengths.push_back(2 * dimensions["D"]);
        lengths.push_back(std::numbers::pi / 4 * (p + h));
        lengths.push_back(std::numbers::pi / 4 * (s + h));

        areas.push_back(q * p);
        areas.push_back(q * h);
        areas.push_back(s * q);
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i] / 2;
            c2 += lengths[i] / pow(areas[i], 2) / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceUr : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(dimensions["F"]) == 0)) {
                windingWindowWidth = dimensions["A"] - dimensions["C"] -
                                     dimensions["H"];
            }
            else {
                windingWindowWidth = dimensions["A"] - dimensions["F"] -
                                     dimensions["H"];
            }
        }
        else {
            windingWindowWidth = dimensions["E"];
        }

        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width(windingWindowWidth);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({(dimensions["A"] - windingWindowWidth) / 2 + windingWindowWidth / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();

        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(dimensions["F"]) == 0)) {
                windingWindowWidth = dimensions["A"] - dimensions["C"] -
                                     dimensions["H"];
            }
            else {
                windingWindowWidth = dimensions["A"] - dimensions["F"] -
                                     dimensions["H"];
            }
        }
        else {
            windingWindowWidth = dimensions["E"];
        }

        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
        if (familySubtype == "1" || familySubtype == "2" || familySubtype == "4") {
            mainColumn.set_width(roundFloat<6>(dimensions["C"]));
            mainColumn.set_depth(roundFloat<6>(dimensions["C"]));
        }
        else {
            mainColumn.set_width(roundFloat<6>(dimensions["F"]));
            mainColumn.set_depth(roundFloat<6>(dimensions["F"]));
        }
        mainColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        if (familySubtype == "1" || familySubtype == "3") {
            lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
            lateralColumn.set_width(roundFloat<6>(dimensions["H"]));
            lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
            lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        }
        else {
            lateralColumn.set_shape(OpenMagnetics::ColumnShape::ROUND);
            lateralColumn.set_width(roundFloat<6>(dimensions["H"]));
            lateralColumn.set_depth(roundFloat<6>(dimensions["H"]));
            lateralColumn.set_area(roundFloat<6>(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        }
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_coordinates({roundFloat<6>((dimensions["A"] + windingWindowWidth) / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths;
        std::vector<double> areas;
        double pi = std::numbers::pi;

        double h = dimensions["B"] - dimensions["D"];
        double a1;
        double a3;
        double l4;
        double l5;
        double e;

        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(dimensions["F"]) == 0)) {
                e = dimensions["A"] - dimensions["C"] -
                    dimensions["H"];
            }
            else {
                e = dimensions["A"] - dimensions["F"] -
                    dimensions["H"];
            }
        }
        else {
            e = dimensions["E"];
        }

        if (familySubtype == "1") {
            a1 = dimensions["C"] * dimensions["H"];
            a3 = pi * pow(dimensions["C"] / 2, 2);
            l4 = std::numbers::pi / 4 * (dimensions["H"] + h);
            l5 = std::numbers::pi / 4 * (dimensions["C"] + h);
        }
        else if (familySubtype == "2") {
            a1 = pi * pow(dimensions["C"] / 2, 2);
            a3 = pi * pow(dimensions["C"] / 2, 2);
            l4 = std::numbers::pi / 4 * (dimensions["C"] + h);
            l5 = std::numbers::pi / 4 * (dimensions["C"] + h);
        }
        else if (familySubtype == "3") {
            a1 = dimensions["C"] * dimensions["H"];
            a3 = pi * pow(dimensions["F"] / 2, 2);
            l4 = std::numbers::pi / 4 * (dimensions["H"] + h);
            l5 = std::numbers::pi / 4 * (dimensions["F"] + h);
        }
        else if (familySubtype == "4") {
            a1 =
                pi * pow(dimensions["F"] / 2, 2) - pi * pow(dimensions["G"] / 2, 2);
            a3 =
                pi * pow(dimensions["F"] / 2, 2) - pi * pow(dimensions["G"] / 2, 2);
            l4 = std::numbers::pi / 4 * (dimensions["C"] + h);
            l5 = std::numbers::pi / 4 * (dimensions["C"] + h);
        }

        lengths.push_back(2 * dimensions["D"]);
        lengths.push_back(2 * e);
        lengths.push_back(2 * dimensions["D"]);
        lengths.push_back(l4);
        lengths.push_back(l5);

        areas.push_back(a1);
        areas.push_back(dimensions["C"] * h);
        areas.push_back(a3);
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i] / 2;
            c2 += lengths[i] / pow(areas[i], 2) / 2;
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceUt : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width(dimensions["E"]);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({(dimensions["A"] - dimensions["E"]) / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        mainColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(dimensions["H"]) == 0)) {
            mainColumn.set_width(roundFloat<6>((dimensions["A"] - dimensions["E"]) / 2));
        }
        else {
            mainColumn.set_width(roundFloat<6>(dimensions["H"]));
        }
        mainColumn.set_depth(roundFloat<6>(dimensions["C"]));
        mainColumn.set_height(roundFloat<6>(dimensions["D"]));
        mainColumn.set_area(roundFloat<6>(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(OpenMagnetics::ColumnType::LATERAL);
        lateralColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        lateralColumn.set_width(mainColumn.get_width());
        lateralColumn.set_depth(roundFloat<6>(dimensions["C"]));
        lateralColumn.set_height(roundFloat<6>(dimensions["D"]));
        lateralColumn.set_area(roundFloat<6>(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat<6>((dimensions["A"] + dimensions["E"]) / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = (dimensions["B"] - dimensions["D"]) / 2;
        double q = dimensions["C"];
        double s;
        double p;
        s = dimensions["A"] - dimensions["E"] - dimensions["F"];
        p = dimensions["F"];

        lengths.push_back(dimensions["D"]);
        lengths.push_back(2 * dimensions["E"]);
        lengths.push_back(dimensions["D"]);
        lengths.push_back(std::numbers::pi / 4 * (p + h));
        lengths.push_back(std::numbers::pi / 4 * (s + h));

        areas.push_back(q * p);
        areas.push_back(q * h);
        areas.push_back(s * q);
        areas.push_back((areas[0] + areas[1]) / 2);
        areas.push_back((areas[1] + areas[2]) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

class CorePieceT : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["A"]);
        set_depth(dimensions["C"]);
    }

    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        windingWindow.set_radial_height(dimensions["B"] / 2);
        windingWindow.set_angle(360);
        windingWindow.set_area(std::numbers::pi * pow(dimensions["B"] / 2, 2));
        windingWindow.set_coordinates(std::vector<double>({(dimensions["A"] - dimensions["B"]) / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        ColumnElement lateralColumn;
        double columnWidth = (dimensions["A"] - dimensions["B"]) / 2;
        mainColumn.set_type(OpenMagnetics::ColumnType::CENTRAL);
        mainColumn.set_shape(OpenMagnetics::ColumnShape::RECTANGULAR);
        mainColumn.set_width(columnWidth);
        mainColumn.set_depth(roundFloat<6>(dimensions["C"]));
        mainColumn.set_height(2 * std::numbers::pi * (dimensions["B"] / 2 + columnWidth / 2));
        mainColumn.set_area(roundFloat<6>(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<double> lengths;
        std::vector<double> areas;
        double columnWidth = (dimensions["A"] - dimensions["B"]) / 2;

        lengths.push_back(2 * std::numbers::pi * (dimensions["B"] / 2 + columnWidth / 2));

        areas.push_back(columnWidth * dimensions["C"]);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

std::shared_ptr<CorePiece> CorePiece::factory(CoreShape shape, bool process) {
    auto family = shape.get_family();
    if (family == CoreShapeFamily::E) {

        auto piece = std::make_shared<CorePieceE>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EC) {
        auto piece = std::make_shared<CorePieceEc>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EFD) {
        auto piece = std::make_shared<CorePieceEfd>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EL) {
        auto piece = std::make_shared<CorePieceEl>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EP) {
        auto piece = std::make_shared<CorePieceEp>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPX) {
        auto piece = std::make_shared<CorePieceEpx>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::LP) {
        auto piece = std::make_shared<CorePieceLp>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EQ) {
        auto piece = std::make_shared<CorePieceEq>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::ER) {
        auto piece = std::make_shared<CorePieceEr>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::ETD) {
        auto piece = std::make_shared<CorePieceEtd>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::P) {
        auto piece = std::make_shared<CorePieceP>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_E) {
        auto piece = std::make_shared<CorePiecePlanarE>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_EL) {
        auto piece = std::make_shared<CorePiecePlanarEl>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_ER) {
        auto piece = std::make_shared<CorePiecePlanarEr>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PM) {
        auto piece = std::make_shared<CorePiecePm>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PQ) {
        auto piece = std::make_shared<CorePiecePq>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::RM) {
        auto piece = std::make_shared<CorePieceRm>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::U) {
        auto piece = std::make_shared<CorePieceU>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::UR) {
        auto piece = std::make_shared<CorePieceUr>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::UT) {
        auto piece = std::make_shared<CorePieceUt>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::T) {
        auto piece = std::make_shared<CorePieceT>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else
        throw std::runtime_error("Unknown shape family: " + std::string{magic_enum::enum_name(family)} + ", available options are: {E, EC, EFD, EL, EP, EPX, LP, EQ, ER, "
                                 "ETD, P, PLANAR_E, PLANAR_EL, PLANAR_ER, PM, PQ, RM, U, UR, UT, T}");
}

inline void from_json(const json& j, OpenMagnetics::CorePiece& x) {
    x.set_columns(j.at("columns").get<std::vector<OpenMagnetics::ColumnElement>>());
    x.set_depth(j.at("depth").get<double>());
    x.set_height(j.at("height").get<double>());
    x.set_width(j.at("width").get<double>());
    x.set_shape(j.at("shape").get<OpenMagnetics::CoreShape>());
    x.set_winding_window(j.at("winding_window").get<OpenMagnetics::WindingWindowElement>());
    x.set_partial_effective_parameters(j.at("partial_effective_parameters").get<OpenMagnetics::EffectiveParameters>());
}

inline void to_json(json& j, const OpenMagnetics::CorePiece& x) {
    j = json::object();
    j["columns"] = x.get_columns();
    j["depth"] = x.get_depth();
    j["height"] = x.get_height();
    j["width"] = x.get_width();
    j["shape"] = x.get_winding_window();
    j["winding_window"] = x.get_shape();
    j["partial_effective_parameters"] = x.get_partial_effective_parameters();
}

std::optional<std::vector<CoreGeometricalDescriptionElement>> CoreWrapper::create_geometrical_description() {
    std::vector<CoreGeometricalDescriptionElement> geometricalDescription;
    auto numberStacks = *(get_functional_description().get_number_stacks());
    auto gapping = get_functional_description().get_gapping();

    auto corePiece =
        OpenMagnetics::CorePiece::factory(std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()));
    auto corePieceHeight = corePiece->get_height();
    auto corePieceDepth = corePiece->get_depth();

    CoreGeometricalDescriptionElement spacer;
    std::vector<Machining> machining;
    CoreGeometricalDescriptionElement piece;
    double currentDepth = roundFloat<6>((-corePieceDepth * (numberStacks - 1)) / 2);
    double spacerThickness = 0;

    for (auto& gap : gapping) {
        if (gap.get_type() == OpenMagnetics::GapType::ADDITIVE) {
            spacerThickness = gap.get_length();
        }
        else if (gap.get_type() == OpenMagnetics::GapType::SUBTRACTIVE) {
            Machining aux;
            aux.set_length(gap.get_length());
            if (!gap.get_coordinates()) {
                aux.set_coordinates({0, 0, 0});
            }
            else {
                aux.set_coordinates(gap.get_coordinates().value());
            }
            // else {  // Gap coordinates are centered in the gap, while machining coordinates start at the base of the
            // central column surface, therefore, we must correct them
            //     aux["coordinates"][1] = aux["coordinates"].at(1).get<double>() - gap.get_length() / 2;
            // }
            machining.push_back(aux);
        }
    }

    piece.set_material(get_functional_description().get_material());
    piece.set_shape(std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()));
    auto topPiece = piece;
    auto bottomPiece = piece;
    switch (get_functional_description().get_type()) {
        case OpenMagnetics::CoreType::TOROIDAL:
            piece.set_type(OpenMagnetics::CoreGeometricalDescriptionElementType::TOROIDAL);
            for (auto i = 0; i < numberStacks; ++i) {
                std::vector<double> coordinates = {0, 0, currentDepth};
                piece.set_coordinates(coordinates);
                piece.set_rotation(std::vector<double>({std::numbers::pi / 2, std::numbers::pi / 2, 0}));

                geometricalDescription.push_back(CoreGeometricalDescriptionElement(piece));

                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }
            break;
        case OpenMagnetics::CoreType::CLOSED_SHAPE:
            piece.set_type(OpenMagnetics::CoreGeometricalDescriptionElementType::CLOSED);
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat<6>(corePieceHeight);
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                piece.set_coordinates(coordinates);
                piece.set_rotation(std::vector<double>({0, 0, 0}));
                piece.set_machining(std::nullopt);
                geometricalDescription.push_back(CoreGeometricalDescriptionElement(piece));
                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }
            break;
        case OpenMagnetics::CoreType::TWO_PIECE_SET:
            topPiece.set_type(OpenMagnetics::CoreGeometricalDescriptionElementType::HALF_SET);
            bottomPiece.set_type(OpenMagnetics::CoreGeometricalDescriptionElementType::HALF_SET);
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat<6>(spacerThickness / 2);
                std::vector<Machining> topHalfMachining;
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                topPiece.set_coordinates(coordinates);
                topPiece.set_rotation(std::vector<double>({std::numbers::pi, std::numbers::pi, 0}));
                for (auto& operating : machining) {
                    if (operating.get_coordinates()[1] >= 0 &&
                        operating.get_coordinates()[1] < operating.get_length() / 2) {
                        Machining brokenDownOperating;
                        brokenDownOperating.set_coordinates(operating.get_coordinates());
                        brokenDownOperating.set_length(operating.get_length() / 2 + operating.get_coordinates()[1]);
                        brokenDownOperating.get_mutable_coordinates()[1] = brokenDownOperating.get_length() / 2;
                        topHalfMachining.push_back(brokenDownOperating);
                    }
                    else if (operating.get_coordinates()[1] < 0 &&
                             (operating.get_coordinates()[1] + operating.get_length() / 2) >
                                 0) {
                        Machining brokenDownOperating;
                        brokenDownOperating.set_coordinates(operating.get_coordinates());
                        brokenDownOperating.set_length(operating.get_length() / 2 + operating.get_coordinates()[1]);
                        brokenDownOperating.get_mutable_coordinates()[1] = brokenDownOperating.get_length() / 2;
                        topHalfMachining.push_back(brokenDownOperating);
                    }
                    else if (operating.get_coordinates()[1] > 0) {
                        topHalfMachining.push_back(operating);
                    }
                }

                if (topHalfMachining.size() > 0) {
                    topPiece.set_machining(topHalfMachining);
                }
                geometricalDescription.push_back(topPiece);

                std::vector<Machining> bottomHalfMachining;

                for (auto& operating : machining) {
                    if (operating.get_coordinates()[1] <= 0 &&
                        (-operating.get_coordinates()[1] < operating.get_length() / 2)) {
                        Machining brokenDownOperating;
                        brokenDownOperating.set_coordinates(operating.get_coordinates());
                        brokenDownOperating.set_length(operating.get_length() / 2 - operating.get_coordinates()[1]);
                        brokenDownOperating.get_mutable_coordinates()[1] = -brokenDownOperating.get_length() / 2;
                        bottomHalfMachining.push_back(brokenDownOperating);
                    }
                    else if (operating.get_coordinates()[1] > 0 &&
                             (operating.get_coordinates()[1] - operating.get_length() / 2) <
                                 0) {
                        Machining brokenDownOperating;
                        brokenDownOperating.set_coordinates(operating.get_coordinates());
                        brokenDownOperating.set_length(operating.get_length() / 2 - operating.get_coordinates()[1]);
                        brokenDownOperating.get_mutable_coordinates()[1] = -brokenDownOperating.get_length() / 2;
                        bottomHalfMachining.push_back(brokenDownOperating);
                    }
                    else if (operating.get_coordinates()[1] < 0) {
                        bottomHalfMachining.push_back(operating);
                    }
                }

                if (std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::UR ||
                    std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::U) {
                    bottomPiece.set_rotation(std::vector<double>({0, std::numbers::pi, 0}));
                }
                else {
                    bottomPiece.set_rotation(std::vector<double>({0, 0, 0}));
                }

                if (bottomHalfMachining.size() > 0) {
                    bottomPiece.set_machining(bottomHalfMachining);
                }
                currentHeight = -currentHeight;
                coordinates = {0, currentHeight, currentDepth};
                bottomPiece.set_coordinates(coordinates);
                geometricalDescription.push_back(bottomPiece);

                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }

            if (spacerThickness > 0) {
                for (auto& column : corePiece->get_columns()) {
                    auto shape_data = std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape());
                    if (column.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                        spacer.set_type(OpenMagnetics::CoreGeometricalDescriptionElementType::SPACER);
                        spacer.set_material("plastic");
                        // We cannot use directly column.get_width()
                        auto dimensions = flatten_dimensions(shape_data.get_dimensions().value());
                        double windingWindowWidth;
                        if (dimensions.find("E") == dimensions.end() ||
                            (roundFloat<6>(dimensions["E"]) == 0)) {
                            if (dimensions.find("F") == dimensions.end() ||
                                (roundFloat<6>(dimensions["F"]) == 0)) {
                                windingWindowWidth = dimensions["A"] -
                                                     dimensions["C"] -
                                                     dimensions["H"];
                            }
                            else {
                                windingWindowWidth = dimensions["A"] -
                                                     dimensions["F"] -
                                                     dimensions["H"];
                            }
                        }
                        else {
                            windingWindowWidth = dimensions["E"];
                        }
                        double minimum_column_width;
                        double minimum_column_depth;
                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_width = dimensions["A"];
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::U ||
                                 shape_data.get_family() == CoreShapeFamily::UR) {
                            if (dimensions.find("H") == dimensions.end() ||
                                (roundFloat<6>(dimensions["H"]) == 0)) {
                                minimum_column_width = (dimensions["A"] - windingWindowWidth) / 2;
                            }
                            else {
                                minimum_column_width = dimensions["H"];
                            }
                        }
                        else {
                            minimum_column_width = (dimensions["A"] - windingWindowWidth) / 2;
                        }

                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_depth = column.get_depth();
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::P ||
                                 shape_data.get_family() == CoreShapeFamily::PM) {
                            minimum_column_depth = dimensions["F"];
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::RM) {
                            if (dimensions.find("J") != dimensions.end() &&
                                (roundFloat<6>(dimensions["J"]) != 0)) {
                                minimum_column_depth =
                                    sqrt(2) * dimensions["J"] - dimensions["A"];
                            }
                            else if (dimensions.find("H") != dimensions.end() &&
                                     (roundFloat<6>(dimensions["H"]) != 0)) {
                                minimum_column_depth = dimensions["H"];
                            }
                            else {
                                minimum_column_depth = dimensions["F"];
                            }
                        }
                        else {
                            minimum_column_depth =
                                std::min(dimensions["C"], column.get_depth()) * numberStacks;
                        }
                        minimum_column_width *= (1 + constants.spacerProtudingPercentage);
                        minimum_column_depth *= (1 + constants.spacerProtudingPercentage);
                        double protuding_width = minimum_column_width * constants.spacerProtudingPercentage;
                        double protuding_depth = minimum_column_depth * constants.spacerProtudingPercentage;
                        spacer.set_dimensions(std::vector<double>({minimum_column_width, spacerThickness, minimum_column_depth}));
                        spacer.set_rotation(std::vector<double>({0, 0, 0}));
                        if (column.get_coordinates()[0] == 0) {
                            spacer.set_coordinates({0, column.get_coordinates()[1],
                                                         -dimensions["C"] / 2 +
                                                             minimum_column_depth / 2 - protuding_depth});
                        }
                        else if (column.get_coordinates()[0] < 0) {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR) {
                                spacer.set_coordinates(std::vector<double>({column.get_coordinates()[0] - column.get_width() / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                            else {
                                spacer.set_coordinates(std::vector<double>({-dimensions["A"] / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                        }
                        else {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR) {
                                spacer.set_coordinates(std::vector<double>({column.get_coordinates()[0] + column.get_width() / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                            else {
                                spacer.set_coordinates(std::vector<double>({dimensions["A"] / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                        }
                        geometricalDescription.push_back(spacer);
                    }
                }
            }
            break;
        case OpenMagnetics::CoreType::PIECE_AND_PLATE:
            // TODO add for toroPIECE_AND_PLATE
            break;
        default:
            throw std::runtime_error(
                "Unknown type of core, options are {TOROIDAL, TWO_PIECE_SET, PIECE_AND_PLATE, CLOSED_SHAPE}");
    }

    return geometricalDescription;
}

std::vector<ColumnElement> CoreWrapper::find_columns_by_type(ColumnType columnType) {
    std::vector<ColumnElement> foundColumns;
    auto processedDescription = get_processed_description().value();
    for (auto& column : processedDescription.get_columns()) {
        if (column.get_type() == columnType) {
            foundColumns.push_back(column);
        }
    }
    return foundColumns;
}

int CoreWrapper::find_closest_column_index_by_coordinates(std::vector<double> coordinates) {
    double closestDistance = std::numeric_limits<double>::infinity();
    int closestColumnIndex = -1;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (size_t index = 0; index < columns.size(); ++index) {
        double distance = 0;
        auto column_coordinates = columns[index].get_coordinates();
        for (size_t i = 0; i < column_coordinates.size(); ++i) {
            if (i != 1) { // We don't care about how high in the column the gap is, just about its projection, with are
                          // axis X and Z
                distance += fabs(column_coordinates[i] - coordinates[i]);
            }
        }
        if (distance < closestDistance) {
            closestColumnIndex = index;
            closestDistance = distance;
        }
    }
    return closestColumnIndex;
}

int CoreWrapper::find_exact_column_index_by_coordinates(std::vector<double> coordinates) {
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (size_t index = 0; index < columns.size(); ++index) {
        double distance = 0;
        auto column_coordinates = columns[index].get_coordinates();
        for (size_t i = 0; i < column_coordinates.size(); ++i) {
            if (i != 1) { // We don't care about how high in the column the gap is, just about its projection, with are
                          // axis X and Z
                distance += fabs(column_coordinates[i] - coordinates[i]);
            }
        }
        if (distance == 0) {
            return index;
        }
    }
    return -1;
}

ColumnElement CoreWrapper::find_closest_column_by_coordinates(std::vector<double> coordinates) {
    double closestDistance = std::numeric_limits<double>::infinity();
    ColumnElement closestColumn;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (auto& column : columns) {
        double distance = 0;
        auto column_coordinates = column.get_coordinates();
        for (size_t i = 0; i < column_coordinates.size(); ++i) {
            distance += fabs(column_coordinates[i] - coordinates[i]);
        }
        if (distance < closestDistance) {
            closestColumn = column;
            closestDistance = distance;
        }
    }
    return closestColumn;
}

std::vector<CoreGap> CoreWrapper::find_gaps_by_type(GapType gappingType) {
    std::vector<CoreGap> foundGaps;
    for (auto& gap : get_functional_description().get_gapping()) {
        if (gap.get_type() == gappingType) {
            foundGaps.push_back(gap);
        }
    }
    return foundGaps;
}

std::vector<CoreGap> CoreWrapper::find_gaps_by_column(ColumnElement column) {
    std::vector<CoreGap> foundGaps;
    double columnLeftLimit = column.get_coordinates()[0] - column.get_width() / 2;
    double columnRightLimit = column.get_coordinates()[0] + column.get_width() / 2;
    double columnFrontLimit = column.get_coordinates()[2] + column.get_depth() / 2;
    double columnBackLimit = column.get_coordinates()[2] - column.get_depth() / 2;

    bool gapIsComplete = true;
    for (auto& gap : get_functional_description().get_gapping()) {
        if (!gap.get_coordinates()) {
            gapIsComplete = false;
        }
    }

    if (!gapIsComplete) {
        process_gap();
    }

    for (auto& gap : get_functional_description().get_gapping()) {
        if (gap.get_coordinates().value()[0] >= columnLeftLimit &&
            gap.get_coordinates().value()[0] <= columnRightLimit &&
            gap.get_coordinates().value()[2] <= columnFrontLimit &&
            gap.get_coordinates().value()[2] >= columnBackLimit) {
            foundGaps.push_back(gap);
        }
    }
    return foundGaps;
}

void CoreWrapper::scale_to_stacks(int64_t numberStacks) {
    auto processedDescription = get_processed_description().value();
    processedDescription.get_mutable_effective_parameters().set_effective_area(
        processedDescription.get_effective_parameters().get_effective_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_minimum_area(
        processedDescription.get_effective_parameters().get_minimum_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_effective_volume(
        processedDescription.get_effective_parameters().get_effective_volume() * numberStacks);
    processedDescription.set_depth(processedDescription.get_depth() * numberStacks);
    for (auto& column : processedDescription.get_mutable_columns()) {
        column.set_area(column.get_area() * numberStacks);
        column.set_depth(column.get_depth() * numberStacks);
    }
    set_processed_description(processedDescription);

    auto functionalDescription = get_functional_description();
    auto gapping = functionalDescription.get_gapping();
    std::vector<CoreGap> scaledGapping;
    for (auto& gapInfo : gapping) {
        if (gapInfo.get_section_dimensions()) {
            auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
            gapSectionDimensions[1] *= numberStacks;
            gapInfo.set_section_dimensions(gapSectionDimensions);
        }
        if (gapInfo.get_area()) {
            gapInfo.set_area(gapInfo.get_area().value() * numberStacks);
        }
        scaledGapping.push_back(gapInfo);
    }
    functionalDescription.set_gapping(scaledGapping);
    set_functional_description(functionalDescription);
}

bool CoreWrapper::distribute_and_process_gap() {
    auto constants = Constants();
    std::vector<CoreGap> newGapping;
    auto gapping = get_functional_description().get_gapping();
    double centralColumnGapsHeightOffset;
    double distanceClosestNormalSurface;
    double coreChunkSizePlusGap = 0;
    auto nonResidualGaps = find_gaps_by_type(GapType::SUBTRACTIVE);
    auto additiveGaps = find_gaps_by_type(GapType::ADDITIVE);
    nonResidualGaps.insert(nonResidualGaps.end(), additiveGaps.begin(), additiveGaps.end());
    auto residualGaps = find_gaps_by_type(GapType::RESIDUAL);
    int numberNonResidualGaps = nonResidualGaps.size();
    int numberResidualGaps = residualGaps.size();
    int numberGaps = numberNonResidualGaps + numberResidualGaps;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    int numberColumns = columns.size();

    // Check if current information is valid
    // if (numberNonResidualGaps + numberResidualGaps < numberColumns && (numberNonResidualGaps + numberResidualGaps) >
    // 0) {
    //     numberNonResidualGaps = 0;
    //     numberResidualGaps = 0;
    //     // if (!(get_functional_description().get_type() == OpenMagnetics::CoreType::TOROIDAL && numberColumns == 1))
    //     {
    //     //     throw std::runtime_error("A TWO_PIECE_SET core cannot have less gaps than columns");
    //     // }
    // }

    if (numberNonResidualGaps == 0 && numberResidualGaps > numberColumns) {
        gapping = std::vector<CoreGap>(gapping.begin(), gapping.end() - (numberResidualGaps - numberColumns));
        get_mutable_functional_description().set_gapping(gapping);
        residualGaps = find_gaps_by_type(GapType::RESIDUAL);
        numberResidualGaps = residualGaps.size();
        numberGaps = numberNonResidualGaps + numberResidualGaps;
    }

    if (numberNonResidualGaps + numberResidualGaps == 0) {
        for (size_t i = 0; i < columns.size(); ++i) {
            CoreGap gap;
            gap.set_type(GapType::RESIDUAL);
            gap.set_length(constants.residualGap);
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - constants.residualGap / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - constants.residualGap / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else if (numberNonResidualGaps + numberResidualGaps < numberColumns) {
        for (size_t i = 0; i < columns.size(); ++i) {
            size_t gapIndex = i;
            if (i >= gapping.size()) {
                gapIndex = gapping.size() - 1;
            }
            CoreGap gap;
            gap.set_type(gapping[gapIndex].get_type());
            gap.set_length(gapping[gapIndex].get_length());
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - gapping[gapIndex].get_length() / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - gapping[gapIndex].get_length() / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else if ((numberResidualGaps == numberColumns || numberNonResidualGaps == numberColumns) &&
             (numberGaps == numberColumns)) {
        for (size_t i = 0; i < columns.size(); ++i) {
            CoreGap gap;
            gap.set_type(gapping[i].get_type());
            gap.set_length(gapping[i].get_length());
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - gapping[i].get_length() / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - gapping[i].get_length() / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else {
        auto lateralColumns = find_columns_by_type(ColumnType::LATERAL);
        auto centralColumns = find_columns_by_type(ColumnType::CENTRAL);

        ColumnElement windingColumn;
        std::vector<ColumnElement> returnColumns;
        if (centralColumns.size() == 0) {
            windingColumn = lateralColumns[0];
            returnColumns = std::vector<ColumnElement>(lateralColumns.begin() + 1, lateralColumns.end());
        }
        else {
            windingColumn = centralColumns[0];
            returnColumns = lateralColumns;
        }

        if (numberGaps == numberColumns) {
            if (windingColumn.get_height() > nonResidualGaps[0].get_length()) {
                centralColumnGapsHeightOffset = roundFloat<6>(nonResidualGaps[0].get_length() / 2);
            }
            else {
                centralColumnGapsHeightOffset = 0;
            }
            distanceClosestNormalSurface = roundFloat<6>(windingColumn.get_height() / 2 - nonResidualGaps[0].get_length() / 2);
        }
        else {
            coreChunkSizePlusGap = roundFloat<6>(windingColumn.get_height() / (nonResidualGaps.size() + 1));
            centralColumnGapsHeightOffset = roundFloat<6>(-coreChunkSizePlusGap * (nonResidualGaps.size() - 1) / 2);
            distanceClosestNormalSurface = roundFloat<6>(coreChunkSizePlusGap - nonResidualGaps[0].get_length() / 2);
        }

        for (size_t i = 0; i < nonResidualGaps.size(); ++i) {
            CoreGap gap;
            gap.set_type(nonResidualGaps[i].get_type());
            gap.set_length(nonResidualGaps[i].get_length());
            gap.set_coordinates(std::vector<double>({windingColumn.get_coordinates()[0],
                                      windingColumn.get_coordinates()[1] + centralColumnGapsHeightOffset,
                                      windingColumn.get_coordinates()[2]}));
            gap.set_shape(windingColumn.get_shape());
            if (distanceClosestNormalSurface < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", non residual gap of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(distanceClosestNormalSurface);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(windingColumn.get_area());
            gap.set_section_dimensions(std::vector<double>({windingColumn.get_width(), windingColumn.get_depth()}));
            newGapping.push_back(gap);

            centralColumnGapsHeightOffset += roundFloat<6>(windingColumn.get_height() / (nonResidualGaps.size() + 1));
            if (i < nonResidualGaps.size() / 2. - 1) {
                distanceClosestNormalSurface = roundFloat<6>(distanceClosestNormalSurface + coreChunkSizePlusGap);
            }
            else if (i > nonResidualGaps.size() / 2. - 1) {
                distanceClosestNormalSurface = roundFloat<6>(distanceClosestNormalSurface - coreChunkSizePlusGap);
            }
        }

        if (residualGaps.size() < returnColumns.size()) {
            for (size_t i = 0; i < returnColumns.size(); ++i) {
                CoreGap gap;
                gap.set_type(GapType::RESIDUAL);
                gap.set_length(constants.residualGap);
                gap.set_coordinates(returnColumns[i].get_coordinates());
                gap.set_shape(returnColumns[i].get_shape());
                if (returnColumns[i].get_height() / 2 - constants.residualGap / 2 < 0) {
                    return false;
                    // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", return column of index: " + std::to_string(i));

                }
                gap.set_distance_closest_normal_surface(returnColumns[i].get_height() / 2 - constants.residualGap / 2);
                gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
                gap.set_area(returnColumns[i].get_area());
                gap.set_section_dimensions(std::vector<double>({returnColumns[i].get_width(), returnColumns[i].get_depth()}));
                newGapping.push_back(gap);
            }
        }
        else {
            for (size_t i = 0; i < returnColumns.size(); ++i) {
                CoreGap gap;
                gap.set_type(residualGaps[i].get_type());
                gap.set_length(residualGaps[i].get_length());
                gap.set_coordinates(returnColumns[i].get_coordinates());
                gap.set_shape(returnColumns[i].get_shape());
                if (returnColumns[i].get_height() / 2 < 0) {
                    return false;
                    // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", return column of index: " + std::to_string(i));

                }
                gap.set_distance_closest_normal_surface(returnColumns[i].get_height() / 2);
                gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
                gap.set_area(returnColumns[i].get_area());
                gap.set_section_dimensions(std::vector<double>({returnColumns[i].get_width(), returnColumns[i].get_depth()}));
                newGapping.push_back(gap);
            }
        }
    }

    get_mutable_functional_description().set_gapping(newGapping);
    return true;
}

bool CoreWrapper::is_gapping_missaligned() {
    auto gapping = get_functional_description().get_gapping();
    for (size_t i = 0; i < gapping.size(); ++i) {
        if (!gapping[i].get_coordinates()) {
            return true;
        }
        auto columnIndex = find_exact_column_index_by_coordinates(*gapping[i].get_coordinates());
        if (columnIndex == -1) {
            return true;
        }
    }
    return false;
}

bool CoreWrapper::is_gap_processed() {
    for (auto& gap : get_gapping()) {
        if (!gap.get_coordinates()) {
            return false;
        }
    }
    return true;
}

bool CoreWrapper::process_gap() {
    std::vector<CoreGap> newGapping;
    auto gapping = get_functional_description().get_gapping();
    auto family = std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family();
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();

    if (family == CoreShapeFamily::T && gapping.size() > 0 ) {
        throw std::runtime_error("Toroids cannot be gapped: " + std::to_string(gapping[0].get_length()));
    }

    if (family != CoreShapeFamily::T) {
        if (gapping.size() == 0 || !gapping[0].get_coordinates() || is_gapping_missaligned()) {
            return distribute_and_process_gap();
        }

        for (size_t i = 0; i < gapping.size(); ++i) {
            CoreGap gap;
            auto columnIndex = find_closest_column_index_by_coordinates(*gapping[i].get_coordinates());

            gap.set_type(gapping[i].get_type());
            gap.set_length(gapping[i].get_length());
            gap.set_coordinates(gapping[i].get_coordinates());
            gap.set_shape(columns[columnIndex].get_shape());
            if (roundFloat<6>(columns[columnIndex].get_height() / 2 - fabs((*gapping[i].get_coordinates())[1]) - gapping[i].get_length() / 2) < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", gap of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(roundFloat<6>(columns[columnIndex].get_height() / 2 - fabs((*gapping[i].get_coordinates())[1]) - gapping[i].get_length() / 2));
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[columnIndex].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[columnIndex].get_width(), columns[columnIndex].get_depth()}));
            newGapping.push_back(gap);
        }
    }

    get_mutable_functional_description().set_gapping(newGapping);
    return true;
}

CoreMaterial CoreWrapper::resolve_material() {
    return resolve_material(get_functional_description().get_material());
}

CoreMaterial CoreWrapper::resolve_material(CoreMaterialDataOrNameUnion coreMaterial) {
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(coreMaterial)) {
        auto coreMaterialData = OpenMagnetics::find_core_material_by_name(std::get<std::string>(coreMaterial));
        coreMaterial = coreMaterialData;
        return coreMaterialData;
    }
    else {
        return std::get<CoreMaterial>(coreMaterial);
    }
}

CoreShape CoreWrapper::resolve_shape() {
    return resolve_shape(get_functional_description().get_shape());
}

CoreShape CoreWrapper::resolve_shape(CoreShapeDataOrNameUnion coreShape) {
    // If the shape is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(coreShape)) {
        auto coreShapeData = OpenMagnetics::find_core_shape_by_name(std::get<std::string>(coreShape));
        coreShape = coreShapeData;
        return coreShapeData;
    }
    else {
        return std::get<CoreShape>(coreShape);
    }
}

void CoreWrapper::process_data() {
    // If the shape is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(get_functional_description().get_shape())) {

        auto shape_data = OpenMagnetics::find_core_shape_by_name(std::get<std::string>(get_functional_description().get_shape()));
        shape_data.set_name(std::get<std::string>(get_functional_description().get_shape()));
        get_mutable_functional_description().set_shape(shape_data);
    }

    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid
    // long loading operatings)
    if (_includeMaterialData && std::holds_alternative<std::string>(get_functional_description().get_material()) &&
        std::get<std::string>(get_functional_description().get_material()) != "dummy") {
        auto material_data = OpenMagnetics::find_core_material_by_name(
            std::get<std::string>(get_functional_description().get_material()));
        get_mutable_functional_description().set_material(material_data);
    }

    auto corePiece =
        OpenMagnetics::CorePiece::factory(std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()));
    CoreProcessedDescription processedDescription;
    auto coreColumns = corePiece->get_columns();
    auto coreWindingWindow = corePiece->get_winding_window();
    auto coreEffectiveParameters = corePiece->get_partial_effective_parameters();
    switch (get_functional_description().get_type()) {
        case OpenMagnetics::CoreType::TOROIDAL:
        case OpenMagnetics::CoreType::CLOSED_SHAPE:
            processedDescription.set_columns(coreColumns);
            processedDescription.set_effective_parameters(corePiece->get_partial_effective_parameters());
            processedDescription.get_mutable_winding_windows().push_back(corePiece->get_winding_window());
            processedDescription.set_depth(corePiece->get_depth());
            processedDescription.set_height(corePiece->get_height());
            processedDescription.set_width(corePiece->get_width());
            break;

        case OpenMagnetics::CoreType::TWO_PIECE_SET:
            for (auto& column : coreColumns) {
                column.set_height(2 * column.get_height());
            }
            processedDescription.set_columns(coreColumns);

            coreEffectiveParameters.set_effective_length(2 * coreEffectiveParameters.get_effective_length());
            coreEffectiveParameters.set_effective_area(coreEffectiveParameters.get_effective_area());
            coreEffectiveParameters.set_effective_volume(2 * coreEffectiveParameters.get_effective_volume());
            coreEffectiveParameters.set_minimum_area(coreEffectiveParameters.get_minimum_area());
            processedDescription.set_effective_parameters(coreEffectiveParameters);

            coreWindingWindow.set_area(2 * coreWindingWindow.get_area().value());
            coreWindingWindow.set_height(2 * coreWindingWindow.get_height().value());
            processedDescription.get_mutable_winding_windows().push_back(coreWindingWindow);
            processedDescription.set_depth(corePiece->get_depth());
            processedDescription.set_height(corePiece->get_height() * 2);
            processedDescription.set_width(corePiece->get_width());
            break;
        default:
            throw std::runtime_error("Unknown type of core, available options are {TOROIDAL, TWO_PIECE_SET}");
    }
    set_processed_description(processedDescription);
    scale_to_stacks(*(get_functional_description().get_number_stacks()));
}

double CoreWrapper::get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, double temperature, bool proportion) {
    auto defaults = Defaults();
    auto saturationData = coreMaterial.get_saturation();
    std::vector<std::pair<double, double>> data;

    if (saturationData.size() == 0) {
        throw std::runtime_error("Missing saturation data in core material");
    }
    for (auto& datum : saturationData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_flux_density()});
    }

    double saturationMagneticFluxDensity = interp(data, temperature);

    if (proportion)
        return defaults.maximumProportionMagneticFluxDensitySaturation * saturationMagneticFluxDensity;
    else
        return saturationMagneticFluxDensity;
}

double CoreWrapper::get_magnetic_flux_density_saturation(double temperature, bool proportion) {
    auto coreMaterial = resolve_material();
    return get_magnetic_flux_density_saturation(coreMaterial, temperature, proportion);
}

double CoreWrapper::get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, bool proportion) {
    return get_magnetic_flux_density_saturation(coreMaterial, 25, proportion);
}
double CoreWrapper::get_magnetic_flux_density_saturation(bool proportion) {
    return get_magnetic_flux_density_saturation(25, proportion);
}

double CoreWrapper::get_magnetic_field_strength_saturation(CoreMaterial coreMaterial, double temperature) {
    auto saturationData = coreMaterial.get_saturation();
    std::vector<std::pair<double, double>> data;

    if (saturationData.size() == 0) {
        throw std::runtime_error("Missing saturation data in core material");
    }
    for (auto& datum : saturationData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_field()});
    }

    double saturationMagneticFieldStrength = interp(data, temperature);

    return saturationMagneticFieldStrength;
}

double CoreWrapper::get_magnetic_field_strength_saturation(double temperature) {
    auto coreMaterial = resolve_material();
    return get_magnetic_field_strength_saturation(coreMaterial, temperature);
}

double CoreWrapper::get_remanence(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing) {
    if (!coreMaterial.get_remanence()) {
        if (returnZeroIfMissing) {
            return 0;
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto remanenceData = coreMaterial.get_remanence().value();
    std::vector<std::pair<double, double>> data;

    if (remanenceData.size() == 0) {
        throw std::runtime_error("Missing remanence data in core material");
    }
    for (auto& datum : remanenceData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_flux_density()});
    }

    double remanence = interp(data, temperature);

    return remanence;
}

double CoreWrapper::get_remanence(double temperature) {
    auto coreMaterial = resolve_material();
    return get_remanence(coreMaterial, temperature);
}

double CoreWrapper::get_curie_temperature(CoreMaterial coreMaterial) {
    if (!coreMaterial.get_curie_temperature()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return coreMaterial.get_curie_temperature().value();
}

double CoreWrapper::get_curie_temperature() {
    auto coreMaterial = resolve_material();
    return get_curie_temperature(coreMaterial);
}

double CoreWrapper::get_coercive_force(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing) {
    if (!coreMaterial.get_coercive_force()) {
        if (returnZeroIfMissing) {
            return 0;
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto coerciveForceData = coreMaterial.get_coercive_force().value();
    std::vector<std::pair<double, double>> data;

    if (coerciveForceData.size() == 0) {
        throw std::runtime_error("Missing coercive force data in core material");
    }
    for (auto& datum : coerciveForceData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_field()});
    }

    double coerciveForce = interp(data, temperature);

    return coerciveForce;
}

double CoreWrapper::get_coercive_force(double temperature) {
    auto coreMaterial = resolve_material();
    return get_coercive_force(coreMaterial, temperature);
}


double CoreWrapper::get_initial_permeability(double temperature){
    auto coreMaterial = resolve_material();
    return get_initial_permeability(coreMaterial, temperature);
}

double CoreWrapper::get_initial_permeability(CoreMaterial coreMaterial, double temperature){
    InitialPermeability initialPermeability;
    auto initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);
    return initialPermeabilityValue;
}

double CoreWrapper::get_effective_permeability(double temperature){
    auto constants = Constants();
    auto reluctance = get_reluctance(temperature);
    if (!get_processed_description()) {
        process_data();
    }
    auto effectiveLength = get_processed_description().value().get_effective_parameters().get_effective_length();
    auto effectiveArea = get_processed_description().value().get_effective_parameters().get_effective_area();
    double effectivePermeability = effectiveLength / (reluctance * effectiveArea) / constants.vacuumPermeability;
    return effectivePermeability;
}

double CoreWrapper::get_reluctance(double temperature){
    auto coreMaterial = resolve_material();
    InitialPermeability initialPermeability;
    auto initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();

    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(*this, initialPermeabilityValue);
    double calculatedReluctance = magnetizingInductanceOutput.get_core_reluctance();
    return calculatedReluctance;
}

double CoreWrapper::get_resistivity(double temperature){
    auto coreMaterial = resolve_material();
    return get_resistivity(coreMaterial, temperature);
}

double CoreWrapper::get_resistivity(CoreMaterial coreMaterial, double temperature){
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::CORE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(coreMaterial, temperature);
    return resistivity;
}

double CoreWrapper::get_density(CoreMaterial coreMaterial) {
    if (!coreMaterial.get_density()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return coreMaterial.get_density().value();
}

double CoreWrapper::get_density() {
    auto coreMaterial = resolve_material();
    return get_density(coreMaterial);
}

std::vector<ColumnElement> CoreWrapper::get_columns() {
    if (get_processed_description()) {
        return get_processed_description().value().get_columns();
    }
    else {
        return std::vector<ColumnElement>();
    }
}

std::vector<WindingWindowElement> CoreWrapper::get_winding_windows() {
    if (get_processed_description()) {
        return get_processed_description().value().get_winding_windows();
    }
    else {
        return std::vector<WindingWindowElement>();
    }
}

CoreShapeFamily CoreWrapper::get_shape_family() {
    return resolve_shape().get_family();
}

std::string CoreWrapper::get_shape_name() {
    if (resolve_shape().get_name())
        return resolve_shape().get_name().value();
    else
        return "Custom";
}

int64_t CoreWrapper::get_number_stacks() {
    if (get_functional_description().get_number_stacks()) {
        return get_functional_description().get_number_stacks().value();
    }
    else {
        return 1;
    }
}

std::string CoreWrapper::get_material_name() {
    return resolve_material().get_name();
}

std::vector<CoreLossesMethodType> CoreWrapper::get_available_core_losses_methods(){
    std::vector<CoreLossesMethodType> methods;
    auto volumetricLossesMethodsVariants = resolve_material().get_volumetric_losses();
    for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
        auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
        for (auto& volumetricLossesMethod : volumetricLossesMethods) {
            if (std::holds_alternative<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod)) {
                auto methodData = std::get<OpenMagnetics::CoreLossesMethodData>(volumetricLossesMethod);
                if (std::find(methods.begin(), methods.end(), methodData.get_method()) == methods.end()) {
                    methods.push_back(methodData.get_method());
                }
            }
        }
    }
    return methods;
}

OpenMagnetics::CoreType CoreWrapper::get_type() {
    return get_functional_description().get_type();
}

bool fits_one_dimension(CoreProcessedDescription coreProcessedDescription, double dimension) {
    if (coreProcessedDescription.get_depth() <= dimension || coreProcessedDescription.get_height() <= dimension || coreProcessedDescription.get_width() <= dimension) {
        return true;
    }
    return false;
}

bool fits_two_dimensions(CoreProcessedDescription coreProcessedDescription, double firstDimension, double secondDimension) {
    if ((coreProcessedDescription.get_depth() <= firstDimension && (coreProcessedDescription.get_height() <= secondDimension || coreProcessedDescription.get_width() <= secondDimension)) ||
        (coreProcessedDescription.get_height() <= firstDimension && (coreProcessedDescription.get_depth() <= secondDimension || coreProcessedDescription.get_width() <= secondDimension)) ||
        (coreProcessedDescription.get_width() <= firstDimension && (coreProcessedDescription.get_height() <= secondDimension || coreProcessedDescription.get_depth() <= secondDimension))) {
        return true;
    }
    return false;
}

bool fits_three_dimensions(CoreProcessedDescription coreProcessedDescription, double firstDimension, double secondDimension, double thirdDimension) {
    if ((coreProcessedDescription.get_depth() <= firstDimension && coreProcessedDescription.get_height() <= secondDimension && coreProcessedDescription.get_width() <= thirdDimension) ||
        (coreProcessedDescription.get_depth() <= firstDimension && coreProcessedDescription.get_height() <= thirdDimension && coreProcessedDescription.get_width() <= secondDimension) ||
        (coreProcessedDescription.get_depth() <= secondDimension && coreProcessedDescription.get_height() <= firstDimension && coreProcessedDescription.get_width() <= thirdDimension) ||
        (coreProcessedDescription.get_depth() <= secondDimension && coreProcessedDescription.get_height() <= thirdDimension && coreProcessedDescription.get_width() <= firstDimension) ||
        (coreProcessedDescription.get_depth() <= thirdDimension && coreProcessedDescription.get_height() <= firstDimension && coreProcessedDescription.get_width() <= secondDimension) ||
        (coreProcessedDescription.get_depth() <= thirdDimension && coreProcessedDescription.get_height() <= secondDimension && coreProcessedDescription.get_width() <= firstDimension)){
        return true;
    }
    return false;
}

bool CoreWrapper::fits(MaximumDimensions maximumDimensions, bool allowRotation) {
    if (!get_processed_description()) {
        process_data();
    }

    auto coreProcessedDescription = get_processed_description().value();

    if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        return true;
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, depth);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, height);
        }
        else {
            return coreProcessedDescription.get_height() <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, width);
        }
        else {
            return coreProcessedDescription.get_width() <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, depth, height);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_height() <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, width, height);
        }
        else {
            return coreProcessedDescription.get_width() <= width && coreProcessedDescription.get_height() <= height;
        }
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, width, depth);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_width() <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_three_dimensions(coreProcessedDescription, depth, height, width);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_height() <= height && coreProcessedDescription.get_width() <= width;
        }
    }
    else {
        throw std::runtime_error("Not sure how this happened");
    }
}
} // namespace OpenMagnetics
