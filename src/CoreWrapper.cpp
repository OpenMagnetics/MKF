#include "CoreWrapper.h"
#include "WireWrapper.h"
#include <libInterpolate/Interpolate.hpp>

#include "Defaults.h"
#include "Constants.h"
#include "json.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {

double roundFloat(double value, size_t decimals) {
    return round(value * pow(10, decimals)) / pow(10, decimals);
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
        double value = resolve_dimensional_values<OpenMagnetics::DimensionalValues::NOMINAL>(dimension.second);
        flattenedDimensions[dimension.first] = value;
    }
    flattenedShape.set_dimensions(flattenedDimensions);
    return flattenedShape;
}

class CorePieceE : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_winding_window() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(jsonMainColumn["width"].get<double>() * jsonMainColumn["depth"].get<double>());
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] =
            roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 +
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 -
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double q = std::get<double>(dimensions["C"]);
        double s = std::get<double>(dimensions["F"]) / 2;
        double p = (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2;

        lengths.push_back(std::get<double>(dimensions["D"]));
        lengths.push_back((std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2);
        lengths.push_back(std::get<double>(dimensions["D"]));
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
        auto dimensions = get_shape().get_dimensions().value();
        double tetha;
        double aperture;
        if ((dimensions.find("G") == dimensions.end()) || (std::get<double>(dimensions["G"]) == 0)) {
            tetha = asin(std::get<double>(dimensions["C"]) / std::get<double>(dimensions["E"]));
            aperture = std::get<double>(dimensions["E"]) / 2 * cos(tetha);
        }
        else {
            if (std::get<double>(dimensions["G"]) > 0) {
                aperture = std::get<double>(dimensions["G"]) / 2;
                tetha = acos(aperture / (std::get<double>(dimensions["E"]) / 2));
            }
            else {
                tetha = asin(std::get<double>(dimensions["C"]) / std::get<double>(dimensions["E"]));
                aperture = std::get<double>(dimensions["E"]) / 2 * cos(tetha);
            }
        }
        double segmentArea = pow(std::get<double>(dimensions["E"]) / 2, 2) / 2 * (2 * tetha - sin(2 * tetha));
        double area =
            std::get<double>(dimensions["C"]) * (std::get<double>(dimensions["A"]) / 2 - aperture) - segmentArea;
        return area;
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>(std::get<double>(dimensions["A"]) / 2 - std::get<double>(dimensions["E"]) / 2);
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double q = std::get<double>(dimensions["C"]);
        double s = std::get<double>(dimensions["F"]) / 2;
        double s1 = 0.5959 * s;
        double p = get_lateral_leg_area() / std::get<double>(dimensions["C"]);

        lengths.push_back(std::get<double>(dimensions["D"]));
        lengths.push_back((std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2);
        lengths.push_back(std::get<double>(dimensions["D"]));
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::OBLONG;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F2"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2) +
                          (std::get<double>(dimensions["F2"]) - std::get<double>(dimensions["F"])) *
                              std::get<double>(dimensions["F"]));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] =
            roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 +
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 -
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double a = std::get<double>(dimensions["A"]);
        double b = std::get<double>(dimensions["B"]);
        double c = std::get<double>(dimensions["C"]);
        double d = std::get<double>(dimensions["D"]);
        double e = std::get<double>(dimensions["E"]);
        double f = std::get<double>(dimensions["F"]);
        double f2 = std::get<double>(dimensions["F2"]);
        double r = 0;
        if ((dimensions.find("R") != dimensions.end()) && (std::get<double>(dimensions["R"]) != 0)) {
            r = std::get<double>(dimensions["R"]);
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F2"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(jsonMainColumn["width"].get<double>() * jsonMainColumn["depth"].get<double>());
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] =
            roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 +
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 -
                          (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 4),
            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]) + std::max(0., std::get<double>(dimensions["K"])));
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double a = std::get<double>(dimensions["A"]);
        double b = std::get<double>(dimensions["B"]);
        double c = std::get<double>(dimensions["C"]);
        double d = std::get<double>(dimensions["D"]);
        double e = std::get<double>(dimensions["E"]);
        double f = std::get<double>(dimensions["F"]);
        double f2 = std::get<double>(dimensions["F2"]);
        double k = std::get<double>(dimensions["K"]);
        double q = std::get<double>(dimensions["q"]);

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
        auto dimensions = get_shape().get_dimensions().value();
        double tetha = asin(std::get<double>(dimensions["C"]) / std::get<double>(dimensions["E"]));
        double aperture = std::get<double>(dimensions["E"]) / 2 * cos(tetha);
        double segmentArea = pow(std::get<double>(dimensions["E"]) / 2, 2) / 2 * (2 * tetha - sin(2 * tetha));
        double clipHoleArea = std::numbers::pi * pow(std::get<double>(dimensions["s"]), 2) / 2;
        double area = std::get<double>(dimensions["C"]) * (std::get<double>(dimensions["A"]) / 2 - aperture) -
                      segmentArea - clipHoleArea;
        return area;
    }
};

class CorePieceEq : public CorePieceEtd {};

class CorePieceEp : public CorePieceE {
  public:
    double get_lateral_leg_area() {
        auto dimensions = get_shape().get_dimensions().value();

        double baseArea;
        double windingArea;
        double apertureArea;
        double k;
        if ((dimensions.find("K") == dimensions.end()) || (std::get<double>(dimensions["K"]) == 0)) {
            k = std::get<double>(dimensions["F"]) / 2;
        }
        else {
            k = std::get<double>(dimensions["K"]);
        }
        if ((dimensions.find("G") == dimensions.end()) || (std::get<double>(dimensions["G"]) == 0)) {
            baseArea = std::get<double>(dimensions["A"]) * std::get<double>(dimensions["C"]);
            windingArea = k * std::get<double>(dimensions["E"]) +
                          1 / 2 * std::numbers::pi * pow(std::get<double>(dimensions["E"]) / 2, 2);
            apertureArea = 0;
        }
        else {
            double aperture = std::get<double>(dimensions["G"]) / 2;
            double tetha = asin(aperture / (std::get<double>(dimensions["E"]) / 2));
            double segmentArea = (pow(std::get<double>(dimensions["E"]) / 2, 2) / 2 * (2 * tetha - sin(2 * tetha))) / 2;
            double apertureMaximumDepth =
                std::get<double>(dimensions["C"]) - k - std::get<double>(dimensions["E"]) / 2 * cos(tetha);
            apertureArea = aperture * apertureMaximumDepth - segmentArea;
            baseArea = std::get<double>(dimensions["A"]) / 2 * std::get<double>(dimensions["C"]);
            windingArea = k * std::get<double>(dimensions["E"]) / 2 +
                          1 / 4 * std::numbers::pi * pow(std::get<double>(dimensions["E"]) / 2, 2);
        }
        double area = baseArea - windingArea - apertureArea;
        return area;
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        if ((dimensions.find("G") == dimensions.end()) || (std::get<double>(dimensions["G"]) == 0)) {
            jsonLateralColumn["depth"] =
                roundFloat<6>(std::get<double>(dimensions["C"]) - std::get<double>(dimensions["E"]) / 2 -
                              std::get<double>(dimensions["K"]));
            jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
            jsonLateralColumn["width"] =
                roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["depth"].get<double>());
            jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
            jsonLateralColumn["coordinates"] = {
                0, 0,
                roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["depth"].get<double>() / 2)};
            jsonWindingWindows.push_back(jsonLateralColumn);
        }
        else {
            jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
            jsonLateralColumn["width"] =
                roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
            jsonLateralColumn["depth"] =
                roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
            jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
            jsonLateralColumn["coordinates"] = {
                roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0,
                0};
            jsonWindingWindows.push_back(jsonLateralColumn);
            jsonLateralColumn["coordinates"] = {
                roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0,
                0};
            jsonWindingWindows.push_back(jsonLateralColumn);
        }
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double h1 = 2 * std::get<double>(dimensions["B"]);
        double h2 = 2 * std::get<double>(dimensions["D"]);
        double d1 = std::get<double>(dimensions["E"]);
        double d2 = std::get<double>(dimensions["F"]);
        double a = std::get<double>(dimensions["A"]);
        double b = std::get<double>(dimensions["C"]);
        double k;
        if ((dimensions.find("K") == dimensions.end()) || (std::get<double>(dimensions["K"]) == 0)) {
            k = std::get<double>(dimensions["F"]) / 2;
        }
        else {
            k = std::get<double>(dimensions["K"]);
        }
        double pi = std::numbers::pi;
        double a1 = a * b - pi * pow(d1, 2) / 8 - d1 * k;
        double a3 = pi * pow(d2 / 2, 2) + (k - d2 / 2) * d2;
        double alpha = atan(std::get<double>(dimensions["E"]) / 2 / k);
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["depth"] =
            roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }
};

class CorePieceEpx : public CorePieceEp {
  public:
    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::OBLONG;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] =
            roundFloat<6>(std::get<double>(dimensions["F"])) / 2 + roundFloat<6>(std::get<double>(dimensions["K"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2) +
                          (std::get<double>(dimensions["K"]) - std::get<double>(dimensions["F"]) / 2) *
                              std::get<double>(dimensions["F"]));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        if ((dimensions.find("G") == dimensions.end()) || (std::get<double>(dimensions["G"]) == 0)) {
            jsonLateralColumn["depth"] =
                roundFloat<6>(std::get<double>(dimensions["C"]) - std::get<double>(dimensions["E"]) / 2 -
                              std::get<double>(dimensions["K"]));
            jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
            jsonLateralColumn["width"] =
                roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["depth"].get<double>());
            jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
            jsonLateralColumn["coordinates"] = {
                0, 0,
                roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["depth"].get<double>() / 2 -
                              (std::get<double>(dimensions["K"]) - std::get<double>(dimensions["F"]) / 2) / 2)};
            jsonWindingWindows.push_back(jsonLateralColumn);
        }
        else {
            jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
            jsonLateralColumn["width"] =
                roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
            jsonLateralColumn["depth"] =
                roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
            jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
            jsonLateralColumn["coordinates"] = {
                roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0,
                0};
            jsonWindingWindows.push_back(jsonLateralColumn);
            jsonLateralColumn["coordinates"] = {
                roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0,
                0};
            jsonWindingWindows.push_back(jsonLateralColumn);
        }
        set_columns(jsonWindingWindows);
    }
};

class CorePieceRm : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["E"]));
    }

    double get_lateral_leg_area() {
        auto dimensions = get_shape().get_dimensions().value();

        double d2 = std::get<double>(dimensions["E"]);
        double a = std::get<double>(dimensions["J"]);
        double e = std::get<double>(dimensions["G"]);
        double p = sqrt(2) * std::get<double>(dimensions["J"]) - std::get<double>(dimensions["A"]);
        double pi = std::numbers::pi;
        double alpha = pi / 2;
        double beta = alpha - asin(e / d2);

        double a1 = 1. / 2 * pow(a, 2) * (1 + tan(beta - pi / 4)) - beta / 2 * pow(d2, 2) - 1. / 2 * pow(p, 2);
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["depth"] =
            roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double d2 = std::get<double>(dimensions["E"]);
        double d3 = std::get<double>(dimensions["F"]);
        double d4 = std::get<double>(dimensions["H"]);
        double a = std::get<double>(dimensions["J"]);
        double c = std::get<double>(dimensions["C"]);
        double e = std::get<double>(dimensions["G"]);
        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double p = sqrt(2) * std::get<double>(dimensions["J"]) - std::get<double>(dimensions["A"]);
        double b = 0;
        double pi = std::numbers::pi;
        double alpha = pi / 2;
        double gamma = pi / 2;
        double beta = alpha - asin(e / d2);
        double lmin = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
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

        double l1 = 2 * std::get<double>(dimensions["D"]);
        double a1 = 1. / 2 * pow(a, 2) * (1 + tan(beta - pi / 4)) - beta / 2 * pow(d2, 2) - 1. / 2 * pow(p, 2);

        double l3 = 2 * std::get<double>(dimensions["D"]);
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
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_winding_window() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    double get_lateral_leg_area() {
        auto dimensions = get_shape().get_dimensions().value();

        double A = std::get<double>(dimensions["A"]);
        double C = std::get<double>(dimensions["C"]);
        double E = std::get<double>(dimensions["E"]);
        double G = std::get<double>(dimensions["G"]);

        double beta = acos(G / E);
        double I = E * sin(beta);

        double a1 = C * (A - G) - beta * pow(E, 2) / 2 + 1. / 2 * G * I;
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["depth"] = std::get<double>(dimensions["C"]);
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["width"] =
            roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;

        double A = std::get<double>(dimensions["A"]);
        double B = std::get<double>(dimensions["B"]);
        double C = std::get<double>(dimensions["C"]);
        double D = std::get<double>(dimensions["D"]);
        double E = std::get<double>(dimensions["E"]);
        double F = std::get<double>(dimensions["F"]);
        double G = std::get<double>(dimensions["G"]);
        double J;
        double L;
        if ((dimensions.find("J") == dimensions.end()) || (std::get<double>(dimensions["J"]) == 0)) {
            J = std::get<double>(dimensions["F"]) / 2; // Totally made up base on drawings
            L = F + (C - F) / 3; // Totally made up base on drawings
        }
        else {
            J = std::get<double>(dimensions["J"]);
            L = std::get<double>(dimensions["L"]);
        }

        double pi = std::numbers::pi;
        double beta = acos(G / E);
        double alpha = atan(L / J);
        double I = E * sin(beta);
        double a7 = 1. / 8 * (beta * pow(E, 2) - alpha * pow(F, 2) + G * L - J * I);
        double a8 = pi / 16 * (pow(E, 2) - pow(F, 2));
        double a9 = 2 * alpha * F * (B - D);
        double a10 = 2 * beta * E * (B - D);
        double lmin = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["E"]));
    }

    double get_lateral_leg_area() {
        auto dimensions = get_shape().get_dimensions().value();

        double d1 = std::get<double>(dimensions["A"]);
        double d2 = std::get<double>(dimensions["E"]);
        double f = std::get<double>(dimensions["G"]);
        double b = std::get<double>(dimensions["b"]);
        double t = std::get<double>(dimensions["t"]);
        double pi = std::numbers::pi;

        double alpha = pi / 2;
        double beta = alpha - asin(f / d2);

        double a1 = beta / 2 * (pow(d1, 2) - pow(d2, 2)) - 2 * b * t;
        double area = a1 / 2;
        return area;
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["depth"] =
            roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
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

        double d1 = std::get<double>(dimensions["A"]);
        double h1 = 2 * std::get<double>(dimensions["B"]);
        double h2 = 2 * std::get<double>(dimensions["D"]);
        double d2 = std::get<double>(dimensions["E"]);
        double d3 = std::get<double>(dimensions["F"]);
        double f = std::get<double>(dimensions["G"]);
        double d4 = std::get<double>(dimensions["H"]);
        double gamma = std::get<double>(dimensions["alpha"]) / 180 * pi;
        double b = std::get<double>(dimensions["b"]);
        double t = std::get<double>(dimensions["t"]);

        double alpha = pi / 2;
        double beta = alpha - asin(f / d2);
        double lmin = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = (std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"])) / 2;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {std::get<double>(dimensions["F"]) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["A"]));
    }

    double get_lateral_leg_area() {
        auto dimensions = get_shape().get_dimensions().value();
        auto familySubtype = *get_shape().get_family_subtype();
        double pi = std::numbers::pi;
        double d1 = std::get<double>(dimensions["A"]);
        double d2 = std::get<double>(dimensions["E"]);
        double b = std::get<double>(dimensions["G"]);
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::IRREGULAR;
        jsonLateralColumn["width"] =
            roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        jsonLateralColumn["area"] = roundFloat<6>(get_lateral_leg_area());
        jsonLateralColumn["depth"] =
            roundFloat<6>(jsonLateralColumn["area"].get<double>() / jsonLateralColumn["width"].get<double>());
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(std::get<double>(dimensions["E"]) / 2 + jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>(-std::get<double>(dimensions["E"]) / 2 - jsonLateralColumn["width"].get<double>() / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths_areas;
        std::vector<double> lengths_areas_2;
        std::vector<double> areas;
        double pi = std::numbers::pi;

        double r4 = std::get<double>(dimensions["A"]) / 2;
        double r3 = std::get<double>(dimensions["E"]) / 2;
        double r2 = std::get<double>(dimensions["F"]) / 2;
        double r1 = std::get<double>(dimensions["H"]) / 2;
        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double h2 = 2 * std::get<double>(dimensions["D"]);
        double b = std::get<double>(dimensions["G"]);

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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["E"])) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["F"])) == 0)) {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["C"]) -
                                     std::get<double>(dimensions["H"]);
            }
            else {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["F"]) -
                                     std::get<double>(dimensions["H"]);
            }
        }
        else {
            windingWindowWidth = std::get<double>(dimensions["E"]);
        }

        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = windingWindowWidth;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {(std::get<double>(dimensions["A"]) - windingWindowWidth) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["H"])) == 0)) {
            jsonMainColumn["width"] =
                roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        }
        else {
            jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["H"]));
        }
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(jsonMainColumn["width"].get<double>() * jsonMainColumn["depth"].get<double>());
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonLateralColumn["width"] = jsonMainColumn["width"].get<double>();
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] =
            roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>((std::get<double>(dimensions["A"]) + std::get<double>(dimensions["E"])) / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double q = std::get<double>(dimensions["C"]);
        double s;
        double p;
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["H"])) == 0)) {
            s = (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2;
            p = (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2;
        }
        else {
            s = std::get<double>(dimensions["H"]);
            p = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"]) -
                std::get<double>(dimensions["H"]);
        }

        lengths.push_back(2 * std::get<double>(dimensions["D"]));
        lengths.push_back(2 * std::get<double>(dimensions["E"]));
        lengths.push_back(2 * std::get<double>(dimensions["D"]));
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["E"])) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["F"])) == 0)) {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["C"]) -
                                     std::get<double>(dimensions["H"]);
            }
            else {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["F"]) -
                                     std::get<double>(dimensions["H"]);
            }
        }
        else {
            windingWindowWidth = std::get<double>(dimensions["E"]);
        }

        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = windingWindowWidth;
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {(std::get<double>(dimensions["A"]) - windingWindowWidth) / 2, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        auto familySubtype = *get_shape().get_family_subtype();

        double windingWindowWidth;
        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["E"])) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["F"])) == 0)) {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["C"]) -
                                     std::get<double>(dimensions["H"]);
            }
            else {
                windingWindowWidth = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["F"]) -
                                     std::get<double>(dimensions["H"]);
            }
        }
        else {
            windingWindowWidth = std::get<double>(dimensions["E"]);
        }

        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
        if (familySubtype == "1" || familySubtype == "2" || familySubtype == "4") {
            jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["C"]));
            jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        }
        else {
            jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["F"]));
            jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["F"]));
        }
        jsonMainColumn["area"] = roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        if (familySubtype == "1" || familySubtype == "3") {
            jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
            jsonLateralColumn["width"] = roundFloat<6>(std::get<double>(dimensions["H"]));
            jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
            jsonLateralColumn["area"] =
                roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        }
        else {
            jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::ROUND;
            jsonLateralColumn["width"] = roundFloat<6>(std::get<double>(dimensions["H"]));
            jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["H"]));
            jsonLateralColumn["area"] =
                roundFloat<6>(std::numbers::pi * pow(jsonMainColumn["width"].get<double>() / 2, 2));
        }
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["coordinates"] = {roundFloat<6>((std::get<double>(dimensions["A"]) + windingWindowWidth) / 2),
                                            0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        auto familySubtype = *get_shape().get_family_subtype();
        std::vector<double> lengths;
        std::vector<double> areas;
        double pi = std::numbers::pi;

        double h = std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"]);
        double a1;
        double a3;
        double l4;
        double l5;
        double e;

        if (dimensions.find("E") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["E"])) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["F"])) == 0)) {
                e = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["C"]) -
                    std::get<double>(dimensions["H"]);
            }
            else {
                e = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["F"]) -
                    std::get<double>(dimensions["H"]);
            }
        }
        else {
            e = std::get<double>(dimensions["E"]);
        }

        if (familySubtype == "1") {
            a1 = std::get<double>(dimensions["C"]) * std::get<double>(dimensions["H"]);
            a3 = pi * pow(std::get<double>(dimensions["C"]) / 2, 2);
            l4 = std::numbers::pi / 4 * (std::get<double>(dimensions["H"]) + h);
            l5 = std::numbers::pi / 4 * (std::get<double>(dimensions["C"]) + h);
        }
        else if (familySubtype == "2") {
            a1 = pi * pow(std::get<double>(dimensions["C"]) / 2, 2);
            a3 = pi * pow(std::get<double>(dimensions["C"]) / 2, 2);
            l4 = std::numbers::pi / 4 * (std::get<double>(dimensions["C"]) + h);
            l5 = std::numbers::pi / 4 * (std::get<double>(dimensions["C"]) + h);
        }
        else if (familySubtype == "3") {
            a1 = std::get<double>(dimensions["C"]) * std::get<double>(dimensions["H"]);
            a3 = pi * pow(std::get<double>(dimensions["F"]) / 2, 2);
            l4 = std::numbers::pi / 4 * (std::get<double>(dimensions["H"]) + h);
            l5 = std::numbers::pi / 4 * (std::get<double>(dimensions["F"]) + h);
        }
        else if (familySubtype == "4") {
            a1 =
                pi * pow(std::get<double>(dimensions["F"]) / 2, 2) - pi * pow(std::get<double>(dimensions["G"]) / 2, 2);
            a3 =
                pi * pow(std::get<double>(dimensions["F"]) / 2, 2) - pi * pow(std::get<double>(dimensions["G"]) / 2, 2);
            l4 = std::numbers::pi / 4 * (std::get<double>(dimensions["C"]) + h);
            l5 = std::numbers::pi / 4 * (std::get<double>(dimensions["C"]) + h);
        }

        lengths.push_back(2 * std::get<double>(dimensions["D"]));
        lengths.push_back(2 * e);
        lengths.push_back(2 * std::get<double>(dimensions["D"]));
        lengths.push_back(l4);
        lengths.push_back(l5);

        areas.push_back(a1);
        areas.push_back(std::get<double>(dimensions["C"]) * h);
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
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["height"] = std::get<double>(dimensions["D"]);
        jsonWindingWindow["width"] = std::get<double>(dimensions["E"]);
        jsonWindingWindow["area"] =
            jsonWindingWindow["height"].get<double>() * jsonWindingWindow["width"].get<double>();
        jsonWindingWindow["coordinates"] = {(std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2,
                                            0};
        set_winding_window(jsonWindingWindow);
    }

    void process_extra_data() {
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["B"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        if (dimensions.find("H") == dimensions.end() || (roundFloat<6>(std::get<double>(dimensions["H"])) == 0)) {
            jsonMainColumn["width"] =
                roundFloat<6>((std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"])) / 2);
        }
        else {
            jsonMainColumn["width"] = roundFloat<6>(std::get<double>(dimensions["H"]));
        }
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonMainColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonMainColumn["area"] =
            roundFloat<6>(jsonMainColumn["width"].get<double>() * jsonMainColumn["depth"].get<double>());
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        jsonLateralColumn["type"] = OpenMagnetics::ColumnType::LATERAL;
        jsonLateralColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonLateralColumn["width"] = jsonMainColumn["width"].get<double>();
        jsonLateralColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonLateralColumn["height"] = roundFloat<6>(std::get<double>(dimensions["D"]));
        jsonLateralColumn["area"] =
            roundFloat<6>(jsonLateralColumn["width"].get<double>() * jsonLateralColumn["depth"].get<double>());
        jsonLateralColumn["coordinates"] = {
            roundFloat<6>((std::get<double>(dimensions["A"]) + std::get<double>(dimensions["E"])) / 2), 0, 0};
        jsonWindingWindows.push_back(jsonLateralColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;

        double h = (std::get<double>(dimensions["B"]) - std::get<double>(dimensions["D"])) / 2;
        double q = std::get<double>(dimensions["C"]);
        double s;
        double p;
        s = std::get<double>(dimensions["A"]) - std::get<double>(dimensions["E"]) - std::get<double>(dimensions["F"]);
        p = std::get<double>(dimensions["F"]);

        lengths.push_back(std::get<double>(dimensions["D"]));
        lengths.push_back(2 * std::get<double>(dimensions["E"]));
        lengths.push_back(std::get<double>(dimensions["D"]));
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
        auto dimensions = get_shape().get_dimensions().value();
        set_width(std::get<double>(dimensions["A"]));
        set_height(std::get<double>(dimensions["A"]));
        set_depth(std::get<double>(dimensions["C"]));
    }

    void process_winding_window() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindow;
        jsonWindingWindow["radialHeight"] = std::get<double>(dimensions["B"]) / 2;
        jsonWindingWindow["angle"] = 2 * std::numbers::pi;
        jsonWindingWindow["area"] = std::numbers::pi * pow(std::get<double>(dimensions["B"]) / 2, 2);
        jsonWindingWindow["coordinates"] = {(std::get<double>(dimensions["A"]) - std::get<double>(dimensions["B"])) / 4, 0};
        set_winding_window(jsonWindingWindow);
    }

    void process_columns() {
        auto dimensions = get_shape().get_dimensions().value();
        json jsonWindingWindows = json::array();
        json jsonMainColumn;
        json jsonLateralColumn;
        double columnWidth = (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["B"])) / 2;
        jsonMainColumn["type"] = OpenMagnetics::ColumnType::CENTRAL;
        jsonMainColumn["shape"] = OpenMagnetics::ColumnShape::RECTANGULAR;
        jsonMainColumn["width"] = columnWidth;
        jsonMainColumn["depth"] = roundFloat<6>(std::get<double>(dimensions["C"]));
        jsonMainColumn["height"] = 2 * std::numbers::pi * (std::get<double>(dimensions["B"]) / 2 + columnWidth / 2);
        jsonMainColumn["area"] =
            roundFloat<6>(jsonMainColumn["width"].get<double>() * jsonMainColumn["depth"].get<double>());
        jsonMainColumn["coordinates"] = {0, 0, 0};
        jsonWindingWindows.push_back(jsonMainColumn);
        set_columns(jsonWindingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = get_shape().get_dimensions().value();
        std::vector<double> lengths;
        std::vector<double> areas;
        double columnWidth = (std::get<double>(dimensions["A"]) - std::get<double>(dimensions["B"])) / 2;

        lengths.push_back(2 * std::numbers::pi * (std::get<double>(dimensions["B"]) / 2 + columnWidth / 2));

        areas.push_back(columnWidth * std::get<double>(dimensions["C"]));

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        auto minimumArea = *min_element(areas.begin(), areas.end());

        return {c1, c2, minimumArea};
    }
};

std::shared_ptr<CorePiece> CorePiece::factory(CoreShape shape) {
    auto family = shape.get_family();
    if (family == CoreShapeFamily::E) {
        std::shared_ptr<CorePiece> piece(new CorePieceE);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EC) {
        std::shared_ptr<CorePiece> piece(new CorePieceEc);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EFD) {
        std::shared_ptr<CorePiece> piece(new CorePieceEfd);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EL) {
        std::shared_ptr<CorePiece> piece(new CorePieceEl);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EP) {
        std::shared_ptr<CorePiece> piece(new CorePieceEp);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPX) {
        std::shared_ptr<CorePiece> piece(new CorePieceEpx);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::LP) {
        std::shared_ptr<CorePiece> piece(new CorePieceLp);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EQ) {
        std::shared_ptr<CorePiece> piece(new CorePieceEq);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::ER) {
        std::shared_ptr<CorePiece> piece(new CorePieceEr);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::ETD) {
        std::shared_ptr<CorePiece> piece(new CorePieceEtd);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::P) {
        std::shared_ptr<CorePiece> piece(new CorePieceP);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_E) {
        std::shared_ptr<CorePiece> piece(new CorePiecePlanarE);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_EL) {
        std::shared_ptr<CorePiece> piece(new CorePiecePlanarEl);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PLANAR_ER) {
        std::shared_ptr<CorePiece> piece(new CorePiecePlanarEr);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PM) {
        std::shared_ptr<CorePiece> piece(new CorePiecePm);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PQ) {
        std::shared_ptr<CorePiece> piece(new CorePiecePq);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::RM) {
        std::shared_ptr<CorePiece> piece(new CorePieceRm);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::U) {
        std::shared_ptr<CorePiece> piece(new CorePieceU);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::UR) {
        std::shared_ptr<CorePiece> piece(new CorePieceUr);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::UT) {
        std::shared_ptr<CorePiece> piece(new CorePieceUt);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::T) {
        std::shared_ptr<CorePiece> piece(new CorePieceT);
        piece->set_shape(shape);
        piece->process();
        return piece;
    }
    else
        throw std::runtime_error("Unknown shape family, available options are: {E, EC, EFD, EL, EP, EPX, LP, EQ, ER, "
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

    json jsonSpacer;
    json jsonMachining = json::array();
    json jsonGeometricalDescription;
    double currentDepth = roundFloat<6>((-corePieceDepth * (numberStacks - 1)) / 2);
    double spacerThickness = 0;

    for (auto& gap : gapping) {
        if (gap.get_type() == OpenMagnetics::GapType::ADDITIVE) {
            spacerThickness = gap.get_length();
        }
        else if (gap.get_type() == OpenMagnetics::GapType::SUBTRACTIVE) {
            json aux;
            aux["length"] = gap.get_length();
            aux["coordinates"] = gap.get_coordinates();
            if (aux["coordinates"][0] == nullptr) {
                aux["coordinates"] = {0, 0, 0};
            }
            // else {  // Gap coordinates are centered in the gap, while machining coordinates start at the base of the
            // central column surface, therefore, we must correct them
            //     aux["coordinates"][1] = aux["coordinates"].at(1).get<double>() - gap.get_length() / 2;
            // }
            jsonMachining.push_back(aux);
        }
    }

    jsonGeometricalDescription["material"] = get_functional_description().get_material();
    jsonGeometricalDescription["shape"] = std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape());
    switch (get_functional_description().get_type()) {
        case OpenMagnetics::CoreType::TOROIDAL:
            jsonGeometricalDescription["type"] = OpenMagnetics::CoreGeometricalDescriptionElementType::TOROIDAL;
            for (auto i = 0; i < numberStacks; ++i) {
                std::vector<double> coordinates = {0, 0, currentDepth};
                jsonGeometricalDescription["coordinates"] = coordinates;
                jsonGeometricalDescription["rotation"] = {std::numbers::pi / 2, std::numbers::pi / 2, 0};

                geometricalDescription.push_back(CoreGeometricalDescriptionElement(jsonGeometricalDescription));

                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }
            break;
        case OpenMagnetics::CoreType::CLOSED_SHAPE:
            jsonGeometricalDescription["type"] = OpenMagnetics::CoreGeometricalDescriptionElementType::CLOSED;
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat<6>(corePieceHeight);
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                jsonGeometricalDescription["coordinates"] = coordinates;
                jsonGeometricalDescription["rotation"] = {0, 0, 0};
                if (jsonMachining.size() > 0) {
                    jsonGeometricalDescription["machining"] = jsonMachining;
                }
                geometricalDescription.push_back(CoreGeometricalDescriptionElement(jsonGeometricalDescription));

                if (jsonGeometricalDescription.find("machining") != jsonGeometricalDescription.end()) {
                    jsonGeometricalDescription.erase("machining");
                }
                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }
            break;
        case OpenMagnetics::CoreType::TWO_PIECE_SET:
            jsonGeometricalDescription["type"] = OpenMagnetics::CoreGeometricalDescriptionElementType::HALF_SET;
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat<6>(spacerThickness / 2);
                json topHalfMachining = json::array();
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                jsonGeometricalDescription["coordinates"] = coordinates;
                jsonGeometricalDescription["rotation"] = {std::numbers::pi, std::numbers::pi, 0};
                for (auto& operation : jsonMachining) {
                    if (operation["coordinates"].at(1).get<double>() >= 0 &&
                        operation["coordinates"].at(1).get<double>() < operation["length"].get<double>() / 2) {
                        json brokenDownOperation;
                        brokenDownOperation["coordinates"] = operation["coordinates"];
                        brokenDownOperation["length"] =
                            operation["length"].get<double>() / 2 + operation["coordinates"].at(1).get<double>();
                        brokenDownOperation["coordinates"][1] = brokenDownOperation["length"].get<double>() / 2;
                        topHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operation["coordinates"].at(1).get<double>() < 0 &&
                             (operation["coordinates"].at(1).get<double>() + operation["length"].get<double>() / 2) >
                                 0) {
                        json brokenDownOperation;
                        brokenDownOperation["coordinates"] = operation["coordinates"];
                        brokenDownOperation["length"] =
                            operation["length"].get<double>() / 2 + operation["coordinates"].at(1).get<double>();
                        brokenDownOperation["coordinates"][1] = brokenDownOperation["length"].get<double>() / 2;
                        topHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operation["coordinates"].at(1).get<double>() > 0) {
                        topHalfMachining.push_back(operation);
                    }
                }

                if (topHalfMachining.size() > 0) {
                    jsonGeometricalDescription["machining"] = topHalfMachining;
                }
                geometricalDescription.push_back(CoreGeometricalDescriptionElement(jsonGeometricalDescription));

                json bottomHalfMachining = json::array();

                for (auto& operation : jsonMachining) {
                    if (operation["coordinates"].at(1).get<double>() <= 0 &&
                        (-operation["coordinates"].at(1).get<double>() < operation["length"].get<double>() / 2)) {
                        json brokenDownOperation;
                        brokenDownOperation["coordinates"] = operation["coordinates"];
                        brokenDownOperation["length"] =
                            operation["length"].get<double>() / 2 - operation["coordinates"].at(1).get<double>();
                        brokenDownOperation["coordinates"][1] = -brokenDownOperation["length"].get<double>() / 2;
                        bottomHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operation["coordinates"].at(1).get<double>() > 0 &&
                             (operation["coordinates"].at(1).get<double>() - operation["length"].get<double>() / 2) <
                                 0) {
                        json brokenDownOperation;
                        brokenDownOperation["coordinates"] = operation["coordinates"];
                        brokenDownOperation["length"] =
                            operation["length"].get<double>() / 2 - operation["coordinates"].at(1).get<double>();
                        brokenDownOperation["coordinates"][1] = -brokenDownOperation["length"].get<double>() / 2;
                        bottomHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operation["coordinates"].at(1).get<double>() < 0) {
                        bottomHalfMachining.push_back(operation);
                    }
                }

                if (std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::UR ||
                    std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::U) {
                    jsonGeometricalDescription["rotation"] = {0, std::numbers::pi, 0};
                }
                else {
                    jsonGeometricalDescription["rotation"] = {0, 0, 0};
                }

                if (bottomHalfMachining.size() > 0) {
                    jsonGeometricalDescription["machining"] = bottomHalfMachining;
                }
                currentHeight = -currentHeight;
                coordinates = {0, currentHeight, currentDepth};
                jsonGeometricalDescription["coordinates"] = coordinates;
                geometricalDescription.push_back(CoreGeometricalDescriptionElement(jsonGeometricalDescription));

                currentDepth = roundFloat<6>(currentDepth + corePieceDepth);
            }

            if (spacerThickness > 0) {
                for (auto& column : corePiece->get_columns()) {
                    auto shape_data = flatten_dimensions(
                        std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()));
                    if (column.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                        jsonSpacer["type"] = OpenMagnetics::CoreGeometricalDescriptionElementType::SPACER;
                        jsonSpacer["material"] = "plastic";
                        // We cannot use directly column.get_width()
                        auto dimensions = shape_data.get_dimensions().value();
                        double windingWindowWidth;
                        if (dimensions.find("E") == dimensions.end() ||
                            (roundFloat<6>(std::get<double>(dimensions["E"])) == 0)) {
                            if (dimensions.find("F") == dimensions.end() ||
                                (roundFloat<6>(std::get<double>(dimensions["F"])) == 0)) {
                                windingWindowWidth = std::get<double>(dimensions["A"]) -
                                                     std::get<double>(dimensions["C"]) -
                                                     std::get<double>(dimensions["H"]);
                            }
                            else {
                                windingWindowWidth = std::get<double>(dimensions["A"]) -
                                                     std::get<double>(dimensions["F"]) -
                                                     std::get<double>(dimensions["H"]);
                            }
                        }
                        else {
                            windingWindowWidth = std::get<double>(dimensions["E"]);
                        }
                        double minimum_column_width;
                        double minimum_column_depth;
                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_width = std::get<double>(dimensions["A"]);
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::U ||
                                 shape_data.get_family() == CoreShapeFamily::UR) {
                            if (dimensions.find("H") == dimensions.end() ||
                                (roundFloat<6>(std::get<double>(dimensions["H"])) == 0)) {
                                minimum_column_width = (std::get<double>(dimensions["A"]) - windingWindowWidth) / 2;
                            }
                            else {
                                minimum_column_width = std::get<double>(dimensions["H"]);
                            }
                        }
                        else {
                            minimum_column_width = (std::get<double>(dimensions["A"]) - windingWindowWidth) / 2;
                        }

                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_depth = column.get_depth();
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::P ||
                                 shape_data.get_family() == CoreShapeFamily::PM) {
                            minimum_column_depth = std::get<double>(dimensions["F"]);
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::RM) {
                            if (dimensions.find("J") != dimensions.end() &&
                                (roundFloat<6>(std::get<double>(dimensions["J"])) != 0)) {
                                minimum_column_depth =
                                    sqrt(2) * std::get<double>(dimensions["J"]) - std::get<double>(dimensions["A"]);
                            }
                            else if (dimensions.find("H") != dimensions.end() &&
                                     (roundFloat<6>(std::get<double>(dimensions["H"])) != 0)) {
                                minimum_column_depth = std::get<double>(dimensions["H"]);
                            }
                            else {
                                minimum_column_depth = std::get<double>(dimensions["F"]);
                            }
                        }
                        else {
                            minimum_column_depth =
                                std::min(std::get<double>(dimensions["C"]), column.get_depth()) * numberStacks;
                        }
                        minimum_column_width *= (1 + constants.spacer_protuding_percentage);
                        minimum_column_depth *= (1 + constants.spacer_protuding_percentage);
                        double protuding_width = minimum_column_width * constants.spacer_protuding_percentage;
                        double protuding_depth = minimum_column_depth * constants.spacer_protuding_percentage;
                        jsonSpacer["dimensions"] = {minimum_column_width, spacerThickness, minimum_column_depth};
                        jsonSpacer["rotation"] = {0, 0, 0};
                        if (column.get_coordinates()[0] == 0) {
                            jsonSpacer["coordinates"] = {0, column.get_coordinates()[1],
                                                         -std::get<double>(dimensions["C"]) / 2 +
                                                             minimum_column_depth / 2 - protuding_depth};
                        }
                        else if (column.get_coordinates()[0] < 0) {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR) {
                                jsonSpacer["coordinates"] = {column.get_coordinates()[0] - column.get_width() / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]};
                            }
                            else {
                                jsonSpacer["coordinates"] = {-std::get<double>(dimensions["A"]) / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]};
                            }
                        }
                        else {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR) {
                                jsonSpacer["coordinates"] = {column.get_coordinates()[0] + column.get_width() / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]};
                            }
                            else {
                                jsonSpacer["coordinates"] = {std::get<double>(dimensions["A"]) / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]};
                            }
                        }
                        geometricalDescription.push_back(CoreGeometricalDescriptionElement(jsonSpacer));
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

void CoreWrapper::scale_to_stacks(int64_t numberStacks) {
    auto processedDescription = get_processed_description().value();
    processedDescription.get_mutable_effective_parameters().set_effective_area(
        processedDescription.get_effective_parameters().get_effective_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_minimum_area(
        processedDescription.get_effective_parameters().get_minimum_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_effective_volume(
        processedDescription.get_effective_parameters().get_effective_volume() * numberStacks);
    processedDescription.set_depth(processedDescription.get_depth() * numberStacks);
    for (auto& windingWindow : processedDescription.get_mutable_columns()) {
        windingWindow.set_area(windingWindow.get_area() * numberStacks);
        windingWindow.set_depth(windingWindow.get_depth() * numberStacks);
    }
    set_processed_description(processedDescription);
}

void CoreWrapper::distribute_and_process_gap() {
    auto constants = Constants();
    json jsonGapping = json::array();
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
            json jsonGap;
            jsonGap["type"] = GapType::RESIDUAL;
            jsonGap["length"] = constants.residualGap;
            jsonGap["coordinates"] = columns[i].get_coordinates();
            jsonGap["shape"] = columns[i].get_shape();
            jsonGap["distanceClosestNormalSurface"] = columns[i].get_height() / 2 - constants.residualGap / 2;
            jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
            jsonGap["area"] = columns[i].get_area();
            jsonGap["sectionDimensions"] = {columns[i].get_width(), columns[i].get_depth()};
            jsonGapping.push_back(jsonGap);
        }
    }
    else if (numberNonResidualGaps + numberResidualGaps < numberColumns) {
        for (size_t i = 0; i < columns.size(); ++i) {
            size_t gapIndex = i;
            if (i >= gapping.size()) {
                gapIndex = gapping.size() - 1;
            }
            json jsonGap;
            jsonGap["type"] = gapping[gapIndex].get_type();
            jsonGap["length"] = gapping[gapIndex].get_length();
            jsonGap["coordinates"] = columns[i].get_coordinates();
            jsonGap["shape"] = columns[i].get_shape();
            jsonGap["distanceClosestNormalSurface"] = columns[i].get_height() / 2 - gapping[gapIndex].get_length() / 2;
            jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
            jsonGap["area"] = columns[i].get_area();
            jsonGap["sectionDimensions"] = {columns[i].get_width(), columns[i].get_depth()};
            jsonGapping.push_back(jsonGap);
        }
    }
    else if ((numberResidualGaps == numberColumns || numberNonResidualGaps == numberColumns) &&
             (numberGaps == numberColumns)) {
        for (size_t i = 0; i < columns.size(); ++i) {
            json jsonGap;
            jsonGap["type"] = gapping[i].get_type();
            jsonGap["length"] = gapping[i].get_length();
            jsonGap["coordinates"] = columns[i].get_coordinates();
            jsonGap["shape"] = columns[i].get_shape();
            jsonGap["distanceClosestNormalSurface"] = columns[i].get_height() / 2 - gapping[i].get_length() / 2;
            jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
            jsonGap["area"] = columns[i].get_area();
            jsonGap["sectionDimensions"] = {columns[i].get_width(), columns[i].get_depth()};
            jsonGapping.push_back(jsonGap);
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
            centralColumnGapsHeightOffset = roundFloat<6>(nonResidualGaps[0].get_length() / 2);
            distanceClosestNormalSurface =
                roundFloat<6>(windingColumn.get_height() / 2 - nonResidualGaps[0].get_length() / 2);
        }
        else {
            coreChunkSizePlusGap = roundFloat<6>(windingColumn.get_height() / (nonResidualGaps.size() + 1));
            centralColumnGapsHeightOffset = roundFloat<6>(-coreChunkSizePlusGap * (nonResidualGaps.size() - 1) / 2);
            distanceClosestNormalSurface = roundFloat<6>(coreChunkSizePlusGap - nonResidualGaps[0].get_length() / 2);
        }

        for (size_t i = 0; i < nonResidualGaps.size(); ++i) {
            json jsonGap;
            jsonGap["type"] = nonResidualGaps[i].get_type();
            jsonGap["length"] = nonResidualGaps[i].get_length();
            jsonGap["coordinates"] = {windingColumn.get_coordinates()[0],
                                      windingColumn.get_coordinates()[0] + centralColumnGapsHeightOffset,
                                      windingColumn.get_coordinates()[2]};
            jsonGap["shape"] = windingColumn.get_shape();
            jsonGap["distanceClosestNormalSurface"] = distanceClosestNormalSurface;
            jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
            jsonGap["area"] = windingColumn.get_area();
            jsonGap["sectionDimensions"] = {windingColumn.get_width(), windingColumn.get_depth()};
            jsonGapping.push_back(jsonGap);

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
                json jsonGap;
                jsonGap["type"] = GapType::RESIDUAL;
                jsonGap["length"] = constants.residualGap;
                jsonGap["coordinates"] = returnColumns[i].get_coordinates();
                jsonGap["shape"] = returnColumns[i].get_shape();
                jsonGap["distanceClosestNormalSurface"] = returnColumns[i].get_height() / 2 - constants.residualGap / 2;
                jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
                jsonGap["area"] = returnColumns[i].get_area();
                jsonGap["sectionDimensions"] = {returnColumns[i].get_width(), returnColumns[i].get_depth()};
                jsonGapping.push_back(jsonGap);
            }
        }
        else {
            for (size_t i = 0; i < returnColumns.size(); ++i) {
                json jsonGap;
                jsonGap["type"] = residualGaps[i].get_type();
                jsonGap["length"] = residualGaps[i].get_length();
                jsonGap["coordinates"] = returnColumns[i].get_coordinates();
                jsonGap["shape"] = returnColumns[i].get_shape();
                jsonGap["distanceClosestNormalSurface"] = returnColumns[i].get_height() / 2;
                jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
                jsonGap["area"] = returnColumns[i].get_area();
                jsonGap["sectionDimensions"] = {returnColumns[i].get_width(), returnColumns[i].get_depth()};
                jsonGapping.push_back(jsonGap);
            }
        }
    }

    get_mutable_functional_description().set_gapping(jsonGapping);
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

void CoreWrapper::process_gap() {
    json jsonGap;
    json jsonGapping = json::array();
    auto gapping = get_functional_description().get_gapping();
    auto family = std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()).get_family();
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();

    if (family == CoreShapeFamily::T && gapping.size() > 0 ) {
        throw std::runtime_error("Toroids cannot be gapped");
    }

    if (family != CoreShapeFamily::T) {
        if (gapping.size() == 0 || !gapping[0].get_coordinates() || is_gapping_missaligned()) {
            return distribute_and_process_gap();
        }

        std::vector<int> numberGapsPerColumn;
        for (size_t i = 0; i < columns.size(); ++i) {
            numberGapsPerColumn.push_back(0);
        }

        for (size_t i = 0; i < gapping.size(); ++i) {
            auto columnIndex = find_closest_column_index_by_coordinates(*gapping[i].get_coordinates());
            numberGapsPerColumn[columnIndex]++;
        }

        for (size_t i = 0; i < gapping.size(); ++i) {
            auto columnIndex = find_closest_column_index_by_coordinates(*gapping[i].get_coordinates());

            jsonGap["type"] = gapping[i].get_type();
            jsonGap["length"] = gapping[i].get_length();
            jsonGap["coordinates"] = gapping[i].get_coordinates();
            jsonGap["shape"] = columns[columnIndex].get_shape();
            jsonGap["distanceClosestNormalSurface"] =
                roundFloat<6>(columns[columnIndex].get_height() / 2 - fabs((*gapping[i].get_coordinates())[1]) -
                              gapping[i].get_length() / 2);
            jsonGap["distanceClosestParallelSurface"] = processedDescription.get_winding_windows()[0].get_width();
            jsonGap["area"] = columns[columnIndex].get_area();
            jsonGap["sectionDimensions"] = {columns[columnIndex].get_width(), columns[columnIndex].get_depth()};
            jsonGapping.push_back(jsonGap);
        }
    }

    get_mutable_functional_description().set_gapping(jsonGapping);
}

CoreMaterial CoreWrapper::get_material() {
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(get_functional_description().get_material())) {
        auto material_data = OpenMagnetics::find_core_material_by_name(
            std::get<std::string>(get_functional_description().get_material()));
        return material_data;
    }
    else {
        return std::get<CoreMaterial>(get_functional_description().get_material());
    }

}

void CoreWrapper::process_data() {
    // If the shape is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(get_functional_description().get_shape())) {
        auto shape_data = OpenMagnetics::find_core_shape_by_name(
            std::get<std::string>(get_functional_description().get_shape()));
        shape_data.set_name(std::get<std::string>(get_functional_description().get_shape()));
        get_mutable_functional_description().set_shape(shape_data);
    }

    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid
    // long loading operations)
    if (_includeMaterialData && std::holds_alternative<std::string>(get_functional_description().get_material()) &&
        std::get<std::string>(get_functional_description().get_material()) != "dummy") {
        auto material_data = OpenMagnetics::find_core_material_by_name(
            std::get<std::string>(get_functional_description().get_material()));
        get_mutable_functional_description().set_material(material_data);
    }

    auto corePiece =
        OpenMagnetics::CorePiece::factory(std::get<OpenMagnetics::CoreShape>(get_functional_description().get_shape()));
    CoreProcessedDescription processedDescription;
    json coreEffectiveParameters;
    json coreWindingWindow;
    auto coreColumns = corePiece->get_columns();
    switch (get_functional_description().get_type()) {
        case OpenMagnetics::CoreType::TOROIDAL:
        case OpenMagnetics::CoreType::CLOSED_SHAPE:
            processedDescription.set_columns(coreColumns);
            to_json(coreEffectiveParameters, corePiece->get_partial_effective_parameters());
            processedDescription.set_effective_parameters(coreEffectiveParameters);
            to_json(coreWindingWindow, corePiece->get_winding_window());
            processedDescription.get_mutable_winding_windows().push_back(coreWindingWindow);
            processedDescription.set_depth(corePiece->get_depth());
            processedDescription.set_height(corePiece->get_height());
            processedDescription.set_width(corePiece->get_width());
            break;

        case OpenMagnetics::CoreType::TWO_PIECE_SET:
            for (auto& column : coreColumns) {
                column.set_height(2 * column.get_height());
            }
            processedDescription.set_columns(coreColumns);

            to_json(coreEffectiveParameters, corePiece->get_partial_effective_parameters());
            coreEffectiveParameters["effectiveLength"] =
                2 * coreEffectiveParameters.at("effectiveLength").get<double>();
            coreEffectiveParameters["effectiveArea"] = coreEffectiveParameters.at("effectiveArea").get<double>();
            coreEffectiveParameters["effectiveVolume"] =
                2 * coreEffectiveParameters.at("effectiveVolume").get<double>();
            coreEffectiveParameters["minimumArea"] = coreEffectiveParameters.at("minimumArea").get<double>();
            processedDescription.set_effective_parameters(coreEffectiveParameters);

            to_json(coreWindingWindow, corePiece->get_winding_window());
            coreWindingWindow.at("area") = 2 * coreWindingWindow.at("area").get<double>();
            coreWindingWindow.at("height") = 2 * coreWindingWindow.at("height").get<double>();
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

double CoreWrapper::get_magnetic_flux_density_saturation(double temperature, bool proportion) {
    auto defaults = Defaults();
    auto coreMaterial =  get_material();
    auto saturationData = coreMaterial.get_saturation();
    double saturationMagneticFluxDensity;

    if (saturationData.size() == 0) {
        throw std::runtime_error("Missing saturation data in core material");
    }
    else if (saturationData.size() == 1) {
        saturationMagneticFluxDensity = saturationData[0].get_magnetic_flux_density();
    }
    else if (saturationData.size() == 2) {
        saturationMagneticFluxDensity = saturationData[0].get_magnetic_flux_density() -
        (saturationData[0].get_temperature() - temperature) *
            (saturationData[0].get_magnetic_flux_density() - saturationData[1].get_magnetic_flux_density()) /
            (saturationData[0].get_temperature() - saturationData[1].get_temperature());
    }
    else {
        size_t n = saturationData.size();
        _1D::CubicSplineInterpolator<double> interp;
        _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

        for (size_t i = 0; i < n; i++) {
            xx(i) = saturationData[i].get_temperature();
            yy(i) = saturationData[i].get_magnetic_flux_density();
        }
        interp.setData(xx, yy);
        saturationMagneticFluxDensity = interp(temperature);
    }

    if (proportion)
        return defaults.maximumProportionMagneticFluxDensitySaturation * saturationMagneticFluxDensity;
    else
        return saturationMagneticFluxDensity;
}

double CoreWrapper::get_magnetic_fielda_strength_saturation(double temperature) {
    auto coreMaterial =  get_material();
    auto saturationData = coreMaterial.get_saturation();
    double saturationMagneticFieldStrength;

    if (saturationData.size() == 0) {
        throw std::runtime_error("Missing saturation data in core material");
    }
    else if (saturationData.size() == 1) {
        saturationMagneticFieldStrength = saturationData[0].get_magnetic_field();
    }
    else if (saturationData.size() == 2) {
        saturationMagneticFieldStrength = saturationData[0].get_magnetic_field() -
        (saturationData[0].get_temperature() - temperature) *
            (saturationData[0].get_magnetic_field() - saturationData[1].get_magnetic_field()) /
            (saturationData[0].get_temperature() - saturationData[1].get_temperature());
    }
    else {
        size_t n = saturationData.size();
        _1D::CubicSplineInterpolator<double> interp;
        _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

        for (size_t i = 0; i < n; i++) {
            xx(i) = saturationData[i].get_temperature();
            yy(i) = saturationData[i].get_magnetic_field();
        }
        interp.setData(xx, yy);
        saturationMagneticFieldStrength = interp(temperature);
    }

    return saturationMagneticFieldStrength;
}

double CoreWrapper::get_remanence(double temperature) {
    auto coreMaterial =  get_material();
    auto remanenceData = coreMaterial.get_remanence().value();
    double remanence;

    if (remanenceData.size() == 0) {
        throw std::runtime_error("Missing remanence data in core material");
    }
    else if (remanenceData.size() == 1) {
        remanence = remanenceData[0].get_magnetic_flux_density();
    }
    else if (remanenceData.size() == 2) {
        remanence = remanenceData[0].get_magnetic_flux_density() -
        (remanenceData[0].get_temperature() - temperature) *
            (remanenceData[0].get_magnetic_flux_density() - remanenceData[1].get_magnetic_flux_density()) /
            (remanenceData[0].get_temperature() - remanenceData[1].get_temperature());
    }
    else {
        size_t n = remanenceData.size();
        _1D::CubicSplineInterpolator<double> interp;
        _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

        for (size_t i = 0; i < n; i++) {
            xx(i) = remanenceData[i].get_temperature();
            yy(i) = remanenceData[i].get_magnetic_flux_density();
        }
        interp.setData(xx, yy);
        remanence = interp(temperature);
    }

    return remanence;
}

double CoreWrapper::get_coercive_force(double temperature) {
    auto coreMaterial =  get_material();
    auto coerciveForceData = coreMaterial.get_coercive_force().value();
    double coerciveForce;

    if (coerciveForceData.size() == 0) {
        throw std::runtime_error("Missing coerciveForce data in core material");
    }
    else if (coerciveForceData.size() == 1) {
        coerciveForce = coerciveForceData[0].get_magnetic_field();
    }
    else if (coerciveForceData.size() == 2) {
        coerciveForce = coerciveForceData[0].get_magnetic_field() -
        (coerciveForceData[0].get_temperature() - temperature) *
            (coerciveForceData[0].get_magnetic_field() - coerciveForceData[1].get_magnetic_field()) /
            (coerciveForceData[0].get_temperature() - coerciveForceData[1].get_temperature());
    }
    else {
        size_t n = coerciveForceData.size();
        _1D::CubicSplineInterpolator<double> interp;
        _1D::CubicSplineInterpolator<double>::VectorType xx(n), yy(n);

        for (size_t i = 0; i < n; i++) {
            xx(i) = coerciveForceData[i].get_temperature();
            yy(i) = coerciveForceData[i].get_magnetic_field();
        }
        interp.setData(xx, yy);
        coerciveForce = interp(temperature);
    }

    return coerciveForce;
}
} // namespace OpenMagnetics
