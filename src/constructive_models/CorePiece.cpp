#include "constructive_models/CorePiece.h"

#include "json.hpp"
#include "support/Utils.h"

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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 -
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_minimum_width(roundFloat(dimensions["A"] / 2 - dimensions["E"] / 2));
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_width(roundFloat(lateralColumn.get_area() / lateralColumn.get_depth()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::OBLONG);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F2"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2) )+
                          (dimensions["F2"] - dimensions["F"]) *
                              dimensions["F"]);
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 -
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::IRREGULAR);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F2"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 +
                          (dimensions["A"] - dimensions["E"]) / 4),
            0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 -
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            lateralColumn.set_depth(roundFloat(dimensions["C"] - dimensions["E"] / 2 ) -
                              dimensions["K"]);
            lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
            lateralColumn.set_minimum_width(roundFloat(dimensions["A"] / 2 - dimensions["E"] / 2));
            lateralColumn.set_width(roundFloat(lateralColumn.get_area() / lateralColumn.get_depth()));
            lateralColumn.set_height(roundFloat(dimensions["D"]));
            lateralColumn.set_coordinates({0, 0, roundFloat(-dimensions["E"] / 2 - lateralColumn.get_depth() / 2)});
            windingWindows.push_back(lateralColumn);
        }
        else {
            lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
            lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
            lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
            lateralColumn.set_height(roundFloat(dimensions["D"]));
            lateralColumn.set_coordinates({
                roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
            lateralColumn.set_coordinates({
                roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0,
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::OBLONG);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]) / 2 + roundFloat(dimensions["K"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2) )+
                          (dimensions["K"] - dimensions["F"] / 2) *
                              dimensions["F"]);
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        if ((dimensions.find("G") == dimensions.end()) || (dimensions["G"] == 0)) {
            lateralColumn.set_depth(roundFloat(dimensions["C"] - dimensions["E"] / 2 )-
                              dimensions["K"]);
            lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
            lateralColumn.set_minimum_width(roundFloat(dimensions["A"] / 2 - dimensions["E"] / 2));
            lateralColumn.set_width(roundFloat(lateralColumn.get_area() / lateralColumn.get_depth()));
            lateralColumn.set_height(roundFloat(dimensions["D"]));
            lateralColumn.set_coordinates({
                0, 0,
                roundFloat(-dimensions["E"] / 2 - lateralColumn.get_depth() / 2 -
                              (dimensions["K"] - dimensions["F"] / 2) / 2)});
            windingWindows.push_back(lateralColumn);
        }
        else {
            lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
            lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
            lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
            lateralColumn.set_height(roundFloat(dimensions["D"]));
            lateralColumn.set_coordinates({
                roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0,
                0});
            windingWindows.push_back(lateralColumn);
            lateralColumn.set_coordinates({
                roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0,
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_depth(dimensions["C"]);
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_minimum_width(roundFloat(dimensions["A"] / 2 - dimensions["E"] / 2));
        lateralColumn.set_width(roundFloat(lateralColumn.get_area() / lateralColumn.get_depth()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_area(roundFloat(get_lateral_leg_area()));
        lateralColumn.set_depth(roundFloat(lateralColumn.get_area() / lateralColumn.get_width()));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({
            roundFloat(dimensions["E"] / 2 + lateralColumn.get_width() / 2), 0, 0});
        windingWindows.push_back(lateralColumn);
        lateralColumn.set_coordinates({
            roundFloat(-dimensions["E"] / 2 - lateralColumn.get_width() / 2), 0, 0});
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
        if (dimensions.find("E") == dimensions.end() || (roundFloat(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat(dimensions["F"]) == 0)) {
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        if (dimensions.find("H") == dimensions.end() || (roundFloat(dimensions["H"]) == 0)) {
            mainColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        }
        else {
            mainColumn.set_width(roundFloat(dimensions["H"]));
        }
        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(mainColumn.get_width());
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat((dimensions["A"] + dimensions["E"]) / 2), 0, 0});
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
        if (dimensions.find("H") == dimensions.end() || (roundFloat(dimensions["H"]) == 0)) {
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
        if (dimensions.find("E") == dimensions.end() || (roundFloat(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat(dimensions["F"]) == 0)) {
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
        if (dimensions.find("E") == dimensions.end() || (roundFloat(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat(dimensions["F"]) == 0)) {
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        if (familySubtype == "1" || familySubtype == "2" || familySubtype == "4") {
            mainColumn.set_width(roundFloat(dimensions["C"]));
            mainColumn.set_depth(roundFloat(dimensions["C"]));
        }
        else {
            mainColumn.set_width(roundFloat(dimensions["F"]));
            mainColumn.set_depth(roundFloat(dimensions["F"]));
        }
        mainColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        if (familySubtype == "1" || familySubtype == "3") {
            lateralColumn.set_shape(ColumnShape::RECTANGULAR);
            lateralColumn.set_width(roundFloat(dimensions["H"]));
            lateralColumn.set_depth(roundFloat(dimensions["C"]));
            lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        }
        else {
            lateralColumn.set_shape(ColumnShape::ROUND);
            lateralColumn.set_width(roundFloat(dimensions["H"]));
            lateralColumn.set_depth(roundFloat(dimensions["H"]));
            lateralColumn.set_area(roundFloat(std::numbers::pi * pow(mainColumn.get_width() / 2, 2)));
        }
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_coordinates({roundFloat((dimensions["A"] + windingWindowWidth) / 2), 0, 0});
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

        if (dimensions.find("E") == dimensions.end() || (roundFloat(dimensions["E"]) == 0)) {
            if (dimensions.find("F") == dimensions.end() || (roundFloat(dimensions["F"]) == 0)) {
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
        mainColumn.set_type(ColumnType::LATERAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        if (dimensions.find("H") == dimensions.end() || (roundFloat(dimensions["H"]) == 0)) {
            mainColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        }
        else {
            mainColumn.set_width(roundFloat(dimensions["H"]));
        }
        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(mainColumn.get_width());
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({
            roundFloat((dimensions["A"] + dimensions["E"]) / 2), 0, 0});
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        mainColumn.set_width(columnWidth);
        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(2 * std::numbers::pi * (dimensions["B"] / 2 + columnWidth / 2));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
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

class CorePieceC : public CorePiece {
  public:
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;

        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width(dimensions["E"]);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        windingWindow.set_coordinates(std::vector<double>({(dimensions["A"] - dimensions["E"]) / 2 + dimensions["E"] / 2, 0}));
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
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        mainColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));

        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(mainColumn.get_width() * mainColumn.get_depth()));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(mainColumn.get_width());
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({roundFloat((dimensions["A"] + dimensions["E"]) / 2), 0, 0});
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

        s = (dimensions["A"] - dimensions["E"]) / 2;
        p = (dimensions["A"] - dimensions["E"]) / 2;

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
    else if (family == CoreShapeFamily::C) {
        auto piece = std::make_shared<CorePieceC>();
        piece->set_shape(shape);
        if (process)
            piece->process();
        return piece;
    }
    else
        throw std::runtime_error("Unknown shape family: " + std::string{magic_enum::enum_name(family)} + ", available options are: {E, EC, EFD, EL, EP, EPX, LP, EQ, ER, "
                                 "ETD, P, PLANAR_E, PLANAR_EL, PLANAR_ER, PM, PQ, RM, U, UR, UT, T, C}");
}

inline void from_json(const json& j, CorePiece& x) {
    x.set_columns(j.at("columns").get<std::vector<ColumnElement>>());
    x.set_depth(j.at("depth").get<double>());
    x.set_height(j.at("height").get<double>());
    x.set_width(j.at("width").get<double>());
    x.set_shape(j.at("shape").get<CoreShape>());
    x.set_winding_window(j.at("winding_window").get<WindingWindowElement>());
    x.set_partial_effective_parameters(j.at("partial_effective_parameters").get<EffectiveParameters>());
}

inline void to_json(json& j, const CorePiece& x) {
    j = json::object();
    j["columns"] = x.get_columns();
    j["depth"] = x.get_depth();
    j["height"] = x.get_height();
    j["width"] = x.get_width();
    j["shape"] = x.get_winding_window();
    j["winding_window"] = x.get_shape();
    j["partial_effective_parameters"] = x.get_partial_effective_parameters();
}

} // namespace OpenMagnetics
