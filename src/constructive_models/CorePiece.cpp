#include "constructive_models/CorePiece.h"
#include "support/Settings.h"

#include <magic_enum.hpp>

#include "json.hpp"
#include "support/Utils.h"

#include <algorithm>
#include <map>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

using json = nlohmann::json;

namespace OpenMagnetics {

void CorePiece::process() {

    process_winding_window();
    process_columns();
    process_extra_data();

    auto& settings = Settings::GetInstance();
    auto standard = settings.get_effective_parameter_standard();

    if (standard == EffectiveParameterStandard::IEC_63182) {
        // IEC 63182: get_shape_constants_iec63182() returns (le, Ae, Amin) directly.
        auto [le, Ae, minimumArea] = get_shape_constants_iec63182();
        if (le <= 0 || Ae <= 0 || minimumArea <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "IEC 63182 effective parameters cannot be negative or 0");
        }
        EffectiveParameters ep;
        ep.set_effective_length(le);
        ep.set_effective_area(Ae);
        ep.set_effective_volume(le * Ae);
        ep.set_minimum_area(minimumArea);
        set_partial_effective_parameters(ep);
    }
    else {
        // IEC 60205 (DEFAULT): shape-constant integral method.
        auto [c1, c2, minimumArea] = get_shape_constants();
        if (c1 <= 0 || c2 <= 0 || minimumArea <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "Shape constants cannot be negative or 0");
        }
        EffectiveParameters ep;
        ep.set_effective_length(pow(c1, 2) / c2);
        ep.set_effective_area(c1 / c2);
        ep.set_effective_volume(pow(c1, 3) / pow(c2, 2));
        ep.set_minimum_area(minimumArea);
        set_partial_effective_parameters(ep);
    }
}


// Default IEC 63182 fallback: derives le,Ae from IEC 60205 shape constants.
std::tuple<double, double, double> CorePiece::get_shape_constants_iec63182() {
    auto [c1, c2, minimumArea] = get_shape_constants();
    double le = pow(c1, 2) / c2;
    double Ae = c1 / c2;
    return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double h = dimensions["B"] - dimensions["D"];
        double q = dimensions["C"];
        double s = dimensions["F"] / 2;
        double p = (dimensions["A"] - dimensions["E"]) / 2;
        std::vector<double> L = { dimensions["D"], (dimensions["E"]-dimensions["F"])/2,
            dimensions["D"], std::numbers::pi/8*(p+h), std::numbers::pi/8*(s+h) };
        std::vector<double> A = { 2*q*p, 2*q*h, 2*s*q, (2*q*p+2*q*h)/2, (2*q*h+2*s*q)/2 };
        double le=0, sla=0;
        for (size_t i=0;i<L.size();++i){ le+=L[i]; sla+=L[i]/A[i]; }
        return {le, le/sla, *min_element(A.begin(),A.end())};
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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

        // FIX L-CP-1: Division by 2 in loop and *2 for minimumArea accounts for
        // half-core shape constants per IEC 60205 — areas[] describe one half-piece,
        // so effective parameters need correction for the full magnetic circuit.
        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i] / 2;
            c2 += lengths[i] / (2 * pow(areas[i], 2)) / 2;
        }
        auto minimumArea = 2 * (*min_element(areas.begin(), areas.end()));

        return {c1, c2, minimumArea};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        if (dimensions.find("q") == dimensions.end() || dimensions["q"] == 0) { // FIX M-CP-1: Validate "q" dimension exists
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Missing 'q' dimension for EFD shape");
        }
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
                          1. / 2 * std::numbers::pi * pow(dimensions["E"] / 2, 2);
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
                          1. / 4 * std::numbers::pi * pow(dimensions["E"] / 2, 2);
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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
            // DOCUMENTED APPROXIMATION: 24 catalog PQ shapes lack J/L. These
            // proportions were derived from PQ drawings and are load-bearing for
            // their published Ae/le; replace only with datasheet values.
            J = dimensions["F"] / 2;
            L = F + (C - F) / 3;
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE (innerEdge F/2 + width/2,
        // width = (E-F)/2). Matches the U/Ur/C pieces, the bobbin processors and consumers.
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
    }
};

// Slab cores: DS, HS and RS (ABT #263). Magnetics' "modified pot cores with the sides removed"
// (Magnetics 2022 Ferrite Catalog, pp. 56-57; the shape-code legend is in the 2013 edition):
//   DS = two slab halves with a SOLID centre post        (symmetric pair, no H)
//   HS = two slab halves with a DRILLED centre post      (symmetric pair, H = bore)
//   RS = one slab half + one plain pot-core round        (ASYMMETRIC pair)
//
// The letters are the pot core's (the catalogue uses the same A B 2B C D 2D E F G H table), with
// C meaning what the slab drawing dimensions it as: the distance ACROSS THE TWO FLATS that cut
// the whole half, flange and wall alike. G, the wire opening, is the chord where a flat breaks
// the wall and is DERIVED from C and E, so the geometry does not read it.
//
// Why not CorePieceP: its wire-opening term k1 = n*b*(r4 - r3) is linear in the slot width,
// right for a pot core's narrow notch and wrong when the "opening" is a removed side: it
// over-predicts Ae by 12-34% on these rows. Here the same IEC 60205 pot clause is kept but every
// annular section is weighted by the angle the flats leave, theta(r) = 2 pi - 4 acos(C / 2r)
// for r > C/2, which is exact for this outline: wall area, flange spreading, corner areas.
//
// ACCURACY -- READ BEFORE TRUSTING A NUMBER FROM THIS CLASS. Against the catalogue's own
// le/Ae/Ve for the 19 rows shipped in MAS (Test_Slab_Cores_Match_Catalogue_Within_Accepted_Band):
//   le:  mean |error| 4.1%, worst 7.7%   (the model runs short)
//   Ae:  mean |error| 10.5%, worst +23.4% (DS 26/16; the model runs high on every row)
// The owner accepted ~15-20% for this family rather than wait for better data (2026-09-23). The
// residual is not the flats: the catalogue gives most slab D/E/F only as a min OR a max, never a
// nominal, and resolve_dimensional_values then has to take that bound as the value, so the model
// is fed a post at its largest and a window at its smallest. Nominal vendor drawings would fix it.
//
// NOT SHIPPED: DS/HS 23/11 and DS/HS 23/18. Their catalogue Amin (37.8 and 40.7 mm^2) is far
// below even the solid post the drawing gives (pi/4 * 9.9^2 = 77 mm^2), so the drawing and the
// magnetic data describe different parts; this model lands +27-49% on Ae there. The RS rows of
// the same sizes carry their own dimensions and fit (+3-6%), so they are shipped.
class CorePieceSlab : public CorePiece {
  protected:
    // RS pairs one slab with a plain pot-core round; DS/HS pair two slabs.
    virtual bool paired_with_round() const { return false; }

    struct SlabDimensions {
        double outerRadius;       // A/2, r4
        double windowOuterRadius; // E/2, r3
        double postRadius;        // F/2, r2
        double boreRadius;        // H/2, r1 (0 for a solid post)
        double halfAcrossFlats;   // C/2
        double halfHeight;        // B
        double windowHalfHeight;  // D
    };

    SlabDimensions read_dimensions() const {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        for (auto letter : {"A", "B", "C", "D", "E", "F"}) {
            if (!dimensions.count(letter) || !(dimensions[letter] > 0)) {
                throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                    std::string("slab core: dimension ") + letter + " is missing or non-positive");
            }
        }
        SlabDimensions slab;
        slab.outerRadius = dimensions["A"] / 2;
        slab.windowOuterRadius = dimensions["E"] / 2;
        slab.postRadius = dimensions["F"] / 2;
        slab.boreRadius = dimensions.count("H") ? dimensions["H"] / 2 : 0.0;
        slab.halfAcrossFlats = dimensions["C"] / 2;
        slab.halfHeight = dimensions["B"];
        slab.windowHalfHeight = dimensions["D"];
        if (slab.boreRadius < 0 || slab.boreRadius >= slab.postRadius) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "slab core: bore H must be non-negative and smaller than the post F");
        }
        if (slab.postRadius >= slab.windowOuterRadius || slab.windowOuterRadius >= slab.outerRadius) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "slab core: needs F < E < A (post inside the window inside the outline)");
        }
        if (slab.halfAcrossFlats <= slab.postRadius || slab.halfAcrossFlats > slab.outerRadius) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "slab core: the flats C must clear the post F and lie within the outline A");
        }
        if (slab.windowHalfHeight >= slab.halfHeight) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "slab core: window height D must be smaller than the half height B");
        }
        return slab;
    }

    // Angle of the annulus at radius r that survives the two flats (2 pi when r is inside them).
    static double surviving_angle(double radius, double halfAcrossFlats) {
        if (radius <= halfAcrossFlats) {
            return 2 * std::numbers::pi;
        }
        return 2 * std::numbers::pi - 4 * acos(halfAcrossFlats / radius);
    }

    // Integral of theta(r) r dr from innerRadius to outerRadius, closed form:
    // for r > c, d/dr[pi r^2 - 4 ((r^2/2) acos(c/r) - (c/2) sqrt(r^2 - c^2))] = theta(r) r.
    static double surviving_area(double innerRadius, double outerRadius, double halfAcrossFlats) {
        double pi = std::numbers::pi;
        double c = halfAcrossFlats;
        auto primitive = [c, pi](double r) {
            return pi * r * r - 4 * (r * r / 2 * acos(c / r) - c / 2 * sqrt(r * r - c * c));
        };
        double area = 0;
        if (innerRadius < c) {
            area += pi * (pow(std::min(c, outerRadius), 2) - pow(innerRadius, 2));
        }
        double lower = std::max(innerRadius, c);
        if (outerRadius > lower) {
            area += primitive(outerRadius) - primitive(lower);
        }
        return area;
    }

    // Integral of (theta(r) r)^-power dr over [innerRadius, outerRadius]: the flange's radial
    // spreading term. No closed form once theta varies, so composite Simpson, split at the flats
    // where theta has a kink. 2000 panels per side is far below any tolerance that matters here.
    static double radial_spreading_integral(double innerRadius, double outerRadius,
                                            double halfAcrossFlats, int power) {
        auto integrand = [halfAcrossFlats, power](double r) {
            return pow(surviving_angle(r, halfAcrossFlats) * r, -power);
        };
        auto simpson = [&integrand](double from, double to) {
            const int panels = 2000;
            double step = (to - from) / panels;
            double sum = integrand(from) + integrand(to);
            for (int index = 1; index < panels; ++index) {
                sum += integrand(from + index * step) * ((index % 2 == 1) ? 4 : 2);
            }
            return sum * step / 3;
        };
        if (halfAcrossFlats > innerRadius && halfAcrossFlats < outerRadius) {
            return simpson(innerRadius, halfAcrossFlats) + simpson(halfAcrossFlats, outerRadius);
        }
        return simpson(innerRadius, outerRadius);
    }

    struct SetConstants {
        double c1;
        double c2;
        double minimumArea;
        double wallArea;
    };

    // IEC 60205 pot clause for a SYMMETRIC SET of two identical halves, exactly as
    // CorePieceP::get_shape_constants lays it out (wall, flange, post, two corners), with each
    // annular area weighted by the angle the flats leave. With flats=false it is CorePieceP
    // without wire slots (n = 0), which is what the plain pot round of an RS set is.
    static SetConstants symmetric_set_constants(const SlabDimensions& slab, bool flats) {
        double pi = std::numbers::pi;
        double c = flats ? slab.halfAcrossFlats : std::numeric_limits<double>::infinity();
        double r4 = slab.outerRadius;
        double r3 = slab.windowOuterRadius;
        double r2 = slab.postRadius;
        double r1 = slab.boreRadius;
        double h = slab.halfHeight - slab.windowHalfHeight;
        double h2 = 2 * slab.windowHalfHeight;

        double wall = surviving_area(r3, r4, c);
        double a1 = wall;
        double l1 = h2;
        double a3 = pi * (r2 * r2 - r1 * r1);
        double l3 = h2;
        double s1 = r2 - sqrt((r1 * r1 + r2 * r2) / 2);
        double s2 = sqrt((r3 * r3 + r4 * r4) / 2) - r3;
        double l4 = pi / 4 * (2 * s2 + h);
        double a4 = (wall + surviving_angle(r3, c) * r3 * h) / 2;
        double l5 = pi / 4 * (2 * s1 + h);
        double a5 = pi / 2 * (r2 * r2 - r1 * r1 + 2 * r2 * h);
        // Two flanges, one per half.
        double radial1 = 2 / h * radial_spreading_integral(r2, r3, c, 1);
        double radial2 = 2 / (h * h) * radial_spreading_integral(r2, r3, c, 2);

        SetConstants set;
        set.c1 = l1 / a1 + radial1 + l3 / a3 + l4 / a4 + l5 / a5;
        set.c2 = l1 / (a1 * a1) + radial2 + l3 / (a3 * a3) + l4 / (a4 * a4) + l5 / (a5 * a5);
        set.minimumArea = std::min({a1, a3, a4, a5});
        set.wallArea = wall;
        return set;
    }

    // The whole set's constants. For RS each half contributes half of its own symmetric set.
    SetConstants set_constants() const {
        auto slab = read_dimensions();
        auto slabSet = symmetric_set_constants(slab, true);
        if (!paired_with_round()) {
            return slabSet;
        }
        auto roundSet = symmetric_set_constants(slab, false);
        SetConstants set;
        set.c1 = (slabSet.c1 + roundSet.c1) / 2;
        set.c2 = (slabSet.c2 + roundSet.c2) / 2;
        set.minimumArea = std::min(slabSet.minimumArea, roundSet.minimumArea);
        set.wallArea = (slabSet.wallArea + roundSet.wallArea) / 2;
        return set;
    }

  public:
    void process_winding_window() {
        auto slab = read_dimensions();
        WindingWindowElement windingWindow;
        double width = slab.windowOuterRadius - slab.postRadius;
        windingWindow.set_height(slab.windowHalfHeight);
        windingWindow.set_width(width);
        windingWindow.set_area(slab.windowHalfHeight * width);
        windingWindow.set_coordinates(std::vector<double>({slab.postRadius + width / 2, 0}));
        set_winding_window(windingWindow);
    }

    void process_extra_data() {
        auto slab = read_dimensions();
        // Seen through the walls (the flats face the other way): full diameter across, the
        // flats' width deep. An RS set's round half is not cut, so its envelope is the full A.
        set_width(2 * slab.outerRadius);
        set_height(slab.halfHeight);
        set_depth(paired_with_round() ? 2 * slab.outerRadius : 2 * slab.halfAcrossFlats);
    }

    void process_columns() {
        auto slab = read_dimensions();
        auto set = set_constants();
        std::vector<ColumnElement> columns;
        ColumnElement mainColumn;
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(2 * slab.postRadius));
        mainColumn.set_depth(roundFloat(2 * slab.postRadius));
        mainColumn.set_height(roundFloat(slab.windowHalfHeight));
        mainColumn.set_area(roundFloat(std::numbers::pi * (pow(slab.postRadius, 2) - pow(slab.boreRadius, 2))));
        mainColumn.set_coordinates({0, 0, 0});
        columns.push_back(mainColumn);

        // The two surviving wall arcs, each carrying half the wall's area.
        ColumnElement lateralColumn;
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::IRREGULAR);
        double lateralWidth = slab.outerRadius - slab.windowOuterRadius;
        lateralColumn.set_width(roundFloat(lateralWidth));
        lateralColumn.set_area(roundFloat(set.wallArea / 2));
        lateralColumn.set_depth(roundFloat(set.wallArea / 2 / lateralWidth));
        lateralColumn.set_height(roundFloat(slab.windowHalfHeight));
        lateralColumn.set_coordinates({roundFloat(slab.windowOuterRadius + lateralWidth / 2), 0, 0});
        columns.push_back(lateralColumn);
        lateralColumn.set_coordinates({roundFloat(-slab.windowOuterRadius - lateralWidth / 2), 0, 0});
        columns.push_back(lateralColumn);
        set_columns(columns);
    }

    // Per piece, half the set, as every two-piece family here reports it.
    std::tuple<double, double, double> get_shape_constants() {
        auto set = set_constants();
        return {set.c1 / 2, set.c2 / 2, set.minimumArea};
    }
};

class CorePieceDs : public CorePieceSlab {};
class CorePieceHs : public CorePieceSlab {};
class CorePieceRs : public CorePieceSlab {
  protected:
    bool paired_with_round() const override { return true; }
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
            lateralColumn.set_area(roundFloat(std::numbers::pi * pow(lateralColumn.get_width() / 2, 2)));
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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
        // ABT #107: coordinates[0] is the winding-window CENTRE = innerEdge (A-E)/2 + width/2
        // (width = E). Previously stored only the inner edge; the U/Ur/C pieces already centre.
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        // IEC 63182: Calculation method for effective parameters of ring-cores
        // le = π × (A - B) / ln(A/B)  (logarithmic mean magnetic path length)
        // Ag = ((A - B) / 2) × C - (4 - π) × r₀²  (geometric cross-sectional area)
        // Ve = le × Ag  (effective volume)
        // Note: Ag is taken as the effective cross-sectional area (Ae)
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double A = dimensions["A"];   // Outer diameter (OD)
        double B = dimensions["B"];   // Inner diameter (ID)
        double C = dimensions["C"];   // Height (h)
        
        // Effective magnetic path length (logarithmic mean)
        double le = std::numbers::pi * (A - B) / std::log(A / B);
        
        // Geometric cross-sectional area
        double Ag = ((A - B) / 2.0) * C;
        
        // Apply chamfer correction if rounding radius r₀ is available
        if (dimensions.find("r0") != dimensions.end() || dimensions.find("R") != dimensions.end()) {
            double r0 = 0.0;
            if (dimensions.find("r0") != dimensions.end()) {
                r0 = dimensions["r0"];
            } else if (dimensions.find("R") != dimensions.end()) {
                r0 = dimensions["R"];
            }
            if (r0 > 0) {
                Ag -= (4.0 - std::numbers::pi) * r0 * r0;
            }
        }
        
        return {le, Ag, Ag};  // For toroid, minimum area = geometric area = effective area
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

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        double le = pow(c1, 2) / c2;
        double Ae = c1 / c2;
        return {le, Ae, minimumArea};
    }
};

// ---- Additional families (2026-07, ABT #263-#272): geometry variants ----
// Each is modelled as sharing the A-F dimension convention and magnetic-path model
// of an existing piece. The mapping is validated per-family by comparing the
// MKF-computed Ae/Le against the vendor datasheet's stated Ae/Le for a
// representative size before any cores of that family are emitted (see the
// validation notes in the accompanying commit / ABT tickets).
// EER: E-type with a round centre leg -> identical model to ETD/ER.
class CorePieceEer : public CorePieceEtd {};
// EF: flat/economy E core -> E geometry.
class CorePieceEf : public CorePieceE {};
// EPC (TDK's low-profile EPC series, ABT #270). An EPC is NOT an EP. Seen on its mating face it is
// an E core lying on its side: a back plate, two outer legs, and a centre pole whose section is a
// STADIUM (a rectangle with a half-disc on each end), long along the width A and short across the
// thin dimension. It is not EP's round pole in a round bore. Forcing it through CorePieceEp gave
// -20 % le / +11 % Ae on EPC13, so it has its own piece.
//
// Letters. TDK's drawing (catalogue "EPC series", 2014-03 ferrite_mz_sw_epc_en.pdf and 2026-04
// ferrite_mz_rf-power_epc_en.pdf, Fig. 1) uses its own labels. The shape record uses the house E
// labels plus the oblong-column pair F/F2 that CorePieceEl/CorePieceMolded already use:
//     house  TDK     meaning
//     A      A       overall width
//     B      D       height of one half, mating face to back
//     C      F       depth (the low-profile thickness)
//     D      H       height of the winding window in one half
//     E      B min.  window width between the outer legs
//     F      C1      pole length along A
//     F2     C2      pole thickness along C; the stadium's end discs have diameter F2
//     G      E min.  opening between the ears at the top of the window (recorded, not modelled)
// The stadium area F*F2 - (4 - pi)/4 * F2^2 reproduces TDK's published centre-pole area Acp on
// every size: EPC13 10.58 vs 10.6 mm2, EPC17 19.88 vs 19.9, EPC25 42.57 vs 42.6, EPC27 48.57 vs 48.6.
//
// THIS MODEL IS APPROXIMATE, AND THE OWNER ACCEPTED IT AS SUCH (ABT #270, 2026-09-23). The
// magnetic path is the house E-core IEC 60205 walk (CorePieceE) with the true stadium pole area.
// Against TDK's published effective parameters, for the three records MAS ships (EPC13/EPC17 on
// the 2026 dimensions, EPC25 on the 2014 ones):
//     EPC 13   le -7.4 %   Ae +0.4 %   Ve  -7.0 %
//     EPC 17   le -5.3 %   Ae -6.7 %   Ve -11.6 %
//     EPC 25   le -1.3 %   Ae +2.8 %   Ve  +1.3 %
// For comparison, from the 2014 table on sizes not shipped: EPC27 le -0.8 % / Ae +0.3 %, EPC30
// le -0.4 % / Ae +2.3 %. The EPC25B variant (TDK Fig. 2/3) is off by +14 % in le and is NOT
// covered by this piece. The small sizes come out short on le. The drawing does not dimension the ears, the curved
// inner walls of the outer legs or the pole's offset toward the bottom face. The ferrite mass says
// there is LESS material than this walk assumes (EPC17: TDK 4.5 g/set against about 4.8 g
// modelled), so the shortfall is not simply missing ears. Other variants were tried and did no
// better: EL-style pole-perimeter spreading gave EPC13 le -10.8 %; taking the window width from
// TDK's published winding area gave -7.6 %; adding the ears to the legs moved Ae but not le.
// Nothing is fitted to close the gap. Do not read EPC13/EPC17 results as better than about 8 % in
// le, 7 % in Ae and 12 % in Ve (and so in core loss).
//
// The pole is a RECTANGULAR column F x F2 with a corner radius of F2 / 2, which is exactly the
// stadium. It is NOT ColumnShape::OBLONG, because MKF's OBLONG assumes the long axis lies along the
// depth (width F < depth F2; bobbin corner radius = depth / 2; turn length 2*pi*r + 4*(depth -
// width)). The EPC pole's long axis faces the window, and that formula would silently cut its
// turn length by about 19 %. Under real winding geometry the rectangular frame charges the corner
// radius, and the turn around a rectangle whose corner radius is half its short side IS the
// stadium turn. The classic (non-real-winding) path treats every rectangular column as sharp, so
// there an EPC turn is charged the bounding rectangle: (4 - pi) times the short side of the bobbin
// bore longer per turn than the stadium. That is the same sharp-corner approximation every other
// rectangular column gets in that mode. On EPC 13 with a quick bobbin, the first layer comes out at
// 19.63 mm under real winding (the stadium exactly) and 22.62 mm classic (+15 %).
class CorePieceEpc : public CorePiece {
  public:
    // The dimensions this piece reads, checked once, loudly: an EPC record whose pole is not a
    // stadium lying along A cannot be modelled by the walk below.
    std::map<std::string, double> checked_dimensions() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        for (const auto& key : {"A", "B", "C", "D", "E", "F", "F2"}) {
            if (dimensions.find(key) == dimensions.end() || dimensions[key] <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                    "EPC shape '" + get_shape().get_name().value_or("<unnamed>") + "' has no dimension " + key);
            }
        }
        if (dimensions["F2"] > dimensions["F"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "EPC shape '" + get_shape().get_name().value_or("<unnamed>") + "': the pole thickness F2 (TDK C2) "
                "is larger than its length F (TDK C1); an EPC pole is a stadium long along A");
        }
        if (dimensions["F2"] > dimensions["C"] || dimensions["F"] >= dimensions["E"] ||
            dimensions["E"] >= dimensions["A"] || dimensions["D"] >= dimensions["B"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "EPC shape '" + get_shape().get_name().value_or("<unnamed>") + "' is not a closed E-like "
                "section: need F2 <= C, F < E < A and D < B");
        }
        return dimensions;
    }

    static double stadium_area(double length, double thickness) {
        return length * thickness - (4 - std::numbers::pi) / 4 * pow(thickness, 2);
    }

    void process_extra_data() {
        auto dimensions = checked_dimensions();
        set_width(dimensions["A"]);
        set_height(dimensions["B"]);
        set_depth(dimensions["C"]);
    }

    void process_winding_window() {
        auto dimensions = checked_dimensions();
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["D"]);
        windingWindow.set_width((dimensions["E"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        // Centre of the window, as for CorePieceE (ABT #107).
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = checked_dimensions();
        std::vector<ColumnElement> columns;
        ColumnElement mainColumn;
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        mainColumn.set_width(roundFloat(dimensions["F"]));
        mainColumn.set_depth(roundFloat(dimensions["F2"]));
        mainColumn.set_corner_radius(roundFloat(dimensions["F2"] / 2));
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_area(roundFloat(stadium_area(dimensions["F"], dimensions["F2"])));
        mainColumn.set_coordinates({0, 0, 0});
        columns.push_back(mainColumn);

        ColumnElement lateralColumn;
        lateralColumn.set_type(ColumnType::LATERAL);
        lateralColumn.set_shape(ColumnShape::RECTANGULAR);
        lateralColumn.set_width(roundFloat((dimensions["A"] - dimensions["E"]) / 2));
        lateralColumn.set_depth(roundFloat(dimensions["C"]));
        lateralColumn.set_height(roundFloat(dimensions["D"]));
        lateralColumn.set_area(roundFloat(lateralColumn.get_width() * lateralColumn.get_depth()));
        lateralColumn.set_coordinates({roundFloat(dimensions["E"] / 2 + (dimensions["A"] - dimensions["E"]) / 4), 0, 0});
        columns.push_back(lateralColumn);
        lateralColumn.set_coordinates({roundFloat(-dimensions["E"] / 2 - (dimensions["A"] - dimensions["E"]) / 4), 0, 0});
        columns.push_back(lateralColumn);
        set_columns(columns);
    }

    // CorePieceE's IEC 60205 walk (same sections, same corner rule, same per-piece convention),
    // with the centre section's area replaced by the stadium's. s is half the pole's extent
    // along the path (F / 2), as CorePieceE's s is half its rectangular leg.
    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = checked_dimensions();
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
        areas.push_back(stadium_area(dimensions["F"], dimensions["F2"]));
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

// A PIECE-AND-PLATE core: a shaped piece (U, PQ, E...) closed by a flat I plate rather than by a
// mirrored second half. IEC 60205:2016 has no clause for this combination -- every clause in 5.x is
// "Pair of X-cores" -- so the section list below applies the standard's GENERAL method instead:
// C1 = sum(l_i/A_i), C2 = sum(l_i/A_i^2), le = C1^2/C2, Ae = C1/C2 (clause 4.6 and following), with
// corner sections taken as "the mean circular path joining the centres of area of the two adjacent
// uniform sections", their area "the average area of the two adjacent uniform sections" (clause 4.6).
//
// Two consequences of the plate, versus the pair the parent class models:
//   * the window is ONE piece tall (D), not two (2D), so each leg contributes D and not 2D;
//   * one of the two yokes is the plate, of thickness B2 instead of the piece's own B - D.
// IEC's l4 / l5 each already lump TWO physical corners ("l4 = l4' + l4'' = pi/4 (p + h)"), so a
// single corner is pi/8 of the same sum; a piece-and-plate has four corners of two different
// thicknesses, hence four pi/8 terms rather than two pi/4 ones.
//
// These constants are the WHOLE assembled core, not a half-piece. CoreType::PIECE_AND_PLATE must
// therefore NOT double them the way TWO_PIECE_SET doubles a mirrored pair.
//
// VALIDATED against the Magnetics 2022 Ferrite Catalog (printed pages 38-39). Note the catalogue
// lists a U+I combination's parameters on the I-core row -- an I bar has no closed magnetic path of
// its own, so its published le/Ae are those of the combination:
//     I 93/28/16   computed le 257.0 mm / Ae 451.5 mm2   vs published 257 / 450
//     I 25/6/6     computed le  64.1 mm / Ae  40.8 mm2   vs published 64.3 / 40.3
// The same section list with 2D legs and two piece yokes reproduces the U-pair row exactly
// (U 93/76/16: computed 353.0 / 452.3 vs published 353 / 452), which is the cross-check that the
// method and not just the fit is right.
static std::tuple<double, double, double> piece_and_plate_shape_constants(
        double legHeight, double yokeLength, double depth,
        double legWidthOne, double legWidthTwo, double pieceYokeThickness, double plateThickness) {
    double areaLegOne = depth * legWidthOne;
    double areaLegTwo = depth * legWidthTwo;
    double areaPieceYoke = depth * pieceYokeThickness;
    double areaPlate = depth * plateThickness;

    std::vector<double> lengths;
    std::vector<double> areas;
    // Straight sections: two legs, the piece's own yoke, and the plate closing the circuit.
    lengths.push_back(legHeight);      areas.push_back(areaLegOne);
    lengths.push_back(legHeight);      areas.push_back(areaLegTwo);
    lengths.push_back(yokeLength);     areas.push_back(areaPieceYoke);
    lengths.push_back(yokeLength);     areas.push_back(areaPlate);
    // Four corners, IEC 60205 clause 4.6.
    const double quarterCircle = std::numbers::pi / 8;
    lengths.push_back(quarterCircle * (legWidthOne + pieceYokeThickness));
    areas.push_back((areaLegOne + areaPieceYoke) / 2);
    lengths.push_back(quarterCircle * (legWidthTwo + pieceYokeThickness));
    areas.push_back((areaLegTwo + areaPieceYoke) / 2);
    lengths.push_back(quarterCircle * (legWidthOne + plateThickness));
    areas.push_back((areaLegOne + areaPlate) / 2);
    lengths.push_back(quarterCircle * (legWidthTwo + plateThickness));
    areas.push_back((areaLegTwo + areaPlate) / 2);

    double c1 = 0, c2 = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        c1 += lengths[i] / areas[i];
        c2 += lengths[i] / pow(areas[i], 2);
    }
    return {c1, c2, *min_element(areas.begin(), areas.end())};
}

class CorePieceUi : public CorePieceU {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        // The assembled height is the U plus the plate laid on top of it.
        set_height(dimensions["B"] + dimensions["B2"]);
        set_depth(dimensions["C"]);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double legWidthOne, legWidthTwo;
        if (dimensions.find("H") == dimensions.end() || (roundFloat(dimensions["H"]) == 0)) {
            legWidthOne = (dimensions["A"] - dimensions["E"]) / 2;
            legWidthTwo = (dimensions["A"] - dimensions["E"]) / 2;
        }
        else {
            legWidthOne = dimensions["H"];
            legWidthTwo = dimensions["A"] - dimensions["E"] - dimensions["H"];
        }
        return piece_and_plate_shape_constants(dimensions["D"], dimensions["E"], dimensions["C"],
                                               legWidthTwo, legWidthOne,
                                               dimensions["B"] - dimensions["D"], dimensions["B2"]);
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

// An E-PIECE-AND-PLATE core (ABT #625): an E-shaped piece closed by a flat I plate rather than by a
// mirrored second E-half. This is NOT the two-leg piece_and_plate_shape_constants above: E's own
// single-piece formula (CorePieceE::get_shape_constants) treats the piece as TWO parallel,
// geometrically-identical sub-loops sharing a half-width centre post (areas carry a x2 factor,
// lengths use the piece's own D once per leg -- see that method's comments), and implicitly assumes
// the loop closes at the leg tips with zero extra length -- exactly right when mirrored by an
// identical E-half (TWO_PIECE_SET doubles le/Ve externally in Core.cpp) but wrong for a flat plate:
// the plate is itself a second closing element, of its OWN thickness B2, at the leg tips -- it is
// NOT another D-tall leg. So: keep the piece's own leg terms (D height, once each, not doubled) and
// its own yoke term (thickness B-D), and ADD a second "yoke" of thickness B2 standing in for the
// plate, with its own pair of corner terms to the same legs.
//
// VALIDATED against the Magnetics 2022 Ferrite Catalog (p.32-35, planar E/I dimensions and magnetic
// data), which publishes the E-piece-ALONE (mirrored-pair) le/Ae under the E-code row and the E+I
// ASSEMBLY's le/Ae under the mated I-code row (same convention as the U/I validation above -- "an I
// bar has no closed magnetic path of its own"):
//     E 40/8/10 (0_44008EC, A=40.65 B=8.51 C=10.7 D=4.06 E=30.45 F=10.15): CorePieceE's own formula
//         doubled (TWO_PIECE_SET) gives le=51.4mm/Ae=100.9mm2 vs published 51.9/101.
//     I 40/4/10 (0_44008IC, B2=4.45, symmetric plate/yoke): this method gives le=43.3mm/Ae=99.6mm2
//         vs published 43.8/99.5.
//     I 43/4/28 (0_44308IC, B2=4.1 vs the piece's own yoke B-D=4.32 -- an ASYMMETRIC plate/yoke
//         case): this method gives le=47.9mm/Ae=235mm2 vs published 48.6/227.
static std::tuple<double, double, double> e_piece_and_plate_shape_constants(
        double legHeight, double yokeSpan, double depth,
        double outerLegWidth, double halfPostWidth,
        double pieceYokeThickness, double plateThickness) {
    double areaOuterLeg = 2 * depth * outerLegWidth;
    double areaPost = 2 * depth * halfPostWidth;
    double areaPieceYoke = 2 * depth * pieceYokeThickness;
    double areaPlateYoke = 2 * depth * plateThickness;

    std::vector<double> lengths;
    std::vector<double> areas;
    lengths.push_back(legHeight);  areas.push_back(areaOuterLeg);
    lengths.push_back(legHeight);  areas.push_back(areaPost);
    lengths.push_back(yokeSpan);   areas.push_back(areaPieceYoke);
    lengths.push_back(yokeSpan);   areas.push_back(areaPlateYoke);
    const double quarterCircle = std::numbers::pi / 8;
    lengths.push_back(quarterCircle * (outerLegWidth + pieceYokeThickness));
    areas.push_back((areaOuterLeg + areaPieceYoke) / 2);
    lengths.push_back(quarterCircle * (halfPostWidth + pieceYokeThickness));
    areas.push_back((areaPost + areaPieceYoke) / 2);
    lengths.push_back(quarterCircle * (outerLegWidth + plateThickness));
    areas.push_back((areaOuterLeg + areaPlateYoke) / 2);
    lengths.push_back(quarterCircle * (halfPostWidth + plateThickness));
    areas.push_back((areaPost + areaPlateYoke) / 2);

    double c1 = 0, c2 = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        c1 += lengths[i] / areas[i];
        c2 += lengths[i] / pow(areas[i], 2);
    }
    return {c1, c2, *min_element(areas.begin(), areas.end())};
}

class CorePieceEi : public CorePieceE {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        // The assembled height is the E piece plus the flat plate closing it.
        set_height(dimensions["B"] + dimensions["B2"]);
        set_depth(dimensions["C"]);
    }
    // process_columns() and process_winding_window() are inherited unmodified from CorePieceE: a
    // flat plate has no legs of its own, so the column height and window height stay the single
    // piece's own D, exactly like CorePieceUi inherits CorePieceU's (see CoreType::PIECE_AND_PLATE
    // in Core.cpp -- "the column height and the winding window stay those of the single piece").

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pieceYokeThickness = dimensions["B"] - dimensions["D"];
        double depth = dimensions["C"];
        double halfPostWidth = dimensions["F"] / 2;
        double outerLegWidth = (dimensions["A"] - dimensions["E"]) / 2;
        double yokeSpan = (dimensions["E"] - dimensions["F"]) / 2;
        return e_piece_and_plate_shape_constants(dimensions["D"], yokeSpan, depth, outerLegWidth,
                                                 halfPostWidth, pieceYokeThickness, dimensions["B2"]);
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};


// DRUM (bobbin / dumbbell) core: a round centre post between two flange discs, wound in the
// groove; the magnetic circuit CLOSES THROUGH THE SURROUNDING AIR (magneticCircuit = open,
// CoreType::OPEN_SHAPE). Dimension convention (defined with the MAS records, ABT #331):
//   A flange OD, B total height, C post OD, D top flange thickness, E winding groove height,
//   F bottom flange thickness, H bore (optional). Every catalogued row satisfies D+E+F == B.
// Sources: TDG DRH datasheet (letter key printed in TDG's ordering system) and Fair-Rite
// bobbins whose letter mapping was weight-verified (volume x density == listed grams).
//
// IMPORTANT: get_shape_constants below describes the FERRITE INTERNAL PATH ONLY (post +
// two spreading flanges + corners, IEC 60205 general method). It is honest as the piece's
// partial parameters — but it is NOT a closed magnetic circuit, and feeding it to the
// closed-circuit reluctance path would silently drop the dominant air-return reluctance.
// Open-core magnetizing inductance goes through the dedicated model in
// MagnetizingInductance.cpp (demagnetising-factor bracket, validated against published
// Fair-Rite AL), which routes on CoreShapeFamily::DRUM.
class CorePieceDrum : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double largerFlange = dimensions["A"];
        if (dimensions.find("A2") != dimensions.end() && dimensions["A2"] > 0) {
            largerFlange = std::max(dimensions["A"], dimensions["A2"]);
        }
        set_width(largerFlange);
        set_height(dimensions["B"]);
        set_depth(largerFlange);
    }

    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        // A2 = the SECOND flange OD for asymmetric drums (TDG DRC A1/A2 convention; e.g. the
        // WE-TI radial chokes have a larger base flange). The usable groove is bounded by the
        // SMALLER flange — wire beyond it is not retained.
        double smallerFlange = dimensions["A"];
        if (dimensions.find("A2") != dimensions.end() && dimensions["A2"] > 0) {
            smallerFlange = std::min(dimensions["A"], dimensions["A2"]);
        }
        WindingWindowElement windingWindow;
        windingWindow.set_height(dimensions["E"]);
        windingWindow.set_width((smallerFlange - dimensions["C"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        // ABT #107 convention: coordinates[0] = window CENTRE (post edge + half width).
        windingWindow.set_coordinates(std::vector<double>({dimensions["C"] / 2 + (smallerFlange - dimensions["C"]) / 4, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double bore = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] : 0.0;
        std::vector<ColumnElement> windingWindows;
        ColumnElement mainColumn;
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(dimensions["C"]));
        mainColumn.set_depth(roundFloat(dimensions["C"]));
        mainColumn.set_height(roundFloat(dimensions["E"]));
        // Bore subtracted: the mounting hole carries no flux.
        mainColumn.set_area(roundFloat(std::numbers::pi / 4 * (pow(dimensions["C"], 2) - pow(bore, 2))));
        mainColumn.set_coordinates({0, 0, 0});
        windingWindows.push_back(mainColumn);
        set_columns(windingWindows);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pi = std::numbers::pi;
        double flangeRadius = dimensions["A"] / 2;
        double postRadius = dimensions["C"] / 2;
        double boreRadius = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] / 2 : 0.0;
        double grooveHeight = dimensions["E"];

        double postArea = pi * (pow(postRadius, 2) - pow(boreRadius, 2));
        std::vector<double> areas;
        areas.push_back(postArea);
        double c1 = grooveHeight / postArea;
        double c2 = grooveHeight / pow(postArea, 2);
        // Two flanges: radial spreading discs from post edge to flange edge, plus the
        // post-to-flange corners per IEC 60205 clause 4.6.
        for (auto flangeThickness : {dimensions["D"], dimensions["F"]}) {
            c1 += 1.0 / (2 * pi * flangeThickness) * log(flangeRadius / postRadius);
            c2 += 1.0 / (2 * pow(pi * flangeThickness, 2)) * (flangeRadius - postRadius) / (flangeRadius * postRadius);
            areas.push_back(2 * pi * postRadius * flangeThickness);
            double s1 = postRadius - sqrt((pow(boreRadius, 2) + pow(postRadius, 2)) / 2);
            double cornerLength = pi / 4 * (2 * s1 + flangeThickness);
            double cornerArea = 0.5 * (postArea + 2 * pi * postRadius * flangeThickness);
            areas.push_back(cornerArea);
            c1 += cornerLength / cornerArea;
            c2 += cornerLength / pow(cornerArea, 2);
        }
        return {c1, c2, *min_element(areas.begin(), areas.end())};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

// ROD core (ABT #933): a bare cylinder — no flanges, no return limb. The magnetic circuit closes
// entirely through the surrounding air, so this is the most open shape MKF carries
// (CoreType::OPEN_SHAPE). Letters, following the drum convention with the flanges removed:
//   A rod diameter, B rod length, H axial bore (optional, carries no flux).
// C is deliberately NOT read: on a drum C is the post OD, and a rod IS its own post, so a rod
// record that also states C would be describing two different diameters for one cylinder.
//
// IMPORTANT, same caveat as CorePieceDrum: get_shape_constants below is the FERRITE INTERNAL
// PATH ONLY (the cylinder, straight through). It is honest as the piece's partial parameters,
// but on its own it is not a magnetic circuit at all — there is no return path in the ferrite.
// Feeding it to the closed-circuit reluctance path would report an inductance many times too
// high. Rod magnetizing inductance routes to the demagnetising-factor model in
// MagnetizingInductance.cpp, which is what makes the open circuit's air return appear.
class CorePieceRod : public CorePiece {
  public:
    static double rod_diameter(std::map<std::string, double>& dimensions) {
        if (!dimensions.count("A") || dimensions["A"] <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "rod: dimension A (rod diameter) is missing or non-positive");
        }
        return dimensions["A"];
    }

    static double rod_length(std::map<std::string, double>& dimensions) {
        if (!dimensions.count("B") || dimensions["B"] <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "rod: dimension B (rod length) is missing or non-positive");
        }
        return dimensions["B"];
    }

    static double rod_bore(std::map<std::string, double>& dimensions) {
        double bore = dimensions.count("H") ? dimensions["H"] : 0.0;
        if (bore < 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "rod: bore H cannot be negative");
        }
        if (bore >= rod_diameter(dimensions)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "rod: bore H (" + std::to_string(bore) + ") is not smaller than the rod diameter A (" +
                std::to_string(rod_diameter(dimensions)) + "); nothing would be left to carry flux");
        }
        return bore;
    }

    static double rod_area(std::map<std::string, double>& dimensions) {
        double radius = rod_diameter(dimensions) / 2;
        double boreRadius = rod_bore(dimensions) / 2;
        return std::numbers::pi * (pow(radius, 2) - pow(boreRadius, 2));
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(rod_diameter(dimensions));
        set_depth(rod_diameter(dimensions));
        set_height(rod_length(dimensions));
    }

    // A rod has no groove, so it constrains the winding only in length: the wire is laid along
    // the cylinder over at most its full length B. The window WIDTH is the build height the
    // winding may occupy radially; a bare rod imposes no outer bound, so the natural, honest
    // statement of "unbounded" would be infinite. Instead the window is one rod radius deep,
    // which is the largest build that keeps the coil's outer diameter within 2A — beyond that
    // the solenoid is no longer a rod-core inductor in any useful sense. A caller who knows the
    // real build volume (a rod potted in a body) states it with an explicit bobbin, which is
    // what Coil::resolve_bobbin honours and what the inductance model reads its length from.
    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double diameter = rod_diameter(dimensions);
        WindingWindowElement windingWindow;
        windingWindow.set_height(rod_length(dimensions));
        windingWindow.set_width(diameter / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        // ABT #107 convention: coordinates[0] = window CENTRE (rod surface + half the build).
        windingWindow.set_coordinates(std::vector<double>({diameter / 2 + diameter / 4, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> columns;
        ColumnElement mainColumn;
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::ROUND);
        mainColumn.set_width(roundFloat(rod_diameter(dimensions)));
        mainColumn.set_depth(roundFloat(rod_diameter(dimensions)));
        mainColumn.set_height(roundFloat(rod_length(dimensions)));
        // Bore subtracted: the mounting hole carries no flux.
        mainColumn.set_area(roundFloat(rod_area(dimensions)));
        mainColumn.set_coordinates({0, 0, 0});
        columns.push_back(mainColumn);
        set_columns(columns);
    }

    // The cylinder, straight through: one uniform section, so c1 = l/A and c2 = l/A^2 and the
    // IEC 60205 reduction gives back le = l and Ae = A exactly.
    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double area = rod_area(dimensions);
        double length = rod_length(dimensions);
        return {length / area, length / pow(area, 2), area};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

// DRUM + RING ("shielded drum", ABT #366): a drum closed by a concentric shield ring,
// CoreType::PIECE_AND_PLATE. Letters: the drum convention above (A, A2, B, C, D, E, F, H)
// plus J ring OD, K ring ID, L ring height (defined with the MAS drumRing records, sourced
// from the ACME DR + SRI catalogue families).
//
// The magnetic circuit is CLOSED, but through two annular RADIAL clearance gaps between
// flange rim and ring bore — (K - A)/2 each — which dominate the reluctance and give these
// parts their soft saturation. Core::process_gap synthesizes them as GapType::RESIDUAL;
// they are deliberately NOT part of the ferrite constants here. get_shape_constants
// returns the WHOLE ferrite assembly (post + flanges + corners + ring + ring corners),
// following the piece-and-plate convention that the piece class reports the full circuit
// and Core::process_data doubles nothing.
class CorePieceDrumRing : public CorePieceDrum {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        if (dimensions["K"] <= dimensions["A"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumRing: ring inner diameter K (" + std::to_string(dimensions["K"]) +
                ") must exceed the flange OD A (" + std::to_string(dimensions["A"]) + ")");
        }
        // The envelope includes the ring: J >= A always (the ring wraps the flanges).
        set_width(dimensions["J"]);
        set_depth(dimensions["J"]);
        set_height(std::max(dimensions["B"], dimensions["L"]));
    }

    // ABT #576: the drum and its closing ring are routinely DIFFERENT grades (e.g. WE-DPC-6040 is
    // an ACME P47 MnZn drum, mu_i 3000, inside an ACME B45 NiZn ring, mu_i 450). The two sections
    // are therefore kept apart here so MagnetizingInductance can apply a permeability to each;
    // get_shape_constants() just sums them, so every single-material caller is unaffected.
    std::tuple<double, double, double, double, double> get_split_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pi = std::numbers::pi;
        double flangeRadius = dimensions["A"] / 2;
        double postRadius = dimensions["C"] / 2;
        double boreRadius = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] / 2 : 0.0;
        double grooveHeight = dimensions["E"];
        double ringOuterRadius = dimensions["J"] / 2;
        double ringInnerRadius = dimensions["K"] / 2;

        // Post + two spreading flanges + post corners: identical to the bare drum. DRUM SIDE.
        double postArea = pi * (pow(postRadius, 2) - pow(boreRadius, 2));
        std::vector<double> areas;
        areas.push_back(postArea);
        double drumC1 = grooveHeight / postArea;
        double drumC2 = grooveHeight / pow(postArea, 2);
        for (auto flangeThickness : {dimensions["D"], dimensions["F"]}) {
            drumC1 += 1.0 / (2 * pi * flangeThickness) * log(flangeRadius / postRadius);
            drumC2 += 1.0 / (2 * pow(pi * flangeThickness, 2)) * (flangeRadius - postRadius) / (flangeRadius * postRadius);
            areas.push_back(2 * pi * postRadius * flangeThickness);
            double s1 = postRadius - sqrt((pow(boreRadius, 2) + pow(postRadius, 2)) / 2);
            double cornerLength = pi / 4 * (2 * s1 + flangeThickness);
            double cornerArea = 0.5 * (postArea + 2 * pi * postRadius * flangeThickness);
            areas.push_back(cornerArea);
            drumC1 += cornerLength / cornerArea;
            drumC2 += cornerLength / pow(cornerArea, 2);
        }

        // Ring: axial annulus running between the two flange mid-planes (bounded by the
        // ring's own height when the ring is shorter than the drum). RING SIDE.
        double ringArea = pi * (pow(ringOuterRadius, 2) - pow(ringInnerRadius, 2));
        double ringLength = std::min(dimensions["L"],
                                     dimensions["B"] - (dimensions["D"] + dimensions["F"]) / 2);
        double ringC1 = ringLength / ringArea;
        double ringC2 = ringLength / pow(ringArea, 2);
        areas.push_back(ringArea);

        // Two radial->axial turns inside the ring wall, IEC 60205 clause 4.6 idiom: mean of
        // the entry band (ring bore over one flange thickness, where the gap flux lands) and
        // the axial annulus. These sit in the RING material, so they belong to the ring side.
        for (auto flangeThickness : {dimensions["D"], dimensions["F"]}) {
            double entryArea = 2 * pi * ringInnerRadius * flangeThickness;
            double ringWall = ringOuterRadius - ringInnerRadius;
            double cornerLength = pi / 4 * (ringWall / 2 + flangeThickness / 2);
            double cornerArea = 0.5 * (entryArea + ringArea);
            areas.push_back(cornerArea);
            ringC1 += cornerLength / cornerArea;
            ringC2 += cornerLength / pow(cornerArea, 2);
        }
        return {drumC1, drumC2, ringC1, ringC2, *min_element(areas.begin(), areas.end())};
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto [drumC1, drumC2, ringC1, ringC2, minimumArea] = get_split_constants();
        return {drumC1 + ringC1, drumC2 + ringC2, minimumArea};
    }

    std::optional<std::array<double, 4>> get_mixed_material_constants() override {
        auto [drumC1, drumC2, ringC1, ringC2, minimumArea] = get_split_constants();
        return std::array<double, 4>{drumC1, drumC2, ringC1, ringC2};
    }
};

// SEMI-SHIELDED DRUM (ABT #362): a wound drum whose winding is overcoated with MAGNETIC
// EPOXY (polymer binder + magnetic powder, mu ~3-15) acting as a low-permeability return
// shell — the third drum variant after bare (open circuit) and drumRing (ferrite ring).
// CoreType::PIECE_AND_PLATE. Letters: drum convention (A, A2, B, C, D, E, F, H) plus the
// shell OUTER ENVELOPE J width, K depth, L height (square LQS-class bodies have J = K).
// The glue is cast in contact with the flange rims — NO clearance gaps.
//
// get_shape_constants returns the GEOMETRIC single-material constants (vendor Ae/le/Ve
// convention: geometry alone). Because the circuit crosses two materials, the piece also
// exposes get_mixed_material_c1() — {c1 through the FERRITE drum, c1 through the SHELL} —
// which MagnetizingInductance uses with the drum material's mu and the magneticEpoxy
// coating material's mu (per-section reluctance). The shell material rides
// functionalDescription.coating {type: magneticEpoxy, material: <core material name>}.
class CorePieceDrumSemishielded : public CorePieceDrum {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        if (dimensions["J"] < dimensions["A"] || dimensions["K"] < dimensions["A"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumSemishielded: shell envelope J x K (" + std::to_string(dimensions["J"]) + " x " +
                std::to_string(dimensions["K"]) + ") cannot be smaller than the flange OD A (" +
                std::to_string(dimensions["A"]) + ")");
        }
        if (dimensions["L"] < dimensions["B"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumSemishielded: shell envelope height L (" + std::to_string(dimensions["L"]) +
                ") cannot be smaller than the drum height B (" + std::to_string(dimensions["B"]) +
                ") — the letters describe the FINISHED body");
        }
        set_width(dimensions["J"]);
        set_depth(dimensions["K"]);
        set_height(dimensions["L"]);
    }

    // The two halves of the circuit, computed once: {c1_ferrite, c2_ferrite, c1_shell,
    // c2_shell, minimumArea}. Ferrite = post + spreading flanges + post corners (the bare
    // drum's internal path). Shell = axial glue annulus outside the flange discs between the
    // flange mid-planes, plus the two rim (radial->axial) corners, clause-4.6 idiom.
    std::tuple<double, double, double, double, double> get_split_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pi = std::numbers::pi;
        double flangeRadius = dimensions["A"] / 2;
        double postRadius = dimensions["C"] / 2;
        double boreRadius = (dimensions.find("H") != dimensions.end()) ? dimensions["H"] / 2 : 0.0;
        double grooveHeight = dimensions["E"];

        double postArea = pi * (pow(postRadius, 2) - pow(boreRadius, 2));
        std::vector<double> areas;
        areas.push_back(postArea);
        double ferriteC1 = grooveHeight / postArea;
        double ferriteC2 = grooveHeight / pow(postArea, 2);
        for (auto flangeThickness : {dimensions["D"], dimensions["F"]}) {
            ferriteC1 += 1.0 / (2 * pi * flangeThickness) * log(flangeRadius / postRadius);
            ferriteC2 += 1.0 / (2 * pow(pi * flangeThickness, 2)) * (flangeRadius - postRadius) / (flangeRadius * postRadius);
            areas.push_back(2 * pi * postRadius * flangeThickness);
            double s1 = postRadius - sqrt((pow(boreRadius, 2) + pow(postRadius, 2)) / 2);
            double cornerLength = pi / 4 * (2 * s1 + flangeThickness);
            double cornerArea = 0.5 * (postArea + 2 * pi * postRadius * flangeThickness);
            areas.push_back(cornerArea);
            ferriteC1 += cornerLength / cornerArea;
            ferriteC2 += cornerLength / pow(cornerArea, 2);
        }

        // Shell: the glue outside the flange discs. Axial cross-section = envelope minus the
        // flange circle (conservative: only the glue radially beyond the flange OD conducts
        // the axial return; corner regions of a square envelope are included by construction).
        double shellArea = dimensions["J"] * dimensions["K"] - pi * pow(flangeRadius, 2);
        if (shellArea <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumSemishielded: shell envelope leaves no return cross-section outside the flanges");
        }
        double shellLength = dimensions["B"] - (dimensions["D"] + dimensions["F"]) / 2;
        double shellC1 = shellLength / shellArea;
        double shellC2 = shellLength / pow(shellArea, 2);
        areas.push_back(shellArea);
        double equivalentOuterRadius = sqrt(dimensions["J"] * dimensions["K"] / pi);
        double shellWall = equivalentOuterRadius - flangeRadius;
        for (auto flangeThickness : {dimensions["D"], dimensions["F"]}) {
            // Rim corner: flux leaves the flange rim (contact, no gap) and turns axial in
            // the glue. Assigned to the SHELL side: the turning region is glue.
            double entryArea = 2 * pi * flangeRadius * flangeThickness;
            double cornerLength = pi / 4 * (shellWall / 2 + flangeThickness / 2);
            double cornerArea = 0.5 * (entryArea + shellArea);
            areas.push_back(cornerArea);
            shellC1 += cornerLength / cornerArea;
            shellC2 += cornerLength / pow(cornerArea, 2);
        }
        double minimumArea = *min_element(areas.begin(), areas.end());
        return {ferriteC1, ferriteC2, shellC1, shellC2, minimumArea};
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto [ferriteC1, ferriteC2, shellC1, shellC2, minimumArea] = get_split_constants();
        return {ferriteC1 + shellC1, ferriteC2 + shellC2, minimumArea};
    }

    std::optional<std::array<double, 4>> get_mixed_material_constants() override {
        auto [ferriteC1, ferriteC2, shellC1, shellC2, minimumArea] = get_split_constants();
        return std::array<double, 4>{ferriteC1, ferriteC2, shellC1, shellC2};
    }
};

// MOLDED (ABT #357): a metal-composite molded inductor body (WE-MAPI / IHLP / XAL class) —
// a coil compression-molded inside a homogeneous low-permeability SMC block. The distributed
// gap lives in the MATERIAL (mu_eff ~15-40), so the piece is a single solid CLOSED circuit
// (CoreType::CLOSED_SHAPE) with NO discrete gaps: magnetically a pot core with a rectangular
// outer boundary.
//
// CLOSED_SHAPE IS LOAD-BEARING, and getting it wrong is silent (ABT #1002). A moulded core
// written as CoreType::TWO_PIECE_SET is mirrored: process_data doubles le and Ve, so a 4020 body
// reports le 8.0 mm instead of 4.0 and its inductance comes out HALF. Nothing throws, and the
// number stays plausible — so a permeability FITTED against that geometry silently absorbs the
// factor and comes out at twice the powder's. Four Würth grades did exactly that and were wrong
// by 2.00 each until an independent measurement caught them. If you fit a material through this
// piece, check the core type first; the geometry is not self-evidently right just because the
// inductance agrees with the part it was fitted on. Letters follow the pot-core (P/PM) convention: A body width, B body height
// (coil axis), C body depth, D coil-cavity height (internal height), E coil-cavity outer
// diameter, F coil-cavity inner diameter (the composite post under the coil bore).
//
// Sections, IEC 60205 general method: post (pi/4 D^2 over F), two plates ((B - F)/2 thick,
// radial spreading D/2 -> E/2), outer shell (block cross-section minus the cavity circle,
// axial mid-plate to mid-plate), plus clause-4.6 corner terms at both transitions. The
// rectangular outer is handled through the exact shell AREA and an area-equivalent outer
// radius for the corner lengths — the same first-order treatment RM/PQ pieces use for their
// non-round outlines. Validated forward against the REDEXPERT measured L(I) data in the
// ABT #357 phase-2 fit; this class only owns the geometry.
// DRUM_PLATE (ABT #996): the core a wire-wound chip common-mode choke is built on -- a RECTANGULAR
// drum, two end flanges with a post between them carrying the winding, closed by a flat ferrite lid
// bonded across the top. Phonon Meiwa's "Ferrite Core for Chip Common Mode Choke Coil" catalogue
// dimensions it L/W/H1/H2 and draws exactly that; the WE works documents name the two pieces
// outright, "dr core" and "i core", with glue "used to bond the core and cover".
//
// NOT derived from CorePieceDrum, though it is the same idea: that class computes a RADIAL winding
// window from a flange OD and a round post (width = (A - C)/2), and rectangular is the one thing
// that does not carry over. Everything else follows its structure closely on purpose.
//
// Letters, per the MAS schema commit that added the family:
//   A overall length along the post axis   D flange thickness, each end   J lid thickness
//   B drum height, lid excluded            E winding window length, A-2D  K lid length (dflt A)
//   C overall depth                        F post height                  L lid depth  (dflt C)
//                                          G post depth (dflt C)
//
// Flux path, and so the IEC 60205 sections: along the post between the flanges, up through one
// flange, across the lid, down the other flange. Corner terms follow the pi/4 (p + h) form the
// drum and E pieces already use.
class CorePieceDrumPlate : public CorePiece {
  public:
    static double lidLength(std::map<std::string, double>& dimensions) {
        auto it = dimensions.find("K");
        return (it != dimensions.end() && it->second > 0) ? it->second : dimensions["A"];
    }

    static double lidDepth(std::map<std::string, double>& dimensions) {
        auto it = dimensions.find("L");
        return (it != dimensions.end() && it->second > 0) ? it->second : dimensions["C"];
    }

    static double postDepth(std::map<std::string, double>& dimensions) {
        auto it = dimensions.find("G");
        return (it != dimensions.end() && it->second > 0) ? it->second : dimensions["C"];
    }

    static double windowLength(std::map<std::string, double>& dimensions) {
        auto it = dimensions.find("E");
        if (it != dimensions.end() && it->second > 0) {
            return it->second;
        }
        return dimensions["A"] - 2 * dimensions["D"];
    }

    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        // Validated rather than trusted, in the style drumSemishielded already uses: a window that
        // does not exist, or a lid larger than the body, is a data error and must not become a
        // plausible-looking core.
        if (windowLength(dimensions) <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumPlate: the two flanges (D = " + std::to_string(dimensions["D"]) +
                ") leave no winding window in the length A = " + std::to_string(dimensions["A"]));
        }
        if (dimensions["F"] <= 0 || dimensions["F"] >= dimensions["B"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumPlate: the post height F (" + std::to_string(dimensions["F"]) +
                ") must be positive and smaller than the drum height B (" +
                std::to_string(dimensions["B"]) + ")");
        }
        if (postDepth(dimensions) > dimensions["C"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumPlate: the post depth G (" + std::to_string(postDepth(dimensions)) +
                ") cannot exceed the overall depth C (" + std::to_string(dimensions["C"]) + ")");
        }
        if (dimensions["J"] <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "drumPlate: the lid thickness J must be positive; a drumPlate without a lid is a "
                "drum, and its magnetic circuit does not close");
        }
        // The FINISHED body, drumSemishielded's convention: the lid is part of what you hold.
        set_width(std::max(dimensions["A"], lidLength(dimensions)));
        set_depth(std::max(dimensions["C"], lidDepth(dimensions)));
        set_height(dimensions["B"] + dimensions["J"]);
    }

    void process_winding_window() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        WindingWindowElement windingWindow;
        // Mirrors the drum: the window's HEIGHT is its extent along the post, and its WIDTH is the
        // room the turns have off the post before they reach the body -- radial there, vertical
        // here, because the lid closes over the top.
        windingWindow.set_height(windowLength(dimensions));
        windingWindow.set_width((dimensions["B"] - dimensions["F"]) / 2);
        windingWindow.set_area(windingWindow.get_height().value() * windingWindow.get_width().value());
        // ABT #107 convention: coordinates[0] is the window CENTRE, not its inner edge.
        windingWindow.set_coordinates(std::vector<double>(
            {dimensions["F"] / 2 + (dimensions["B"] - dimensions["F"]) / 4, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> columns;
        ColumnElement mainColumn;
        // The post is what the winding encircles, so it is what sets the mean turn length and, in
        // turn, the bobbin MKF derives as half this column plus the bobbin thickness.
        mainColumn.set_type(ColumnType::CENTRAL);
        mainColumn.set_shape(ColumnShape::RECTANGULAR);
        mainColumn.set_width(roundFloat(postDepth(dimensions)));
        mainColumn.set_depth(roundFloat(dimensions["F"]));
        mainColumn.set_height(roundFloat(windowLength(dimensions)));
        mainColumn.set_area(roundFloat(postDepth(dimensions) * dimensions["F"]));
        mainColumn.set_coordinates({0, 0, 0});
        columns.push_back(mainColumn);
        set_columns(columns);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pi = std::numbers::pi;
        double window = windowLength(dimensions);
        double post = postDepth(dimensions);
        double postArea = post * dimensions["F"];
        double flangeArea = dimensions["C"] * dimensions["D"];
        double lidArea = lidDepth(dimensions) * dimensions["J"];
        // From the post's mid-height up to the middle of the lid: the length the flux travels in
        // each flange before it turns into the lid.
        double flangeRun = dimensions["B"] - dimensions["F"] / 2 + dimensions["J"] / 2;

        std::vector<double> areas = {postArea, flangeArea, lidArea};
        double c1 = window / postArea;
        double c2 = window / pow(postArea, 2);
        // Two flanges and the lid, each with the corner that joins it to what came before.
        for (auto [sectionLength, sectionArea, previousArea] :
             {std::make_tuple(flangeRun, flangeArea, postArea),
              std::make_tuple(flangeRun, flangeArea, postArea),
              std::make_tuple(window, lidArea, flangeArea)}) {
            c1 += sectionLength / sectionArea;
            c2 += sectionLength / pow(sectionArea, 2);
            double cornerLength = pi / 4 * (sqrt(previousArea) + sqrt(sectionArea)) / 2;
            double cornerArea = 0.5 * (previousArea + sectionArea);
            areas.push_back(cornerArea);
            c1 += cornerLength / cornerArea;
            c2 += cornerLength / pow(cornerArea, 2);
        }
        return {c1, c2, *min_element(areas.begin(), areas.end())};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

class CorePieceMolded : public CorePiece {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        if (dimensions["E"] <= dimensions["F"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "molded: cavity outer diameter E (" + std::to_string(dimensions["E"]) +
                ") must exceed the cavity inner diameter F (" + std::to_string(dimensions["F"]) + ")");
        }
        if (dimensions["D"] >= dimensions["B"]) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "molded: cavity height D (" + std::to_string(dimensions["D"]) +
                ") must be smaller than the body height B (" + std::to_string(dimensions["B"]) + ")");
        }
        if (dimensions["E"] >= std::min(dimensions["A"], dimensions["C"])) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                "molded: cavity outer diameter E (" + std::to_string(dimensions["E"]) +
                ") must fit inside the body footprint A x C");
        }
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
        // ABT #107 convention: coordinates[0] = window CENTRE (post edge + half width).
        windingWindow.set_coordinates(std::vector<double>({dimensions["F"] / 2 + (dimensions["E"] - dimensions["F"]) / 4, 0}));
        set_winding_window(windingWindow);
    }

    void process_columns() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        std::vector<ColumnElement> columns;
        ColumnElement mainColumn;
        mainColumn.set_type(ColumnType::CENTRAL);
        // An OVAL post, when the shape states one. Several Wurth moulded families are wound on an
        // obround pole and their sheets give both axes ('1.5/1.0'); read as one diameter both the
        // pole area and the winding window (E - F)/2 come out wrong. F is the smaller axis and F2
        // the larger -- the nomenclature the oblong-column families already use (F width, F2
        // depth) -- and the area is theirs too: a stadium is a rectangle with a half-disc on each
        // end. Absent F2, or F2 <= F, nothing changes and the post stays round.
        auto f2 = dimensions.find("F2");
        bool oblongPost = f2 != dimensions.end() && f2->second > dimensions["F"];
        if (oblongPost) {
            mainColumn.set_shape(ColumnShape::OBLONG);
            mainColumn.set_width(roundFloat(dimensions["F"]));
            mainColumn.set_depth(roundFloat(f2->second));
            mainColumn.set_area(roundFloat(std::numbers::pi / 4 * pow(dimensions["F"], 2) +
                                           (f2->second - dimensions["F"]) * dimensions["F"]));
        }
        else {
            mainColumn.set_shape(ColumnShape::ROUND);
            mainColumn.set_width(roundFloat(dimensions["F"]));
            mainColumn.set_depth(roundFloat(dimensions["F"]));
            mainColumn.set_area(roundFloat(std::numbers::pi / 4 * pow(dimensions["F"], 2)));
        }
        mainColumn.set_height(roundFloat(dimensions["D"]));
        mainColumn.set_coordinates({0, 0, 0});
        columns.push_back(mainColumn);
        // The return shell: one annular-equivalent lateral column wrapping the cavity.
        ColumnElement shellColumn;
        shellColumn.set_type(ColumnType::LATERAL);
        shellColumn.set_shape(ColumnShape::IRREGULAR);
        double shellArea = dimensions["A"] * dimensions["C"] - std::numbers::pi / 4 * pow(dimensions["E"], 2);
        double equivalentOuterRadius = sqrt(dimensions["A"] * dimensions["C"] / std::numbers::pi);
        double shellWall = equivalentOuterRadius - dimensions["E"] / 2;
        shellColumn.set_width(roundFloat(shellWall));
        shellColumn.set_area(roundFloat(shellArea));
        shellColumn.set_depth(roundFloat(shellArea / shellWall));
        shellColumn.set_height(roundFloat(dimensions["D"]));
        shellColumn.set_coordinates({roundFloat(dimensions["E"] / 2 + shellWall / 2), 0, 0});
        columns.push_back(shellColumn);
        set_columns(columns);
    }

    // The IEC 60205 walk, section by section, each assigned to the REGION of the body that is
    // pressed around it (ABT #1002). The WE molded families are made in up to three pressings --
    // SUB (base plate), COR (the post the coil sits on) and COV (the cover over the coil) in the
    // MXGI list of parts, Inner (post) / Outer (the rest) in the MAPI/MAIA ones -- and each
    // pressing can be a different powder. Regions, in the order the magnetic circuit crosses
    // them from the post upward, which is the order functionalDescription.material lists them:
    //   post   the axial run of the composite post through the cavity height D
    //   cover  the top plate, its post corner, the outer shell down to the top face of the base
    //          plate, and the top rim corner
    //   base   the bottom plate, its post corner, the base plate's share of the shell (from its
    //          mid-plane, where the IEC walk turns, up to its top face) and the bottom rim corner
    // A body pressed from one powder sums the regions back into get_shape_constants() below, so
    // the single-material and the per-region views can never disagree.
    std::optional<std::vector<RegionShapeConstants>> get_region_shape_constants() override {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double pi = std::numbers::pi;
        double postRadius = dimensions["F"] / 2;
        double cavityRadius = dimensions["E"] / 2;
        double cavityHeight = dimensions["D"];
        double plateThickness = (dimensions["B"] - dimensions["D"]) / 2;
        double shellArea = dimensions["A"] * dimensions["C"] - pi * pow(cavityRadius, 2);
        double equivalentOuterRadius = sqrt(dimensions["A"] * dimensions["C"] / pi);
        double shellWall = equivalentOuterRadius - cavityRadius;
        double postArea = pi * pow(postRadius, 2);

        RegionShapeConstants post{"post", cavityHeight / postArea, cavityHeight / pow(postArea, 2), postArea};

        // One plate: radial spreading disc post edge -> cavity edge, the post corner (clause
        // 4.6, same idiom as the drum flanges) and the radial->axial rim corner into the shell.
        auto plate = [&](std::string name) {
            RegionShapeConstants region{name, 0, 0, std::numeric_limits<double>::max()};
            region.c1 += 1.0 / (2 * pi * plateThickness) * log(cavityRadius / postRadius);
            region.c2 += 1.0 / (2 * pow(pi * plateThickness, 2)) * (cavityRadius - postRadius) / (cavityRadius * postRadius);
            double plateEntryArea = 2 * pi * postRadius * plateThickness;
            region.minimumArea = std::min(region.minimumArea, plateEntryArea);
            double s1 = postRadius - sqrt(pow(postRadius, 2) / 2);
            double postCornerLength = pi / 4 * (2 * s1 + plateThickness);
            double postCornerArea = 0.5 * (postArea + plateEntryArea);
            region.minimumArea = std::min(region.minimumArea, postCornerArea);
            region.c1 += postCornerLength / postCornerArea;
            region.c2 += postCornerLength / pow(postCornerArea, 2);
            double rimEntryArea = 2 * pi * cavityRadius * plateThickness;
            double rimCornerLength = pi / 4 * (shellWall / 2 + plateThickness / 2);
            double rimCornerArea = 0.5 * (rimEntryArea + shellArea);
            region.minimumArea = std::min(region.minimumArea, rimCornerArea);
            region.c1 += rimCornerLength / rimCornerArea;
            region.c2 += rimCornerLength / pow(rimCornerArea, 2);
            return region;
        };
        RegionShapeConstants cover = plate("cover");
        RegionShapeConstants base = plate("base");

        // Outer shell: the axial run between the two plate mid-planes, B - t. The base plate
        // owns the half-thickness of it that lies inside the base pressing; the cover owns the
        // rest, which is the side wall of the cover pressing.
        double baseShellLength = plateThickness / 2;
        double coverShellLength = dimensions["B"] - plateThickness - baseShellLength;
        base.c1 += baseShellLength / shellArea;
        base.c2 += baseShellLength / pow(shellArea, 2);
        base.minimumArea = std::min(base.minimumArea, shellArea);
        cover.c1 += coverShellLength / shellArea;
        cover.c2 += coverShellLength / pow(shellArea, 2);
        cover.minimumArea = std::min(cover.minimumArea, shellArea);

        return std::vector<RegionShapeConstants>{post, cover, base};
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto regions = get_region_shape_constants().value();
        double c1 = 0;
        double c2 = 0;
        double minimumArea = std::numeric_limits<double>::max();
        for (auto& region : regions) {
            c1 += region.c1;
            c2 += region.c2;
            minimumArea = std::min(minimumArea, region.minimumArea);
        }
        return {c1, c2, minimumArea};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

// PQI: a PQ half closed by a flat I plate.
//
// IEC 60205:2016 DOES cover this case, contrary to a first reading of its clause titles. Every
// clause is headed "Pair of X-cores", but clause 5.12 (PQ) carries TWO figures -- Figure 15
// "PQ-cores" and Figure 16 "PLT(plate)-cores" -- and the same PLT pairing appears in 5.10 (EL),
// 5.11 (ER) and 5.14 (E planar). The plate is treated as its own piece, so a PQ+plate is the PQ
// pair's section list with the second half's yoke replaced by the plate.
//
// Concretely that means ONE change from CorePiecePq: the window is one piece tall, so the legs
// contribute D rather than 2D. The yoke and corner sections are split between the PQ's own
// thickness (B - D) and the plate's (B2); since a2, a9 and a10 are all proportional to the
// thickness, the plate's share scales by B2/(B - D). When B2 == B - D the split collapses back to
// CorePiecePq's own single terms exactly, which is the case for every PQI record MAS holds.
//
// Crucially the yoke keeps CorePiecePq's LOGARITHMIC radial-spreading section (a2/l2) rather than
// a straight run: a PQ carries flux outward from a round centre post to the lateral legs, and a
// straight yoke of length (E-F)/2 is not that path. An earlier attempt that used a straight yoke
// got Ae to 1.1% but le 69% long -- Ae alone is not a sufficient test of this geometry.
//
// VALIDATED against TDK's published planar PQI data:
//     PQI 16/7.8   computed le 19.48 mm / Ae 41.73 mm2 / Ve 812.9 mm3
//                  published Ae 41.8 mm2, Ve 815 mm3 (so le = Ve/Ae = 19.50 mm)   -> within 0.3%
// Further published rows exist for cross-checking if the record set grows: ACME PQI35F/29
// (le 75.20, Ae 181.20), PQI35.2 (68.71, 183.66) and PQI40B/28/14.6/5 (51.29, 169.59) -- note the
// last of those is internally inconsistent in the source (le*Ae = 8698 against a printed Ve of
// 8871.51), so treat it with suspicion.
class CorePiecePqi : public CorePiecePq {
  public:
    void process_extra_data() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        set_width(dimensions["A"]);
        set_height(dimensions["B"] + dimensions["B2"]);
        set_depth(dimensions["C"]);
    }

    std::tuple<double, double, double> get_shape_constants() {
        auto dimensions = flatten_dimensions(get_shape().get_dimensions().value());
        double A = dimensions["A"];
        double B = dimensions["B"];
        double C = dimensions["C"];
        double D = dimensions["D"];
        double E = dimensions["E"];
        double F = dimensions["F"];
        double G = dimensions["G"];
        double plateThickness = dimensions["B2"];
        double J, L;
        if ((dimensions.find("J") == dimensions.end()) || (dimensions["J"] == 0)) {
            // Same documented approximation as CorePiecePq for shapes lacking J/L.
            J = F / 2;
            L = F + (C - F) / 3;
        }
        else {
            J = dimensions["J"];
            L = dimensions["L"];
        }

        double pi = std::numbers::pi;
        double yokeThickness = B - D;
        double beta = acos(G / E);
        double alpha = atan(L / J);
        double I = E * sin(beta);
        double a7 = 1. / 8 * (beta * pow(E, 2) - alpha * pow(F, 2) + G * L - J * I);
        double a8 = pi / 16 * (pow(E, 2) - pow(F, 2));
        double K = a7 / a8;
        double lmin = (E - F) / 2;
        double lmax = sqrt(pow(E, 2) + pow(F, 2) - 2 * E * F * cos(alpha - beta)) / 2;
        double f = (lmin + lmax) / (2 * lmin);

        // Lateral legs and round centre post, each spanning ONE window height.
        double a1 = C * (A - G) - beta * pow(E, 2) / 2 + 1. / 2 * G * I;
        double a3 = pi / 4 * pow(F, 2);
        // Radial-spreading yoke, split between the PQ half and the plate.
        double l2 = f * E * F / (E - F) * pow(log(E / F), 2);
        double a2Piece = pi * K * E * F * yokeThickness / (E - F) * log(E / F);
        double a2Plate = pi * K * E * F * plateThickness / (E - F) * log(E / F);
        // Corner sections, likewise split; a9/a10 scale with the thickness they belong to.
        double a9Piece = 2 * alpha * F * yokeThickness;
        double a10Piece = 2 * beta * E * yokeThickness;
        double a9Plate = 2 * alpha * F * plateThickness;
        double a10Plate = 2 * beta * E * plateThickness;

        std::vector<double> lengths;
        std::vector<double> areas;
        lengths.push_back(D);                areas.push_back(a1);
        lengths.push_back(D);                areas.push_back(a3);
        lengths.push_back(l2 / 2);           areas.push_back(a2Piece);
        lengths.push_back(l2 / 2);           areas.push_back(a2Plate);
        lengths.push_back(pi / 8 * (yokeThickness + A / 2 - E / 2));
        areas.push_back((a1 + a10Piece) / 2);
        lengths.push_back(pi / 8 * (plateThickness + A / 2 - E / 2));
        areas.push_back((a1 + a10Plate) / 2);
        lengths.push_back(pi / 8 * (yokeThickness + (1 - 1. / sqrt(2)) * F));
        areas.push_back((a3 + a9Piece) / 2);
        lengths.push_back(pi / 8 * (plateThickness + (1 - 1. / sqrt(2)) * F));
        areas.push_back((a3 + a9Plate) / 2);

        double c1 = 0, c2 = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            c1 += lengths[i] / areas[i];
            c2 += lengths[i] / pow(areas[i], 2);
        }
        return {c1, c2, *min_element(areas.begin(), areas.end())};
    }

    std::tuple<double, double, double> get_shape_constants_iec63182() override {
        auto [c1, c2, minimumArea] = get_shape_constants();
        return {pow(c1, 2) / c2, c1 / c2, minimumArea};
    }
};

class CorePieceEpq : public CorePieceEp {};
class CorePieceEpw : public CorePieceEp {};
class CorePieceEpt : public CorePieceEp {};
class CorePieceLep : public CorePieceEp {};

// The shape families CorePiece::factory below actually constructs. Kept as data so the
// "available options" in its error can be GENERATED rather than hand-typed: the literal
// that used to live there had gone stale twice over — it stopped at "T, C" while the
// factory had long handled UI/EI/EER/EF/EPC/EPQ/EPW/EPT/LEP/PQI, and it never listed the
// drum families at all (DRUM/DRUM_RING/DRUM_SEMISHIELDED/ROD/MOLDED, added in 52d99c35).
// A caller reading that message reasonably concluded those cores were unimplemented
// upstream and went looking for a feature gap that did not exist (ABT #975).
//
// Add a family here in the same edit that adds its branch below, and the message can
// never disagree with the code again.
static constexpr CoreShapeFamily kSupportedShapeFamilies[] = {
    CoreShapeFamily::E,          CoreShapeFamily::EC,         CoreShapeFamily::EFD,
    CoreShapeFamily::EL,         CoreShapeFamily::EP,         CoreShapeFamily::EPX,
    CoreShapeFamily::LP,         CoreShapeFamily::EQ,         CoreShapeFamily::ER,
    CoreShapeFamily::ETD,        CoreShapeFamily::P,          CoreShapeFamily::PLANAR_E,
    CoreShapeFamily::PLANAR_EL,  CoreShapeFamily::PLANAR_ER,  CoreShapeFamily::PM,
    CoreShapeFamily::PQ,         CoreShapeFamily::RM,         CoreShapeFamily::U,
    CoreShapeFamily::UR,         CoreShapeFamily::UT,         CoreShapeFamily::T,
    CoreShapeFamily::C,          CoreShapeFamily::EER,        CoreShapeFamily::EF,
    CoreShapeFamily::EPC,        CoreShapeFamily::UI,         CoreShapeFamily::EI,
    CoreShapeFamily::DRUM,       CoreShapeFamily::DRUM_RING,  CoreShapeFamily::DRUM_SEMISHIELDED,
    CoreShapeFamily::DRUM_PLATE, CoreShapeFamily::H,
    CoreShapeFamily::ROD,        CoreShapeFamily::MOLDED,     CoreShapeFamily::PQI,
    CoreShapeFamily::EPQ,        CoreShapeFamily::EPW,        CoreShapeFamily::EPT,
    CoreShapeFamily::LEP,
    CoreShapeFamily::DS,         CoreShapeFamily::HS,         CoreShapeFamily::RS,
};

// The dimensions each family's geometry actually READS, keyed by family rather than
// discovered from whatever shapes happen to be published. Every entry was extracted from
// the CorePiece subclass above (including the reads it inherits and does not override),
// and covers the REQUIRED reads only: a key behind a `dimensions.find(...)` guard — drum's
// A2, a toroid's R/r0, PM's alpha — is optional geometry and reaches the caller through the
// catalogue union in get_shape_family_dimensions instead, so declaring it here would put an
// input field in front of every user for a dimension almost no part has.
//
// This exists because the question "what does this family need?" was being answered by
// scanning the shape database, which returns NOTHING for a family that ships no bare-core
// record (ABT #1007). Eleven buildable families were selectable in the builder and then
// offered no dimension fields at all, so a custom shape in them could not be defined.
//
// Add a family here in the same edit that adds its CorePiece subclass. The guard is
// Test_Family_Dimensions_Are_Declared_For_Every_Supported_Family.
static const std::map<CoreShapeFamily, std::vector<std::string>> kFamilyRequiredDimensions = {
    {CoreShapeFamily::E,                   {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EC,                  {"A", "B", "C", "D", "E", "F", "s"}},
    {CoreShapeFamily::EFD,                 {"A", "B", "C", "D", "E", "F", "F2", "K"}},
    {CoreShapeFamily::EL,                  {"A", "B", "C", "D", "E", "F", "F2"}},
    {CoreShapeFamily::EP,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EPX,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::LP,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EQ,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::ER,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::ETD,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::P,                   {"A", "B", "D", "E", "F", "G", "H"}},
    {CoreShapeFamily::PLANAR_E,            {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::PLANAR_EL,           {"A", "B", "C", "D", "E", "F", "F2"}},
    {CoreShapeFamily::PLANAR_ER,           {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::PM,                  {"A", "B", "D", "E", "F", "G", "H", "b", "t"}},
    {CoreShapeFamily::PQ,                  {"A", "B", "C", "D", "E", "F", "G", "L"}},
    {CoreShapeFamily::RM,                  {"A", "B", "C", "D", "E", "F", "G", "H", "J"}},
    {CoreShapeFamily::U,                   {"A", "B", "C", "D"}},
    {CoreShapeFamily::UR,                  {"A", "B", "C", "D", "G", "H"}},
    {CoreShapeFamily::UT,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::T,                   {"A", "B", "C"}},
    {CoreShapeFamily::C,                   {"A", "B", "C", "D", "E"}},
    {CoreShapeFamily::EER,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EF,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EPC,                 {"A", "B", "C", "D", "E", "F", "F2"}},
    {CoreShapeFamily::UI,                  {"A", "B", "B2", "C", "D"}},
    {CoreShapeFamily::EI,                  {"A", "B", "B2", "C", "D", "E", "F"}},
    {CoreShapeFamily::DRUM,                {"A", "B", "C", "D", "E", "F"}},
    // H is an alias of DRUM (ABT #277): same piece, so the same dimensions, by construction.
    {CoreShapeFamily::H,                   {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::DRUM_RING,           {"A", "B", "C", "D", "E", "F", "J", "K", "L"}},
    {CoreShapeFamily::DRUM_SEMISHIELDED,   {"A", "B", "C", "D", "E", "F", "J", "K", "L"}},
    {CoreShapeFamily::DRUM_PLATE,          {"A", "B", "C", "D", "F", "G"}},
    {CoreShapeFamily::ROD,                 {"A", "B", "H"}},
    {CoreShapeFamily::MOLDED,              {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::PQI,                 {"A", "B", "B2", "C", "D", "E", "F", "G", "L"}},
    {CoreShapeFamily::EPQ,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EPW,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::EPT,                 {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::LEP,                 {"A", "B", "C", "D", "E", "F"}},
    // Slab cores (ABT #263): H (the bore) is optional -- DS has a solid post -- and G is derived.
    {CoreShapeFamily::DS,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::HS,                  {"A", "B", "C", "D", "E", "F"}},
    {CoreShapeFamily::RS,                  {"A", "B", "C", "D", "E", "F"}},
};

std::vector<std::string> get_core_shape_family_required_dimensions(CoreShapeFamily family) {
    auto entry = kFamilyRequiredDimensions.find(family);
    if (entry == kFamilyRequiredDimensions.end()) {
        throw std::runtime_error("No dimensions declared for shape family: " +
                                 std::string{magic_enum::enum_name(family)});
    }
    return entry->second;
}

std::vector<CoreShapeFamily> get_supported_core_shape_families() {
    return {std::begin(kSupportedShapeFamilies), std::end(kSupportedShapeFamilies)};
}

// Answered from the same array, so a family can never be "supported" for the error message
// and unsupported for the guard. This used to be a 37-case switch that mirrored the array
// by hand: a third copy of the list that had to be edited in lockstep with the array and
// with factory()'s dispatch, and the array's own comment already records what happens when
// one copy is forgotten (ABT #975).
bool CorePiece::is_family_supported(CoreShapeFamily family) {
    return std::find(std::begin(kSupportedShapeFamilies), std::end(kSupportedShapeFamilies),
                     family) != std::end(kSupportedShapeFamilies);
}

CoreShapeFamily canonical_core_shape_family(CoreShapeFamily family) {
    if (family == CoreShapeFamily::H) {
        return CoreShapeFamily::DRUM;
    }
    return family;
}

CoreShape canonicalize_core_shape_family(CoreShape shape) {
    shape.set_family(canonical_core_shape_family(shape.get_family()));
    return shape;
}

static std::string supported_shape_family_list() {
    std::string list = "{";
    for (size_t index = 0; index < std::size(kSupportedShapeFamilies); ++index) {
        if (index > 0) {
            list += ", ";
        }
        list += magic_enum::enum_name(kSupportedShapeFamilies[index]);
    }
    return list + "}";
}

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
    else if (family == CoreShapeFamily::EER) {
        auto piece = std::make_shared<CorePieceEer>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EF) {
        auto piece = std::make_shared<CorePieceEf>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPC) {
        auto piece = std::make_shared<CorePieceEpc>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::UI) {
        auto piece = std::make_shared<CorePieceUi>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EI) {
        auto piece = std::make_shared<CorePieceEi>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::DS) {
        auto piece = std::make_shared<CorePieceDs>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::HS) {
        auto piece = std::make_shared<CorePieceHs>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::RS) {
        auto piece = std::make_shared<CorePieceRs>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::DRUM || family == CoreShapeFamily::H) {
        // H is the drum under another name (ABT #277); see canonical_core_shape_family.
        auto piece = std::make_shared<CorePieceDrum>();
        piece->set_shape(canonicalize_core_shape_family(shape));
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::DRUM_RING) {
        auto piece = std::make_shared<CorePieceDrumRing>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::DRUM_PLATE) {
        auto piece = std::make_shared<CorePieceDrumPlate>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::DRUM_SEMISHIELDED) {
        auto piece = std::make_shared<CorePieceDrumSemishielded>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::ROD) {
        auto piece = std::make_shared<CorePieceRod>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::MOLDED) {
        auto piece = std::make_shared<CorePieceMolded>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::PQI) {
        auto piece = std::make_shared<CorePiecePqi>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPQ) {
        auto piece = std::make_shared<CorePieceEpq>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPW) {
        auto piece = std::make_shared<CorePieceEpw>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::EPT) {
        auto piece = std::make_shared<CorePieceEpt>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else if (family == CoreShapeFamily::LEP) {
        auto piece = std::make_shared<CorePieceLep>();
        piece->set_shape(shape);
        if (process) piece->process();
        return piece;
    }
    else
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Unknown shape family: " + to_string(family) + " for shape '" + shape.get_name().value_or("<unnamed>") + "', available options are: " + supported_shape_family_list());
}

// ============================================================================
// Thermal Surface Area Calculations
// ============================================================================

double CorePiece::get_column_right_face_area(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }
    
    auto& column = columns[columnIndex];
    auto shape = column.get_shape();
    double columnHeight = column.get_height();
    double columnDepth = column.get_depth();
    double columnWidth = column.get_width();
    
    if (shape == ColumnShape::RECTANGULAR || shape == ColumnShape::IRREGULAR) {
        // For rectangular: height * depth
        return columnHeight * columnDepth;
    }
    else if (shape == ColumnShape::ROUND || shape == ColumnShape::OBLONG) {
        // For round: curved surface area = π * d * h (for half cylinder facing winding)
        double diameter = columnWidth;
        return std::numbers::pi * diameter * columnHeight / 2.0;  // Half circumference * height
    }
    else {
        // Default to rectangular approximation
        return columnHeight * columnDepth;
    }
}

double CorePiece::get_column_top_face_area(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }
    
    auto& column = columns[columnIndex];
    auto shape = column.get_shape();
    double columnWidth = column.get_width();
    double columnDepth = column.get_depth();
    
    if (shape == ColumnShape::RECTANGULAR || shape == ColumnShape::IRREGULAR) {
        // For rectangular: width * depth
        return columnWidth * columnDepth;
    }
    else if (shape == ColumnShape::ROUND) {
        // For round: circular area = π * r²
        double radius = columnWidth / 2.0;
        return std::numbers::pi * radius * radius;
    }
    else if (shape == ColumnShape::OBLONG) {
        // For oblong: rectangle + two half circles
        // Approximate as width * depth (conservative)
        return columnWidth * columnDepth;
    }
    else {
        return columnWidth * columnDepth;
    }
}

double CorePiece::get_column_bottom_face_area(size_t columnIndex) {
    // Same as top face for symmetrical cores
    return get_column_top_face_area(columnIndex);
}

double CorePiece::get_yoke_interior_face_area(bool isTopYoke) {
    // The interior yoke face is the face facing the winding window
    // For E-shaped cores, this is the bottom of top yoke or top of bottom yoke
    
    double windingWindowWidth = get_winding_window_width();
    double depth = get_depth();
    
    // For most E-core families, the yoke interior face spans the winding window
    // Area = winding window width * depth (for one side of the core piece)
    return windingWindowWidth * depth;
}

double CorePiece::get_yoke_exterior_face_area(bool isTopYoke) {
    // The exterior yoke face faces away from the winding window
    // This is typically the outer edge of the core
    
    double totalWidth = get_width();
    double centralColumnWidth = (columns.size() > 0) ? columns[0].get_width() : 0;
    double depth = get_depth();
    
    // For E-cores: total width minus central column gives the yoke span
    // We approximate the exterior face as (totalWidth - centralColumnWidth) / 2 * depth
    // This accounts for one side of the yoke
    double yokeSpan = (totalWidth - centralColumnWidth) / 2.0;
    return yokeSpan * depth;
}

double CorePiece::get_yoke_right_face_area(bool isTopYoke) {
    // The right face of the yoke is the vertical face at the end of the yoke
    // For E-cores, this connects to the lateral column
    
    double yokeHeight = (height - windingWindow.get_height().value()) / 2.0;  // Top/bottom yoke thickness
    double depth = get_depth();
    
    // Area = yoke thickness * depth
    return yokeHeight * depth;
}

double CorePiece::get_winding_window_height() {
    return windingWindow.get_height().value();
}

double CorePiece::get_winding_window_width() {
    return windingWindow.get_width().value();
}

double CorePiece::get_column_width(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }
    return columns[columnIndex].get_width();
}

double CorePiece::get_column_depth(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }
    return columns[columnIndex].get_depth();
}

ColumnShape CorePiece::get_column_shape(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }
    return columns[columnIndex].get_shape();
}


// ============================================================================
// Volume-Proportional Core Loss Distribution
// ============================================================================

double CorePiece::calculate_column_cross_section(size_t columnIndex) {
    if (columnIndex >= columns.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Column index out of range");
    }

    const auto& col = columns[columnIndex];

    // Prefer the pre-computed area from the schema when available.
    double area = col.get_area();
    if (area > 0) return area;

    double w = col.get_width();
    double d = col.get_depth();

    switch (col.get_shape()) {
        case ColumnShape::ROUND:
            return std::numbers::pi / 4.0 * w * d;

        case ColumnShape::OBLONG: {
            double r = std::min(w, d) / 2.0;
            return w * d - (4.0 - std::numbers::pi) * r * r;
        }

        case ColumnShape::RECTANGULAR:
        case ColumnShape::IRREGULAR:
        default:
            return w * d;
    }
}

CorePartVolumes CorePiece::calculate_core_part_volumes() {
    CorePartVolumes v;

    if (columns.empty()) {
        return v;
    }

    // =========================================================================
    // FULL-CORE volumes.
    //
    // The thermal model meshes the RIGHT HALF of the core at HALF DEPTH.
    // The depth symmetry (Z = 0 plane) means the meshed front-half
    // thermally represents the FULL right half.  So this is a HALF model,
    // not a quarter model.
    //
    // We return full-core volumes here; the caller (calculate_core_loss_fractions)
    // applies the correct /2 factors to get the right-half fractions.
    //
    // The modeled right half contains:
    //   - Half of the central column    (right portion)
    //   - One lateral column            (the right one, 1 of 2)
    //   - Right half of top yoke        (center to right edge)
    //   - Right half of bottom yoke     (center to right edge)
    //
    // Total modeled = V_full / 2  -->  losses should sum to 50%.
    // =========================================================================

    // --- Central column (ONE, full depth, full cross-section) ---
    double mainArea   = calculate_column_cross_section(0);
    double mainHeight = columns[0].get_height();
    v.centralColumn   = mainArea * mainHeight;

    // --- Lateral columns: store ONE lateral's full volume ---
    // The full core has TWO (left + right); we model only the right one.
    if (columns.size() > 1) {
        double latArea   = calculate_column_cross_section(1);
        double latHeight = columns[1].get_height();
        v.lateralColumn  = latArea * latHeight;
    }

    // --- Yokes (top / bottom, each spans full width x full depth) ---
    // Yoke thickness matches Temperature.cpp: mainColumn.get_width() / 2
    double yokeThickness = columns[0].get_width() / 2.0;
    double yokeVolume    = width * yokeThickness * depth;

    v.topYoke    = yokeVolume;
    v.bottomYoke = yokeVolume;

    return v;
}

CoreLossFractions CorePiece::calculate_core_loss_fractions() {
    CoreLossFractions f;

    CorePartVolumes v = calculate_core_part_volumes();

    // Full-core total volume.
    // lateralColumn = ONE lateral's volume; the full core has TWO.
    double fullCoreVolume = v.centralColumn
                          + 2.0 * v.lateralColumn
                          + v.topYoke
                          + v.bottomYoke;

    if (fullCoreVolume < 1e-18) {
        // A zero/degenerate core volume means the upstream geometry is broken;
        // substituting made-up loss fractions would mask that bug.
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "Degenerate core volume (" + std::to_string(fullCoreVolume) + " m^3) while computing core loss fractions");
    }

    // =========================================================================
    // FULL-core loss fractions (ABT #461).
    //
    // The thermal model builds its core nodes at HALF depth (symmetry), but at the
    // end of the convection build it DOUBLES every core surface back to the full
    // component (Temperature.cpp, "Half-core symmetry correction"), so each node
    // stands for BOTH symmetric halves of its part. The losses must follow the same
    // convention or energy is not conserved: with the old right-half fractions
    // (sum = 0.5) every concentric component shed only Pcore/2 + Pcu — a core with
    // full-size cooling surfaces but half its heat, running ~2x too cold in the
    // core-loss-driven part of its temperature rise. Measured on
    // concentric_transformer (Pcore 1.26 W dominant): network input was 0.83 W of
    // a real 1.46 W, and that missing ~43% of the heat had been silently
    // compensating the missing radiation-to-ambient path this ticket fixed.
    //
    //   Central column:  the one full column                    -> V_central
    //   Lateral node:    BOTH lateral legs (node represents 2)  -> 2 * V_one_lateral
    //   Top yoke:        full yoke                              -> V_top_yoke
    //   Bottom yoke:     full yoke                              -> V_bot_yoke
    //
    // Sum of all fractions = 1.0 (full component, matching the doubled surfaces).
    // =========================================================================

    f.centralColumn =        v.centralColumn / fullCoreVolume;
    f.lateralColumn = (2.0 * v.lateralColumn) / fullCoreVolume;
    f.topYoke       =        v.topYoke       / fullCoreVolume;
    f.bottomYoke    =        v.bottomYoke    / fullCoreVolume;

    return f;
}

inline void from_json(const json& j, CorePiece& x) {
    x.set_columns(j.at("columns").get<std::vector<ColumnElement>>());
    x.set_depth(j.at("depth").get<double>());
    x.set_height(j.at("height").get<double>());
    x.set_width(j.at("width").get<double>());
    x.set_shape(j.at("shape").get<CoreShape>());
    x.set_winding_window(j.at("windingWindow").get<WindingWindowElement>());
    x.set_partial_effective_parameters(j.at("partialEffectiveParameters").get<EffectiveParameters>());
}

inline void to_json(json& j, const CorePiece& x) {
    j = json::object();
    j["columns"] = x.get_columns();
    j["depth"] = x.get_depth();
    j["height"] = x.get_height();
    j["width"] = x.get_width();
    j["shape"] = x.get_shape();
    j["windingWindow"] = x.get_winding_window();
    j["partialEffectiveParameters"] = x.get_partial_effective_parameters();
}

} // namespace OpenMagnetics
