#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11_json/pybind11_json.hpp"
#include "json.hpp"
#include <MAS.hpp>
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Painter.h"
#include "support/Settings.h"

using namespace MAS;
using json = nlohmann::json;
namespace py = pybind11;

std::string plot_turns(json magneticJson, std::string outFile) {
    try {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        return outFile;
    }
    catch(const std::runtime_error& re) {
        return re.what();
    }
    catch(const std::exception& ex) {
        return ex.what();
    }
    catch(...) {
        return "Unknown failure occurred. Possible memory corruption";
    }
}

std::string plot_field(json magneticJson, json operatingPointJson, std::string outFile) {
    try {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        OperatingPoint operatingPoint(operatingPointJson);
        OpenMagnetics::Painter painter(outFile, true);
        painter.paint_magnetic_field(operatingPoint, magnetic);
        painter.paint_core(magnetic);
        painter.paint_bobbin(magnetic);
        painter.paint_coil_turns(magnetic);
        painter.export_svg();
        return outFile;
    }
    catch(const std::runtime_error& re) {
        return re.what();
    }
    catch(const std::exception& ex) {
        return ex.what();
    }
    catch(...) {
        return "Unknown failure occurred. Possible memory corruption";
    }
}

std::string plot_wire(json wireJson, std::string outFile, std::string cci_coords_path) {
    try {
        OpenMagnetics::Wire wire(wireJson);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_wire(wire);
        painter.export_svg();
        return outFile;
    }
    catch(const std::runtime_error& re) {
        return re.what();
    }
    catch(const std::exception& ex) {
        return ex.what();
    }
    catch(...) {
        return "Unknown failure occurred. Possible memory corruption";
    }
}

std::string plot_current_density(json wireJson, json operatingPointJson, std::string outFile) {
    try {
        OpenMagnetics::Wire wire(wireJson);
        OperatingPoint operatingPoint(operatingPointJson);
        OpenMagnetics::Painter painter(outFile);
        painter.paint_wire_with_current_density(wire, operatingPoint);
        painter.export_svg();
        return outFile;
    }
    catch(const std::runtime_error& re) {
        return re.what();
    }
    catch(const std::exception& ex) {
        return ex.what();
    }
    catch(...) {
        return "Unknown failure occurred. Possible memory corruption";
    }
}

json get_settings() {
    try {
        json settingsJson;
        settingsJson["coilAllowMarginTape"] = OpenMagnetics::settings.get_coil_allow_margin_tape();
        settingsJson["coilAllowInsulatedWire"] = OpenMagnetics::settings.get_coil_allow_insulated_wire();
        settingsJson["coilFillSectionsWithMarginTape"] = OpenMagnetics::settings.get_coil_fill_sections_with_margin_tape();
        settingsJson["coilWindEvenIfNotFit"] = OpenMagnetics::settings.get_coil_wind_even_if_not_fit();
        settingsJson["coilDelimitAndCompact"] = OpenMagnetics::settings.get_coil_delimit_and_compact();
        settingsJson["coilTryRewind"] = OpenMagnetics::settings.get_coil_try_rewind();
        settingsJson["coilOnlyOneTurnPerLayerInContiguousRectangular"] = OpenMagnetics::settings.get_coil_only_one_turn_per_layer_in_contiguous_rectangular();
        settingsJson["useOnlyCoresInStock"] = OpenMagnetics::settings.get_use_only_cores_in_stock();
        settingsJson["coilMaximumLayersPlanar"] = OpenMagnetics::settings.get_coil_maximum_layers_planar();

        settingsJson["painterNumberPointsX"] = OpenMagnetics::settings.get_painter_number_points_x();
        settingsJson["painterNumberPointsY"] = OpenMagnetics::settings.get_painter_number_points_y();
        settingsJson["painterMode"] = OpenMagnetics::settings.get_painter_mode();
        settingsJson["painterLogarithmicScale"] = OpenMagnetics::settings.get_painter_logarithmic_scale();
        settingsJson["painterIncludeFringing"] = OpenMagnetics::settings.get_painter_include_fringing();
        if (OpenMagnetics::settings.get_painter_maximum_value_colorbar()) {
            settingsJson["painterMaximumValueColorbar"] = OpenMagnetics::settings.get_painter_maximum_value_colorbar().value();
        }
        if (OpenMagnetics::settings.get_painter_minimum_value_colorbar()) {
            settingsJson["painterMinimumValueColorbar"] = OpenMagnetics::settings.get_painter_minimum_value_colorbar().value();
        }
        settingsJson["painterColorLines"] = OpenMagnetics::settings.get_painter_color_lines();
        settingsJson["painterColorText"] = OpenMagnetics::settings.get_painter_color_text();
        settingsJson["painterColorFerrite"] = OpenMagnetics::settings.get_painter_color_ferrite();
        settingsJson["painterColorBobbin"] = OpenMagnetics::settings.get_painter_color_bobbin();
        settingsJson["painterColorCopper"] = OpenMagnetics::settings.get_painter_color_copper();
        settingsJson["painterColorInsulation"] = OpenMagnetics::settings.get_painter_color_insulation();
        settingsJson["painterColorMargin"] = OpenMagnetics::settings.get_painter_color_margin();
        settingsJson["painterMirroringDimension"] = OpenMagnetics::settings.get_painter_mirroring_dimension();
        settingsJson["magneticFieldNumberPointsX"] = OpenMagnetics::settings.get_magnetic_field_number_points_x();
        settingsJson["magneticFieldNumberPointsY"] = OpenMagnetics::settings.get_magnetic_field_number_points_y();
        settingsJson["magneticFieldIncludeFringing"] = OpenMagnetics::settings.get_magnetic_field_include_fringing();
        settingsJson["magneticFieldMirroringDimension"] = OpenMagnetics::settings.get_magnetic_field_mirroring_dimension();
        settingsJson["painterSimpleLitz"] = OpenMagnetics::settings.get_painter_simple_litz();
        settingsJson["painterAdvancedLitz"] = OpenMagnetics::settings.get_painter_advanced_litz();
        settingsJson["painterCciCoordinatesPath"] = OpenMagnetics::settings.get_painter_cci_coordinates_path();
        return settingsJson;
    }
    catch (const std::exception &exc) {
        json exception;
        exception["data"] = "Exception: " + std::string{exc.what()};
        return exception;
    }
}

void set_settings(json settingsJson) {
    OpenMagnetics::settings.set_coil_allow_margin_tape(settingsJson["coilAllowMarginTape"]);
    OpenMagnetics::settings.set_coil_allow_insulated_wire(settingsJson["coilAllowInsulatedWire"]);
    OpenMagnetics::settings.set_coil_fill_sections_with_margin_tape(settingsJson["coilFillSectionsWithMarginTape"]);
    OpenMagnetics::settings.set_coil_wind_even_if_not_fit(settingsJson["coilWindEvenIfNotFit"]);
    OpenMagnetics::settings.set_coil_delimit_and_compact(settingsJson["coilDelimitAndCompact"]);
    OpenMagnetics::settings.set_coil_try_rewind(settingsJson["coilTryRewind"]);
    OpenMagnetics::settings.set_coil_only_one_turn_per_layer_in_contiguous_rectangular(settingsJson["coilOnlyOneTurnPerLayerInContiguousRectangular"]);
    OpenMagnetics::settings.set_coil_maximum_layers_planar(settingsJson["coilMaximumLayersPlanar"]);

    OpenMagnetics::settings.set_painter_mode(settingsJson["painterMode"]);
    OpenMagnetics::settings.set_use_only_cores_in_stock(settingsJson["useOnlyCoresInStock"]);
    OpenMagnetics::settings.set_painter_number_points_x(settingsJson["painterNumberPointsX"]);
    OpenMagnetics::settings.set_painter_number_points_y(settingsJson["painterNumberPointsY"]);
    OpenMagnetics::settings.set_painter_logarithmic_scale(settingsJson["painterLogarithmicScale"]);
    OpenMagnetics::settings.set_painter_include_fringing(settingsJson["painterIncludeFringing"]);
    if (settingsJson.contains("painterMaximumValueColorbar")) {
        OpenMagnetics::settings.set_painter_maximum_value_colorbar(settingsJson["painterMaximumValueColorbar"]);
    }
    if (settingsJson.contains("painterMinimumValueColorbar")) {
        OpenMagnetics::settings.set_painter_minimum_value_colorbar(settingsJson["painterMinimumValueColorbar"]);
    }
    OpenMagnetics::settings.set_painter_color_lines(settingsJson["painterColorLines"]);
    OpenMagnetics::settings.set_painter_color_text(settingsJson["painterColorText"]);
    OpenMagnetics::settings.set_painter_color_ferrite(settingsJson["painterColorFerrite"]);
    OpenMagnetics::settings.set_painter_color_bobbin(settingsJson["painterColorBobbin"]);
    OpenMagnetics::settings.set_painter_color_copper(settingsJson["painterColorCopper"]);
    OpenMagnetics::settings.set_painter_color_insulation(settingsJson["painterColorInsulation"]);
    OpenMagnetics::settings.set_painter_color_margin(settingsJson["painterColorMargin"]);
    OpenMagnetics::settings.set_painter_mirroring_dimension(settingsJson["painterMirroringDimension"]);
    OpenMagnetics::settings.set_magnetic_field_number_points_x(settingsJson["magneticFieldNumberPointsX"]);
    OpenMagnetics::settings.set_magnetic_field_number_points_y(settingsJson["magneticFieldNumberPointsY"]);
    OpenMagnetics::settings.set_magnetic_field_include_fringing(settingsJson["magneticFieldIncludeFringing"]);
    OpenMagnetics::settings.set_magnetic_field_mirroring_dimension(settingsJson["magneticFieldMirroringDimension"]);
    OpenMagnetics::settings.set_painter_simple_litz(settingsJson["painterSimpleLitz"]);
    OpenMagnetics::settings.set_painter_advanced_litz(settingsJson["painterAdvancedLitz"]);
    OpenMagnetics::settings.set_painter_cci_coordinates_path(settingsJson["painterCciCoordinatesPath"]);
}

void reset_settings() {
    OpenMagnetics::settings.reset();
}

PYBIND11_MODULE(PyMKF, m) {
    m.doc() = "OpenMagnetics Python bindings (legacy PyMKF compatibility)";
    m.def("plot_turns", &plot_turns, "Generate SVG plot of individual turns");
    m.def("plot_field", &plot_field, "Generate SVG plot of magnetic field");
    m.def("plot_wire", &plot_wire, "Generate SVG plot of wire cross-section");
    m.def("plot_current_density", &plot_current_density, "Generate SVG plot of current density distribution");
    m.def("set_settings", &set_settings, "Configure simulation settings");
    m.def("get_settings", &get_settings, "Retrieve current simulation settings");
    m.def("reset_settings", &reset_settings, "Reset settings to defaults");
}
