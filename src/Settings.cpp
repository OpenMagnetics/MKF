#include "Settings.h"

namespace OpenMagnetics {
    Settings* Settings::settings_= nullptr;

    Settings *Settings::GetInstance()
    {
        if(settings_==nullptr){
            settings_ = new Settings();
        }
        return settings_;
    }

    void Settings::reset() {
        _allowMarginTape = true;
        _allowInsulatedWire = true;
        _fillCoilSectionsWithMarginTape = false;
        _windEvenIfNotFit = false;
        _useOnlyCoresInStock = true;

        _painterNumberPointsX = 25;
        _painterNumberPointsY = 50;
        _painterMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _painterMode;
        _painterLogarithmicScale = false;
        _painterIncludeFringing = true;
        _painterMaximumValueColorbar = std::nullopt;
        _painterMinimumValueColorbar = std::nullopt;
        _painterColorFerrite = "0x8F7b7c7d";
        _painterColorBobbin = "0x8F1b1b1b";
        _painterColorCopper = "0x8Fb87333";
        _painterColorInsulation = "0x18E37E00";
        _painterColorMargin = "0x8Ffff05b";

        _magneticFieldNumberPointsX = 25;
        _magneticFieldNumberPointsY = 50;
        _magneticFieldMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _magneticFieldIncludeFringing = true;
    } 

    bool Settings::get_allow_margin_tape() const {
        return _allowMarginTape;
    } 
    void Settings::set_allow_margin_tape(bool value) {
        _allowMarginTape = value;
    }

    bool Settings::get_allow_insulated_wire() const {
        return _allowInsulatedWire;
    } 
    void Settings::set_allow_insulated_wire(bool value) {
        _allowInsulatedWire = value;
    }

    bool Settings::get_fill_coil_sections_with_margin_tape() const {
        return _fillCoilSectionsWithMarginTape;
    } 
    void Settings::set_fill_coil_sections_with_margin_tape(bool value) {
        _fillCoilSectionsWithMarginTape = value;
    } 

    bool Settings::get_wind_even_if_not_fit() const {
        return _windEvenIfNotFit;
    } 
    void Settings::set_wind_even_if_not_fit(bool value) {
        _windEvenIfNotFit = value;
    } 

    bool Settings::get_use_only_cores_in_stock() const {
        return _useOnlyCoresInStock;
    } 
    void Settings::set_use_only_cores_in_stock(bool value) {
        _useOnlyCoresInStock = value;
    }

    size_t Settings::get_painter_number_points_x() const {
        return _painterNumberPointsX;
    }
    void Settings::set_painter_number_points_x(size_t value) {
        _painterNumberPointsX = value;
    }

    size_t Settings::get_painter_number_points_y() const {
        return _painterNumberPointsY;
    }
    void Settings::set_painter_number_points_y(size_t value) {
        _painterNumberPointsY = value;
    }

    Painter::PainterModes Settings::get_painter_mode() const {
        return _painterMode;
    }
    void Settings::set_painter_mode(Painter::PainterModes value) {
        _painterMode = value;
    }

    bool Settings::get_painter_logarithmic_scale() const {
        return _painterLogarithmicScale;
    }
    void Settings::set_painter_logarithmic_scale(bool value) {
        _painterLogarithmicScale = value;
    }

    bool Settings::get_painter_include_fringing() const {
        return _painterIncludeFringing;
    }
    void Settings::set_painter_include_fringing(bool value) {
        _painterIncludeFringing = value;
    }

    std::optional<double> Settings::get_painter_maximum_value_colorbar() const {
        return _painterMaximumValueColorbar;
    }
    void Settings::set_painter_maximum_value_colorbar(std::optional<double> value) {
        _painterMaximumValueColorbar = value;
    }

    std::optional<double> Settings::get_painter_minimum_value_colorbar() const {
        return _painterMinimumValueColorbar;
    }
    void Settings::set_painter_minimum_value_colorbar(std::optional<double> value) {
        _painterMinimumValueColorbar = value;
    }

    std::string Settings::get_painter_color_ferrite() const {
        return _painterColorFerrite;
    }
    void Settings::set_painter_color_ferrite(std::string value) {
        _painterColorFerrite = value;
    }

    std::string Settings::get_painter_color_bobbin() const {
        return _painterColorBobbin;
    }
    void Settings::set_painter_color_bobbin(std::string value) {
        _painterColorBobbin = value;
    }

    std::string Settings::get_painter_color_copper() const {
        return _painterColorCopper;
    }
    void Settings::set_painter_color_copper(std::string value) {
        _painterColorCopper = value;
    }

    std::string Settings::get_painter_color_insulation() const {
        return _painterColorInsulation;
    }
    void Settings::set_painter_color_insulation(std::string value) {
        _painterColorInsulation = value;
    }

    std::string Settings::get_painter_color_margin() const {
        return _painterColorMargin;
    }
    void Settings::set_painter_color_margin(std::string value) {
        _painterColorMargin = value;
    }

    int Settings::get_painter_mirroring_dimension() const {
        return _painterMirroringDimension;
    }
    void Settings::set_painter_mirroring_dimension(int value) {
        _painterMirroringDimension = value;
    }

    size_t Settings::get_magnetic_field_number_points_x() const {
        return _magneticFieldNumberPointsX;
    }
    void Settings::set_magnetic_field_number_points_x(size_t value) {
        _magneticFieldNumberPointsX = value;
    }

    size_t Settings::get_magnetic_field_number_points_y() const {
        return _magneticFieldNumberPointsY;
    }
    void Settings::set_magnetic_field_number_points_y(size_t value) {
        _magneticFieldNumberPointsY = value;
    }

    int Settings::get_magnetic_field_mirroring_dimension() const {
        return _magneticFieldMirroringDimension;
    }
    void Settings::set_magnetic_field_mirroring_dimension(int value) {
        _magneticFieldMirroringDimension = value;
    }

    bool Settings::get_magnetic_field_include_fringing() const {
        return _magneticFieldIncludeFringing;
    }
    void Settings::set_magnetic_field_include_fringing(bool value) {
        _magneticFieldIncludeFringing = value;
    }

} // namespace OpenMagnetics