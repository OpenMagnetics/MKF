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

        _useToroidalCores = true;
        _useConcentricCores = true;

        _inputsTrimHarmonics = true;

        _inputsNumberPointsSampledWaveforms = Constants().numberPointsSampledWaveforms;

        _magnetizingInductanceIncludeAirInductance = false;

        _coilAllowMarginTape = true;
        _coilAllowInsulatedWire = true;
        _coilFillSectionsWithMarginTape = false;
        _coilWindEvenIfNotFit = false;
        _coilDelimitAndCompact = true;
        _coilTryRewind = true;
        _coilIncludeAdditionalCoordinates = true;
        _coilEqualizeMargins = true;

        _useOnlyCoresInStock = true;

        _painterNumberPointsX = 25;
        _painterNumberPointsY = 50;
        _painterMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _painterMode = Painter::PainterModes::CONTOUR;
        _painterLogarithmicScale = false;
        _painterIncludeFringing = true;
        _painterDrawSpacer = true;
        _painterMaximumValueColorbar = std::nullopt;
        _painterMinimumValueColorbar = std::nullopt;
        _painterColorFerrite = "0x007b7c7d";
        _painterColorBobbin = "0x8f1b1b1b";
        _painterColorCopper = "0x00b87333";
        _painterColorInsulation = "0x18539796";
        _painterColorMargin = "0x00fff05b";
        _painterColorSpacer = "0x003b3b3b";

        _magneticFieldNumberPointsX = 25;
        _magneticFieldNumberPointsY = 50;
        _magneticFieldMirroringDimension = Defaults().magneticFieldMirroringDimension;
        _magneticFieldIncludeFringing = true;
        
        _coilAdviserMaximumNumberWires = 100;
        _coreIncludeStacks = true;
        _coreIncludeDistributedGaps = true;

        _verbose = false;

    }

    bool Settings::get_verbose() const {
        return _verbose;
    }
    void Settings::set_verbose(bool value) {
        _verbose = value;
    }

    bool Settings::get_use_toroidal_cores() const {
        return _useToroidalCores;
    }
    void Settings::set_use_toroidal_cores(bool value) {
        _useToroidalCores = value;
    }

    bool Settings::get_use_concentric_cores() const {
        return _useConcentricCores;
    }
    void Settings::set_use_concentric_cores(bool value) {
        _useConcentricCores = value;
    }

    bool Settings::get_inputs_trim_harmonics() const {
        return _inputsTrimHarmonics;
    }
    void Settings::set_inputs_trim_harmonics(bool value) {
        _inputsTrimHarmonics = value;
    }

    size_t Settings::get_inputs_number_points_sampled_waveforms() const {
        return _inputsNumberPointsSampledWaveforms;
    }
    void Settings::set_inputs_number_points_sampled_waveforms(size_t value) {
        _inputsNumberPointsSampledWaveforms = value;
    }

    bool Settings::get_magnetizing_inductance_include_air_inductance() const {
        return _magnetizingInductanceIncludeAirInductance;
    }
    void Settings::set_magnetizing_inductance_include_air_inductance(bool value) {
        _magnetizingInductanceIncludeAirInductance = value;
    }

    bool Settings::get_coil_allow_margin_tape() const {
        return _coilAllowMarginTape;
    }
    void Settings::set_coil_allow_margin_tape(bool value) {
        _coilAllowMarginTape = value;
    }

    bool Settings::get_coil_allow_insulated_wire() const {
        return _coilAllowInsulatedWire;
    }
    void Settings::set_coil_allow_insulated_wire(bool value) {
        _coilAllowInsulatedWire = value;
    }

    bool Settings::get_coil_fill_sections_with_margin_tape() const {
        return _coilFillSectionsWithMarginTape;
    }
    void Settings::set_coil_fill_sections_with_margin_tape(bool value) {
        _coilFillSectionsWithMarginTape = value;
    }

    bool Settings::get_coil_wind_even_if_not_fit() const {
        return _coilWindEvenIfNotFit;
    }
    void Settings::set_coil_wind_even_if_not_fit(bool value) {
        _coilWindEvenIfNotFit = value;
    }

    bool Settings::get_coil_delimit_and_compact() const {
        return _coilDelimitAndCompact;
    }
    void Settings::set_coil_delimit_and_compact(bool value) {
        _coilDelimitAndCompact = value;
    }

    bool Settings::get_coil_try_rewind() const {
        return _coilTryRewind;
    }
    void Settings::set_coil_try_rewind(bool value) {
        _coilTryRewind = value;
    }

    bool Settings::get_coil_include_additional_coordinates() const {
        return _coilIncludeAdditionalCoordinates;
    }
    void Settings::set_coil_include_additional_coordinates(bool value) {
        _coilIncludeAdditionalCoordinates = value;
    }

    bool Settings::get_coil_equalize_margins() const {
        return _coilEqualizeMargins;
    }
    void Settings::set_coil_equalize_margins(bool value) {
        _coilEqualizeMargins = value;
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

    bool Settings::get_painter_painter_draw_spacer() const {
        return _painterDrawSpacer;
    }
    void Settings::set_painter_painter_draw_spacer(bool value) {
        _painterDrawSpacer = value;
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

    std::string Settings::get_painter_color_spacer() const {
        return _painterColorSpacer;
    }
    void Settings::set_painter_color_spacer(std::string value) {
        _painterColorSpacer = value;
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

    size_t Settings::get_coil_adviser_maximum_number_wires() const {
        return _coilAdviserMaximumNumberWires;
    }
    void Settings::set_coil_adviser_maximum_number_wires(size_t value) {
        _coilAdviserMaximumNumberWires = value;
    }

    bool Settings::get_core_include_stacks() const {
        return _coreIncludeStacks;
    }
    void Settings::set_core_include_stacks(bool value) {
        _coreIncludeStacks = value;
    }

    bool Settings::get_core_include_distributed_gaps() const {
        return _coreIncludeDistributedGaps;
    }
    void Settings::set_core_include_distributed_gaps(bool value) {
        _coreIncludeDistributedGaps = value;
    }

} // namespace OpenMagnetics