#pragma once
#include "MAS.hpp"
#include "Defaults.h"
#include "support/Painter.h"

namespace OpenMagnetics {

class Settings
{
    protected:
        Settings() {}

        static Settings* settings_;
        std::string selfFilePath = __FILE__;

        bool _useToroidalCores = true;
        bool _useConcentricCores = true;

        bool _inputsTrimHarmonics = true;
        size_t _inputsNumberPointsSampledWaveforms = Constants().numberPointsSampledWaveforms;

        bool _magnetizingInductanceIncludeAirInductance = false;

        bool _coilAllowMarginTape = true;
        bool _coilAllowInsulatedWire = true;
        bool _coilFillSectionsWithMarginTape = false;
        bool _coilWindEvenIfNotFit = false;
        bool _coilDelimitAndCompact = true;
        bool _coilTryRewind = true;
        bool _coilIncludeAdditionalCoordinates = true;
        bool _coilEqualizeMargins = true;
        size_t _coilMaximumLayersPlanar = 4;

        bool _useOnlyCoresInStock = true;

        size_t _painterNumberPointsX = 25;
        size_t _painterNumberPointsY = 50;
        int _painterMirroringDimension = Defaults().magneticFieldMirroringDimension;
        Painter::PainterModes _painterMode = Painter::PainterModes::CONTOUR;
        bool _painterLogarithmicScale = false;
        bool _painterIncludeFringing = true;
        bool _painterDrawSpacer = true;
        bool _painterSimpleLitz = true;
        bool _painterAdvancedLitz = false;
        std::optional<double> _painterMaximumValueColorbar = std::nullopt;
        std::optional<double> _painterMinimumValueColorbar = std::nullopt;
        std::string _painterColorFerrite = "0x7b7c7d";
        std::string _painterColorBobbin = "0x539796";
        std::string _painterColorCopper = "0xb87333";
        std::string _painterColorInsulation = "0xfff05b";
        std::string _painterColorEnamel = "0xc63032";
        std::string _painterColorFEP = "0x252525";
        std::string _painterColorETFE = "0xb42811";
        std::string _painterColorTCA = "0x696969";
        std::string _painterColorPFA = "0xedbe1c";
        std::string _painterColorSilk = "0xe7e7e8";
        std::string _painterColorMargin = "0xfff05b";
        std::string _painterColorSpacer = "0x3b3b3b";
        std::string _painterColorLines = "0x010000";
        std::string _painterColorText = "0x000000";
        std::string _painterColorCurrentDensity = "0x0892D0";
        std::string _painterCciCoordinatesPath = std::string{selfFilePath}.substr(0, std::string{selfFilePath}.rfind("/")).append("/../cci_coords/coordinates/");

        size_t _magneticFieldNumberPointsX = 25;
        size_t _magneticFieldNumberPointsY = 50;
        int _magneticFieldMirroringDimension = Defaults().magneticFieldMirroringDimension;
        bool _magneticFieldIncludeFringing = true;

        size_t _coilAdviserMaximumNumberWires = 100;
        bool _coreAdviserIncludeStacks = true;
        bool _coreAdviserIncludeDistributedGaps = true;
        bool _coreAdviserIncludeMargin = false;


        bool _wireAdviserIncludePlanar = false;
        bool _wireAdviserIncludeFoil = false;
        bool _wireAdviserIncludeRectangular = true;
        bool _wireAdviserIncludeLitz = true;
        bool _wireAdviserIncludeRound = true;
        bool _wireAdviserAllowRectangularInToroidalCores = false;

        bool _harmonicAmplitudeThresholdQuickMode = true;
        double _harmonicAmplitudeThreshold = Defaults().harmonicAmplitudeThreshold;


        std::vector<CoreLossesModels> _coreLossesModelNames = {Defaults().coreLossesModelDefault, CoreLossesModels::PROPRIETARY, CoreLossesModels::LOSS_FACTOR, CoreLossesModels::STEINMETZ, CoreLossesModels::ROSHEN};

        bool _verbose = false;

    public:
        bool _debug = false;
        Settings(Settings &other) = delete;
        void operator=(const Settings &) = delete;

        static Settings *GetInstance();

        void reset();

        bool get_verbose() const;
        void set_verbose(bool value);

        bool get_use_toroidal_cores() const;
        void set_use_toroidal_cores(bool value);

        bool get_use_concentric_cores() const;
        void set_use_concentric_cores(bool value);

        bool get_inputs_trim_harmonics() const;
        void set_inputs_trim_harmonics(bool value);

        size_t get_inputs_number_points_sampled_waveforms() const;
        void set_inputs_number_points_sampled_waveforms(size_t value);

        bool get_magnetizing_inductance_include_air_inductance() const;
        void set_magnetizing_inductance_include_air_inductance(bool value);

        bool get_coil_allow_margin_tape() const;
        void set_coil_allow_margin_tape(bool value);

        bool get_coil_allow_insulated_wire() const;
        void set_coil_allow_insulated_wire(bool value);

        bool get_coil_fill_sections_with_margin_tape() const;
        void set_coil_fill_sections_with_margin_tape(bool value);

        bool get_coil_wind_even_if_not_fit() const;
        void set_coil_wind_even_if_not_fit(bool value);

        bool get_coil_delimit_and_compact() const;
        void set_coil_delimit_and_compact(bool value);

        bool get_coil_try_rewind() const;
        void set_coil_try_rewind(bool value);

        bool get_coil_include_additional_coordinates() const;
        void set_coil_include_additional_coordinates(bool value);

        bool get_coil_equalize_margins() const;
        void set_coil_equalize_margins(bool value);

        size_t get_coil_maximum_layers_planar() const;
        void set_coil_maximum_layers_planar(size_t value);

        bool get_use_only_cores_in_stock() const;
        void set_use_only_cores_in_stock(bool value);

        size_t get_painter_number_points_x() const;
        void set_painter_number_points_x(size_t value);

        size_t get_painter_number_points_y() const;
        void set_painter_number_points_y(size_t value);

        Painter::PainterModes get_painter_mode() const;
        void set_painter_mode(Painter::PainterModes value);

        bool get_painter_logarithmic_scale() const;
        void set_painter_logarithmic_scale(bool value);

        bool get_painter_include_fringing() const;
        void set_painter_include_fringing(bool value);

        bool get_painter_draw_spacer() const;
        void set_painter_draw_spacer(bool value);

        bool get_painter_simple_litz() const;
        void set_painter_simple_litz(bool value);

        bool get_painter_advanced_litz() const;
        void set_painter_advanced_litz(bool value);

        std::optional<double> get_painter_maximum_value_colorbar() const;
        void set_painter_maximum_value_colorbar(std::optional<double> value);

        std::optional<double> get_painter_minimum_value_colorbar() const;
        void set_painter_minimum_value_colorbar(std::optional<double> value);

        std::string get_painter_color_ferrite() const;
        void set_painter_color_ferrite(std::string value);

        std::string get_painter_color_bobbin() const;
        void set_painter_color_bobbin(std::string value);

        std::string get_painter_color_copper() const;
        void set_painter_color_copper(std::string value);

        std::string get_painter_color_insulation() const;
        void set_painter_color_insulation(std::string value);

        std::string get_painter_color_margin() const;
        void set_painter_color_margin(std::string value);

        std::string get_painter_color_spacer() const;
        void set_painter_color_spacer(std::string value);

        std::string get_painter_color_lines() const;
        void set_painter_color_lines(std::string value);

        std::string get_painter_color_text() const;
        void set_painter_color_text(std::string value);

        std::string get_painter_cci_coordinates_path() const;
        void set_painter_cci_coordinates_path(std::string value);

        std::string get_painter_color_enamel() const;
        void set_painter_color_enamel(std::string value);

        std::string get_painter_color_fep() const;
        void set_painter_color_fep(std::string value);

        std::string get_painter_color_etfe() const;
        void set_painter_color_etfe(std::string value);

        std::string get_painter_color_tca() const;
        void set_painter_color_tca(std::string value);

        std::string get_painter_color_pfa() const;
        void set_painter_color_pfa(std::string value);

        std::string get_painter_color_silk() const;
        void set_painter_color_silk(std::string value);

        std::string get_painter_color_current_density() const;
        void set_painter_color_current_density(std::string value);

        int get_painter_mirroring_dimension() const;
        void set_painter_mirroring_dimension(int value);

        size_t get_magnetic_field_number_points_x() const;
        void set_magnetic_field_number_points_x(size_t value);

        size_t get_magnetic_field_number_points_y() const;
        void set_magnetic_field_number_points_y(size_t value);

        bool get_magnetic_field_include_fringing() const;
        void set_magnetic_field_include_fringing(bool value);

        int get_magnetic_field_mirroring_dimension() const;
        void set_magnetic_field_mirroring_dimension(int value);

        size_t get_coil_adviser_maximum_number_wires() const;
        void set_coil_adviser_maximum_number_wires(size_t value);

        bool get_core_adviser_include_stacks() const;
        void set_core_adviser_include_stacks(bool value);

        bool get_core_adviser_include_distributed_gaps() const;
        void set_core_adviser_include_distributed_gaps(bool value);

        bool get_core_adviser_include_margin() const;
        void set_core_adviser_include_margin(bool value);

        bool get_wire_adviser_include_planar() const;
        void set_wire_adviser_include_planar(bool value);

        bool get_wire_adviser_include_foil() const;
        void set_wire_adviser_include_foil(bool value);

        bool get_wire_adviser_include_rectangular() const;
        void set_wire_adviser_include_rectangular(bool value);

        bool get_wire_adviser_include_litz() const;
        void set_wire_adviser_include_litz(bool value);

        bool get_wire_adviser_include_round() const;
        void set_wire_adviser_include_round(bool value);

        bool get_wire_adviser_allow_rectangular_in_toroidal_cores() const;
        void set_wire_adviser_allow_rectangular_in_toroidal_cores(bool value);

        bool get_harmonic_amplitude_threshold_quick_mode() const;
        void set_harmonic_amplitude_threshold_quick_mode(bool value);

        double get_harmonic_amplitude_threshold() const;
        void set_harmonic_amplitude_threshold(double value);

        std::vector<CoreLossesModels> get_core_losses_model_names() const;
        void set_core_losses_preferred_model_name(CoreLossesModels value);

    };
} // namespace OpenMagnetics