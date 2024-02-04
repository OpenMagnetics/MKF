#pragma once
#include "MAS.hpp"
#include "Defaults.h"
#include "Painter.h"

namespace OpenMagnetics {

class Settings
{
    protected:
        Settings() {}

        static Settings* settings_;

        bool _coilAllowMarginTape = true;
        bool _coilAllowInsulatedWire = true;
        bool _coilFillSectionsWithMarginTape = false;
        bool _coilWindEvenIfNotFit = false;
        bool _coilDelimitAndCompact = true;
        bool _coilTryRewind = true;

        bool _useOnlyCoresInStock = true;

        size_t _painterNumberPointsX = 25;
        size_t _painterNumberPointsY = 50;
        int _painterMirroringDimension = Defaults().magneticFieldMirroringDimension;
        Painter::PainterModes _painterMode = Painter::PainterModes::CONTOUR;
        bool _painterLogarithmicScale = false;
        bool _painterIncludeFringing = true;
        std::optional<double> _painterMaximumValueColorbar = std::nullopt;
        std::optional<double> _painterMinimumValueColorbar = std::nullopt;
        std::string _painterColorFerrite = "0x8F7b7c7d";
        std::string _painterColorBobbin = "0x8F1b1b1b";
        std::string _painterColorCopper = "0x8Fb87333";
        std::string _painterColorInsulation = "0x18E37E00";
        std::string _painterColorMargin = "0x8Ffff05b";

        size_t _magneticFieldNumberPointsX = 25;
        size_t _magneticFieldNumberPointsY = 50;
        int _magneticFieldMirroringDimension = Defaults().magneticFieldMirroringDimension;
        bool _magneticFieldIncludeFringing = true;
        size_t _coilAdviserMaximumNumberWires = 100;
        bool _verbose = false;

    public:
        Settings(Settings &other) = delete;
        void operator=(const Settings &) = delete;

        static Settings *GetInstance();

        void reset();

        bool get_verbose() const;
        void set_verbose(bool value);

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

    };
} // namespace OpenMagnetics