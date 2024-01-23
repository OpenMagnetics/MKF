#pragma once
#include "MAS.hpp"
#include "Defaults.h"

namespace OpenMagnetics {

class Settings
{
    protected:
        Settings() {}

        static Settings* settings_;

        bool _allowMarginTape = true;
        bool _allowInsulatedWire = true;
        bool _fillCoilSectionsWithMarginTape = false;
        bool _windEvenIfNotFit = false;
        bool _useOnlyCoresInStock = true;

    public:
        Settings(Settings &other) = delete;
        void operator=(const Settings &) = delete;

        static Settings *GetInstance();

        void reset();

        bool get_allow_margin_tape() const;
        void set_allow_margin_tape(bool value);

        bool get_allow_insulated_wire() const;
        void set_allow_insulated_wire(bool value);

        bool get_fill_coil_sections_with_margin_tape() const;
        void set_fill_coil_sections_with_margin_tape(bool value);

        bool get_wind_even_if_not_fit() const;
        void set_wind_even_if_not_fit(bool value);

        bool get_use_only_cores_in_stock() const;
        void set_use_only_cores_in_stock(bool value);
    };
} // namespace OpenMagnetics