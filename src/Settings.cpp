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
} // namespace OpenMagnetics