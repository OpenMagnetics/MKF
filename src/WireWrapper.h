#pragma once

#include "Utils.h"
#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using json = nlohmann::json;

namespace OpenMagnetics {

class WireWrapper : public WireS {
    public:
        WireWrapper(WireS wire) {
            if (wire.get_coating())
                set_coating(wire.get_coating().value());
            if (wire.get_conducting_diameter())
                set_conducting_diameter(wire.get_conducting_diameter().value());
            if (wire.get_conducting_height())
                set_conducting_height(wire.get_conducting_height().value());
            if (wire.get_conducting_width())
                set_conducting_width(wire.get_conducting_width().value());
            if (wire.get_manufacturer_info())
                set_manufacturer_info(wire.get_manufacturer_info().value());
            if (wire.get_material())
                set_material(wire.get_material().value());
            if (wire.get_name())
                set_name(wire.get_name().value());
            if (wire.get_number_conductors())
                set_number_conductors(wire.get_number_conductors().value());
            if (wire.get_outer_diameter())
                set_outer_diameter(wire.get_outer_diameter().value());
            if (wire.get_outer_height())
                set_outer_height(wire.get_outer_height().value());
            if (wire.get_outer_width())
                set_outer_width(wire.get_outer_width().value());
            if (wire.get_standard())
                set_standard(wire.get_standard().value());
            if (wire.get_standard_name())
                set_standard_name(wire.get_standard_name().value());
            if (wire.get_type())
                set_type(wire.get_type().value());
            if (wire.get_strand())
                set_strand(wire.get_strand().value());
        }

        WireWrapper(WireSolid wire) {
            if (wire.get_coating())
                set_coating(wire.get_coating().value());
            if (wire.get_conducting_diameter())
                set_conducting_diameter(wire.get_conducting_diameter().value());
            if (wire.get_conducting_height())
                set_conducting_height(wire.get_conducting_height().value());
            if (wire.get_conducting_width())
                set_conducting_width(wire.get_conducting_width().value());
            if (wire.get_manufacturer_info())
                set_manufacturer_info(wire.get_manufacturer_info().value());
            if (wire.get_material())
                set_material(wire.get_material().value());
            if (wire.get_name())
                set_name(wire.get_name().value());
            if (wire.get_number_conductors())
                set_number_conductors(wire.get_number_conductors().value());
            if (wire.get_outer_diameter())
                set_outer_diameter(wire.get_outer_diameter().value());
            if (wire.get_outer_height())
                set_outer_height(wire.get_outer_height().value());
            if (wire.get_outer_width())
                set_outer_width(wire.get_outer_width().value());
            if (wire.get_standard())
                set_standard(wire.get_standard().value());
            if (wire.get_standard_name())
                set_standard_name(wire.get_standard_name().value());
            if (wire.get_type())
                set_type(wire.get_type().value());
        }
        WireWrapper() = default;
        virtual ~WireWrapper() = default;

        static std::optional<InsulationWireCoating> get_coating(WireS wire);
        static WireWrapper get_strand(WireS wire);
        static double get_filling_factor(double conductingDiameter, int grade = 1, WireStandard standard = WireStandard::IEC_60317, bool includeAirInCell = false);
        void set_nominal_value_conducting_diameter(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_diameter(aux);
        }

        void set_nominal_value_conducting_height(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_height(aux);
        }

        void set_nominal_value_conducting_width(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_conducting_width(aux);
        }

        void set_nominal_value_outer_diameter(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_diameter(aux);
        }

        void set_nominal_value_outer_height(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_height(aux);
        }

        void set_nominal_value_outer_width(double value) {
            DimensionWithTolerance aux;
            aux.set_nominal(value);
            set_outer_width(aux);
        }

        int get_equivalent_insulation_layers(double voltageToInsulate);
};
} // namespace OpenMagnetics
