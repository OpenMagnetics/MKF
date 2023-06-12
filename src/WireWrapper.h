#pragma once

#include "Utils.h"
#include "json.hpp"

#include <MAS.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {

class WireWrapper : public WireS {
    public:
        WireWrapper(WireS wire) {
            set_coating(wire.get_coating());
            set_conducting_diameter(wire.get_conducting_diameter());
            set_conducting_height(wire.get_conducting_height());
            set_conducting_width(wire.get_conducting_width());
            set_manufacturer_info(wire.get_manufacturer_info());
            set_material(wire.get_material());
            set_name(wire.get_name());
            set_number_conductors(wire.get_number_conductors());
            set_outer_diameter(wire.get_outer_diameter());
            set_outer_height(wire.get_outer_height());
            set_outer_width(wire.get_outer_width());
            set_standard(wire.get_standard());
            set_standard_name(wire.get_standard_name());
            set_type(wire.get_type());
            set_strand(wire.get_strand());
        }
        WireWrapper() = default;
        virtual ~WireWrapper() = default;

        static InsulationWireCoating get_coating(WireS wire);
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
