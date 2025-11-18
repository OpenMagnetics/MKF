#pragma once

#include <MAS.hpp>
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Bobbin.h"

using namespace MAS;

namespace OpenMagnetics {

class Magnetic : public MAS::Magnetic {
    private:
        Core core;
        std::optional<std::vector<DistributorInfo>> distributors_info;
        std::optional<ManufacturerInfo> manufacturer_info;
        std::optional<std::vector<double>> _maximumDimensions;
        Coil coil;
    public:
        Magnetic() = default;
        virtual ~Magnetic() = default;

        /**
         * Data describing the coil
         */
        const Coil & get_coil() const { return coil; }
        Coil & get_mutable_coil() { return coil; }
        void set_coil(const Coil & value) { this->coil = value; }

        /**
         * Data describing the magnetic core.
         */
        const Core & get_core() const { return core; }
        Core & get_mutable_core() { return core; }
        void set_core(const Core & value) { this->core = value; }

        Magnetic(const MAS::Magnetic magnetic) {
            set_core(magnetic.get_core());
            set_coil(magnetic.get_coil());

            if (magnetic.get_distributors_info()) {
                set_distributors_info(magnetic.get_distributors_info());
            }
            if (magnetic.get_manufacturer_info()) {
                set_manufacturer_info(magnetic.get_manufacturer_info());
            }
        }


        Bobbin get_bobbin();
        std::vector<Wire> get_wires();
        std::vector<double> get_turns_ratios();
        static std::vector<double> get_turns_ratios(MAS::Magnetic magnetic);
        Wire get_wire(size_t windingIndex=0);
        std::string get_reference();
        std::vector<double> get_maximum_dimensions();
        bool fits(MaximumDimensions maximumDimensions, bool allowRotation);

        double calculate_saturation_current(double temperature = Defaults().ambientTemperature);

};

bool operator==(Magnetic lhs, Magnetic rhs);

inline bool operator==(Magnetic lhs, Magnetic rhs) {
    bool isEqual = lhs.get_reference() == rhs.get_reference() && 
                   lhs.get_mutable_core().get_shape_name() == rhs.get_mutable_core().get_shape_name() && 
                   lhs.get_mutable_core().get_material_name() == rhs.get_mutable_core().get_material_name() && 
                   lhs.get_mutable_core().get_number_stacks() == rhs.get_mutable_core().get_number_stacks() && 
                   lhs.get_coil().get_functional_description().size() == rhs.get_coil().get_functional_description().size();
    if (isEqual) {
        for (size_t i = 0; i < lhs.get_coil().get_functional_description().size(); ++i) {
            auto lhsWinding = lhs.get_coil().get_functional_description()[i];
            auto rhsWinding = rhs.get_coil().get_functional_description()[i];
            auto lhsWire = lhsWinding.resolve_wire();
            auto rhsWire = rhsWinding.resolve_wire();
            isEqual &= lhsWinding.get_number_turns() == rhsWinding.get_number_turns() &&
                       lhsWinding.get_number_parallels() == rhsWinding.get_number_parallels() &&
                       lhsWire.get_type() == rhsWire.get_type();
        }
    }

    return isEqual;
}

void from_json(const json & j, Magnetic & x);
void to_json(json & j, const Magnetic & x);
void from_json(const json& j, std::vector<Magnetic>& v);
void to_json(json& j, const std::vector<Magnetic>& v);
void from_file(std::filesystem::path filepath, Magnetic & x);
void to_file(std::filesystem::path filepath, const Magnetic & x);

inline void from_file(std::filesystem::path filepath, Magnetic & x) {
    std::ifstream f(filepath);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    auto masJson = json::parse(data);
    auto magneticJson = masJson["magnetic"];

    x = OpenMagnetics::Magnetic(magneticJson);
}

inline void from_json(const json & j, Magnetic& x) {
    x.set_coil(j.at("coil").get<Coil>());
    x.set_core(j.at("core").get<Core>());
    x.set_distributors_info(get_stack_optional<std::vector<DistributorInfo>>(j, "distributorsInfo"));
    x.set_manufacturer_info(get_stack_optional<MagneticManufacturerInfo>(j, "manufacturerInfo"));
    x.set_rotation(get_stack_optional<std::vector<double>>(j, "rotation"));
}

inline void to_json(json & j, const Magnetic & x) {
    j = json::object();
    j["coil"] = x.get_coil();
    j["core"] = x.get_core();
    j["distributorsInfo"] = x.get_distributors_info();
    j["manufacturerInfo"] = x.get_manufacturer_info();
    j["rotation"] = x.get_rotation();
}

inline void from_json(const json& j, std::vector<Magnetic>& v) {
    for (auto e : j) {
        Magnetic x;
        x.set_coil(e.at("coil").get<Coil>());
        x.set_core(e.at("core").get<Core>());
        x.set_distributors_info(get_stack_optional<std::vector<DistributorInfo>>(e, "distributorsInfo"));
        x.set_manufacturer_info(get_stack_optional<MagneticManufacturerInfo>(e, "manufacturerInfo"));
        x.set_rotation(get_stack_optional<std::vector<double>>(e, "rotation"));
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<Magnetic>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        e["coil"] = x.get_coil();
        e["core"] = x.get_core();
        e["distributorsInfo"] = x.get_distributors_info();
        e["manufacturerInfo"] = x.get_manufacturer_info();
        e["rotation"] = x.get_rotation();
        j.push_back(e);
    }
}
} // namespace OpenMagnetics
