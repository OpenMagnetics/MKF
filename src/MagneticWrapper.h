#pragma once

#include <MAS.hpp>
#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include "WireWrapper.h"
#include "BobbinWrapper.h"

namespace OpenMagnetics {

class MagneticWrapper : public Magnetic {
    private:
        CoreWrapper core;
        std::optional<std::vector<DistributorInfo>> distributors_info;
        std::optional<ManufacturerInfo> manufacturer_info;
        CoilWrapper coil;
    public:
        MagneticWrapper() = default;
        virtual ~MagneticWrapper() = default;

        /**
         * Data describing the coil
         */
        const CoilWrapper & get_coil() const { return coil; }
        CoilWrapper & get_mutable_coil() { return coil; }
        void set_coil(const CoilWrapper & value) { this->coil = value; }

        /**
         * Data describing the magnetic core.
         */
        const CoreWrapper & get_core() const { return core; }
        CoreWrapper & get_mutable_core() { return core; }
        void set_core(const CoreWrapper & value) { this->core = value; }

        MagneticWrapper(const Magnetic magnetic) {
            set_core(magnetic.get_core());
            set_coil(magnetic.get_coil());

            if (magnetic.get_distributors_info()) {
                set_distributors_info(magnetic.get_distributors_info());
            }
            if (magnetic.get_manufacturer_info()) {
                set_manufacturer_info(magnetic.get_manufacturer_info());
            }
        }


        BobbinWrapper get_bobbin();
        std::vector<WireWrapper> get_wires();
        std::vector<double> get_turns_ratios();
        WireWrapper get_wire(size_t windingIndex=0);
        std::string get_reference();
        std::vector<double> get_maximum_dimensions();
        bool fits(MaximumDimensions maximumDimensions, bool allowRotation);

        double calculate_saturation_current(double temperature = Defaults().ambientTemperature);

};

void from_json(const json & j, MagneticWrapper & x);
void to_json(json & j, const MagneticWrapper & x);

inline void from_json(const json & j, MagneticWrapper& x) {
    x.set_coil(j.at("coil").get<CoilWrapper>());
    x.set_core(j.at("core").get<CoreWrapper>());
    x.set_distributors_info(get_stack_optional<std::vector<DistributorInfo>>(j, "distributorsInfo"));
    x.set_manufacturer_info(get_stack_optional<MagneticManufacturerInfo>(j, "manufacturerInfo"));
    x.set_rotation(get_stack_optional<std::vector<double>>(j, "rotation"));
}

inline void to_json(json & j, const MagneticWrapper & x) {
    j = json::object();
    j["coil"] = x.get_coil();
    j["core"] = x.get_core();
    j["distributorsInfo"] = x.get_distributors_info();
    j["manufacturerInfo"] = x.get_manufacturer_info();
    j["rotation"] = x.get_rotation();
}

} // namespace OpenMagnetics
