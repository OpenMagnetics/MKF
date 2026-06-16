#pragma once

#include "constructive_models/MasMigration.h"
#include <MAS.hpp>
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Bobbin.h"

using namespace MAS;

namespace OpenMagnetics {

class Magnetic : public MAS::Magnetic {
    private:
        // core/coil are optional because, as of the MAS chip_beads_datasheet_schema
        // change, a Magnetic may describe a component (e.g. a chip bead specified
        // purely by its datasheet impedance/model) that has no physical core or
        // coil construction. The accessors below throw when the value is absent
        // rather than returning a fabricated default — callers that need a core or
        // coil must gate on has_core()/has_coil() first.
        std::optional<Core> core;
        std::optional<std::vector<DistributorInfo>> distributors_info;
        std::optional<ManufacturerInfo> manufacturer_info;
        std::optional<std::vector<double>> _maximumDimensions;
        std::optional<Coil> coil;
    public:
        Magnetic() = default;
        virtual ~Magnetic() = default;

        bool has_coil() const { return coil.has_value(); }
        bool has_core() const { return core.has_value(); }

        /**
         * Data describing the coil
         */
        const Coil & get_coil() const {
            if (!coil) throw std::runtime_error("Magnetic has no coil (e.g. a chip-bead datasheet entry with no winding construction)");
            return coil.value();
        }
        Coil & get_mutable_coil() {
            if (!coil) throw std::runtime_error("Magnetic has no coil (e.g. a chip-bead datasheet entry with no winding construction)");
            return coil.value();
        }
        void set_coil(const Coil & value) { this->coil = value; }

        /**
         * Data describing the magnetic core.
         */
        const Core & get_core() const {
            if (!core) throw std::runtime_error("Magnetic has no core (e.g. a chip-bead datasheet entry with no core construction)");
            return core.value();
        }
        Core & get_mutable_core() {
            if (!core) throw std::runtime_error("Magnetic has no core (e.g. a chip-bead datasheet entry with no core construction)");
            return core.value();
        }
        void set_core(const Core & value) { this->core = value; }

        Magnetic(const MAS::Magnetic magnetic) {
            if (magnetic.get_core()) {
                set_core(magnetic.get_core().value());
            }
            if (magnetic.get_coil()) {
                set_coil(magnetic.get_coil().value());
            }

            if (magnetic.get_distributors_info()) {
                set_distributors_info(magnetic.get_distributors_info());
            }
            if (magnetic.get_manufacturer_info()) {
                set_manufacturer_info(magnetic.get_manufacturer_info());
            }
        }


        Bobbin get_bobbin();
        std::vector<Wire> get_wires();
        std::vector<double> get_turns_ratios() const;
        static std::vector<double> get_turns_ratios(MAS::Magnetic magnetic);
        Wire get_wire(size_t windingIndex=0);
        std::string get_reference() const;
        std::vector<double> get_maximum_dimensions();
        bool fits(MaximumDimensions maximumDimensions, bool allowRotation);

        // I_sat = B_sat(T)·N·A_e/L. By default B_sat is the PROPORTION-derated
        // (0.7·B_sat_raw) usable value the realism gate expects. Pass
        // proportion=false for the RAW saturation current — the CoreAdviser
        // saturation gate uses raw B_sat so its derating is exactly
        // (hot temperature × margin), with no flux-proportion factor stacked on
        // top of the margin (ABT #13).
        double calculate_saturation_current(double temperature = Defaults().ambientTemperature,
                                            bool proportion = true);

        /**
         * @brief Saturation current evaluated at a specific operating point.
         *
         * Same I_sat = B_sat(T)·N·A_e/L identity as the no-arg overload,
         * but L is the *operating* inductance (with DC-bias rolloff of
         * the core's permeability under the operating point's primary
         * current). For heavily-gapped cores the gap dominates and L
         * barely moves with bias; for low-gap / ungapped cores the
         * operating L can be 30–50 % lower than the μ_init L and the
         * I_sat values differ accordingly. Callers comparing I_sat
         * against I_peak (a realism-gate idiom) should use this overload
         * so both numbers are evaluated at the same operating conditions.
         */
        double calculate_saturation_current(OperatingPoint& operatingPoint,
                                            double temperature = Defaults().ambientTemperature);

        /**
         * @brief Rated (heat-limited) current: the steady DC current that raises the
         * component temperature by @p temperatureRise degrees above @p ambientTemperature.
         *
         * This is the datasheet "Irms" / rated-current figure. By industry convention
         * (Coilcraft, Würth, TDK, ...) it is a *DC* current and the heating is the
         * winding's ohmic loss I²·R_dc(T) only — core loss is excluded because a pure DC
         * bias produces no flux swing. AC operation must be derated against this number.
         * The standard rise is 40 K for power inductors (smaller, e.g. 15-20 K, for chip
         * inductors); the default here is 40 K.
         *
         * The current is taken to flow in the primary winding (index 0). The per-turn R_dc
         * ohmic loss comes from WindingOhmicLosses and the temperature rise from the full
         * Temperature thermal-network model, so the value stays consistent with the rest of
         * MKF. The copper-temperature dependence of R_dc is resolved by an inner fixed-point;
         * the rated current itself is found by bisection on the rise.
         */
        double calculate_rated_current(double temperatureRise = 40,
                                       double ambientTemperature = Defaults().ambientTemperature);

};

bool operator==(Magnetic lhs, Magnetic rhs);

inline bool operator==(Magnetic lhs, Magnetic rhs) {
    // A magnetic without a core/coil (e.g. a chip-bead datasheet entry) can only
    // be equal to another that is missing the same parts; guard before touching
    // the accessors, which throw when the value is absent.
    if (lhs.has_core() != rhs.has_core() || lhs.has_coil() != rhs.has_coil()) {
        return false;
    }

    bool isEqual = lhs.get_reference() == rhs.get_reference();

    if (lhs.has_core()) {
        isEqual &= lhs.get_mutable_core().get_shape_name() == rhs.get_mutable_core().get_shape_name() &&
                   lhs.get_mutable_core().get_material_name() == rhs.get_mutable_core().get_material_name() &&
                   lhs.get_mutable_core().get_number_stacks() == rhs.get_mutable_core().get_number_stacks();
    }

    if (lhs.has_coil()) {
        isEqual &= lhs.get_coil().get_functional_description().size() == rhs.get_coil().get_functional_description().size();
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
    OpenMagnetics::compat::migrate_pre_1_0(masJson);
    auto magneticJson = masJson["magnetic"];

    x = OpenMagnetics::Magnetic(magneticJson);
}

inline void from_json(const json & j, Magnetic& x) {
    if (j.contains("coil") && !j.at("coil").is_null()) {
        x.set_coil(j.at("coil").get<Coil>());
    }
    if (j.contains("core") && !j.at("core").is_null()) {
        x.set_core(j.at("core").get<Core>());
    }
    x.set_distributors_info(get_stack_optional<std::vector<DistributorInfo>>(j, "distributorsInfo"));
    x.set_manufacturer_info(get_stack_optional<MagneticManufacturerInfo>(j, "manufacturerInfo"));
    x.set_rotation(get_stack_optional<std::vector<double>>(j, "rotation"));
}

inline void to_json(json & j, const Magnetic & x) {
    j = json::object();
    if (x.has_coil()) {
        j["coil"] = x.get_coil();
    }
    if (x.has_core()) {
        j["core"] = x.get_core();
    }
    j["distributorsInfo"] = x.get_distributors_info();
    j["manufacturerInfo"] = x.get_manufacturer_info();
    j["rotation"] = x.get_rotation();
}

inline void from_json(const json& j, std::vector<Magnetic>& v) {
    for (auto e : j) {
        Magnetic x;
        if (e.contains("coil") && !e.at("coil").is_null()) {
            x.set_coil(e.at("coil").get<Coil>());
        }
        if (e.contains("core") && !e.at("core").is_null()) {
            x.set_core(e.at("core").get<Core>());
        }
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
        if (x.has_coil()) {
            e["coil"] = x.get_coil();
        }
        if (x.has_core()) {
            e["core"] = x.get_core();
        }
        e["distributorsInfo"] = x.get_distributors_info();
        e["manufacturerInfo"] = x.get_manufacturer_info();
        e["rotation"] = x.get_rotation();
        j.push_back(e);
    }
}
} // namespace OpenMagnetics
