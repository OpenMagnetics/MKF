// The core thermal properties MKF hands the thermal solver and the FEM exporters: a material's own
// MAS heatConductivity / heatCapacity, else the approved fallback (Alf, 2026-10-11) -- the median
// over the MAS materials of the same materialComposition that carry the value -- else a throw.
#include "physical_models/ThermalResistance.h"
#include "support/Utils.h"
#include "support/Exceptions.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <map>
#include <string>

using namespace OpenMagnetics;

namespace {

DimensionWithTolerance nominal(double value) {
    DimensionWithTolerance d;
    d.set_nominal(value);
    return d;
}

// A real record to clone, so every required field is present; read before any test swaps the
// database out.
const CoreMaterial& template_material() {
    static const CoreMaterial m = find_core_material_by_name("3C95");
    return m;
}

// A material of the given composition and thermal values, cloned from the template.
CoreMaterial material(const std::string& name, std::optional<MaterialComposition> composition,
                      std::optional<double> k, std::optional<double> cp) {
    CoreMaterial m = template_material();
    m.set_name(name);
    m.set_material_composition(composition);
    m.set_heat_conductivity(k ? std::optional<DimensionWithTolerance>(nominal(*k)) : std::nullopt);
    m.set_heat_capacity(cp ? std::optional<DimensionWithTolerance>(nominal(*cp)) : std::nullopt);
    return m;
}

// Swaps the global material database for a controlled one and restores it on scope exit.
struct ScopedMaterialDatabase {
    std::map<std::string, MAS::CoreMaterial> saved;
    explicit ScopedMaterialDatabase(const std::vector<CoreMaterial>& materials) {
        template_material();
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        saved = coreMaterialDatabase;
        coreMaterialDatabase.clear();
        for (const auto& m : materials) {
            coreMaterialDatabase[m.get_name()] = m;
        }
    }
    ~ScopedMaterialDatabase() { coreMaterialDatabase = saved; }
};

}  // namespace

TEST_CASE("Core thermal properties: a material's own MAS value wins", "[core-thermal-properties]") {
    // Ferroxcube 3C95: heatConductivity {3.5, 5}, heatCapacity {700, 800} in MAS
    const auto m = find_core_material_by_name("3C95");
    REQUIRE(m.get_heat_conductivity());
    REQUIRE(m.get_heat_capacity());
    CHECK_THAT(ThermalResistance::getCoreMaterialThermalConductivity(m), Catch::Matchers::WithinAbs(4.25, 1e-12));
    CHECK_THAT(ThermalResistance::getCoreMaterialHeatCapacity(m), Catch::Matchers::WithinAbs(750.0, 1e-12));
}

TEST_CASE("Core thermal properties: a material without a value takes the same-composition MEDIAN", "[core-thermal-properties]") {
    // Odd count: median 2 (the mean would be 4.33); the NiZn 100 must not leak into FeSi.
    const ScopedMaterialDatabase db({material("a", MaterialComposition::FE_SI, 1.0, 400.0),
                                     material("b", MaterialComposition::FE_SI, 2.0, 500.0),
                                     material("c", MaterialComposition::FE_SI, 10.0, 900.0),
                                     material("d", MaterialComposition::NI_ZN, 100.0, 9000.0),
                                     material("e", MaterialComposition::FE_SI, std::nullopt, std::nullopt)});
    const auto bare = material("bare", MaterialComposition::FE_SI, std::nullopt, std::nullopt);
    CHECK_THAT(ThermalResistance::getCoreMaterialThermalConductivity(bare), Catch::Matchers::WithinAbs(2.0, 1e-12));
    CHECK_THAT(ThermalResistance::getCoreMaterialHeatCapacity(bare), Catch::Matchers::WithinAbs(500.0, 1e-12));
    // The by-name lookup goes through the same resolver.
    CHECK_THAT(ThermalResistance::getMaterialThermalConductivity("e"), Catch::Matchers::WithinAbs(2.0, 1e-12));
}

TEST_CASE("Core thermal properties: an even count takes the mean of the two middle values", "[core-thermal-properties]") {
    const ScopedMaterialDatabase db({material("a", MaterialComposition::MN_ZN, 3.0, 700.0),
                                     material("b", MaterialComposition::MN_ZN, 4.0, 750.0),
                                     material("c", MaterialComposition::MN_ZN, 5.0, 800.0),
                                     material("d", MaterialComposition::MN_ZN, 9.0, 1100.0)});
    const auto bare = material("bare", MaterialComposition::MN_ZN, std::nullopt, std::nullopt);
    CHECK_THAT(ThermalResistance::getCoreMaterialThermalConductivity(bare), Catch::Matchers::WithinAbs(4.5, 1e-12));
    CHECK_THAT(ThermalResistance::getCoreMaterialHeatCapacity(bare), Catch::Matchers::WithinAbs(775.0, 1e-12));
}

TEST_CASE("Core thermal properties: no value and nothing to pool by throws, never a guess", "[core-thermal-properties]") {
    const ScopedMaterialDatabase db({material("a", MaterialComposition::FE_SI, 1.0, std::nullopt),
                                     material("p", MaterialComposition::PROPRIETARY, 5.0, 500.0)});
    using Catch::Matchers::ContainsSubstring;
    // composition present, but no record of it carries heatCapacity
    CHECK_THROWS_WITH(ThermalResistance::getCoreMaterialHeatCapacity(
                          material("bare", MaterialComposition::FE_SI, std::nullopt, std::nullopt)),
                      ContainsSubstring("no FE_SI material in the database carries one"));
    // proprietary: no chemistry to pool by, even though a proprietary record has values
    CHECK_THROWS_WITH(ThermalResistance::getCoreMaterialThermalConductivity(
                          material("bare", MaterialComposition::PROPRIETARY, std::nullopt, std::nullopt)),
                      ContainsSubstring("no materialComposition"));
    // no composition at all
    CHECK_THROWS_WITH(ThermalResistance::getCoreMaterialThermalConductivity(
                          material("bare", std::nullopt, std::nullopt, std::nullopt)),
                      ContainsSubstring("no materialComposition"));
}
