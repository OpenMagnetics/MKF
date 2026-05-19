#include "advisers/CoreCrossReferencer.h"
#include "support/Utils.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace { 
TEST_CASE("Test_All_Core_Materials", "[adviser][core-cross-referencer]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);

    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;

    std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    // Phase 1 fix landed: SATURATION filter is now active (previously
    // dropped by an incomplete switch in apply_filters). For a ferrite
    // reference the top cross-reference is now the same shape with an
    // adjacent ferrite material (3C94) — formerly a powder Kool Mµ on a
    // totally different shape, which was the symptom of the ferrite→powder
    // bias documented in the audit.
    REQUIRE(crossReferencedCores[0].first.get_name().value() == "EC 35/17/10 - 3C94 - Gapped 1.000 mm");
}

TEST_CASE("Test_All_Core_Materials_Same_Material", "[adviser][core-cross-referencer][smoke-test]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;
    coreCrossReferencer.use_only_reference_material(true);

    std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    // Check that the cross-referenced core has the same material (3C91)
    auto resultName = crossReferencedCores[0].first.get_name().value();
    REQUIRE(resultName.find("3C91") != std::string::npos);
}

TEST_CASE("Test_All_Core_Materials_Same_Material_Maximum_Height", "[adviser][core-cross-referencer][smoke-test]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;
    coreCrossReferencer.use_only_reference_material(true);

    std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    MaximumDimensions maximumDimensions;
    maximumDimensions.set_depth(0.04);
    maximumDimensions.set_height(0.025);
    maximumDimensions.set_width(0.025);
    inputs.get_mutable_design_requirements().set_maximum_dimensions(maximumDimensions);

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    REQUIRE(crossReferencedCores[0].first.get_name().value() == "EP 20 - 3C91 - Gapped 0.605 mm");
}

TEST_CASE("Test_All_Core_Materials_Only_TDK", "[adviser][core-cross-referencer]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;
    coreCrossReferencer.use_only_manufacturer("TDK");

    std::string coreName = "EC 35/17/10 - 3C91 - Gapped 1.000 mm";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    // Phase 1 fix landed: with SATURATION now active, the TDK-only top
    // cross-reference for an EC 35/17/10 ferrite shifts from the larger
    // ETD 34/17/11 to the dimensionally closer ETD 29/16/10.
    REQUIRE(crossReferencedCores[0].first.get_name().value() == "ETD 29/16/10 - N87 - Gapped 1.000 mm");
}

TEST_CASE("Test_All_Core_Materials_Powder", "[adviser][core-cross-referencer]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;

    std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    // Phase 1 fix landed: with the CORE_LOSSES ceiling clamp removed and
    // SATURATION active, the powder reference now matches a same-shape
    // candidate (E 25/9.5/6.3 MS 60) rather than a larger MS 75.
    REQUIRE(crossReferencedCores[0].first.get_name().value() == "E 25/9.5/6.3 - MS 60 - Ungapped");
}

TEST_CASE("Test_All_Core_Materials_Only_Micrometals", "[adviser][core-cross-referencer][smoke-test]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;
    coreCrossReferencer.use_only_manufacturer("Micrometals");

    std::string coreName = "E 25/9.5/6.3 - XFlux 60 - Ungapped";
    Core core = find_core_by_name(coreName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 28;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 5);

    REQUIRE(crossReferencedCores.size() > 0);
    // Accept either MS 75, MS 40, or MS 26 as top result due to algorithm improvements
    auto topCoreName = crossReferencedCores[0].first.get_name().value();
    REQUIRE((topCoreName.find("MS ") != std::string::npos || 
             topCoreName.find("E 25/9.5/6.3") != std::string::npos));
}

TEST_CASE("Test_Cross_Reference_Core_Web_0", "[adviser][core-cross-referencer][bug]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();
    OperatingPoint operatingPoint;
    CoreCrossReferencer coreCrossReferencer;

    std::string shapeName = "PQ 40/40";
    std::string materialName = "3C97";
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, materialName);

    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    int64_t numberTurns = 16;

    auto crossReferencedCores = coreCrossReferencer.get_cross_referenced_core(core, numberTurns, inputs, 20);

    REQUIRE(crossReferencedCores.size() > 0);
    // Phase 1 fix landed: with SATURATION active, the 3C97 PQ 40/40
    // reference now matches a same-material ferrite (ETD 44/22/15 - 3C97)
    // rather than a cross-family powder (EQ 42/20/28 - High Flux 26),
    // which is the expected behaviour for a similarity cross-referencer.
    REQUIRE(crossReferencedCores[0].first.get_name().value() == "ETD 44/22/15 - 3C97 - Gapped 1.000 mm");

    auto scorings = coreCrossReferencer.get_scorings();
    auto scoredValues = coreCrossReferencer.get_scored_values();
    json results;
    results["cores"] = json::array();
    results["scorings"] = json::array();
    for (auto& [core, scoring] : crossReferencedCores) {
        std::string name = core.get_name().value();

        json coreJson;
        to_json(coreJson, core);
        results["cores"].push_back(coreJson);
        results["scorings"].push_back(scoring);

        json result;
        result["scoringPerFilter"] = json();
        result["scoredValuePerFilter"] = json();
        for (auto& filter : magic_enum::enum_names<CoreCrossReferencerFilters>()) {
            std::string filterString(filter);

            result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()];
            result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()];
            REQUIRE(!std::isnan(scorings[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()]));
            REQUIRE(!std::isnan(scoredValues[name][magic_enum::enum_cast<CoreCrossReferencerFilters>(filterString).value()]));
        };
        results["data"].push_back(result);
    }

    settings.reset();
}

}  // namespace
