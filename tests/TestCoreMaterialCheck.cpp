// check_core_material / compare_core_material_volumetric_losses: a core material judged on its own,
// without a design (ABT #1911). The three refused materials are the embedded copies stored web designs
// carry: DMEGC DMR96 (Steinmetz ct(T) not positive over its data), Ferroxcube 3C94 (ct(100 C) = -218)
// and Fair-Rite 95 (remanence 0.8 T at 100 C, above its 0.39 T saturation flux density).
#include "physical_models/CoreMaterialCheck.h"
#include "physical_models/CoreLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>

using namespace OpenMagnetics;

namespace {

CoreMaterial load_material_fixture(const std::string& fileName) {
    auto path = std::filesystem::path(__FILE__).parent_path() / "testData" / fileName;
    std::ifstream file(path);
    REQUIRE(file.good());
    return CoreMaterial(json::parse(file));
}

std::vector<CoreMaterialProblem> problems_with(const CoreMaterialCheck& check, ErrorCode code, const std::string& property) {
    std::vector<CoreMaterialProblem> found;
    for (const auto& problem : check.problems) {
        if (problem.code == code && problem.property == property) {
            found.push_back(problem);
        }
    }
    return found;
}

bool covers(const std::optional<CoreMaterialSpan>& span, double value) {
    return span && span->minimum <= value && value <= span->maximum;
}

} // namespace

TEST_CASE("Test_Core_Material_Check_Library_Materials_Pass", "[material-check]") {
    for (const std::string name : {"3C94", "DMR96", "95", "N87"}) {
        auto check = check_core_material(find_core_material_by_name(name));
        INFO(name << ": " << json(check).dump(2));
        CHECK(check.valid);
        CHECK(check.problems.empty());
        CHECK(check.coverage.temperature.minimum < check.coverage.temperature.maximum);
        REQUIRE(check.coverage.lossFrequency);
        CHECK(check.coverage.lossFrequency->minimum < check.coverage.lossFrequency->maximum);
    }
}

TEST_CASE("Test_Core_Material_Check_DMR96_Old_Copy_Invalid_Loss_Data", "[material-check]") {
    // MAS 051605e's DMR96: ct2 = -0.446 in the 1 Hz - 500 kHz range, ct(T) > 0 only between -1.07 C and
    // 2.72 C. Its data covers -55 C (initial permeability) to 230 C, its Curie point; it is refused from
    // 2.72 C up to the Curie point, which holds every temperature its saturation and remanence are
    // tabulated at (25 C to 140 C).
    auto material = load_material_fixture("Core_Material_Check_DMR96_Old_Copy.json");
    auto check = check_core_material(material);
    INFO(json(check).dump(2));
    CHECK(!check.valid);
    CHECK(check.coverage.temperature.minimum == -55);
    CHECK(check.coverage.temperature.maximum == 230);
    auto problems = problems_with(check, ErrorCode::INVALID_LOSS_DATA, "volumetricLosses");
    REQUIRE(!problems.empty());
    bool lowRangeRefusedEverywhere = false;
    for (const auto& problem : problems) {
        CHECK_THAT(problem.engineMessage, Catch::Matchers::StartsWith("[INVALID_LOSS_DATA]"));
        CHECK_THAT(problem.message, Catch::Matchers::ContainsSubstring("DMR96"));
        if (problem.frequency && problem.frequency->maximum == 500000 && problem.temperature->minimum < 3 &&
            problem.temperature->maximum >= 140) {
            lowRangeRefusedEverywhere = true;
        }
    }
    CHECK(lowRangeRefusedEverywhere);
}

TEST_CASE("Test_Core_Material_Check_3C94_Old_Copy_Invalid_Loss_Data", "[material-check]") {
    // The 3C94 copy a stored design embeds: ct(T) = -0.0229 T^2 + 0.105 T + 0.0159 in its 1 Hz - 3 MHz
    // range, -218 at 100 C, positive only below 4.74 C.
    auto material = load_material_fixture("Core_Material_Check_3C94_Old_Copy.json");
    REQUIRE(material.get_name() == "3C94");
    auto check = check_core_material(material);
    INFO(json(check).dump(2));
    CHECK(!check.valid);
    auto problems = problems_with(check, ErrorCode::INVALID_LOSS_DATA, "volumetricLosses");
    REQUIRE(!problems.empty());
    bool refusedAt100 = false;
    for (const auto& problem : problems) {
        if (covers(problem.temperature, 100)) {
            refusedAt100 = true;
            REQUIRE(problem.frequency);
            CHECK(problem.frequency->maximum == 3000000);
            // The interval starts at the root of ct(T) (4.74 C) or at the start of the data, whichever is later.
            CHECK(problem.temperature->minimum >= 4.7);
            CHECK(problem.temperature->maximum >= 100);
            CHECK_THAT(problem.message, Catch::Matchers::ContainsSubstring("C to"));
        }
    }
    CHECK(refusedAt100);
}

TEST_CASE("Test_Core_Material_Check_95_Old_Copy_Remanence_Above_Saturation", "[material-check]") {
    // Fair-Rite 95 with remanence 0.8 T at 100 C (the library: 0.08 T) and saturation 0.39 T there.
    auto material = load_material_fixture("Core_Material_Check_95_Old_Copy.json");
    auto check = check_core_material(material);
    INFO(json(check).dump(2));
    CHECK(!check.valid);
    auto problems = problems_with(check, ErrorCode::MATERIAL_INVALID_PROPERTY, "remanence");
    REQUIRE(problems.size() == 1);
    const auto& problem = problems.front();
    CHECK(covers(problem.temperature, 100));
    // At 25 C (0.102 T against 0.5 T) the remanence is accepted: the refusal starts between the two tabulated temperatures.
    CHECK(problem.temperature->minimum > 25);
    CHECK(problem.temperature->minimum < 100);
    CHECK_THAT(problem.engineMessage, Catch::Matchers::StartsWith("[MATERIAL_INVALID_PROPERTY]"));
    // Its loss data is the library's: not refused.
    CHECK(problems_with(check, ErrorCode::INVALID_LOSS_DATA, "volumetricLosses").empty());
}

TEST_CASE("Test_Core_Material_Check_No_Temperature_Data_Throws", "[material-check]") {
    json materialJson;
    to_json(materialJson, find_core_material_by_name("3C94"));
    materialJson["saturation"] = json::array();
    materialJson.erase("remanence");
    materialJson.erase("coerciveForce");
    materialJson["permeability"]["initial"] = json{{"value", 2000}};
    CoreMaterial material(materialJson);
    REQUIRE(CoreLossesModel::get_volumetric_losses_data(material).empty());
    try {
        check_core_material(material);
        FAIL("check_core_material accepted a material with no temperature data");
    }
    catch (const MaterialException& exception) {
        CHECK(exception.code() == ErrorCode::MATERIAL_DATA_MISSING);
    }
}

TEST_CASE("Test_Core_Material_Compare_Library_Copy_With_Itself", "[material-check]") {
    auto material = find_core_material_by_name("3C94");
    auto comparison = compare_core_material_volumetric_losses(material, material);
    INFO(json(comparison).dump(2));
    REQUIRE(!comparison.points.empty());
    CHECK(comparison.pointsRefusedByReference == 0);
    CHECK(comparison.pointsRefusedByOther == 0);
    REQUIRE(comparison.maximumAbsoluteRelativeDifference);
    CHECK(comparison.maximumAbsoluteRelativeDifference.value() == 0);
    // Every point is inside the reference's own data: its saturation temperatures and Steinmetz ranges.
    auto [minimumFrequency, maximumFrequency] = CoreLossesModel::get_steinmetz_fitted_span(material);
    for (const auto& point : comparison.points) {
        CHECK(point.frequency >= minimumFrequency);
        CHECK(point.frequency <= maximumFrequency);
        CHECK((point.temperature == 25 || point.temperature == 100));
        CHECK(point.magneticFluxDensityPeak < Core::get_magnetic_flux_density_saturation(material, point.temperature, false));
    }
}

TEST_CASE("Test_Core_Material_Compare_Library_With_Old_Copies", "[material-check]") {
    SECTION("DMR96: the old copy is refused where its ct is not positive") {
        auto comparison = compare_core_material_volumetric_losses(find_core_material_by_name("DMR96"),
                                                                  load_material_fixture("Core_Material_Check_DMR96_Old_Copy.json"));
        INFO(json(comparison).dump(2));
        CHECK(comparison.pointsRefusedByReference == 0);
        CHECK(comparison.pointsRefusedByOther > 0);
        for (const auto& point : comparison.points) {
            if (!point.other.volumetricLosses) {
                CHECK(point.other.code == ErrorCode::INVALID_LOSS_DATA);
            }
        }
    }
    SECTION("3C94: the old copy is refused at every point") {
        auto comparison = compare_core_material_volumetric_losses(find_core_material_by_name("3C94"),
                                                                  load_material_fixture("Core_Material_Check_3C94_Old_Copy.json"));
        INFO(json(comparison).dump(2));
        CHECK(comparison.pointsRefusedByReference == 0);
        CHECK(comparison.pointsRefusedByOther == comparison.points.size());
        CHECK(!comparison.maximumAbsoluteRelativeDifference);
    }
}

TEST_CASE("Test_Core_Material_Compare_Reference_Without_Loss_Span_Throws", "[material-check]") {
    json materialJson;
    to_json(materialJson, find_core_material_by_name("3C94"));
    json methods = json::array();
    for (const auto& method : materialJson["volumetricLosses"]["default"]) {
        if (!(method.is_object() && method.contains("method") && method["method"] == "steinmetz")) {
            methods.push_back(method);
        }
    }
    materialJson["volumetricLosses"]["default"] = methods;
    CoreMaterial reference(materialJson);
    try {
        compare_core_material_volumetric_losses(reference, find_core_material_by_name("3C94"));
        FAIL("compare_core_material_volumetric_losses compared at points the reference's data does not state");
    }
    catch (const MaterialException& exception) {
        CHECK(exception.code() == ErrorCode::MATERIAL_DATA_MISSING);
    }
}
