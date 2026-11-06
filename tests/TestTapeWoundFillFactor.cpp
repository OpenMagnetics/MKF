// Fill factor of tape-wound cores, looked up from MAS laminationThickness (ribbon thickness) and the makers' published
// values (ABT #2234). One case per priority path of Core::get_tape_wound_fill_factor, plus one that the flat 0.80
// setting is gone. Sources of every number: the comment block above get_tape_wound_fill_factor in Core.cpp.
#include "constructive_models/Core.h"
#include "support/Logger.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string>

using namespace OpenMagnetics;
using json = nlohmann::json;

namespace {

// A ungapped toroid of the given nominal OD / ID / height (m) in the given material.
Core tape_wound_toroid(const std::string& material, double outerDiameter, double innerDiameter, double height) {
    json coreJson;
    coreJson["name"] = "tape-wound fill factor test toroid";
    coreJson["manufacturerInfo"] = nullptr;
    coreJson["functionalDescription"]["type"] = "toroidal";
    coreJson["functionalDescription"]["gapping"] = json::array();
    coreJson["functionalDescription"]["numberStacks"] = 1;
    coreJson["functionalDescription"]["material"] = material;
    coreJson["functionalDescription"]["shape"] = {
        {"type", "custom"}, {"family", "t"}, {"magneticCircuit", "closed"}, {"aliases", json::array()},
        {"name", "T fill factor test"},
        {"dimensions", {{"A", {{"nominal", outerDiameter}}}, {"B", {{"nominal", innerDiameter}}}, {"C", {{"nominal", height}}}}}};
    return Core(coreJson);
}

double effective_area(const Core& core) {
    return core.get_processed_description()->get_effective_parameters().get_effective_area();
}

// The same core in a ferrite: its gross (shape) effective area, no fill factor applied.
double gross_effective_area(Core core) {
    core.get_mutable_functional_description().set_material(std::string("N87"));
    core.process_data();
    return effective_area(core);
}

// The whole core catalogue for the scope of a test, not only the cores in stock; restores the default on exit.
struct AllCataloguedCores {
    AllCataloguedCores() {
        Settings::GetInstance().reset();
        clear_loaded_cores();
        Settings::GetInstance().set_use_only_cores_in_stock(false);
        load_cores();
    }
    ~AllCataloguedCores() {
        clear_loaded_cores();
        Settings::GetInstance().reset();
    }
};

// True when Settings still offers the flat stacking factor that ABT #2234 removed.
template <typename SettingsType>
constexpr bool has_flat_stacking_factor = requires(SettingsType& settings) { settings.get_nanocrystalline_stacking_factor(); };

}  // namespace

TEST_CASE("Test_Tape_Wound_Fill_Factor_Maker_Published_Value", "[constructive-model][core][tape-wound-fill-factor][abt-2234]") {
    AllCataloguedCores allCataloguedCores;  // the catalogue cores below are not all stocked
    SECTION("Metglas 2605SA1 C core: Metglas' lamination factor >=84 % is used as 0.84") {
        auto core = find_core_by_name("C 4 - Metglas 2605SA1 - Ungapped");
        auto material = core.resolve_material();
        REQUIRE(Core::get_tape_wound_fill_factor(material, CoreShapeFamily::C).value() == 0.84);
        CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(0.84, 1e-9));
    }
    SECTION("1K107 toroid: AT&M's 75-84 % over 14-26 um, midpoint 0.795") {
        auto core = find_core_by_name("T 9.8/6.5/4.5 - 1K107 - Ungapped");
        CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(0.795, 1e-9));
    }
    SECTION("Magnetics Supermalloy toroids against the Magnetics catalogue Ae") {
        // Magnetics 2016 Tape Wound Cores Catalog, PDF page 17, core 50026: nominal ID 1.000 in, OD 1.500 in, HT 0.375 in;
        // Effective Core Area 0.299 / 0.454 / 0.514 / 0.545 cm2 at 0.0005 / 0.001 / 0.002 / 0.004 in.
        const double inch = 0.0254;
        std::vector<std::pair<std::string, double>> gaugesAndCatalogueArea = {
            {"Supermalloy 0.5mil", 0.299e-4}, {"Supermalloy 1mil", 0.454e-4},
            {"Supermalloy 2mil", 0.514e-4}, {"Supermalloy 4mil", 0.545e-4}};
        std::vector<double> factors = {0.50, 0.75, 0.85, 0.90};
        for (size_t index = 0; index < factors.size(); ++index) {
            auto core = tape_wound_toroid(gaugesAndCatalogueArea[index].first, 1.5 * inch, 1.0 * inch, 0.375 * inch);
            INFO(gaugesAndCatalogueArea[index].first);
            CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(factors[index], 1e-9));
            // MKF's IEC 60205 gross area of a rectangular ring section is 1.4 % below (OD - ID) / 2 x HT at this
            // OD / ID, so it lands within 3 % of the catalogue, not on it.
            CHECK_THAT(effective_area(core), Catch::Matchers::WithinRel(gaugesAndCatalogueArea[index].second, 0.03));
        }
    }
}

TEST_CASE("Test_Tape_Wound_Fill_Factor_VAC_Guide_Table", "[constructive-model][core][tape-wound-fill-factor][abt-2234]") {
    AllCataloguedCores allCataloguedCores;  // the catalogue cores below are not all stocked
    SECTION("Nanoperm toroid: 17-23 um resolves to 20 um, VAC Tabelle 10.4 between 0.015 (0.70) and 0.025 mm (0.80)") {
        auto core = find_core_by_name("T 50/40/20 - Nanoperm 30000 - Ungapped");
        auto material = core.resolve_material();
        REQUIRE(material.get_lamination_thickness());
        CHECK_THAT(resolve_dimensional_values(material.get_lamination_thickness().value()), Catch::Matchers::WithinRel(20e-6, 1e-9));
        CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(0.75, 1e-9));
    }
    SECTION("FT-3M C core: 14-18 um resolves to 16 um, VAC Tabelle 10.4: 0.70 + 0.1 x (0.016 - 0.015) / 0.010 = 0.71") {
        auto core = find_core_by_name("C 6.3 - FT-3M - Ungapped");
        CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(0.71, 1e-9));
    }
    SECTION("A ribbon thicker than VAC's VITROVAC column (0.025 mm) has no guide value and throws") {
        auto material = Core::resolve_material(std::string("Metglas 2605SA1"));
        material.set_name("Thick amorphous test ribbon");
        MAS::DimensionWithTolerance thickness;
        thickness.set_nominal(50e-6);
        material.set_lamination_thickness(thickness);
        CHECK_THROWS(Core::get_tape_wound_fill_factor(material, CoreShapeFamily::T));
    }
}

TEST_CASE("Test_Tape_Wound_Fill_Factor_Estimated_From_Other_Manufacturers", "[constructive-model][core][tape-wound-fill-factor][abt-2234]") {
    auto& logger = Logger::getInstance();
    logger.drainCollected();
    logger.enableCollector(LogLevel::WARNING);
    // VITROPERM 550 HF has no laminationThickness in MAS and VAC prints no fill factor for it. The other makers'
    // published nanocrystalline values: AT&M 0.795, King Magnetics 0.75, Pourleroi 0.76; mean 0.768333.
    auto core = tape_wound_toroid("VITROPERM 550 HF", 0.045, 0.030, 0.015);
    auto material = core.resolve_material();
    REQUIRE(!material.get_lamination_thickness());
    double expected = (0.795 + 0.75 + 0.76) / 3;
    CHECK_THAT(effective_area(core) / gross_effective_area(core), Catch::Matchers::WithinRel(expected, 1e-9));
    bool warned = false;
    for (const auto& record : logger.drainCollected()) {
        if (record.level == LogLevel::WARNING && record.message.find("VITROPERM 550 HF") != std::string::npos &&
            record.message.find("estimated") != std::string::npos) {
            warned = true;
        }
    }
    logger.disableCollector();
    CHECK(warned);

    // Pourleroi's own 0.76 is left out for a Pourleroi grade: other manufacturers only.
    auto pourleroi = Core::resolve_material(std::string("Pourleroi Nanocrystalline 40000"));
    REQUIRE(!pourleroi.get_lamination_thickness());
    CHECK_THAT(Core::get_tape_wound_fill_factor(pourleroi, CoreShapeFamily::T).value(),
               Catch::Matchers::WithinRel((0.795 + 0.75) / 2, 1e-9));
}

TEST_CASE("Test_Tape_Wound_Fill_Factor_No_Flat_Default", "[constructive-model][core][tape-wound-fill-factor][abt-2234]") {
    // The flat 0.80 stacking-factor setting is gone: nothing can apply one number to every tape-wound core.
    CHECK_FALSE(has_flat_stacking_factor<Settings>);

    // A material the lookup cannot place throws instead of taking a default: electrical steel wound as a toroid
    // with no tape thickness.
    auto steel = Core::resolve_material(std::string("Supermalloy 1mil"));
    steel.set_name("Supermalloy without thickness");
    steel.set_lamination_thickness(std::nullopt);
    CHECK_THROWS(Core::get_tape_wound_fill_factor(steel, CoreShapeFamily::T));

    // Laminated (E, U, ...) electrical steel stacks are not tape-wound and keep their gross section, as before.
    CHECK_FALSE(Core::get_tape_wound_fill_factor(steel, CoreShapeFamily::E).has_value());
    // Ferrite gets nothing.
    CHECK_FALSE(Core::get_tape_wound_fill_factor(Core::resolve_material(std::string("N87")), CoreShapeFamily::T).has_value());
}

TEST_CASE("Test_Tape_Wound_Fill_Factor_Source", "[constructive-model][core][tape-wound-fill-factor][abt-2234]") {
    // get_tape_wound_fill_factor_with_source gives the same value as get_tape_wound_fill_factor and names where it comes
    // from, one source per priority path, so a UI can show the k_f behind the Ae.
    SECTION("Maker's published value") {
        auto material = Core::resolve_material(std::string("Metglas 2605SA1"));
        auto fillFactor = Core::get_tape_wound_fill_factor_with_source(material, CoreShapeFamily::C);
        REQUIRE(fillFactor);
        CHECK(fillFactor->value == Core::get_tape_wound_fill_factor(material, CoreShapeFamily::C).value());
        CHECK(fillFactor->source.find("Metglas 2605SA1 datasheet") != std::string::npos);
    }
    SECTION("VAC table by ribbon thickness") {
        auto material = Core::resolve_material(std::string("VITROPERM 500F"));
        REQUIRE(material.get_lamination_thickness());
        auto fillFactor = Core::get_tape_wound_fill_factor_with_source(material, CoreShapeFamily::T);
        REQUIRE(fillFactor);
        CHECK(fillFactor->value == Core::get_tape_wound_fill_factor(material, CoreShapeFamily::T).value());
        CHECK(fillFactor->source.find("VAC Tabelle 10.4") != std::string::npos);
        CHECK(fillFactor->source.find("um ribbon") != std::string::npos);
    }
    SECTION("Estimate from other manufacturers") {
        auto material = Core::resolve_material(std::string("Pourleroi Nanocrystalline 40000"));
        REQUIRE(!material.get_lamination_thickness());
        auto fillFactor = Core::get_tape_wound_fill_factor_with_source(material, CoreShapeFamily::T);
        REQUIRE(fillFactor);
        CHECK_THAT(fillFactor->value, Catch::Matchers::WithinRel((0.795 + 0.75) / 2, 1e-9));
        CHECK(fillFactor->source.find("Estimated") != std::string::npos);
    }
    SECTION("Not tape-wound") {
        CHECK_FALSE(Core::get_tape_wound_fill_factor_with_source(Core::resolve_material(std::string("N87")),
                                                                 CoreShapeFamily::T).has_value());
    }
}
