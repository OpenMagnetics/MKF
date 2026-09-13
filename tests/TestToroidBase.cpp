// ABT #1173 (WP4 of docs/2026-09-12_manufacturing_features_proposal.md, MAS-RFC 0014 part B): a
// toroid base is a bobbin of family t. Its pins stand on the seating plane, base selection filters
// on the coated ring, Coil::assign_pins gives winding k the pins of side k, and R16 checks coating,
// barrier, TIW/FIW and the base limits.
//
// The bases below are TEST FIXTURES, not catalogue data: the MAS catalogue carries no toroid base
// yet, because no public drawing reviewed for ABT #1173 (TDK B64291/2/3, Lodestone Pacific,
// CB-Magnetics) states every field bobbin.json requires (the standoff above all). Their numbers are
// chosen to exercise each branch and are named for what they test.
//
// Cores are the MVB++ complete fixtures': buck_inductor_complete (T 10/6/4, High Flux 60, 8 turns x
// 3 parallels of Round 0.63 - Grade 1) and common_mode_choke_complete (T 25.3/14.8/10, N30, two
// windings of 15 turns of Round 1.25 - Grade 1).

#include "advisers/Manufacturability.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <set>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

json base_record(const std::string& name, const std::string& mounting, int pinsPerRow, double pitch, double rowDistance,
                 double standoff, std::optional<double> maximumOuterDiameter, std::optional<double> maximumHeight = std::nullopt,
                 std::optional<double> pocketInnerDiameter = std::nullopt, std::optional<double> boatWidth = std::nullopt) {
    json record;
    record["name"] = name;
    record["manufacturerInfo"] = {{"name", "Test fixture"}};
    auto& functional = record["functionalDescription"];
    functional["type"] = "standard";
    functional["family"] = "t";
    functional["shape"] = "T";
    functional["dimensions"] = json::object();
    functional["pinout"] = {{"numberPins", 2 * pinsPerRow},
                            {"numberRows", 2},
                            {"numberPinsPerRow", {pinsPerRow, pinsPerRow}},
                            {"pitch", pitch},
                            {"rowDistance", rowDistance},
                            {"pinDescription", {{"shape", "round"}, {"type", "tht"}, {"dimensions", {0.0008, 0.0008, 0.004}}}}};
    auto& base = functional["base"];
    base["mounting"] = mounting;
    base["length"] = {{"nominal", 0.03}};
    base["width"] = {{"nominal", 0.02}};
    base["height"] = {{"nominal", 0.004}};
    base["standoff"] = {{"nominal", standoff}};
    if (maximumOuterDiameter) {
        base["maximumCoreOuterDiameter"] = maximumOuterDiameter.value();
    }
    if (maximumHeight) {
        base["maximumCoreHeight"] = maximumHeight.value();
    }
    if (pocketInnerDiameter) {
        base["pocketInnerDiameter"] = {{"minimum", pocketInnerDiameter.value()}, {"maximum", pocketInnerDiameter.value() + 0.0003}};
    }
    if (boatWidth) {
        base["boatWidth"] = {{"nominal", boatWidth.value()}};
    }
    return record;
}

OpenMagnetics::Bobbin base_bobbin(const json& record) {
    return OpenMagnetics::Bobbin(record, false);
}

OpenMagnetics::Core buck_core() {
    return OpenMagneticsTesting::get_quick_core("T 10/6/4", json::parse("[]"), 1, "High Flux 60");
}

OpenMagnetics::Core cmc_core() {
    return OpenMagneticsTesting::get_quick_core("T 25.3/14.8/10", json::parse("[]"), 1, "N30");
}

OpenMagnetics::Coil wind_on(const OpenMagnetics::Bobbin& bobbin, const std::vector<std::tuple<std::string, int64_t, int64_t, std::string>>& windings) {
    json coilJson;
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    coilJson["bobbin"] = bobbinJson;
    coilJson["functionalDescription"] = json::array();
    for (const auto& [name, turns, parallels, wire] : windings) {
        coilJson["functionalDescription"].push_back(
            {{"name", name}, {"numberTurns", turns}, {"numberParallels", parallels}, {"isolationSide", "primary"}, {"wire", wire}});
    }
    OpenMagnetics::Coil coil(coilJson, false);
    REQUIRE(coil.wind());
    return coil;
}

std::string pin_of(const OpenMagnetics::Coil& coil, const std::string& winding, End end, std::optional<int64_t> parallel = std::nullopt) {
    for (const auto& w : coil.get_functional_description()) {
        if (w.get_name() != winding) {
            continue;
        }
        REQUIRE(w.get_connections());
        const auto connections = w.get_connections().value();
        for (const auto& connection : connections) {
            if (connection.get_end() == end && connection.get_parallel() == parallel) {
                REQUIRE(connection.get_pin_name());
                return connection.get_pin_name().value();
            }
        }
    }
    FAIL("no connection " << winding);
    return "";
}

}  // namespace

TEST_CASE("A catalogue toroid base is kept unprocessed and refuses processing until a ring is seated (ABT #1173)",
          "[constructive-model][bobbin][pins][abt1173]") {
    auto record = base_record("fixture base", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117);
    OpenMagnetics::Bobbin base(record, false);
    REQUIRE(base.has_base());
    CHECK(OpenMagnetics::Bobbin::is_toroid_base_record(base.get_functional_description().value()));
    CHECK_THROWS_WITH(OpenMagnetics::Bobbin(record, true), Catch::Matchers::ContainsSubstring("no ring seated"));
}

TEST_CASE("A toroid base missing its height is refused (ABT #1173)", "[constructive-model][bobbin][abt1173]") {
    auto record = base_record("fixture base", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117);
    record["functionalDescription"]["base"].erase("height");
    CHECK_THROWS(OpenMagnetics::Bobbin(record, false));
}

TEST_CASE("A ring seated on a horizontal base gets the base's pins on the seating plane (ABT #1173)",
          "[constructive-model][bobbin][pins][abt1173]") {
    auto core = buck_core();
    const double coating = core.get_coating_thickness();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture base", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117)));
    REQUIRE(bobbin.has_base());
    auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
    CHECK_THAT(dimensions.at("A"), Catch::Matchers::WithinAbs(0.010 + 2 * coating, 1e-9));
    CHECK_THAT(dimensions.at("C"), Catch::Matchers::WithinAbs(0.004 + 2 * coating, 1e-9));
    CHECK(bobbin.get_winding_window_shape() == WindingWindowShape::ROUND);

    auto pins = bobbin.get_processed_description()->get_pins().value();
    REQUIRE(pins.size() == 6);
    // Seating plane C/2 + standoff below the ring centre along -Y, pin centre half its 4 mm further.
    const double pinCentreY = -((0.004 + 2 * coating) / 2 + 0.0015 + 0.002);
    for (const auto& pin : pins) {
        CHECK_THAT(pin.get_coordinates()->at(1), Catch::Matchers::WithinAbs(pinCentreY, 1e-9));
        CHECK(!pin.get_rotation());
    }
    // Counter-clockwise: 1..3 along +X on row 0 (z = -rowDistance/2), 4..6 back along -X on row 1.
    auto pin1 = bobbin.get_pin("1").get_coordinates().value();
    auto pin3 = bobbin.get_pin("3").get_coordinates().value();
    auto pin4 = bobbin.get_pin("4").get_coordinates().value();
    auto pin6 = bobbin.get_pin("6").get_coordinates().value();
    CHECK_THAT(pin1[2], Catch::Matchers::WithinAbs(-0.0051, 1e-9));
    CHECK_THAT(pin4[2], Catch::Matchers::WithinAbs(0.0051, 1e-9));
    CHECK(pin1[0] < pin3[0]);
    CHECK(pin4[0] > pin6[0]);
    CHECK_THAT(pin1[0], Catch::Matchers::WithinAbs(pin6[0], 1e-9));

    // A round trip through json processes to the same pins.
    json bobbinJson;
    to_json(bobbinJson, bobbin);
    OpenMagnetics::Bobbin reloaded(bobbinJson);
    REQUIRE(reloaded.get_processed_description()->get_pins()->size() == 6);
    CHECK_THAT(reloaded.get_pin("4").get_coordinates()->at(1), Catch::Matchers::WithinAbs(pin4[1], 1e-12));
}

TEST_CASE("A ring on a vertical base has its pins along -Z, rows straddling the ring thickness (ABT #1173)",
          "[constructive-model][bobbin][pins][abt1173]") {
    auto core = cmc_core();
    const double coating = core.get_coating_thickness();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture vertical", "vertical", 2, 0.0125, 0.0075, 0.001, 0.03)));
    auto pins = bobbin.get_processed_description()->get_pins().value();
    REQUIRE(pins.size() == 4);
    const double pinCentreZ = -((0.0253 + 2 * coating) / 2 + 0.001 + 0.002);
    for (const auto& pin : pins) {
        CHECK_THAT(pin.get_coordinates()->at(2), Catch::Matchers::WithinAbs(pinCentreZ, 1e-9));
        CHECK_THAT(std::abs(pin.get_coordinates()->at(1)), Catch::Matchers::WithinAbs(0.00375, 1e-9));
        REQUIRE(pin.get_rotation());
        CHECK(pin.get_rotation()->at(0) == 90);
    }
}

TEST_CASE("The pin rail of a base without the ring dimension its mounting needs throws naming it (ABT #1173)",
          "[constructive-model][bobbin][abt1173]") {
    auto record = base_record("fixture base", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117);
    record["functionalDescription"]["dimensions"] = {{"A", {{"nominal", 0.010}}}, {"B", {{"nominal", 0.006}}}};
    OpenMagnetics::Bobbin bobbin(record, false);
    CHECK_THROWS_WITH(OpenMagnetics::Bobbin::get_toroid_base_pin_rail_distance(bobbin.get_functional_description().value()),
                      Catch::Matchers::ContainsSubstring("no ring dimension 'C'"));
}

TEST_CASE("Base selection returns only the bases whose limits hold the coated buck and CMC rings (ABT #1173)",
          "[constructive-model][bobbin][abt1173]") {
    std::vector<OpenMagnetics::Bobbin> candidates = {
        base_bobbin(base_record("max OD 11.7 mm horizontal", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117)),
        base_bobbin(base_record("pocket 12 mm horizontal, no max OD", "horizontal", 3, 0.00254, 0.0102, 0.0015, std::nullopt, std::nullopt, 0.012)),
        base_bobbin(base_record("max OD 11.7 mm but max height 3 mm", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117, 0.003)),
        base_bobbin(base_record("max OD 28 mm vertical, boat 4.8 mm", "vertical", 2, 0.0125, 0.0075, 0.001, 0.028, std::nullopt, std::nullopt, 0.0048)),
        base_bobbin(base_record("max OD 30 mm vertical, boat 12 mm", "vertical", 2, 0.0125, 0.0075, 0.001, 0.030, std::nullopt, std::nullopt, 0.012)),
        base_bobbin(base_record("max OD 30 mm horizontal", "horizontal", 2, 0.0125, 0.0075, 0.001, 0.030)),
        base_bobbin(base_record("no limit at all", "horizontal", 2, 0.0125, 0.0075, 0.001, std::nullopt)),
    };
    auto names = [](const std::vector<OpenMagnetics::Bobbin>& bases) {
        std::set<std::string> result;
        for (const auto& base : bases) {
            result.insert(base.get_name().value());
        }
        return result;
    };
    CHECK(names(find_toroid_bases_for_core(buck_core(), candidates)) ==
          std::set<std::string>{"max OD 11.7 mm horizontal", "pocket 12 mm horizontal, no max OD", "max OD 28 mm vertical, boat 4.8 mm",
                                "max OD 30 mm vertical, boat 12 mm", "max OD 30 mm horizontal"});
    CHECK(names(find_toroid_bases_for_core(cmc_core(), candidates)) ==
          std::set<std::string>{"max OD 30 mm vertical, boat 12 mm", "max OD 30 mm horizontal"});
    CHECK(names(find_toroid_bases_for_core(cmc_core(), candidates, OrientationEnum::VERTICAL)) ==
          std::set<std::string>{"max OD 30 mm vertical, boat 12 mm"});
    CHECK_THROWS(find_toroid_bases_for_core(OpenMagneticsTesting::get_quick_core("PQ 32/30", json::parse("[]"), 1, "N87"), candidates));
}

TEST_CASE("A common-mode choke on a 4-pin base puts winding 1 on pins 1-2 and winding 2 on 4-3 (ABT #1173)",
          "[constructive-model][coil][pins][abt1173]") {
    auto core = cmc_core();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture 4P", "horizontal", 2, 0.0125, 0.0254, 0.001, 0.030)));
    auto coil = wind_on(bobbin, {{"Winding 1", 15, 1, "Round 1.25 - Grade 1"}, {"Winding 2", 15, 1, "Round 1.25 - Grade 1"}});
    auto result = coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK(!result.skipped);
    CHECK(pin_of(coil, "Winding 1", End::START) == "1");
    CHECK(pin_of(coil, "Winding 1", End::FINISH) == "2");
    CHECK(pin_of(coil, "Winding 2", End::START) == "4");
    CHECK(pin_of(coil, "Winding 2", End::FINISH) == "3");
    // Idempotent.
    coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK(pin_of(coil, "Winding 2", End::START) == "4");
}

TEST_CASE("Three thick parallels on a 6-pin base take adjacent pins, starts on one side, finishes on the other (ABT #1173)",
          "[constructive-model][coil][pins][abt1173]") {
    auto core = buck_core();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture 6P", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117)));
    auto coil = wind_on(bobbin, {{"Primary", 8, 3, "Round 0.63 - Grade 1"}});
    coil.assign_pins(coil.resolve_bobbin(), core);
    // Round 0.63 is thick (>= 0.5 mm): at most 2 wires per pin (WP3 rule 4), so 3 strands split.
    CHECK(pin_of(coil, "Primary", End::START, 0) == "1");
    CHECK(pin_of(coil, "Primary", End::START, 1) == "2");
    CHECK(pin_of(coil, "Primary", End::START, 2) == "3");
    // Row 1 walked from -X: 6, 5, 4.
    CHECK(pin_of(coil, "Primary", End::FINISH, 0) == "6");
    CHECK(pin_of(coil, "Primary", End::FINISH, 1) == "5");
    CHECK(pin_of(coil, "Primary", End::FINISH, 2) == "4");
}

TEST_CASE("Three fine parallels share one start pin and one finish pin on the same side (ABT #1173)",
          "[constructive-model][coil][pins][abt1173]") {
    auto core = buck_core();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture 6P", "horizontal", 3, 0.00254, 0.0102, 0.0015, 0.0117)));
    auto coil = wind_on(bobbin, {{"Primary", 8, 3, "Round 0.2 - Grade 1"}});
    coil.assign_pins(coil.resolve_bobbin(), core);
    // The row is walked from the end the winding's first turn lies towards (ABT #1173): a first turn
    // further along +X than the last starts on the +X corner pin 3, otherwise on the -X corner pin 1; the
    // finish takes the inner pin 2 either way. (Which way wind() lays this winding depends on settings an
    // earlier test in the same process may leave behind, so the expectation is read from the turns.)
    std::vector<double> parallel0X;
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        if (turn.get_parallel() == 0) {
            parallel0X.push_back(turn.get_coordinates()[0]);
        }
    }
    REQUIRE(parallel0X.size() == 8);
    REQUIRE(std::abs(parallel0X.front() - parallel0X.back()) > 1e-4);
    INFO("parallel 0 first turn x " << parallel0X.front() * 1e3 << " mm, last " << parallel0X.back() * 1e3 << " mm");
    CHECK(pin_of(coil, "Primary", End::START) == (parallel0X.front() > parallel0X.back() ? "3" : "1"));
    CHECK(pin_of(coil, "Primary", End::FINISH) == "2");
}

TEST_CASE("More windings than base sides, or two windings forced onto one side, throw (ABT #1173)",
          "[constructive-model][coil][pins][abt1173]") {
    auto core = cmc_core();
    auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture 4P", "horizontal", 2, 0.0125, 0.0254, 0.001, 0.030)));
    auto three = wind_on(bobbin, {{"A", 5, 1, "Round 1.25 - Grade 1"}, {"B", 5, 1, "Round 1.25 - Grade 1"}, {"C", 5, 1, "Round 1.25 - Grade 1"}});
    CHECK_THROWS_WITH(three.assign_pins(three.resolve_bobbin(), core), Catch::Matchers::ContainsSubstring("has two sides"));

    auto forced = wind_on(bobbin, {{"Winding 1", 15, 1, "Round 1.25 - Grade 1"}, {"Winding 2", 15, 1, "Round 1.25 - Grade 1"}});
    auto windings = forced.get_functional_description();
    ConnectionElement start;
    start.set_end(End::START);
    start.set_pin_name("1");
    windings[0].set_connections(std::vector<ConnectionElement>{start});
    ConnectionElement otherStart;
    otherStart.set_end(End::START);
    otherStart.set_pin_name("2");
    windings[1].set_connections(std::vector<ConnectionElement>{otherStart});
    forced.set_functional_description(windings);
    CHECK_THROWS_WITH(forced.assign_pins(forced.resolve_bobbin(), core), Catch::Matchers::ContainsSubstring("both have user pins on side 0"));
}

TEST_CASE("R16: toroid base limits and the declared core coating, red and green (ABT #1173)",
          "[adviser][manufacturability][abt1173]") {
    Manufacturability manufacturability;
    auto make = [&](OpenMagnetics::Core core, double maximumOuterDiameter) {
        auto bobbin = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, base_bobbin(base_record("fixture 4P", "horizontal", 2, 0.0125, 0.0254, 0.001, maximumOuterDiameter)));
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(wind_on(bobbin, {{"Winding 1", 15, 1, "Round 1.25 - Grade 1"}, {"Winding 2", 15, 1, "Round 1.25 - Grade 1"}}));
        return magnetic;
    };
    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_nominal(230);
    auto basic = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, Cti::GROUP_I, IsolationClass::BASIC, mainSupplyVoltage,
                                                                   OvervoltageCategory::II, PollutionDegree::PD2,
                                                                   {InsulationStandards::IEC_623681}, 230, 345, 100000);

    SECTION("header: the coated ring exceeds the base's maximum outer diameter") {
        auto magnetic = make(cmc_core(), 0.020);
        auto finding = manufacturability.evaluate_r16_toroid_insulation_and_header(magnetic, basic);
        CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("exceeds the base's maximum"));
    }
    SECTION("coating: an undeclared coating does not count as basic insulation") {
        auto magnetic = make(cmc_core(), 0.030);
        auto finding = manufacturability.evaluate_r16_toroid_insulation_and_header(magnetic, basic);
        CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("declares no coating"));
    }
    SECTION("green: declared coating, ring within the base") {
        auto core = cmc_core();
        auto coreJson = json();
        to_json(coreJson, core);
        coreJson["functionalDescription"]["coating"] = "epoxy";
        OpenMagnetics::Core coated(coreJson);
        auto magnetic = make(coated, 0.030);
        auto finding = manufacturability.evaluate_r16_toroid_insulation_and_header(magnetic, basic);
        CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("outer diameter fits the base"));
    }
    SECTION("not a toroid") {
        auto magnetic = make(cmc_core(), 0.030);
        magnetic.set_core(OpenMagneticsTesting::get_quick_core("PQ 32/30", json::parse("[]"), 1, "N87"));
        CHECK(manufacturability.evaluate_r16_toroid_insulation_and_header(magnetic, basic).get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
    }
}

// ---------------------------------------------------------------------------------------------------
// The REAL catalogue bases (MAS b7cf597): Lodestone Pacific HTM460-6 / HTM600-6 / HTM850-6, horizontal,
// 3+3 rectangular THT tabs 0.762 x 0.3556 x 7.1374 mm at a 3.556 mm pitch, rows 16.20 / 21.53 / 28.70 mm
// apart, standoff 2.032 mm, maximum ring OD 12.7 / 19.05 / 25.4 mm. Unlike the fixtures above these come
// from get_bobbins(), so the test also proves MKF embeds them.
namespace {

const std::vector<std::string> kHtmBases = {"Toroid base 18.28x11.68x2.03 6P horizontal",
                                            "Toroid base 23.79x15.24x2.03 6P horizontal",
                                            "Toroid base 30.68x21.59x2.03 6P horizontal"};

double coated_outer_diameter(const OpenMagnetics::Core& core) {
    auto shape = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    return shape.at("A") + 2 * core.get_coating_thickness();
}

// The base a design takes: the smallest catalogue base that holds the coated ring (tightest maximum
// outer diameter). find_toroid_bases_for_core only filters; this picks.
OpenMagnetics::Bobbin smallest_base_for(const OpenMagnetics::Core& core) {
    auto bases = find_toroid_bases_for_core(core);
    REQUIRE(!bases.empty());
    std::stable_sort(bases.begin(), bases.end(), [](const OpenMagnetics::Bobbin& a, const OpenMagnetics::Bobbin& b) {
        return a.get_base().get_maximum_core_outer_diameter().value() < b.get_base().get_maximum_core_outer_diameter().value();
    });
    return bases.front();
}

// Mean position of a winding's turns in MKF's ring frame (x, y); MVB++ draws (x, 0, y), so y < 0 faces
// row 0 (z = -rowDistance / 2).
std::pair<double, double> winding_centroid(const OpenMagnetics::Coil& coil, const std::string& winding) {
    double sumX = 0, sumY = 0;
    size_t count = 0;
    const auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        if (turn.get_winding() != winding) {
            continue;
        }
        REQUIRE(turn.get_coordinate_system() == CoordinateSystem::CARTESIAN);
        sumX += turn.get_coordinates()[0];
        sumY += turn.get_coordinates()[1];
        ++count;
    }
    REQUIRE(count > 0);
    return {sumX / count, sumY / count};
}

void check_htm_pins(OpenMagnetics::Bobbin seated, const OpenMagnetics::Bobbin& record) {
    const auto base = record.get_base();
    const auto pinout = record.get_functional_description()->get_pinout().value();
    const double pitch = std::get<double>(pinout.get_pitch().value());
    const double rowDistance = pinout.get_row_distance().value();
    const auto ring = flatten_dimensions(seated.get_functional_description()->get_dimensions());
    const double pinLength = pinout.get_pin_description()->get_dimensions()[2];
    const double pinCentreY = -(ring.at("C") / 2 + resolve_dimensional_values(base.get_standoff()) + pinLength / 2);
    auto pins = seated.get_processed_description()->get_pins().value();
    REQUIRE(pins.size() == 6);
    // Counter-clockwise from row 0: 1, 2, 3 along +X at z = -rowDistance/2; 4, 5, 6 back along -X at +rowDistance/2.
    const std::map<std::string, std::pair<double, double>> expected = {
        {"1", {-pitch, -rowDistance / 2}}, {"2", {0, -rowDistance / 2}}, {"3", {pitch, -rowDistance / 2}},
        {"4", {pitch, rowDistance / 2}},   {"5", {0, rowDistance / 2}},  {"6", {-pitch, rowDistance / 2}}};
    for (const auto& [name, xz] : expected) {
        auto pin = seated.get_pin(name);
        auto coordinates = pin.get_coordinates().value();
        CHECK_THAT(coordinates[0], Catch::Matchers::WithinAbs(xz.first, 1e-9));
        CHECK_THAT(coordinates[1], Catch::Matchers::WithinAbs(pinCentreY, 1e-9));
        CHECK_THAT(coordinates[2], Catch::Matchers::WithinAbs(xz.second, 1e-9));
        CHECK(pin.get_shape() == PinShape::RECTANGULAR);
        CHECK_THAT(pin.get_dimensions()[0], Catch::Matchers::WithinAbs(0.000762, 1e-12));
        CHECK_THAT(pin.get_dimensions()[1], Catch::Matchers::WithinAbs(0.0003556, 1e-12));
        CHECK_THAT(pin.get_dimensions()[2], Catch::Matchers::WithinAbs(0.0071374, 1e-12));
    }
}

}  // namespace

TEST_CASE("The buck_inductor_complete ring seats on the smallest Lodestone HTM base that holds it, with its six tabs (ABT #1173)",
          "[constructive-model][bobbin][coil][pins][abt1173]") {
    auto core = buck_core();
    const double coatedOuterDiameter = coated_outer_diameter(core);
    INFO("coated T 10/6/4 outer diameter " << coatedOuterDiameter * 1e3 << " mm");

    std::set<std::string> offered;
    for (const auto& base : find_toroid_bases_for_core(core)) {
        offered.insert(base.get_name().value());
    }
    // All three HTM bases hold a ~10.3 mm ring (max OD 12.7 / 19.05 / 25.4 mm); the smallest is taken.
    for (const auto& name : kHtmBases) {
        CHECK(offered.count(name) == 1);
    }
    auto record = smallest_base_for(core);
    CHECK(record.get_name().value() == kHtmBases[0]);
    CHECK(record.get_manufacturer_info()->get_reference().value() == "HTM460-6");

    auto seated = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, record);
    check_htm_pins(seated, record);

    auto coil = wind_on(seated, {{"Primary", 8, 3, "Round 0.63 - Grade 1"}});
    auto result = coil.assign_pins(coil.resolve_bobbin(), core);
    CHECK(!result.skipped);
    // Round 0.63 takes at most two wires a pin, so the three strands split: starts on side 0 (z < 0),
    // finishes on side 1, each side walked from -X.
    CHECK(pin_of(coil, "Primary", End::START, 0) == "1");
    CHECK(pin_of(coil, "Primary", End::START, 1) == "2");
    CHECK(pin_of(coil, "Primary", End::START, 2) == "3");
    CHECK(pin_of(coil, "Primary", End::FINISH, 0) == "6");
    CHECK(pin_of(coil, "Primary", End::FINISH, 1) == "5");
    CHECK(pin_of(coil, "Primary", End::FINISH, 2) == "4");
}

TEST_CASE("A common-mode choke on a real HTM base: the fixture ring fits none, a T 16/9.6/6.3 takes HTM600-6 and each winding its facing side (ABT #1173)",
          "[constructive-model][bobbin][coil][pins][abt1173]") {
    // common_mode_choke_complete's T 25.3/14.8/10: coated, it exceeds HTM850-6's 25.4 mm, so no HTM base holds it.
    auto fixtureCore = cmc_core();
    INFO("coated T 25.3/14.8/10 outer diameter " << coated_outer_diameter(fixtureCore) * 1e3 << " mm");
    CHECK(coated_outer_diameter(fixtureCore) > 0.0254);
    for (const auto& base : find_toroid_bases_for_core(fixtureCore)) {
        CHECK(std::find(kHtmBases.begin(), kHtmBases.end(), base.get_name().value()) == kHtmBases.end());
    }

    auto core = OpenMagneticsTesting::get_quick_core("T 16/9.6/6.3", json::parse("[]"), 1, "N30");
    const double coatedOuterDiameter = coated_outer_diameter(core);
    INFO("coated T 16/9.6/6.3 outer diameter " << coatedOuterDiameter * 1e3 << " mm");
    REQUIRE(coatedOuterDiameter > 0.0127);   // HTM460-6 is rejected on its 12.7 mm limit
    std::set<std::string> offered;
    for (const auto& base : find_toroid_bases_for_core(core)) {
        offered.insert(base.get_name().value());
    }
    CHECK(offered.count(kHtmBases[0]) == 0);
    CHECK(offered.count(kHtmBases[1]) == 1);
    CHECK(offered.count(kHtmBases[2]) == 1);
    auto record = smallest_base_for(core);
    CHECK(record.get_manufacturer_info()->get_reference().value() == "HTM600-6");

    auto seated = OpenMagnetics::Bobbin::create_toroid_bobbin_on_base(core, record);
    check_htm_pins(seated, record);

    // Wound and pinned the way every consumer gets it (MVB++ included): magnetic_autocomplete, which
    // winds the round window with its own section orientation/alignment and then runs assign_pins.
    json magneticJson;
    json coreJson;
    to_json(coreJson, core);
    magneticJson["core"] = coreJson;
    json bobbinJson;
    to_json(bobbinJson, seated);
    magneticJson["coil"]["bobbin"] = bobbinJson;
    magneticJson["coil"]["functionalDescription"] = json::array();
    for (const auto* name : {"Winding 1", "Winding 2"}) {
        magneticJson["coil"]["functionalDescription"].push_back(
            {{"name", name}, {"numberTurns", 10}, {"numberParallels", 1}, {"isolationSide", "primary"}, {"wire", "Round 0.80 - Grade 1"}});
    }
    OpenMagnetics::Magnetic magnetic(magneticJson);
    // ABT #1237: autocomplete assigns pins only with coil_connect_leads_to_pins on (off by default).
    OpenMagnetics::SettingsGuard<bool> connectLeads(OpenMagnetics::Settings::GetInstance(),
                                                    &OpenMagnetics::Settings::get_coil_connect_leads_to_pins,
                                                    &OpenMagnetics::Settings::set_coil_connect_leads_to_pins, true);
    auto completed = magnetic_autocomplete(magnetic);
    auto coil = completed.get_coil();
    REQUIRE(coil.get_turns_description());
    const auto w1 = winding_centroid(coil, "Winding 1");
    const auto w2 = winding_centroid(coil, "Winding 2");
    INFO("Winding 1 centroid (" << w1.first * 1e3 << ", " << w1.second * 1e3 << ") mm, Winding 2 (" << w2.first * 1e3 << ", "
                                << w2.second * 1e3 << ") mm, pins W1 " << pin_of(coil, "Winding 1", End::START) << "/"
                                << pin_of(coil, "Winding 1", End::FINISH) << " W2 " << pin_of(coil, "Winding 2", End::START) << "/"
                                << pin_of(coil, "Winding 2", End::FINISH));
    // Each winding takes the side its half of the ring faces, so no lead crosses to the far side.
    const auto side_of = [&](const std::string& pinName) { return seated.get_pin(pinName).get_coordinates()->at(2) < 0 ? -1 : 1; };
    CHECK(side_of(pin_of(coil, "Winding 1", End::START)) == (w1.second < 0 ? -1 : 1));
    CHECK(side_of(pin_of(coil, "Winding 1", End::FINISH)) == (w1.second < 0 ? -1 : 1));
    CHECK(side_of(pin_of(coil, "Winding 2", End::START)) == (w2.second < 0 ? -1 : 1));
    CHECK(side_of(pin_of(coil, "Winding 2", End::FINISH)) == (w2.second < 0 ? -1 : 1));
    CHECK((w1.second < 0) != (w2.second < 0));
    // Within its side, the start takes the pin towards the end of the row the winding's first turn lies
    // towards, so the start and finish leads do not cross on their way to the row.
    const auto ends_x = [&](const std::string& winding) {
        std::optional<double> first;
        double last = 0;
        const auto allTurns = coil.get_turns_description().value();
        for (const auto& turn : allTurns) {
            if (turn.get_winding() != winding) {
                continue;
            }
            if (!first) {
                first = turn.get_coordinates()[0];
            }
            last = turn.get_coordinates()[0];
        }
        return std::make_pair(first.value(), last);
    };
    const auto pin_x = [&](const std::string& pinName) { return seated.get_pin(pinName).get_coordinates()->at(0); };
    for (const auto* winding : {"Winding 1", "Winding 2"}) {
        const auto [firstX, lastX] = ends_x(winding);
        INFO(winding << " first turn x " << firstX * 1e3 << " mm, last " << lastX * 1e3 << " mm");
        REQUIRE(std::abs(firstX - lastX) > 1e-3);
        CHECK((pin_x(pin_of(coil, winding, End::START)) < pin_x(pin_of(coil, winding, End::FINISH))) == (firstX < lastX));
    }
}
