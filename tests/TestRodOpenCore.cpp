#include <source_location>
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/InitialPermeability.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Core.h"
#include "constructive_models/CorePiece.h"
#include "constructive_models/Magnetic.h"
#include "TestingUtils.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;
using json = nlohmann::json;

// ABT #933: the ROD family. A bare cylinder is the most open magnetic circuit MKF carries — no
// return limb at all — so it exercises the open-core path in a way a drum does not: for a drum the
// flange envelope dominates and the winding is pinned to the groove, whereas a rod's answer is set
// by the cylinder's own slenderness plus the length of whatever winding the caller puts on it.
//
// MAS carries ZERO catalogued rod records, so every fixture here is an inline shape. That is the
// intended production route as well until rod records exist.

namespace {

// A rod shape stated the way a caller must state it: A = diameter, B = length, H = optional bore.
json rod_shape_json(double diameter, double length, std::optional<double> bore = std::nullopt) {
    json shape;
    shape["name"] = "Test rod " + std::to_string(diameter * 1000) + "x" + std::to_string(length * 1000);
    shape["type"] = "custom";
    shape["family"] = "rod";
    shape["magneticCircuit"] = "open";
    shape["aliases"] = json::array();
    shape["dimensions"] = json();
    shape["dimensions"]["A"] = diameter;
    shape["dimensions"]["B"] = length;
    if (bore) {
        shape["dimensions"]["H"] = *bore;
    }
    return shape;
}

Core rod_core(double diameter, double length, const std::string& materialName,
              std::optional<double> bore = std::nullopt, json gapping = json::array()) {
    json coreJson;
    coreJson["name"] = "Test rod core";
    coreJson["functionalDescription"] = json();
    coreJson["functionalDescription"]["name"] = "Test rod core";
    coreJson["functionalDescription"]["type"] = "openShape";
    coreJson["functionalDescription"]["material"] = materialName;
    coreJson["functionalDescription"]["shape"] = rod_shape_json(diameter, length, bore);
    coreJson["functionalDescription"]["gapping"] = gapping;
    coreJson["functionalDescription"]["numberStacks"] = 1;
    Core core(coreJson);
    core.process_data();
    return core;
}

// Bozorth/Chen axial demagnetising factor of the equivalent prolate spheroid, recomputed here so
// the tests check the model against an INDEPENDENT expression rather than against itself.
double reference_demagnetizing_factor(double lengthOverDiameter) {
    double m = lengthOverDiameter;
    double s = sqrt(m * m - 1);
    return (1 / (m * m - 1)) * ((m / s) * log(m + s) - 1);
}

}  // namespace

// The piece geometry must be the cylinder and nothing else: a rod's IEC reduction is exact, so
// le is the rod length and Ae the rod cross-section to machine precision. Any deviation means the
// shape constants picked up a spurious section.
TEST_CASE("Test_Rod_Effective_Parameters_Are_The_Cylinder", "[physical-model][rod][open-core]") {
    settings.reset();
    clear_databases();

    double diameter = 0.006;
    double length = 0.030;
    auto core = rod_core(diameter, length, "3C97");

    // A rod is a single open piece: nothing may be doubled the way a two-piece set is.
    REQUIRE(core.get_functional_description().get_type() == CoreType::OPEN_SHAPE);

    auto effective = core.get_processed_description()->get_effective_parameters();
    double expectedArea = std::numbers::pi / 4 * pow(diameter, 2);
    CHECK_THAT(effective.get_effective_length(), Catch::Matchers::WithinRel(length, 1e-9));
    CHECK_THAT(effective.get_effective_area(), Catch::Matchers::WithinRel(expectedArea, 1e-9));
    CHECK_THAT(effective.get_effective_volume(), Catch::Matchers::WithinRel(length * expectedArea, 1e-9));
    CHECK_THAT(effective.get_minimum_area(), Catch::Matchers::WithinRel(expectedArea, 1e-9));

    // The column IS the rod, and the winding window runs its whole length.
    REQUIRE(core.get_columns().size() == 1);
    CHECK_THAT(core.get_columns()[0].get_area(), Catch::Matchers::WithinRel(expectedArea, 1e-4));
    CHECK_THAT(core.get_processed_description()->get_winding_windows()[0].get_height().value(), Catch::Matchers::WithinRel(length, 1e-9));

    // A bore carries no flux and must come straight out of the area.
    auto bored = rod_core(diameter, length, "3C97", 0.002);
    double boredArea = std::numbers::pi / 4 * (pow(diameter, 2) - pow(0.002, 2));
    CHECK_THAT(bored.get_processed_description()->get_effective_parameters().get_effective_area(),
               Catch::Matchers::WithinRel(boredArea, 1e-9));
    settings.reset();
}

// THE defining property of an open magnetic circuit, and the reason a rod cannot be modelled by
// the closed-circuit reluctance path: its effective permeability is set by GEOMETRY, not by the
// material. Doubling mu_i barely moves the answer once mu_i >> 1/N_d, and the ceiling it presses
// against is 1/N_d. This is the same physics the Fair-Rite drum test pins (same published AL for
// mu_i 800 and mu_i 2000), stated directly on the shape it comes from.
TEST_CASE("Test_Rod_Effective_Permeability_Is_Geometry_Limited", "[physical-model][rod][open-core]") {
    settings.reset();
    clear_databases();

    double diameter = 0.004;
    double length = 0.040;              // slender rod, l/d = 10
    double demagnetizingFactor = reference_demagnetizing_factor(length / diameter);
    // Sanity on the reference itself before it is used to judge the model.
    CHECK_THAT(demagnetizingFactor, Catch::Matchers::WithinRel(0.0203, 0.02));

    double numberTurns = 100;
    double rodArea = std::numbers::pi / 4 * pow(diameter, 2);
    double vacuumPermeability = 4e-7 * std::numbers::pi;

    // Fully wound: the bracket collapses and the closed form below IS the model.
    auto softInductance = MagnetizingInductance::calculate_rod_magnetizing_inductance(
        rod_core(diameter, length, "3C90"), numberTurns, length, 25);       // mu_i ~ 2300
    auto hardInductance = MagnetizingInductance::calculate_rod_magnetizing_inductance(
        rod_core(diameter, length, "N87"), numberTurns, length, 25);        // mu_i ~ 2200, different grade

    // The geometry ceiling: mu_eff can never exceed 1/N_d however permeable the ferrite is.
    double ceilingInductance = vacuumPermeability / demagnetizingFactor * rodArea / length * pow(numberTurns, 2);
    CHECK(softInductance < ceilingInductance);
    CHECK(hardInductance < ceilingInductance);
    // ...and with mu_i ~ 2000 against 1/N_d ~ 49 it must be pressed hard against it.
    CHECK(softInductance > 0.9 * ceilingInductance);

    // Two different ferrite grades give nearly the same inductance — the material-insensitivity
    // that is the signature of an open circuit. A closed-circuit model cannot do this.
    CHECK_THAT(softInductance, Catch::Matchers::WithinRel(hardInductance, 0.05));

    // A stubby rod is far more demagnetised than a slender one: same material, same turns, the
    // slender rod must win by a wide margin.
    auto stubbyInductance = MagnetizingInductance::calculate_rod_magnetizing_inductance(
        rod_core(0.020, 0.020, "3C90"), numberTurns, 0.020, 25);
    double stubbyPerLength = stubbyInductance / (std::numbers::pi / 4 * pow(0.020, 2) / 0.020);
    double slenderPerLength = softInductance / (rodArea / length);
    CHECK(slenderPerLength > 2 * stubbyPerLength);
    settings.reset();
}

// The winding length is the one coil-side input the model takes, and its behaviour must be the
// textbook one: shortening the coil on a FIXED rod raises inductance per turn (the air-solenoid
// term goes as 1/l), while leaving the demagnetising factor — a property of the rod alone —
// untouched. The bracket between the two physical bounds is what widens to carry the growing
// uncertainty, and it must vanish on a fully wound rod.
TEST_CASE("Test_Rod_Winding_Length_Enters_And_Brackets_Honestly", "[physical-model][rod][open-core]") {
    settings.reset();
    clear_databases();

    double diameter = 0.005;
    double length = 0.025;
    double numberTurns = 50;
    auto core = rod_core(diameter, length, "3C90");

    double fullyWound = MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, length, 25);
    double halfWound = MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, length / 2, 25);
    double quarterWound = MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, length / 4, 25);

    CHECK(halfWound > fullyWound);
    CHECK(quarterWound > halfWound);

    // The rise is bounded and sub-proportional: the geometric mean of a 1/l term and an
    // l-independent term goes as 1/sqrt(l), so halving the coil may not double the inductance.
    CHECK(halfWound < 2 * fullyWound);
    CHECK_THAT(halfWound / fullyWound, Catch::Matchers::WithinRel(sqrt(2.0), 0.02));

    // A rod cannot be gapped: its return path is already air.
    auto gapping = json::array();
    gapping.push_back(json{{"type", "additive"}, {"length", 0.0005}});
    auto gappedCore = rod_core(diameter, length, "3C90", std::nullopt, gapping);
    CHECK_THROWS(MagnetizingInductance::calculate_rod_magnetizing_inductance(gappedCore, numberTurns, length, 25));

    // A winding longer than the rod is not described by this model, and must say so rather than
    // silently extrapolate.
    CHECK_THROWS(MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, length * 1.5, 25));
    CHECK_THROWS(MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, 0, 25));
    settings.reset();
}

// End to end: a rod reaching the ordinary entry point must be ROUTED to the open-core model, not
// to the closed-circuit reluctance machinery. Before ABT #933 this threw in CorePiece::factory;
// the failure mode to guard against now is the quieter one, where a rod loads but is silently
// treated as a closed circuit and reports an inductance many times too high.
TEST_CASE("Test_Rod_Routes_To_Open_Core_Model_End_To_End", "[physical-model][rod][open-core]") {
    settings.reset();
    clear_databases();

    double diameter = 0.005;
    double length = 0.025;
    int64_t numberTurns = 40;
    auto core = rod_core(diameter, length, "3C90");

    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core);
    OpenMagnetics::Coil coil;
    coil.set_bobbin(bobbin);
    MAS::CoilFunctionalDescription winding;
    winding.set_number_turns(numberTurns);
    winding.set_number_parallels(1);
    winding.set_name("primary");
    winding.set_isolation_side(IsolationSide::PRIMARY);
    winding.set_wire("Round 0.25 - Grade 1");
    coil.set_functional_description({winding});

    auto result = MagnetizingInductance().calculate_inductance_and_magnetic_flux_density(core, coil, nullptr);
    REQUIRE(result.first.get_method_used() == "RodOpenCoreDemagnetizingFactor");

    double routed = result.first.get_magnetizing_inductance().get_nominal().value();

    // The winding length the route used must be the COIL's, read from its bobbin window — not the
    // rod's length. A quick bobbin has walls, so its window is genuinely shorter than the rod, and
    // that difference is exactly what proves the coil (not the core) supplied the number.
    double windingLength = coil.resolve_bobbin().get_winding_window_dimensions(0)[1];
    CHECK(windingLength < length);
    double direct = MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, windingLength, 25);
    CHECK_THAT(routed, Catch::Matchers::WithinRel(direct, 1e-9));

    // ...and a shorter coil on the same rod means more inductance, so the routed answer must
    // exceed what a full-length winding would give.
    double fullyWound = MagnetizingInductance::calculate_rod_magnetizing_inductance(core, numberTurns, length, 25);
    CHECK(routed > fullyWound);

    // The reported envelope is the two physical bounds themselves, so it must bracket the answer.
    CHECK(result.first.get_magnetizing_inductance().get_minimum().value() <= routed);
    CHECK(result.first.get_magnetizing_inductance().get_maximum().value() >= routed);

    // The whole point of the open-core route: a closed-circuit reading of the same cylinder
    // (le/Ae over the full material permeability) is enormously larger. Pin the gap so a
    // regression that quietly reverts the routing cannot pass.
    auto effective = core.get_processed_description()->get_effective_parameters();
    double initialPermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), 25.0);
    double closedCircuit = 4e-7 * std::numbers::pi * initialPermeability *
                           effective.get_effective_area() / effective.get_effective_length() * pow(numberTurns, 2);
    CHECK(closedCircuit > 20 * routed);
    settings.reset();
}
