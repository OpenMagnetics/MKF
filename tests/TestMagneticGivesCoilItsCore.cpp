// ABT #1533: a Magnetic that holds both a core and a coil gives the coil that core, whatever way the
// Magnetic was built (JSON, MAS::Magnetic, set_core/set_coil in either order, copies). The ABT #1487
// outside-link rules (real winding, axially stacked sections split in series) need the core outline
// around the winding build; the coil used to get it only from magnetic_autocomplete /
// wind_magnetic_coil_as_described, so simulate() -- which builds Magnetic(json) -- threw
// "[COIL_NOT_PROCESSED] ... the coil was given no core (Coil::set_core_geometry)" on the owner's
// PQ 65/60 design once real winding was on. A coil that really has no core still throws.

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/MagneticSimulator.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <fstream>
#include <source_location>
#include <string>
#include <vector>

using json = nlohmann::json;

namespace {

// The owner's PQ 65/60 file from the web (S,P,S,P,S,P,S; primary 45 turns split 1/3 per section,
// outside series links), as the WebFrontend fixture pq65_stacked_sections_real_winding_unfit.json.
json load_abt1533_design() {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append("abt1533_pq65_stacked_sections_real_winding_unfit.json"));
    REQUIRE(file.good());
    return json::parse(file);
}

void check_coil_has_core_of(const OpenMagnetics::Magnetic& magnetic) {
    REQUIRE(magnetic.has_core());
    REQUIRE(magnetic.has_coil());
    const auto& coreGeometry = magnetic.get_coil().get_core_geometry();
    REQUIRE(coreGeometry);
    CHECK(coreGeometry->get_shape_name() == magnetic.get_core().get_shape_name());
    CHECK(coreGeometry->get_width() == magnetic.get_core().get_width());
    CHECK(coreGeometry->get_depth() == magnetic.get_core().get_depth());
}

}  // namespace

TEST_CASE("A Magnetic gives its coil its core however it is built (ABT #1533)",
          "[constructive-model][magnetic][coil][real-winding][abt1533]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const auto design = load_abt1533_design();
    const auto& magneticJson = design.at("magnetic");

    SECTION("from JSON, the way simulate() builds it") {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        check_coil_has_core_of(magnetic);
    }
    SECTION("from JSON through nlohmann's get<>") {
        auto magnetic = magneticJson.get<OpenMagnetics::Magnetic>();
        check_coil_has_core_of(magnetic);
    }
    SECTION("from a JSON list of magnetics") {
        std::vector<OpenMagnetics::Magnetic> magnetics;
        OpenMagnetics::from_json(json::array({magneticJson, magneticJson}), magnetics);
        REQUIRE(magnetics.size() == 2);
        for (const auto& magnetic : magnetics) {
            check_coil_has_core_of(magnetic);
        }
    }
    SECTION("from a whole MAS document") {
        OpenMagnetics::Mas mas(design);
        check_coil_has_core_of(mas.get_magnetic());
    }
    SECTION("set_coil then set_core, and set_core then set_coil") {
        OpenMagnetics::Core core(magneticJson.at("core"));
        OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
        REQUIRE_FALSE(coil.get_core_geometry());
        {
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_coil(coil);
            magnetic.set_core(core);
            check_coil_has_core_of(magnetic);
        }
        {
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            check_coil_has_core_of(magnetic);
        }
    }
    SECTION("a new core replaces the one the coil was given") {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        auto otherCore = OpenMagnetics::Core::create_quick_core("PQ 50/50", "3C95");
        REQUIRE(otherCore.get_shape_name() != magnetic.get_core().get_shape_name());
        magnetic.set_core(otherCore);
        check_coil_has_core_of(magnetic);
        CHECK(magnetic.get_coil().get_core_geometry()->get_shape_name() == "PQ 50/50");
    }
    SECTION("copies of the magnetic and of its coil keep the core") {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        auto magneticCopy = magnetic;
        check_coil_has_core_of(magneticCopy);
        auto coilCopy = magnetic.get_coil();
        REQUIRE(coilCopy.get_core_geometry());
        CHECK(coilCopy.get_core_geometry()->get_shape_name() == magnetic.get_core().get_shape_name());
    }
    settings.reset();
}

// The throw stays for a coil that really has no core: the outside-link room of this design cannot
// be judged without the core outline, and nothing may assume one.
TEST_CASE("A coil with no core still refuses to judge an outside series link (ABT #1533)",
          "[constructive-model][coil][real-winding][abt1533]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const auto design = load_abt1533_design();
    settings.set_coil_use_real_winding_geometry(true);
    // As stored (wound): not re-wound here, only its connections are laid out.
    OpenMagnetics::Coil coil(design.at("magnetic").at("coil"), false);
    REQUIRE_FALSE(coil.get_core_geometry());
    bool threwNoCore = false;
    try {
        OpenMagnetics::WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(coil);
    }
    catch (const std::exception& exception) {
        const std::string what = exception.what();
        INFO(what);
        threwNoCore = what.find("the coil was given no core") != std::string::npos;
    }
    CHECK(threwNoCore);
    settings.reset();
}

// The web's path: magnetic_autocomplete with real winding on, the autocompleted magnetic serialized
// (the web hands it back as JSON), then Magnetic(json) into the loss and simulation models.
OpenMagnetics::Magnetic autocomplete_and_reload(const json& design, const OpenMagnetics::Inputs& inputs) {
    OpenMagnetics::Settings::GetInstance().set_coil_use_real_winding_geometry(true);
    auto autocompleted = OpenMagnetics::magnetic_autocomplete(OpenMagnetics::Magnetic(design.at("magnetic")), json{}, inputs);
    json autocompletedJson;
    OpenMagnetics::to_json(autocompletedJson, autocompleted);
    return OpenMagnetics::Magnetic(autocompletedJson);
}

// The owner's file as it is: the real winding does not fit (its stack outgrows the window), but the
// winding losses -- the WindingOhmicLosses connection-length path that threw -- are computed.
TEST_CASE("Winding losses of the stacked-section PQ 65 design reloaded from JSON with real winding on (ABT #1533)",
          "[physical-model][winding-losses][real-winding][abt1533]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const auto design = load_abt1533_design();
    OpenMagnetics::Mas mas(design);
    auto magnetic = autocomplete_and_reload(design, mas.get_inputs());

    auto connectionLength = OpenMagnetics::WindingOhmicLosses::calculate_connection_length_per_winding_per_parallel(magnetic.get_coil());
    REQUIRE(connectionLength.size() == 2);
    // The primary is split in series over its stacked sections: its outside links are copper.
    const auto primaryIndex = magnetic.get_mutable_coil().get_winding_index_by_name("Primary");
    CHECK(connectionLength[primaryIndex][0] > 0);

    auto operatingPoint = mas.get_inputs().get_operating_points()[0];
    auto losses = OpenMagnetics::MagneticSimulator().calculate_winding_losses(operatingPoint, magnetic, 25.0);
    CHECK(losses.get_winding_losses() > 0);
    settings.reset();
}

// The same design with a primary that fits (0.5 mm round wire, still split in series over the
// stacked sections, so its outside links still need the core outline): the whole simulate() runs.
TEST_CASE("simulate() runs a fitting stacked-section PQ 65 design reloaded from JSON with real winding on (ABT #1533)",
          "[processor][magnetic-simulator][real-winding][abt1533]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    auto design = load_abt1533_design();
    bool primaryFound = false;
    for (auto& winding : design["magnetic"]["coil"]["functionalDescription"]) {
        if (winding["name"] == "Primary") {
            winding["wire"] = "Round 0.5 - Grade 1";
            primaryFound = true;
        }
    }
    REQUIRE(primaryFound);
    OpenMagnetics::Mas mas(design);
    auto magnetic = autocomplete_and_reload(design, mas.get_inputs());

    auto simulated = OpenMagnetics::MagneticSimulator().simulate(mas.get_inputs(), magnetic);
    REQUIRE(simulated.get_outputs().size() == 1);
    const auto& outputs = simulated.get_outputs()[0];
    REQUIRE(outputs.get_winding_losses());
    CHECK(outputs.get_winding_losses()->get_winding_losses() > 0);
    REQUIRE(outputs.get_core_losses());
    CHECK(outputs.get_core_losses()->get_core_losses() > 0);
    settings.reset();
}

// ABT #1535 (Alf, 2026-09-30, "block the run-in row"): an outside series link whose end needs a stub
// climbs its turn's own column to the run's row and turns the corner there, at the connection
// plane. No station of that column may sit on that row: the turn wound there passes the link's lane
// just before the plane at the run's height (measured in MVB++ on this design: the first link came in
// at y = 2.798 mm, the top row of Primary section 1 layer 1, 0.0378 mm from turn 15's wrap where
// 1.726 mm is needed). MKF blocked that row for the section's outward layers but not for the landing
// column itself.
TEST_CASE("An outside series link's run row is blocked in its own turn's column (ABT #1535)",
          "[constructive-model][coil][real-winding][abt1535]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    const auto design = load_abt1533_design();
    OpenMagnetics::Mas mas(design);
    auto magnetic = autocomplete_and_reload(design, mas.get_inputs());
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    const auto turns = coil.get_turns_description().value();
    auto turnByName = [&](const std::string& name) -> const MAS::Turn& {
        for (const auto& turn : turns) {
            if (turn.get_name() == name) {
                return turn;
            }
        }
        FAIL("no turn named '" << name << "'");
        throw std::logic_error("unreachable");
    };

    const auto layout = coil.get_connection_layout();
    size_t outsideLinks = 0;
    size_t endsWithStub = 0;
    for (const auto& route : layout.routes) {
        if (route.kind != OpenMagnetics::ConnectionKind::EDGE_CONTINUATION || !route.exitSlot) {
            continue;
        }
        ++outsideLinks;
        // {x1,y1} {x1,runA} {linkX,runA} {linkX,runB} {x2,runB} {x2,y2}, with a stub of zero length
        // dropped: the second point and the one before last still stand on the two run rows.
        REQUIRE(route.waypoints.size() >= 4);
        const std::vector<std::pair<std::string, double>> ends = {
            {route.fromTurn, route.waypoints[1][1]},   // exit end: its run row
            {route.toTurn, route.waypoints[route.waypoints.size() - 2][1]},     // landing end: the run-in row
        };
        for (const auto& [turnName, runRow] : ends) {
            const auto& endTurn = turnByName(turnName);
            // The rows are laid one wire apart at the wire's outer height -- the pitch MKF blocks
            // and spreads with (a turn's own dimensions carry a fraction of a micron more).
            const double wireHeight =
                coil.get_wires()[coil.get_winding_index_by_name(endTurn.get_winding())].get_maximum_outer_height();
            if (std::abs(endTurn.get_coordinates()[1] - runRow) <= wireHeight / 2) {
                continue;   // the turn sits on the row itself: no stub, no corner in the column
            }
            ++endsWithStub;
            REQUIRE(endTurn.get_layer());
            for (const auto& turn : turns) {
                if (turn.get_layer() != endTurn.get_layer()) {
                    continue;
                }
                const double clearance = std::abs(turn.get_coordinates()[1] - runRow);
                INFO(route.fromTurn << " -> " << route.toTurn << ": '" << turn.get_name() << "' of "
                                    << endTurn.get_layer().value() << " at y " << turn.get_coordinates()[1] * 1e3
                                    << " mm, run row at " << runRow * 1e3 << " mm");
                // One micron of float bookkeeping between the spread stations and the counted
                // depths; the defect this pins is a station ON the row (0.15 um away, not 1.865 mm).
                CHECK(clearance >= wireHeight - 1e-6);
            }
        }
    }
    CHECK(outsideLinks == 2);    // Primary section 0 -> 1 and 1 -> 2
    CHECK(endsWithStub >= 1);    // the landing of the first link climbs its column
    settings.reset();
}
