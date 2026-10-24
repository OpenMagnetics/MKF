// ABT #1423: a terminal lead whose radial bend cannot fit on the connection face's straight leaves
// TANGENTIALLY, straight off the end of the lateral face straight (Coil.h, TangentDeparture).
//
// Fixture: 18_stacked (E 70/33/32 x2, 30 turns of 5 x 1 mm rectangular wire wound on edge) on its
// derived custom former, whose column corner (9.9925 mm) is the IEC 60317-0-2 flexibility radius
// 12.5425 mm minus the turn's 2.55 mm standoff. The turn's corner is therefore 12.5425 mm, which
// leaves 15.4132 + 2.55 - 12.5425 = 5.4207 mm of connection-face straight either side of the
// crossing, and the lead's own edgewise bend is the same 12.5425 mm: MKF planned a radial exit whose
// bend cannot be drawn anywhere on that face (MVB++: "entrance corner offset 0.012543 m lies outside
// the -Z face straight").
//
// Every expected value below is derived here from the frame and the stations, independently of the
// code under test: the connection face's half straight columnWidth + s - R, the lateral face's
// columnDepth + s - R, the cut-off copper (half straight + quarter corner), the charged turn lengths,
// and the lead's line (collinear with the face straight, climb included).

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "physical_models/WireBend.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "support/Exceptions.h"
#include "processors/ResultsAnalyzer.h"
#include "TestingUtils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <source_location>

using json = nlohmann::json;

namespace {

OpenMagnetics::Coil wind_tangent_fixture() {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append("abt1423_rect_wire_rect_column.json"));
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    // ABT #2012: a rectangular lead is no longer planned at its wire's IEC 60317-0-2 Table 6 EDGEWISE
    // mandrel (a coating test, not a forming limit), so on its own this lead's bend (its coated radius)
    // fits on the connection face and it leaves radially. The tangent departure this test is about is
    // kept by DECLARING that same 12.5425 mm radius as the consumer's minimum lead bend, which is what
    // the reported MVB++ refusal drew it at. Every expected value below is still derived from the frame.
    const OpenMagnetics::Wire fixtureWire(magneticJson.at("coil").at("functionalDescription").at(0).at("wire"));
    settings.set_coil_lead_minimum_bend_radius(
        OpenMagnetics::WireBend::get_flexibility_bend_radius_if_standardised(fixtureWire, OpenMagnetics::BendAxis::EDGEWISE).value());
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    auto wound = enriched.get_coil();
    // The settings stay until the caller has read the layout: get_connection_layout re-plans the
    // leads under the ambient lead settings (the caller resets them).
    REQUIRE(wound.is_real_winding_blocking_applied());
    return wound;
}

}  // namespace

TEST_CASE("Real winding: a rectangular wire's terminal leads leave tangentially off the lateral face (18_stacked)",
          "[constructive-model][coil][real-winding][abt1423]") {
    auto coil = wind_tangent_fixture();
    const auto turns = coil.get_turns_description().value();
    const auto first = turns.front();
    const auto last = turns.back();
    REQUIRE(first.get_name() == "Primary parallel 0 turn 0");
    REQUIRE(last.get_name() == "Primary parallel 0 turn 29_ending");
    const auto frame = coil.get_wound_column_frame_for_section(first.get_section().value());
    const double x = std::abs(first.get_coordinates()[0]);
    const double standoff = x - frame.columnWidth;
    const double bend = frame.cornerRadius + standoff;
    const double connectionHalfStraight = frame.columnWidth + standoff - bend;
    const double lateralHalfStraight = frame.columnDepth + standoff - bend;
    const double removed = connectionHalfStraight + std::numbers::pi / 2 * bend;
    const double fullLength = 4 * frame.columnDepth + 4 * frame.columnWidth + 8 * standoff +
                              (2 * std::numbers::pi - 8) * bend;
    INFO("frame: width " << frame.columnWidth * 1e3 << " mm, depth " << frame.columnDepth * 1e3
                         << " mm, corner " << frame.cornerRadius * 1e3 << " mm; turn corner " << bend * 1e3
                         << " mm, connection-face half straight " << connectionHalfStraight * 1e3 << " mm");
    CHECK(bend == Catch::Approx(0.0125425).margin(1e-9));
    CHECK(connectionHalfStraight == Catch::Approx(0.00542071).margin(1e-8));

    const auto layout = coil.get_connection_layout();
    size_t tangent = 0;
    for (const auto& route : layout.routes) {
        const bool entrance = route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE;
        if (!entrance && route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        INFO((entrance ? "entrance" : "exit") << ": planned bend " << route.plannedBendRadius * 1e3 << " mm");
        // The premise: the radial exit's bend does not fit on the connection-face straight.
        REQUIRE(route.plannedBendRadius > connectionHalfStraight);
        REQUIRE(route.tangentDeparture.has_value());
        ++tangent;
        const auto& departure = route.tangentDeparture.value();
        CHECK(departure.face == (entrance ? -1 : 1));
        CHECK_FALSE(route.exitSlot.has_value());
        CHECK(departure.distanceAlongFace == Catch::Approx(lateralHalfStraight).margin(1e-12));
        CHECK(departure.turnLengthRemoved == Catch::Approx(removed).margin(1e-12));
        REQUIRE(departure.point.size() == 3);
        CHECK(departure.point[0] == Catch::Approx(departure.face * x).margin(1e-12));
        CHECK(departure.point[2] == Catch::Approx(-lateralHalfStraight).margin(1e-12));
        // The pitch-true helix at the lateral face's end: the station moved by the cut-off share of
        // the revolution's advance (into the turn for the entrance, back from it for the exit).
        const auto& station = entrance ? turns[0] : turns[turns.size() - 1];
        const double climb = entrance ? (turns[1].get_coordinates()[1] - turns[0].get_coordinates()[1])
                                      : (turns.back().get_coordinates()[1] - turns[turns.size() - 2].get_coordinates()[1]);
        const double expectedY = station.get_coordinates()[1] + (entrance ? 1.0 : -1.0) * climb / fullLength * removed;
        CHECK(departure.point[1] == Catch::Approx(expectedY).margin(1e-12));
        // The lead IS the lateral face straight, continued: collinear with the wire's travel there.
        const double travelY = climb / fullLength;   // dy per metre of travel
        const double travelZ = entrance ? 1.0 : -1.0;
        const double outwardSign = entrance ? -1.0 : 1.0;   // outward = against the travel at the entrance
        const double norm = std::sqrt(1 + travelY * travelY);
        CHECK(departure.direction[0] == 0.0);
        CHECK(departure.direction[1] == Catch::Approx(outwardSign * travelY / norm).margin(1e-15));
        CHECK(departure.direction[2] == Catch::Approx(outwardSign * travelZ / norm).margin(1e-15));
        CHECK(departure.direction[2] < 0);   // out of the front of the window
        // It runs to the window border moved to the connection face's depth.
        const double border = entrance ? route.waypoints.front()[0] : route.waypoints.back()[0];
        CHECK(border == Catch::Approx(x).margin(1e-12));   // the projection sits at the turn's radius
        CHECK(departure.end[2] < departure.point[2]);
        CHECK(route.routedLength == Catch::Approx(departure.length).margin(1e-9));
        // And it stays inside the window band: its axial extent is within the station and the
        // departure point, both of which are.
        const double yLow = std::min(departure.point[1], departure.end[1]);
        const double yHigh = std::max(departure.point[1], departure.end[1]);
        const double stationY = station.get_coordinates()[1];
        CHECK(((yLow >= stationY - 1e-12) || (yHigh <= stationY + 1e-12)));
        std::cout << "[abt1423] " << (entrance ? "entrance" : "exit") << ": leaves the "
                  << (departure.face < 0 ? "-X" : "+X") << " face at (" << departure.point[0] * 1e3 << ", "
                  << departure.point[1] * 1e3 << ", " << departure.point[2] * 1e3 << ") mm, lead "
                  << departure.length * 1e3 << " mm, revolution '" << departure.shortenedTurn << "' shortened by "
                  << departure.turnLengthRemoved * 1e3 << " mm" << std::endl;
    }
    CHECK(tangent == 2);

    // The charge: the first and last revolutions lose exactly the cut-off copper, the others none.
    size_t full = 0;
    for (size_t t = 1; t < turns.size(); ++t) {
        const bool shortened = t == 1 || t + 1 == turns.size();
        INFO(turns[t].get_name() << " length " << turns[t].get_length() * 1e3 << " mm");
        CHECK(turns[t].get_length() ==
              Catch::Approx(shortened ? fullLength - removed : fullLength).margin(1e-9));
        full += shortened ? 0 : 1;
    }
    CHECK(full == turns.size() - 3);
    OpenMagnetics::Settings::GetInstance().reset();
}

// ABT #2012: the clash check of a tangent lead judges each other turn WHERE its corner meets the
// lead's line, at its own pitch-true height there. Comparing the lead's height at the departure
// point (already climbed slope * removed along the cut-off copper) with the next revolution's
// height at its crossing station put them closer than one coated height in EVERY closely wound
// layer, so every such tangent lead was refused -- round wire too -- on a clash that is not there:
// along the whole straight the next revolution runs one pitch above the lead.
//
// Fixture: round wire on the rectangular column of an E 25/13/7, its terminal leads declared (by
// the consumer's minimum lead bend radius) too wide for the connection-face straight, so they leave
// tangentially. The radius is derived from the wound frame, not pinned.
namespace {

OpenMagnetics::Coil round_wire_coil(int64_t numberTurns) {
    auto wires = std::vector<OpenMagnetics::Wire>({OpenMagnetics::find_wire_by_name("Round 0.4 - Grade 1")});
    return OpenMagneticsTesting::get_quick_coil({numberTurns}, {1}, "E 25/13/7", 1, MAS::WindingOrientation::OVERLAPPING,
                                                MAS::WindingOrientation::OVERLAPPING, MAS::CoilAlignment::CENTERED,
                                                MAS::CoilAlignment::CENTERED, wires);
}

size_t number_conduction_layers(const OpenMagnetics::Coil& coil) {
    const auto layers = coil.get_layers_description().value();   // bound: never range over a temporary
    size_t numberTurnLayers = 0;
    for (const auto& layer : layers) {
        numberTurnLayers += layer.get_type() == MAS::ElectricalType::CONDUCTION ? 1 : 0;
    }
    return numberTurnLayers;
}

// The most turns of round_wire_coil that still wind in one layer (real winding spreads a part-filled
// layer over the window; the layer is closely wound once it is full).
int64_t fullest_single_layer() {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    int64_t numberTurns = 0;
    for (int64_t n = 20; n <= 60; ++n) {
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        auto probe = round_wire_coil(n);
        if (number_conduction_layers(probe) > 1) {
            break;
        }
        numberTurns = n;
    }
    return numberTurns;
}

// The connection face's half straight of the first turn, read off a plain real-winding wind.
double first_turn_face_half_straight(int64_t numberTurns) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    settings.reset();
    settings.set_coil_use_real_winding_geometry(true);
    auto coil = round_wire_coil(numberTurns);
    REQUIRE(coil.wind());
    const auto turns = coil.get_turns_description().value();
    const auto frame = coil.get_wound_column_frame_for_section(turns.front().get_section().value());
    const double x = std::abs(turns.front().get_coordinates()[0]);
    const double standoff = x - frame.columnWidth;
    const double bend = frame.cornerRadius + standoff;
    return frame.columnWidth + standoff - bend;
}

}  // namespace

TEST_CASE("Real winding: a closely wound round-wire layer's tangent lead is not refused by its own next revolution, "
          "while a turn wound across the lead is reported as a collision",
          "[constructive-model][coil][real-winding][tangent-departure][abt-2012]") {
    auto& settings = OpenMagnetics::Settings::GetInstance();

    SECTION("one closely wound layer: both leads leave tangentially") {
        // The fullest single layer: closely wound.
        const int64_t numberTurns = fullest_single_layer();
        REQUIRE(numberTurns >= 20);
        const double faceHalfStraight = first_turn_face_half_straight(numberTurns);
        REQUIRE(faceHalfStraight > 0);
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        settings.set_coil_lead_minimum_bend_radius(1.5 * faceHalfStraight);
        OpenMagnetics::Coil coil;
        bool fits = false;
        // get_quick_coil winds the coil it builds, under these settings
        REQUIRE_NOTHROW(coil = round_wire_coil(numberTurns));
        REQUIRE_NOTHROW(fits = coil.wind());
        REQUIRE(fits);
        const auto turns = coil.get_turns_description().value();
        // The premise: the next revolution's station is closer to the lead's departure height than one
        // coated diameter -- by less than the lead's climb over the cut-off copper -- so judging the
        // station against the departure point (the old check) collides, while along the whole straight
        // the next revolution runs one pitch above the lead.
        const auto frame = coil.get_wound_column_frame_for_section(turns.front().get_section().value());
        const double x = std::abs(turns.front().get_coordinates()[0]);
        const double standoff = x - frame.columnWidth;
        const double bend = frame.cornerRadius + standoff;
        const double removed = frame.columnWidth + standoff - bend + std::numbers::pi / 2 * bend;
        const double fullLength = 4 * frame.columnDepth + 4 * frame.columnWidth + 8 * standoff +
                                  (2 * std::numbers::pi - 8) * bend;
        const double coatedDiameter = turns[0].get_dimensions().value()[1];
        const double pitch = std::abs(turns[1].get_coordinates()[1] - turns[0].get_coordinates()[1]);
        INFO(numberTurns << " turns, pitch " << pitch * 1e3 << " mm, coated diameter " << coatedDiameter * 1e3 << " mm");
        REQUIRE(pitch >= coatedDiameter - 1e-12);
        REQUIRE(pitch - coatedDiameter < pitch / fullLength * removed);
        size_t tangent = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE &&
                route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT) {
                continue;
            }
            REQUIRE(route.plannedBendRadius > faceHalfStraight);   // its bend does not fit on the face
            CHECK(route.tangentDeparture.has_value());
            tangent += route.tangentDeparture ? 1 : 0;
        }
        CHECK(tangent == 2);
    }

    SECTION("a second layer wound back over the first lies across the entrance lead: drawn, and the collision reported") {
        // Alf, 2026-10-24: a lead that cannot be routed without a collision is still drawn along its
        // route -- the customer sees the design -- and analyze_results reports the collision. The clash
        // check is the judge of that finding, not a refusal. Wound in U order: the second layer comes
        // back over the first, so its last turns cross the entrance lead, and the layer link is a
        // tangential U run, not a Z dragback (a tangent lead off a face a dragback lane rides over is a
        // different, unmodelled case).
        // The fullest two layers: the most turns that still wind in two.
        const int64_t singleLayer = fullest_single_layer();
        REQUIRE(singleLayer >= 20);
        int64_t numberTurns = 0;
        for (int64_t n = 2 * singleLayer; n > singleLayer; --n) {
            settings.reset();
            settings.set_coil_use_real_winding_geometry(true);
            if (number_conduction_layers(round_wire_coil(n)) == 2) {
                numberTurns = n;
                break;
            }
        }
        REQUIRE(numberTurns > 0);
        const double faceHalfStraight = first_turn_face_half_straight(numberTurns);
        settings.reset();
        settings.set_coil_use_real_winding_geometry(true);
        OpenMagnetics::Coil coil = round_wire_coil(numberTurns);
        auto bobbin = std::get<OpenMagnetics::Bobbin>(coil.get_bobbin());
        auto processed = bobbin.get_processed_description().value();
        auto windingWindows = processed.get_winding_windows();
        windingWindows[0].set_winding_order(MAS::WindingOrder::U);
        processed.set_winding_windows(windingWindows);
        bobbin.set_processed_description(processed);
        coil.set_bobbin(bobbin);
        settings.set_coil_lead_minimum_bend_radius(1.5 * faceHalfStraight);
        bool fits = false;
        REQUIRE_NOTHROW(fits = coil.wind());
        REQUIRE(fits);
        REQUIRE(number_conduction_layers(coil) == 2);
        const auto turns = coil.get_turns_description().value();
        REQUIRE(coil.get_winding_order(turns[0].get_section().value()) == MAS::WindingOrder::U);
        const std::string windingName = coil.get_functional_description()[0].get_name();
        // Every turn and both leads are drawn, tangentially; the entrance runs under the second layer.
        std::optional<OpenMagnetics::TangentDeparture> entranceDeparture;
        size_t terminals = 0;
        size_t collidingLeads = 0;
        for (const auto& route : coil.get_connection_layout().routes) {
            if (route.kind != OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE &&
                route.kind != OpenMagnetics::ConnectionKind::TERMINAL_EXIT) {
                continue;
            }
            ++terminals;
            REQUIRE(route.tangentDeparture.has_value());
            if (route.kind == OpenMagnetics::ConnectionKind::TERMINAL_ENTRANCE) {
                entranceDeparture = route.tangentDeparture;
            }
            const auto& departure = route.tangentDeparture.value();
            collidingLeads += departure.collisions.empty() ? 0 : 1;
            // Each collision is an outer turn whose copper really overlaps the lead's.
            const double leadX = std::abs(departure.point[0]);
            for (const auto& collision : departure.collisions) {
                const auto turn = std::find_if(turns.begin(), turns.end(),
                                               [&](const MAS::Turn& candidate) { return candidate.get_name() == collision.turnName; });
                REQUIRE(turn != turns.end());
                CHECK(std::abs(turn->get_coordinates()[0]) > leadX);
                CHECK(collision.overlap > 0);
                CHECK(collision.overlap == Catch::Approx(turns[0].get_dimensions().value()[1] / 2 + turn->get_dimensions().value()[1] / 2 -
                                                         std::abs(collision.leadHeight - collision.turnHeight)).margin(1e-12));
            }
        }
        REQUIRE(terminals == 2);
        REQUIRE(entranceDeparture.has_value());
        REQUIRE_FALSE(entranceDeparture->collisions.empty());

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(OpenMagneticsTesting::get_quick_core("E 25/13/7", OpenMagneticsTesting::get_ground_gap(0.001), 1, "N87"));
        magnetic.set_coil(coil);
        const auto findings = OpenMagnetics::ResultsAnalyzer::analyze_lead_collisions(magnetic);
        REQUIRE(findings.size() == collidingLeads);   // one per colliding lead
        bool entranceReported = false;
        for (const auto& finding : findings) {
            CHECK(finding.code == OpenMagnetics::ResultFindingCode::LEAD_COLLIDES_WITH_TURNS);
            CHECK_FALSE(finding.operatingPointIndex);
            CHECK(finding.message.find("winding '" + windingName + "'") != std::string::npos);
            CHECK(finding.message.find("drawn across turn '") != std::string::npos);
            CHECK(finding.message.find("overlaps by") != std::string::npos);
            CHECK(finding.message.find("ABT") == std::string::npos);
            entranceReported |= finding.message.find("entrance lead of '" + turns[0].get_name() + "'") != std::string::npos;
        }
        CHECK(entranceReported);
    }
    settings.reset();
}
