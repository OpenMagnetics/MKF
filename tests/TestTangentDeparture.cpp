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
#include "support/Settings.h"
#include "support/Utils.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

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
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    auto wound = enriched.get_coil();
    settings.reset();
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
}
