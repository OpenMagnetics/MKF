// ABT #1401: under real winding every wrap is a helix, and two side-by-side wraps whose stations
// sit `s` apart along the column only clear each other by s * cos(alpha), where
// tan(alpha) = (the wrap's realised advance per revolution) / (the turn's own length). So the
// axial spacing of adjacent stations in an overlapping round/litz layer must be at least
// OD / cos(alpha), with alpha the slope the wrap is actually DRAWN with -- the one a consumer
// (MVB++) reads back from the stations MKF publishes:
//   - a turn followed by another turn of its conductor in the same layer climbs to it in one
//     revolution: advance = |y_next - y|;
//   - a conductor's LAST turn entered from another layer is the steep exit landing: one
//     revolution from the arrival height to its own station: advance = |y - y_previous|;
//   - a layer's last turn that leaves through a layer link winds at its layer's own advance, the
//     climb from its previous turn: advance = |y - y_previous|.
// Measured before the fix on the 21_interleaved_flyback fixture: the Secondary's layer 0 stations
// sat 148 nm inside the coated envelope (0.679257 mm against the 0.679405 mm its 2.077 mm/rev
// wrap needs), and its layer 1 steep exit landing 4.30 um inside (0.687192 mm against
// 0.691570 mm at 12.438 mm/rev).

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <functional>
#include <iomanip>
#include <sstream>
#include <numbers>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <tuple>
#include <source_location>

using json = nlohmann::json;

namespace {

json load_sibling_slope_test_data(const std::string& name) {
    std::ifstream file(std::filesystem::path{std::source_location::current().file_name()}
                           .parent_path()
                           .append("testData")
                           .append(name));
    REQUIRE(file.good());
    return json::parse(file);
}

// Wound exactly as MVB++'s magnetic_autocomplete_safe winds a real-winding build: Core and Coil
// from the stored MAS (the coil NOT re-wound by its constructor), the lead bend policy MVB++
// declares (round-corner factor 1.05, no minimum), then magnetic_autocomplete.
OpenMagnetics::Coil wind_real(const std::string& fixture) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const auto masJson = load_sibling_slope_test_data(fixture);
    const auto& magneticJson = masJson.at("magnetic");
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    settings.reset();
    REQUIRE(enriched.get_coil().is_real_winding_blocking_applied());
    return enriched.get_coil();
}

struct SiblingPair {
    std::string layer;
    std::string turnA, turnB;
    double spacing = 0;       // axial station spacing
    double advance = 0;       // steeper of the two wraps' advance per revolution
    double turnLength = 0;    // that wrap's own turn length
    double od = 0;
    double perpendicular = 0; // spacing * cos(alpha)
    bool landingPair = false; // both turns are steep exit landings
};

std::vector<SiblingPair> measure_sibling_pairs(OpenMagnetics::Coil& coil) {
    const auto turns = coil.get_turns_description().value();
    const auto layers = coil.get_layers_description().value();
    auto wires = coil.get_wires();

    std::map<std::pair<std::string, int64_t>, std::vector<size_t>> turnsByConductor;
    for (size_t t = 0; t < turns.size(); ++t) {
        turnsByConductor[{turns[t].get_winding(), turns[t].get_parallel()}].push_back(t);
    }
    // A turns-description entry is a STATION: MVB++ draws one revolution (a wrap) from each entry
    // to the next entry of the same conductor when both sit in the same layer, climbing the
    // stations' axial difference over that revolution, and charged the length of the entry it
    // arrives at (the entry a layer is entered on is a zero-length arrival marker). Per station:
    // the steepest wrap touching it, as tangent = advance / length, and whether it belongs to the
    // steep exit landing (the one revolution a conductor's last layer spends between the arrival
    // and its final station).
    std::vector<double> tangentOf(turns.size(), -1.0);
    std::vector<double> advanceOf(turns.size(), 0.0);
    std::vector<double> lengthOf(turns.size(), 0.0);
    std::vector<bool> steepLanding(turns.size(), false);
    for (const auto& [conductor, sequence] : turnsByConductor) {
        for (size_t s = 0; s + 1 < sequence.size(); ++s) {
            const auto& from = turns[sequence[s]];
            const auto& to = turns[sequence[s + 1]];
            if (from.get_layer().value() != to.get_layer().value()) {
                continue;   // a layer link, not a wrap
            }
            const double advance = std::abs(to.get_coordinates()[1] - from.get_coordinates()[1]);
            const double length = to.get_length();
            REQUIRE(length > 0);
            const double tangent = advance / length;
            const bool landing = s + 2 == sequence.size() && s > 0 &&
                                 turns[sequence[s - 1]].get_layer().value() != from.get_layer().value();
            for (size_t t : {sequence[s], sequence[s + 1]}) {
                if (tangent > tangentOf[t]) {
                    tangentOf[t] = tangent;
                    advanceOf[t] = advance;
                    lengthOf[t] = length;
                }
                steepLanding[t] = steepLanding[t] || landing;
            }
        }
    }

    std::vector<SiblingPair> pairs;
    for (const auto& layer : layers) {
        if (layer.get_type() != MAS::ElectricalType::CONDUCTION ||
            layer.get_orientation() != MAS::WindingOrientation::OVERLAPPING) {
            continue;
        }
        const auto windingIndex =
            coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
        const auto wireType = wires[windingIndex].get_type();
        if (wireType != MAS::WireType::ROUND && wireType != MAS::WireType::LITZ) {
            continue;
        }
        const double od = wires[windingIndex].get_maximum_outer_height();
        std::vector<size_t> layerTurns;
        for (size_t t = 0; t < turns.size(); ++t) {
            if (turns[t].get_layer() && turns[t].get_layer().value() == layer.get_name()) {
                layerTurns.push_back(t);
            }
        }
        std::sort(layerTurns.begin(), layerTurns.end(), [&](size_t a, size_t b) {
            return turns[a].get_coordinates()[1] < turns[b].get_coordinates()[1];
        });
        for (size_t k = 0; k + 1 < layerTurns.size(); ++k) {
            const size_t a = layerTurns[k];
            const size_t b = layerTurns[k + 1];
            SiblingPair pair;
            pair.layer = layer.get_name();
            pair.turnA = turns[a].get_name();
            pair.turnB = turns[b].get_name();
            pair.od = od;
            pair.spacing = turns[b].get_coordinates()[1] - turns[a].get_coordinates()[1];
            // A station no wrap touches has no slope to judge; that must fail loudly.
            REQUIRE(tangentOf[a] >= 0);
            REQUIRE(tangentOf[b] >= 0);
            const size_t steeper = tangentOf[a] >= tangentOf[b] ? a : b;
            const double worstTangent = tangentOf[steeper];
            pair.advance = advanceOf[steeper];
            pair.turnLength = lengthOf[steeper];
            pair.perpendicular = pair.spacing * std::cos(std::atan(worstTangent));
            pair.landingPair = steepLanding[a] && steepLanding[b];
            pairs.push_back(pair);
        }
    }
    return pairs;
}

void require_siblings_clear(OpenMagnetics::Coil& coil, const std::optional<std::string>& landingLayer) {
    const auto pairs = measure_sibling_pairs(coil);
    REQUIRE(pairs.size() > 0);
    size_t landingPairs = 0;
    size_t failures = 0;
    const SiblingPair* tightest = nullptr;
    for (const auto& pair : pairs) {
        if (tightest == nullptr ||
            pair.perpendicular - pair.od < tightest->perpendicular - tightest->od) {
            tightest = &pair;
        }
        if (landingLayer && pair.layer == landingLayer.value() && pair.landingPair) {
            ++landingPairs;
        }
        const double required = pair.od / std::cos(std::atan(pair.advance / pair.turnLength));
        const bool clear = pair.perpendicular >= pair.od - 0.5e-9;
        if (!clear) {
            ++failures;
            std::cout << "[abt1401] " << pair.layer << ": '" << pair.turnA << "' / '" << pair.turnB
                      << "' spacing " << pair.spacing * 1e3 << " mm, advance "
                      << pair.advance * 1e3 << " mm/rev over " << pair.turnLength * 1e3
                      << " mm, required " << required * 1e3 << " mm, perpendicular "
                      << pair.perpendicular * 1e3 << " mm = "
                      << (pair.od - pair.perpendicular) * 1e9 << " nm inside the "
                      << pair.od * 1e3 << " mm envelope" << std::endl;
        }
        CHECK(clear);
    }
    std::cout << "[abt1401] checked " << pairs.size() << " adjacent pairs, " << landingPairs
              << " steep-landing pairs in '" << landingLayer.value_or("-") << "', " << failures
              << " inside the envelope" << std::endl;
    std::cout << "[abt1401] tightest: " << tightest->layer << " '" << tightest->turnA << "' / '"
              << tightest->turnB << "' spacing " << tightest->spacing * 1e3 << " mm at "
              << tightest->advance * 1e3 << " mm/rev over " << tightest->turnLength * 1e3
              << " mm, perpendicular gap - OD = "
              << (tightest->perpendicular - tightest->od) * 1e9 << " nm (OD "
              << tightest->od * 1e3 << " mm)" << std::endl;
    if (landingLayer) {
        CHECK(landingPairs > 0);
    }
}

}  // namespace

TEST_CASE("Real winding: sibling wraps are spaced for the slope they are drawn with (interleaved flyback, ETD39)",
          "[constructive-model][coil][real-winding][abt1401]") {
    auto coil = wind_real("abt1401_interleaved_flyback_etd39.json");
    require_siblings_clear(coil, "Secondary section 0 layer 1");
}



// ABT #1424: the U layer link is one radial step in the connection plane, drawn from the departure
// station to the conductor's OWN station in the next layer. Sibling parallels' links share that
// plane, so their centrelines must stay one coated OD apart as 2D segments -- and they cannot cross.
// Measured before the fix on the same design: the winder handed layer 1's bundle out top-down in
// parallel order (p0 6.904, p1 6.213, p2 5.521 mm) while layer 0 closed bottom-up (p0 5.551,
// p1 6.231, p2 6.910 mm), and the order law that should keep the bundle's stacking ran only on
// layouts with blocked lead slots -- this one has none. p0's link climbed 1.353 mm over the
// 0.679 mm step and p2's fell 1.389 mm: they crossed (distance 0), and MVB++'s gate refused p0/p1
// at 0.578 mm against the 0.63 mm copper.
TEST_CASE("Real winding: sibling layer links keep one coated diameter apart (interleaved flyback, ETD39)",
          "[constructive-model][coil][real-winding][abt1424]") {
    auto coil = wind_real("abt1401_interleaved_flyback_etd39.json");
    const auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    std::map<std::pair<std::string, int64_t>, std::vector<size_t>> turnsByConductor;
    for (size_t t = 0; t < turns.size(); ++t) {
        turnsByConductor[{turns[t].get_winding(), turns[t].get_parallel()}].push_back(t);
    }
    struct Link {
        int64_t parallel;
        double x0, y0, x1, y1;
        std::string name;
    };
    std::map<std::tuple<std::string, std::string, std::string>, std::vector<Link>> links;
    for (const auto& [conductor, sequence] : turnsByConductor) {
        const auto windingIndex = coil.get_winding_index_by_name(conductor.first);
        const double od = wires[windingIndex].get_maximum_outer_height();
        const double parallels = double(coil.get_number_parallels(windingIndex));
        for (size_t k = 0; k + 1 < sequence.size(); ++k) {
            const auto& a = turns[sequence[k]];
            const auto& b = turns[sequence[k + 1]];
            if (a.get_layer().value() == b.get_layer().value() || a.get_section() != b.get_section()) {
                continue;
            }
            const double dx = std::abs(b.get_coordinates()[0] - a.get_coordinates()[0]);
            const double dy = std::abs(b.get_coordinates()[1] - a.get_coordinates()[1]);
            if (dx <= 1e-12 || dy > dx + 2.0 * od * std::max(1.0, parallels)) {
                continue;   // a dragback, not a radial layer step
            }
            links[{conductor.first, a.get_layer().value(), b.get_layer().value()}].push_back(
                {conductor.second, a.get_coordinates()[0], a.get_coordinates()[1], b.get_coordinates()[0],
                 b.get_coordinates()[1], a.get_name() + " -> " + b.get_name()});
        }
    }
    const auto pointToSegment = [](double px, double py, double ax, double ay, double bx, double by) {
        const double vx = bx - ax, vy = by - ay;
        const double l2 = vx * vx + vy * vy;
        const double u = l2 > 0 ? std::clamp(((px - ax) * vx + (py - ay) * vy) / l2, 0.0, 1.0) : 0.0;
        return std::hypot(px - (ax + u * vx), py - (ay + u * vy));
    };
    const auto side = [](double ox, double oy, double ax, double ay, double bx, double by) {
        return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
    };
    size_t pairs = 0;
    for (const auto& [key, group] : links) {
        const double od = wires[coil.get_winding_index_by_name(std::get<0>(key))].get_maximum_outer_height();
        for (size_t i = 0; i < group.size(); ++i) {
            for (size_t j = i + 1; j < group.size(); ++j) {
                const auto& p = group[i];
                const auto& q = group[j];
                const bool crossing = (side(p.x0, p.y0, p.x1, p.y1, q.x0, q.y0) > 0) !=
                                          (side(p.x0, p.y0, p.x1, p.y1, q.x1, q.y1) > 0) &&
                                      (side(q.x0, q.y0, q.x1, q.y1, p.x0, p.y0) > 0) !=
                                          (side(q.x0, q.y0, q.x1, q.y1, p.x1, p.y1) > 0);
                const double distance =
                    crossing ? 0.0
                             : std::min({pointToSegment(p.x0, p.y0, q.x0, q.y0, q.x1, q.y1),
                                         pointToSegment(p.x1, p.y1, q.x0, q.y0, q.x1, q.y1),
                                         pointToSegment(q.x0, q.y0, p.x0, p.y0, p.x1, p.y1),
                                         pointToSegment(q.x1, q.y1, p.x0, p.y0, p.x1, p.y1)});
                ++pairs;
                std::cout << "[abt1424] " << p.name << " | " << q.name << ": " << distance * 1e3
                          << " mm (OD " << od * 1e3 << " mm)" << std::endl;
                INFO(p.name << " | " << q.name << ": " << distance * 1e3 << " mm against OD "
                            << od * 1e3 << " mm");
                CHECK(distance >= od - 0.5e-9);
            }
        }
    }
    // Secondary 7 t x 3 p crosses from layer 0 to layer 1 once per parallel: three sibling pairs.
    CHECK(pairs >= 3);
}

// ABT #1422: a RADIAL STEP -- a U layer link, or an inter-section return laid as one step at the
// crossing (ABT #1360) -- runs in the connection plane from its departure station to the
// conductor's own station in the next layer, so it is TILTED whenever the two layers' grids differ.
// The sibling beside either end, in that end's layer, leaves (or reaches) the crossing on its wrap
// at slope advance / length, perpendicular to that plane. Measured here the way MVB++ draws it: the
// step as the segment between the two stations, the sibling's wrap as the half-line leaving (or
// arriving at) its station at its own slope, both at their layers' radii; the distance is the
// minimum over the step, sampled densely, of the exact point-to-half-line distance. Every sibling
// pair must clear one coated OD.
//
// PSPS E16 flyback 2p/4p before the fix: the Secondary's section 1 stations were 0.542857 mm apart
// -- ABT #1401's OD / cos(alpha) for the 0.18289 landing slope -- while parallel 0's return fell
// 25.85 um over its 0.534 mm step toward parallel 1's wrap: 0.533396 mm, 604 nm inside the 0.534 mm
// envelope (MVB++ #1295 certified 39 pairs, worst 604.192 nm). On the interleaved flyback the same
// law bit the U link into layer 1 once ABT #1424 restored the bundle order.
namespace {

struct StepPair {
    std::string step, neighbour;
    double distance = 0, od = 0;
};

std::vector<StepPair> measure_radial_steps(OpenMagnetics::Coil& coil) {
    const auto turns = coil.get_turns_description().value();
    const auto layers = coil.get_layers_description().value();
    auto wires = coil.get_wires();
    std::map<std::string, double> layerX;
    std::vector<double> conductionX;
    for (const auto& layer : layers) {
        layerX[layer.get_name()] = layer.get_coordinates()[0];
        if (layer.get_type() == MAS::ElectricalType::CONDUCTION) {
            conductionX.push_back(layer.get_coordinates()[0]);
        }
    }
    std::map<std::pair<std::string, int64_t>, std::vector<size_t>> sequenceOf;
    std::vector<std::pair<size_t, size_t>> placeInSequence(turns.size());
    for (size_t t = 0; t < turns.size(); ++t) {
        auto& sequence = sequenceOf[{turns[t].get_winding(), turns[t].get_parallel()}];
        sequence.push_back(t);
    }
    std::map<size_t, std::pair<std::optional<size_t>, std::optional<size_t>>> prevNext;
    for (const auto& [conductor, sequence] : sequenceOf) {
        for (size_t k = 0; k < sequence.size(); ++k) {
            prevNext[sequence[k]] = {k > 0 ? std::optional<size_t>(sequence[k - 1]) : std::nullopt,
                                     k + 1 < sequence.size() ? std::optional<size_t>(sequence[k + 1])
                                                             : std::nullopt};
        }
    }
    const auto layerOf = [&](size_t t) { return turns[t].get_layer().value(); };
    std::vector<StepPair> pairs;
    for (const auto& [conductor, sequence] : sequenceOf) {
        const auto windingIndex = coil.get_winding_index_by_name(conductor.first);
        const double od = wires[windingIndex].get_maximum_outer_height();
        for (size_t k = 0; k + 1 < sequence.size(); ++k) {
            const size_t a = sequence[k], b = sequence[k + 1];
            if (layerOf(a) == layerOf(b)) {
                continue;
            }
            const double ra = turns[a].get_coordinates()[0], ya = turns[a].get_coordinates()[1];
            const double rb = turns[b].get_coordinates()[0], yb = turns[b].get_coordinates()[1];
            if (std::abs(rb - ra) <= 1e-12 || std::abs(yb - ya) > 0.5 * od) {
                continue;   // not one radial step (MKF lays a vertical stub past half a wire)
            }
            bool adjacent = true;
            for (double x : conductionX) {
                if (x > std::min(ra, rb) + 1e-12 && x < std::max(ra, rb) - 1e-12) {
                    adjacent = false;
                }
            }
            if (!adjacent) {
                continue;   // band-routed over another layer, not one step
            }
            for (size_t end : {a, b}) {
                std::vector<size_t> inLayer;
                for (size_t t = 0; t < turns.size(); ++t) {
                    if (layerOf(t) == layerOf(end)) inLayer.push_back(t);
                }
                std::sort(inLayer.begin(), inLayer.end(), [&](size_t p, size_t q) {
                    return turns[p].get_coordinates()[1] < turns[q].get_coordinates()[1];
                });
                const auto at = std::find(inLayer.begin(), inLayer.end(), end) - inLayer.begin();
                for (long side : {-1L, 1L}) {
                    const long j = at + side;
                    if (j < 0 || j >= long(inLayer.size())) continue;
                    const size_t n = inLayer[size_t(j)];
                    if (turns[n].get_winding() != conductor.first || turns[n].get_parallel() == conductor.second) {
                        continue;
                    }
                    const double rn = turns[n].get_coordinates()[0], s = turns[n].get_coordinates()[1];
                    std::vector<std::pair<double, bool>> pieces;   // (slope, leaving)
                    const auto [prev, next] = prevNext.at(n);
                    if (next && layerOf(*next) == layerOf(n)) {
                        REQUIRE(turns[*next].get_length() > 0);
                        pieces.push_back({(turns[*next].get_coordinates()[1] - s) / turns[*next].get_length(), true});
                    }
                    if (prev && layerOf(*prev) == layerOf(n)) {
                        REQUIRE(turns[n].get_length() > 0);
                        pieces.push_back({(s - turns[*prev].get_coordinates()[1]) / turns[n].get_length(), false});
                    }
                    for (const auto& [m, leaving] : pieces) {
                        double best = std::numeric_limits<double>::max();
                        for (int i = 0; i <= 20000; ++i) {
                            const double w = i / 20000.0;
                            const double y = ya + w * (yb - ya), r = ra + w * (rb - ra);
                            const double e = y - s;
                            double x = m * e / (1 + m * m);
                            if (leaving ? x < 0 : x > 0) x = 0;
                            best = std::min(best, std::hypot(x, e - m * x, r - rn));
                        }
                        pairs.push_back({turns[a].get_name() + " -> " + turns[b].get_name(),
                                         turns[n].get_name() + (leaving ? " (leaving)" : " (arriving)"), best, od});
                    }
                }
            }
        }
    }
    return pairs;
}

void require_steps_clear(OpenMagnetics::Coil& coil) {
    const auto pairs = measure_radial_steps(coil);
    REQUIRE(pairs.size() > 0);
    const StepPair* worst = nullptr;
    for (const auto& pair : pairs) {
        if (worst == nullptr || pair.distance - pair.od < worst->distance - worst->od) worst = &pair;
        INFO(pair.step << " vs " << pair.neighbour << ": " << pair.distance * 1e3 << " mm against OD "
                       << pair.od * 1e3 << " mm");
        CHECK(pair.distance >= pair.od - 0.5e-9);
    }
    std::cout << "[abt1422] " << pairs.size() << " step/sibling-wrap pairs; tightest " << worst->step
              << " vs " << worst->neighbour << ": " << worst->distance * 1e3 << " mm, "
              << (worst->distance - worst->od) * 1e9 << " nm off the " << worst->od * 1e3 << " mm envelope"
              << std::endl;
}

}  // namespace

TEST_CASE("Real winding: a tilted radial step clears the sibling's wrap (PSPS E16 flyback, rect column)",
          "[constructive-model][coil][real-winding][abt1422]") {
    auto coil = wind_real("abt1422_psps_e16_flyback_2p4p.json");
    require_steps_clear(coil);
    // ABT #1401's wrap law must still hold where the step law widened the stations.
    require_siblings_clear(coil, std::nullopt);
}

TEST_CASE("Real winding: a tilted radial step clears the sibling's wrap (interleaved flyback, ETD39)",
          "[constructive-model][coil][real-winding][abt1422]") {
    auto coil = wind_real("abt1401_interleaved_flyback_etd39.json");
    require_steps_clear(coil);
    require_siblings_clear(coil, "Secondary section 0 layer 1");
}


// ABT #969 (Alf, 2026-10-24: "fix them exactly, and geometrically exactly, never add a gap just for
// adding it without a geometric reason"). The laws above each judge ONE local model of the drawn
// conductor -- #1401 the wraps unrolled flat, #1422 a radial step against the tangent half-line of a
// sibling's wrap, #1424 two links in the connection plane. What a consumer actually sweeps is the
// whole CURVE: every revolution a pitch-true helix along the turn's own path (a circle round a round
// column; on a rectangular column the racetrack whose corners are the bend radius MKF charges the
// turn for), climbing linearly in arc length from its station to the next one, and every radial
// layer link the straight segment between its two stations at the crossing. Here those curves are
// built from the published layout alone and the TRUE minimum distance between every two pieces of
// copper is measured: dense sampling to find each candidate pair, then a local continuous
// minimisation on the exact parametrisations, so no chord or tangent approximation is left in the
// number. Two coated wires may touch, never interpenetrate: every pair must stay at least the mean
// of the two coated diameters apart (less the half nanometre of the emission grid).
namespace {

struct TrueCurvePiece {
    std::string winding;
    int64_t parallel = 0;
    std::string label;
    bool link = false;
    size_t ordinal = 0;           // position in the conductor's sequence
    double od = 0;
    // Exact parametrisation over [0, 1].
    std::function<std::array<double, 3>(double)> at;
    double length = 0;            // drawn length, for the sampling density
};

struct TrueCurvePair {
    std::string a, b;
    bool aLink = false, bLink = false;
    bool siblings = false;
    double distance = 0, envelope = 0;
};

// One racetrack round a rectangular column: half sides a (x) and c (z), corner radius b, travelled
// from the crossing at the -Z face's midpoint towards -X (Coil.h, TangentDeparture's frame).
std::array<double, 2> racetrack_point(double a, double c, double b, double u) {
    const double sx = a - b, sz = c - b;
    const double pi = std::numbers::pi;
    const double quarter = pi / 2 * b;
    const std::array<double, 9> lengths = {sx, quarter, 2 * sz, quarter, 2 * sx, quarter, 2 * sz, quarter, sx};
    size_t k = 0;
    while (k + 1 < lengths.size() && u > lengths[k]) {
        u -= lengths[k];
        ++k;
    }
    const auto arc = [&](double cx, double cz, double from) {
        const double angle = from - u / b;   // clockwise seen from +Y with X right, Z up
        return std::array<double, 2>{cx + b * std::cos(angle), cz + b * std::sin(angle)};
    };
    switch (k) {
        case 0: return {-u, -c};
        case 1: return arc(-sx, -sz, -pi / 2);
        case 2: return {-a, -sz + u};
        case 3: return arc(-sx, sz, -pi);
        case 4: return {-sx + u, c};
        case 5: return arc(sx, sz, pi / 2);
        case 6: return {a, sz - u};
        case 7: return arc(sx, -sz, 0);
        default: return {sx - u, -c};
    }
}

std::vector<TrueCurvePiece> build_true_curve_pieces(OpenMagnetics::Coil& coil, size_t& skippedRoutes) {
    const auto turns = coil.get_turns_description().value();
    auto wires = coil.get_wires();
    std::map<std::pair<std::string, int64_t>, std::vector<size_t>> sequenceOf;
    for (size_t t = 0; t < turns.size(); ++t) {
        sequenceOf[{turns[t].get_winding(), turns[t].get_parallel()}].push_back(t);
    }
    std::map<std::string, OpenMagnetics::WoundColumnFrame> frameOf;
    const auto frame = [&](const MAS::Turn& turn) -> const OpenMagnetics::WoundColumnFrame& {
        const std::string section = turn.get_section().value();
        if (!frameOf.contains(section)) {
            frameOf.insert({section, coil.get_wound_column_frame_for_section(section)});
        }
        return frameOf.at(section);
    };
    // The crossing point of a station and, for a wrap, the whole revolution through it.
    const auto crossing = [&](const MAS::Turn& turn) -> std::array<double, 3> {
        const auto& f = frame(turn);
        const double radius = f.axisX == 0 ? std::abs(turn.get_coordinates()[0])
                                           : std::abs(std::abs(turn.get_coordinates()[0]) - f.axisX);
        const double y = turn.get_coordinates()[1];
        if (f.shape == MAS::ColumnShape::ROUND) {
            return {radius, y, 0.0};
        }
        if (f.shape == MAS::ColumnShape::RECTANGULAR || f.shape == MAS::ColumnShape::IRREGULAR) {
            return {0.0, y, -(f.columnDepth + radius - f.columnWidth)};
        }
        throw std::runtime_error("true-curve census: column shape not modelled");
    };
    skippedRoutes = 0;
    std::vector<TrueCurvePiece> pieces;
    for (const auto& [conductor, sequence] : sequenceOf) {
        const auto windingIndex = coil.get_winding_index_by_name(conductor.first);
        const double od = wires[windingIndex].get_maximum_outer_height();
        for (size_t k = 0; k + 1 < sequence.size(); ++k) {
            const auto& from = turns[sequence[k]];
            const auto& to = turns[sequence[k + 1]];
            TrueCurvePiece piece;
            piece.winding = conductor.first;
            piece.parallel = conductor.second;
            piece.ordinal = k;
            piece.od = od;
            const double ya = from.get_coordinates()[1], yb = to.get_coordinates()[1];
            if (from.get_layer() && to.get_layer() && from.get_layer().value() == to.get_layer().value()) {
                const auto& f = frame(from);
                const double radius = f.axisX == 0 ? std::abs(from.get_coordinates()[0])
                                                   : std::abs(std::abs(from.get_coordinates()[0]) - f.axisX);
                piece.label = conductor.first + " p" + std::to_string(conductor.second) + " '" + from.get_name() +
                              "' -> '" + to.get_name() + "' (wrap)";
                if (f.shape == MAS::ColumnShape::ROUND) {
                    piece.length = 2 * std::numbers::pi * radius;
                    piece.at = [radius, ya, yb](double t) {
                        const double phi = 2 * std::numbers::pi * t;
                        return std::array<double, 3>{radius * std::cos(phi), ya + (yb - ya) * t, radius * std::sin(phi)};
                    };
                }
                else if (f.shape == MAS::ColumnShape::RECTANGULAR || f.shape == MAS::ColumnShape::IRREGULAR) {
                    const double standoff = radius - f.columnWidth;
                    const double a = f.columnWidth + standoff, c = f.columnDepth + standoff;
                    const double b = coil.get_turn_bend_radius(to.get_name())
                                         ? coil.get_turn_bend_radius(to.get_name()).value()
                                         : f.cornerRadius + standoff;
                    piece.length = 4 * a + 4 * c + (2 * std::numbers::pi - 8) * b;
                    const double length = piece.length;
                    piece.at = [a, c, b, length, ya, yb](double t) {
                        const auto xz = racetrack_point(a, c, b, t * length);
                        return std::array<double, 3>{xz[0], ya + (yb - ya) * t, xz[1]};
                    };
                }
                else {
                    ++skippedRoutes;
                    continue;
                }
            }
            else {
                const double ra = std::abs(from.get_coordinates()[0]), rb = std::abs(to.get_coordinates()[0]);
                // The radial step MKF lays with no vertical stub (Coil::measure_radial_step_clearances).
                if (std::abs(rb - ra) <= 1e-12 || std::abs(yb - ya) > 0.5 * od || from.get_section() != to.get_section()) {
                    ++skippedRoutes;   // a dragback, an L route or a band route: drawn elsewhere
                    continue;
                }
                const auto p = crossing(from), q = crossing(to);
                piece.link = true;
                piece.label = conductor.first + " p" + std::to_string(conductor.second) + " '" + from.get_name() +
                              "' -> '" + to.get_name() + "' (layer link)";
                piece.length = std::hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]);
                piece.at = [p, q](double t) {
                    return std::array<double, 3>{p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t,
                                                 p[2] + (q[2] - p[2]) * t};
                };
            }
            pieces.push_back(piece);
        }
    }
    return pieces;
}

double true_curve_distance(const std::array<double, 3>& p, const std::array<double, 3>& q) {
    return std::hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
}

// Local continuous minimisation of |A(u) - B(v)| from a sampled start: alternating golden-section
// searches over a shrinking bracket. The pieces are smooth (C1 at a racetrack's straight/arc joints),
// so this converges on the true local minimum, not the sampled one.
double refine_true_curve_minimum(const TrueCurvePiece& A, const TrueCurvePiece& B, double u, double v, double h) {
    const double golden = (std::sqrt(5.0) - 1) / 2;
    const auto line = [&](bool overU, double centre, double other, double half) {
        double low = std::max(0.0, centre - half), high = std::min(1.0, centre + half);
        const auto f = [&](double s) {
            return overU ? true_curve_distance(A.at(s), B.at(other)) : true_curve_distance(A.at(other), B.at(s));
        };
        double c = high - golden * (high - low), d = low + golden * (high - low);
        double fc = f(c), fd = f(d);
        for (int it = 0; it < 200 && high - low > 1e-16; ++it) {
            if (fc < fd) { high = d; d = c; fd = fc; c = high - golden * (high - low); fc = f(c); }
            else { low = c; c = d; fc = fd; d = low + golden * (high - low); fd = f(d); }
        }
        const double best = fc < fd ? c : d;
        return f(centre) < std::min(fc, fd) ? centre : best;
    };
    double hu = h / std::max(A.length, 1e-12), hv = h / std::max(B.length, 1e-12);
    for (int pass = 0; pass < 400; ++pass) {
        const double nu = line(true, u, v, 2 * hu);
        const double nv = line(false, v, nu, 2 * hv);
        const bool settled = std::abs(nu - u) < 1e-15 && std::abs(nv - v) < 1e-15;
        u = nu;
        v = nv;
        if (settled) break;
        if (pass % 8 == 7) { hu *= 0.5; hv *= 0.5; }
    }
    return true_curve_distance(A.at(u), B.at(v));
}

std::vector<TrueCurvePair> measure_true_curve_pairs(OpenMagnetics::Coil& coil, size_t& skippedRoutes) {
    const auto pieces = build_true_curve_pieces(coil, skippedRoutes);
    const double step = 20e-6;   // sampling, well inside one coated diameter
    struct Sample { size_t piece; double t; std::array<double, 3> p; };
    std::vector<std::vector<Sample>> samples(pieces.size());
    double maxOd = 0;
    for (size_t i = 0; i < pieces.size(); ++i) {
        maxOd = std::max(maxOd, pieces[i].od);
        const size_t n = std::max<size_t>(2, size_t(std::ceil(pieces[i].length / step)) + 1);
        for (size_t s = 0; s < n; ++s) {
            const double t = double(s) / double(n - 1);
            samples[i].push_back({i, t, pieces[i].at(t)});
        }
    }
    const double cell = 1.5 * maxOd;
    std::map<std::array<int64_t, 3>, std::vector<std::pair<size_t, size_t>>> grid;
    const auto key = [&](const std::array<double, 3>& p) {
        return std::array<int64_t, 3>{int64_t(std::floor(p[0] / cell)), int64_t(std::floor(p[1] / cell)),
                                      int64_t(std::floor(p[2] / cell))};
    };
    for (size_t i = 0; i < pieces.size(); ++i) {
        for (size_t s = 0; s < samples[i].size(); ++s) {
            grid[key(samples[i][s].p)].push_back({i, s});
        }
    }
    // Best sampled approach per unordered piece pair.
    std::map<std::pair<size_t, size_t>, std::tuple<double, double, double>> best;
    for (size_t i = 0; i < pieces.size(); ++i) {
        for (const auto& sample : samples[i]) {
            const auto k = key(sample.p);
            for (int64_t dx = -1; dx <= 1; ++dx) for (int64_t dy = -1; dy <= 1; ++dy) for (int64_t dz = -1; dz <= 1; ++dz) {
                auto found = grid.find({k[0] + dx, k[1] + dy, k[2] + dz});
                if (found == grid.end()) continue;
                for (const auto& [j, s] : found->second) {
                    if (j <= i) continue;
                    const auto& A = pieces[i];
                    const auto& B = pieces[j];
                    if (A.winding == B.winding && A.parallel == B.parallel &&
                        (A.ordinal + 1 == B.ordinal || B.ordinal + 1 == A.ordinal)) {
                        continue;   // consecutive pieces of one conductor are one wire
                    }
                    const double d = true_curve_distance(sample.p, samples[j][s].p);
                    const double envelope = (A.od + B.od) / 2;
                    if (d > envelope + 4 * step) continue;
                    auto it = best.find({i, j});
                    if (it == best.end() || d < std::get<0>(it->second)) {
                        best[{i, j}] = {d, sample.t, samples[j][s].t};
                    }
                }
            }
        }
    }
    std::vector<TrueCurvePair> pairs;
    for (const auto& [ij, start] : best) {
        const auto& A = pieces[ij.first];
        const auto& B = pieces[ij.second];
        TrueCurvePair pair;
        pair.a = A.label;
        pair.b = B.label;
        pair.aLink = A.link;
        pair.bLink = B.link;
        pair.siblings = A.winding == B.winding && A.parallel != B.parallel;
        pair.envelope = (A.od + B.od) / 2;
        pair.distance = refine_true_curve_minimum(A, B, std::get<1>(start), std::get<2>(start), 4 * step);
        pairs.push_back(pair);
    }
    return pairs;
}

OpenMagnetics::Coil wind_real_mas_example(const std::string& example) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    const auto path = std::filesystem::path{std::source_location::current().file_name()}
                          .parent_path()
                          .append("..")
                          .append("MAS")
                          .append("examples")
                          .append(example);
    std::ifstream file(path);
    REQUIRE(file.good());
    const auto masJson = json::parse(file);
    const auto& magneticJson = masJson.at("magnetic");
    settings.set_coil_use_real_winding_geometry(true);
    settings.set_coil_lead_bend_radius_factor(1.05);
    settings.set_coil_lead_minimum_bend_radius(std::nullopt);
    OpenMagnetics::Core core(magneticJson.at("core"));
    OpenMagnetics::Coil coil(magneticJson.at("coil"), false);
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    auto enriched = OpenMagnetics::magnetic_autocomplete(magnetic, json{});
    settings.reset();
    REQUIRE(enriched.get_coil().is_real_winding_blocking_applied());
    return enriched.get_coil();
}

// Every pair closer than the coated envelope (less the half nanometre of the emission grid), worst
// first, for the failure message.
std::string pairs_inside_the_envelope(const std::vector<TrueCurvePair>& pairs) {
    auto sorted = pairs;
    std::sort(sorted.begin(), sorted.end(), [](const auto& p, const auto& q) {
        return p.distance - p.envelope < q.distance - q.envelope;
    });
    std::ostringstream report;
    report << std::setprecision(12);
    for (const auto& pair : sorted) {
        if (!(pair.distance < pair.envelope - 0.5e-9)) break;
        report << pair.a << " | " << pair.b << " | " << pair.distance * 1e3 << " mm against the "
               << pair.envelope * 1e3 << " mm envelope: " << (pair.distance - pair.envelope) * 1e9 << " nm\n";
    }
    return report.str();
}

}  // namespace

// Before the true-curve clearance the tangent law left 14_dab's final landings 306 nm (Primary) and
// 172 nm (Secondary) inside the sibling's helix, and 23_interleaved_llc's 0.84 nm.
TEST_CASE("Real winding: no layer link or wrap of the MAS examples cuts into another conductor's true curve",
          "[constructive-model][coil][real-winding][abt-969]") {
    for (const std::string example : {"02_flyback_efd25_3c95.json", "06_llc_xfmr_eq4128_3c97.json",
                                      "23_interleaved_llc_pq3530_n97.json", "21_interleaved_flyback_etd39_n87.json",
                                      "14_dab_xfmr_pm8770_n97.json"}) {
        INFO(example);
        auto coil = wind_real_mas_example(example);
        size_t skipped = 0;
        const auto pairs = measure_true_curve_pairs(coil, skipped);
        size_t linkPairs = 0;
        for (const auto& pair : pairs) {
            if ((pair.aLink || pair.bLink) && pair.siblings) ++linkPairs;
        }
        // Not vacuous: there are pieces within reach on every design, and the two designs the tangent
        // law failed lay sibling layer links beside a sibling's wrap.
        REQUIRE(!pairs.empty());
        if (example == "14_dab_xfmr_pm8770_n97.json" || example == "23_interleaved_llc_pq3530_n97.json") {
            REQUIRE(linkPairs > 0);
        }
        const std::string inside = pairs_inside_the_envelope(pairs);
        INFO(inside);
        REQUIRE(inside.empty());
    }
}
