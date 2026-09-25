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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
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

