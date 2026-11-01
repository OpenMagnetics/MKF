// =============================================================================
// TestCoilAdviserCharacterisation.cpp
// =============================================================================
// Characterisation tests for CoilAdviser::get_advised_coil.
//
// PURPOSE
//   Lock the CURRENT top-N coil configurations produced for a representative
//   transformer fixture. get_advised_coil returns vector<Mas> without scores
//   (scoring is internal), so this file snapshots the STRUCTURAL output:
//     - number of returned configurations
//     - per-winding wire name (joined as "wire0 || wire1 || ...")
//     - section / layer / turn counts on the top result
//   Any refactor of the 1218-line CoilAdviser that touches the pipeline must
//   reproduce these structural facts exactly.
//
// FIXTURE
//   ETD 59 / 3C91, 0.003 m ground gap, 82:5 turns, IEC 60664-1 / 61558 BASIC
//   insulation, 400 V mains, 25 °C. Its core, turns and insulation are those
//   of TestCoilAdviser.cpp: Test_CoilAdviser_Insulation_No_Margin; since
//   2026-11-01 its operating point is a real 150 W CCM flyback (360-440 V in,
//   20 V out, 175.59 kHz) designed by MKF's converter model, see
//   make_flyback_fixture.
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true (or leave a snapshot empty), run, copy
//   BASELINE lines from stderr into the kSnapshots tables.
//
// BENCHMARKS  (tag: [!benchmark])
//   Time get_advised_coil end-to-end for the fixture. Catch2 defaults: use
//      --benchmark-samples 3 --benchmark-warmup-time 0
// =============================================================================

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include "advisers/CoilAdviser.h"
#include "converter_models/KirchhoffBridge.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "support/Settings.h"

#include "TestingUtils.h"

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

constexpr bool kRegenerateBaselines = false;

// --- Fixture builder --------------------------------------------------------

OpenMagnetics::Mas make_flyback_fixture() {
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<int64_t> numberTurns = {82, 5};

    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 59", gapping, numberTurns, /*stacks*/ 1, "3C91");

    // Operating point replaced with Alf's approval (2026-11-01, ABT #1394): the sinusoidal 20 App point drove
    // 115 A rms through the 5-turn secondary (~420 W), which no coil on this core can carry under the
    // temperature limit. This is a real flyback operating point, designed by MKF's converter model (Kirchhoff's
    // Flyback through KirchhoffBridge, not hand-made waveforms): Vin 360/400/440 V (operating at 400 V),
    // Vout 20 V, 150 W, 175.59 kHz, CCM, turns ratio 82:5 and Lm 1.4927 mH pinned (the inductance of this
    // ETD 59 / 3C91 / 3 mm gap at 82 turns), which gives D = 0.45 at 400 V.
    json spec;
    spec["designRequirements"]["efficiency"] = 1.0;
    spec["designRequirements"]["inputVoltage"] = {{"minimum", 360}, {"nominal", 400}, {"maximum", 440}};
    spec["designRequirements"]["switchingFrequency"]["nominal"] = 175590;
    spec["designRequirements"]["outputs"] = json::array({{{"name", "out"}, {"voltage", {{"nominal", 20}}}}});
    spec["designRequirements"]["turnsRatios"] = json::array({{{"nominal", double(numberTurns[0]) / numberTurns[1]}}});
    spec["designRequirements"]["magnetizingInductance"] = {{"nominal", 1.4927e-3}};
    spec["operatingPoints"] = json::array({{{"inputVoltage", 400}, {"outputs", json::array({{{"power", 150}}})}}});
    spec["config"]["mode"] = "ccm";
    auto inputs = OpenMagnetics::design_requirements_from_tas(OpenMagnetics::design_converter_tas("flyback", spec));

    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_nominal(400);
    auto standards = std::vector<InsulationStandards>{
        InsulationStandards::IEC_606641, InsulationStandards::IEC_623681};
    auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(
        altitude, Cti::GROUP_I, IsolationClass::BASIC, mainSupplyVoltage,
        OvervoltageCategory::IV, PollutionDegree::PD1, standards);
    inputs.get_mutable_design_requirements().set_insulation(insulationRequirements);

    OpenMagnetics::Mas mas;
    inputs.process();
    mas.set_inputs(inputs);
    mas.set_magnetic(magnetic);
    return mas;
}

// --- Snapshot type ----------------------------------------------------------

struct CoilStructural {
    std::string coreShape;
    std::string coreMaterial;
    std::string turnsPerWinding;  // "turns0xparallels0 || turns1xparallels1 || ..."
    std::string wirePerWinding;  // "wire0 || wire1 || ..."
    size_t numberSections;
    size_t numberLayers;
    int64_t totalTurns;
};

std::string join_wire_names(OpenMagnetics::Mas& mas) {
    std::ostringstream oss;
    auto wires = mas.get_mutable_magnetic().get_wires();
    for (size_t i = 0; i < wires.size(); ++i) {
        if (i) oss << " || ";
        oss << wires[i].get_name().value_or("<unnamed>");
    }
    return oss.str();
}

CoilStructural extract_structural(OpenMagnetics::Mas& mas) {
    CoilStructural out;
    out.coreShape = mas.get_magnetic().get_core().get_shape_name();
    out.coreMaterial = mas.get_mutable_magnetic().get_mutable_core().resolve_material().get_name();
    out.wirePerWinding = join_wire_names(mas);
    auto& coil = mas.get_mutable_magnetic().get_mutable_coil();
    {
        std::ostringstream oss;
        const auto& windings = coil.get_functional_description();
        for (size_t i = 0; i < windings.size(); ++i) {
            if (i) oss << " || ";
            oss << windings[i].get_number_turns() << "x" << windings[i].get_number_parallels();
        }
        out.turnsPerWinding = oss.str();
    }
    out.numberSections = coil.get_sections_description().has_value()
                       ? coil.get_sections_description().value().size() : 0;
    out.numberLayers = coil.get_layers_description().has_value()
                     ? coil.get_layers_description().value().size() : 0;
    int64_t total = 0;
    for (auto& w : coil.get_functional_description()) {
        total += w.get_number_turns();
    }
    out.totalTurns = total;
    return out;
}

void check_structural_topN(const std::string& label,
                            std::vector<OpenMagnetics::Mas>& got,
                            const std::vector<CoilStructural>& expectedTop) {
    if (kRegenerateBaselines || expectedTop.empty()) {
        std::cerr << "\nBASELINE COIL_TOPN " << label << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            auto s = extract_structural(got[i]);
            std::cerr << "  [" << i << "] core=\"" << s.coreShape << " / " << s.coreMaterial
                      << "\" turns=\"" << s.turnsPerWinding << "\" wires=\"" << s.wirePerWinding
                      << "\" sections=" << s.numberSections
                      << " layers=" << s.numberLayers
                      << " totalTurns=" << s.totalTurns << "\n";
        }
        REQUIRE(got.size() >= expectedTop.size());
        return;
    }
    INFO("scenario=" << label);
    REQUIRE(got.size() >= expectedTop.size());
    for (size_t i = 0; i < expectedTop.size(); ++i) {
        auto s = extract_structural(got[i]);
        INFO("top[" << i << "] want core=\"" << expectedTop[i].coreShape << " / " << expectedTop[i].coreMaterial
                    << "\" turns=\"" << expectedTop[i].turnsPerWinding << "\" wires=\"" << expectedTop[i].wirePerWinding
                    << "\" got core=\"" << s.coreShape << " / " << s.coreMaterial << "\" turns=\"" << s.turnsPerWinding
                    << "\" wires=\"" << s.wirePerWinding
                    << "\" sections=" << s.numberSections
                    << " layers=" << s.numberLayers
                    << " totalTurns=" << s.totalTurns);
        REQUIRE(s.coreShape == expectedTop[i].coreShape);
        REQUIRE(s.coreMaterial == expectedTop[i].coreMaterial);
        REQUIRE(s.turnsPerWinding == expectedTop[i].turnsPerWinding);
        REQUIRE(s.wirePerWinding == expectedTop[i].wirePerWinding);
        REQUIRE(s.numberSections == expectedTop[i].numberSections);
        REQUIRE(s.numberLayers == expectedTop[i].numberLayers);
        REQUIRE(s.totalTurns == expectedTop[i].totalTurns);
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS. Leave empty to auto-harvest.
// ----------------------------------------------------------------------------
// Re-pinned 2026-11-01 with Alf's approval, when the fixture's sinusoidal operating point (115 A rms on the
// secondary, ~420 W; no coil survived the coil adviser's temperature gate, ABT #1394) was replaced by the real
// 150 W CCM flyback point of make_flyback_fixture. Before pinning, every pick was checked against the gates:
// B peak 0.0578 T on 3C91 (L 1.4927 mH), hot spots 37.9 / 36.5 / 36.4 C against the 130 C limit, losses
// 2.13 / 3.46 / 3.32 W.
const std::vector<CoilStructural> kTopFlyback = {
    // core shape, material, turns x parallels per winding, wires, sections, layers, total turns
    {"ETD 59", "3C91", "82x1 || 5x2", "Litz TXXL19/40TXXX-2(MWXX) || Round T10A01PXXX-3",      4, 35, 87},
    {"ETD 59", "3C91", "82x2 || 5x1", "Round T29A01PXXX-3 || Synthesized Litz 648x0.041mm",   8, 70, 87},
    {"ETD 59", "3C91", "82x2 || 5x1", "Round T29A01PXXX-3 || Synthesized Litz 1458x0.028mm",  8, 70, 87},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOT
// =============================================================================

TEST_CASE("CoilAdviser flyback 2-winding top-3 structural snapshot",
          "[adviser][coil-adviser][characterisation][heavy][flyback]") {
    settings.reset();
    settings.set_coil_allow_margin_tape(false);
    settings.set_coil_allow_insulated_wire(true);
    settings.set_coil_try_rewind(false);
    settings.set_coil_adviser_maximum_number_wires(1000);

    auto mas = make_flyback_fixture();

    CoilAdviser adviser;
    auto results = adviser.get_advised_coil(mas, 3);
    check_structural_topN("FLYBACK_82_5_ETD59_3C91", results, kTopFlyback);
    // Three asked, three returned at the 150 W operating point.
    REQUIRE(results.size() == 3);

    settings.reset();
}

// =============================================================================
// BENCHMARK  (opt-in via [!benchmark])
// =============================================================================
//
// Use:  ./MKF_tests "[benchmark-coil-adviser]" --benchmark-samples 3 --benchmark-warmup-time 0
//

TEST_CASE("Benchmark CoilAdviser flyback 2-winding (top-3)",
          "[!benchmark][benchmark-coil-adviser]") {
    settings.reset();
    settings.set_coil_allow_margin_tape(false);
    settings.set_coil_allow_insulated_wire(true);
    settings.set_coil_try_rewind(false);
    settings.set_coil_adviser_maximum_number_wires(1000);

    auto mas = make_flyback_fixture();

    BENCHMARK("get_advised_coil flyback 2W top-3") {
        CoilAdviser adviser;
        return adviser.get_advised_coil(mas, 3);
    };

    settings.reset();
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                | mean (3 samples) | notes
// -----------+-------------------------+------------------+----------
// 2026-05-19 | flyback 2W ETD 59 top-3 | TBD              | initial
//
// =============================================================================
