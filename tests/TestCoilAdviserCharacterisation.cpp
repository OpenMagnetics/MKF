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
//   ETD 59 / 3C91, 0.003 m ground gap, 82:5 turns, 175.59 kHz sinusoidal,
//   20 Vpp, IEC 60664-1 / 61558 BASIC insulation, 400 V mains, 25 °C.
//   This is the same fixture as TestCoilAdviser.cpp:
//   Test_CoilAdviser_Insulation_No_Margin (the canonical 2-winding flyback
//   test in the existing suite).
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
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "support/Settings.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;

namespace {

constexpr bool kRegenerateBaselines = false;

// --- Fixture builder --------------------------------------------------------

OpenMagnetics::Mas make_flyback_fixture() {
    auto gapping = OpenMagneticsTesting::get_ground_gap(0.003);
    std::vector<int64_t> numberTurns = {82, 5};
    std::vector<double> turnsRatios;
    for (size_t i = 1; i < numberTurns.size(); ++i) {
        turnsRatios.push_back(double(numberTurns[0]) / numberTurns[i]);
    }

    auto magnetic = OpenMagneticsTesting::get_quick_magnetic(
        "ETD 59", gapping, numberTurns, /*stacks*/ 1, "3C91");

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        /*frequency*/         175590,
        /*magInductance*/     10e-6,
        /*temperature*/       25,
        /*waveShape*/         WaveformLabel::SINUSOIDAL,
        /*peakToPeak*/        20,
        /*dutyCycle*/         0.5,
        /*dcCurrent*/         0,
        /*turnsRatios*/       turnsRatios);

    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_nominal(400);
    auto standards = std::vector<InsulationStandards>{
        InsulationStandards::IEC_606641, InsulationStandards::IEC_623681};
    auto insulationRequirements = OpenMagneticsTesting::get_quick_insulation_requirements(
        altitude, Cti::GROUP_I, InsulationType::BASIC, mainSupplyVoltage,
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
    out.wirePerWinding = join_wire_names(mas);
    auto& coil = mas.get_mutable_magnetic().get_mutable_coil();
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
            std::cerr << "  [" << i << "] wires=\"" << s.wirePerWinding
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
        INFO("top[" << i << "] want wires=\"" << expectedTop[i].wirePerWinding
                    << "\" got wires=\"" << s.wirePerWinding
                    << "\" sections=" << s.numberSections
                    << " layers=" << s.numberLayers
                    << " totalTurns=" << s.totalTurns);
        REQUIRE(s.wirePerWinding == expectedTop[i].wirePerWinding);
        REQUIRE(s.numberSections == expectedTop[i].numberSections);
        REQUIRE(s.numberLayers == expectedTop[i].numberLayers);
        REQUIRE(s.totalTurns == expectedTop[i].totalTurns);
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS — captured 2026-05-19. Leave empty to auto-harvest.
// ----------------------------------------------------------------------------
// Note: requested top-3 but pipeline only produces 2 candidates for this
// fixture — locked as the current behaviour. If a future refactor unlocks
// more candidates this snapshot will fail loudly and need regeneration.
// Refreshed 2026-06-16 (ABT #10 / ABT #3): top[1] primary re-ranked from
// Round T20A01TXXX-1 to Round T25A01TXXX-1.5 (12 layers, was 10). Both are
// insulated 3-layer reinforced wires (bdv 4500 → 6000) — the expected wire
// class for this BASIC isolated flyback — so this is a benign within-class
// re-rank from the wire-scoring shift, verified pre-existing in ABT #3, not a
// regression. top[0] unchanged.
const std::vector<CoilStructural> kTopFlyback = {
    // wires                                                              sections layers totalTurns
    {"Round T25A01PXXX-1.5 || Litz SXXL825/44FX-3(MWXX)",                  8,       12,    87},
    {"Round T25A01TXXX-1.5 || Litz SXXL825/44FX-3(MWXX)",                  8,       12,    87},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOT
// =============================================================================

TEST_CASE("CoilAdviser flyback 2-winding top-3 structural snapshot",
          "[adviser][coil-adviser][characterisation][flyback]") {
    settings.reset();
    settings.set_coil_allow_margin_tape(false);
    settings.set_coil_allow_insulated_wire(true);
    settings.set_coil_try_rewind(false);
    settings.set_coil_adviser_maximum_number_wires(1000);

    auto mas = make_flyback_fixture();

    CoilAdviser adviser;
    auto results = adviser.get_advised_coil(mas, 3);
    check_structural_topN("FLYBACK_82_5_ETD59_3C91", results, kTopFlyback);
    // Returns 2 (not 3) by design: wire-advancement yields 6 candidate coils, but
    // EFFECTIVE_CURRENT_DENSITY is a strict filter (CoilAdviser.h:_defaultCustomMagneticFilterFlow),
    // and only 2 of the 6 candidates pass it for this fixture (ETD 59 / 82:5 turns,
    // ~175 kHz flyback). The 4 rejected candidates are under-paralleled wires whose
    // per-strand effective current density exceeds defaults.maximumEffectiveCurrentDensity.
    // See MagneticFilter.cpp:1922 (MagneticFilterEffectiveCurrentDensity). This is correct
    // physical filtering, not a bug — locking the count to detect regressions in either
    // direction.
    REQUIRE(results.size() == 2);

    settings.reset();
}

// =============================================================================
// BENCHMARK  (opt-in via [!benchmark])
// =============================================================================
//
// Use:  ./MKF_tests "[benchmark-coil-adviser]" --benchmark-samples 3 --benchmark-warmup-time 0
//

TEST_CASE("Benchmark CoilAdviser flyback 2-winding (top-3)",
          "[adviser][coil-adviser][!benchmark][benchmark-coil-adviser]") {
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
