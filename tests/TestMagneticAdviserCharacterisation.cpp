// =============================================================================
// TestMagneticAdviserCharacterisation.cpp
// =============================================================================
// End-to-end characterisation test for MagneticAdviser::get_advised_magnetic.
//
// PURPOSE
//   Lock the CURRENT top-N (magnetic-reference, score) output for a
//   representative 3-winding fixture. This exercises the WHOLE pipeline:
//   CoreAdviser → CoilAdviser → MagneticFilter scoring → final ranking.
//
//   Because every internal change (filter weights, core dedupe, wire scoring,
//   etc.) flows through here, this snapshot is the strongest single
//   regression net for Phase 5 restructuring. Conversely, a failure here
//   does NOT tell you which sub-component changed — that's what the
//   per-component characterisation files (TestCoreAdviser*, TestCoilAdviser*,
//   TestWireAdviser*, TestCore*CrossReferencer*Characterisation.cpp) are for.
//
// FIXTURE
//   Same 3-winding fixture as TestMagneticAdviser.cpp:Test_MagneticAdviser:
//   24:78:76 turns, 100 µH, 1 Vpp triangular @ 507.026 kHz, 25 °C, default
//   insulation requirements, IEC 60664-1.
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true (or leave a snapshot empty), run, copy
//   BASELINE lines from stderr.
//
// BENCHMARKS  (tag: [!benchmark])
//   This is the most expensive benchmark in the entire suite. Use
//      --benchmark-samples 3 --benchmark-warmup-time 0
// =============================================================================

#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "advisers/MagneticAdviser.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "support/Settings.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

constexpr bool kRegenerateBaselines = false;
constexpr double kRelTol = 1e-6;

struct MagneticEntry {
    std::string reference;          // magnetic.manufacturerInfo.reference
    std::string wirePerWinding;
    double      score;
};

// --- Fixture builder --------------------------------------------------------

OpenMagnetics::Inputs make_three_winding_inputs() {
    std::vector<int64_t> numberTurns = {24, 78, 76};
    std::vector<double> turnsRatios;
    for (size_t i = 1; i < numberTurns.size(); ++i) {
        turnsRatios.push_back(double(numberTurns[0]) / numberTurns[i]);
    }

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        /*frequency*/         507026,
        /*magInductance*/     100e-6,
        /*temperature*/       25,
        /*waveShape*/         WaveformLabel::TRIANGULAR,
        /*peakToPeak*/        1,
        /*dutyCycle*/         0.5,
        /*dcCurrent*/         0,
        /*turnsRatios*/       turnsRatios);

    auto requirements = inputs.get_mutable_design_requirements();
    auto insulationRequirements = requirements.get_insulation().value();
    auto standards = std::vector<InsulationStandards>{InsulationStandards::IEC_606641};
    insulationRequirements.set_standards(standards);
    requirements.set_insulation(insulationRequirements);
    inputs.set_design_requirements(requirements);
    inputs.process();
    return inputs;
}

std::string extract_reference(OpenMagnetics::Mas& mas) {
    auto& mag = mas.get_mutable_magnetic();
    if (mag.get_manufacturer_info()
        && mag.get_manufacturer_info().value().get_reference()) {
        return mag.get_manufacturer_info().value().get_reference().value();
    }
    return "<no-reference>";
}

std::string join_wire_names(OpenMagnetics::Mas& mas) {
    std::ostringstream oss;
    auto wires = mas.get_mutable_magnetic().get_wires();
    for (size_t i = 0; i < wires.size(); ++i) {
        if (i) oss << " || ";
        oss << wires[i].get_name().value_or("<unnamed>");
    }
    return oss.str();
}

void check_top_n(const std::string& label,
                 std::vector<std::pair<OpenMagnetics::Mas, double>>& got,
                 const std::vector<MagneticEntry>& expectedTop) {
    if (kRegenerateBaselines || expectedTop.empty()) {
        std::cerr << "\nBASELINE MAGNETIC_TOPN " << label
                  << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            std::cerr << "  [" << i << "] ref=\"" << extract_reference(got[i].first)
                      << "\" wires=\"" << join_wire_names(got[i].first)
                      << "\" score=" << std::setprecision(17) << got[i].second
                      << "\n";
        }
        REQUIRE(got.size() >= expectedTop.size());
        for (size_t i = 1; i < got.size(); ++i) {
            REQUIRE(got[i].second <= got[i - 1].second);
        }
        return;
    }
    INFO("scenario=" << label);
    REQUIRE(got.size() >= expectedTop.size());
    for (size_t i = 1; i < got.size(); ++i) {
        REQUIRE(got[i].second <= got[i - 1].second);
    }
    for (size_t i = 0; i < expectedTop.size(); ++i) {
        auto ref = extract_reference(got[i].first);
        auto wires = join_wire_names(got[i].first);
        INFO("top[" << i << "] want ref=\"" << expectedTop[i].reference
                    << "\" got ref=\"" << ref << "\" wires=\"" << wires
                    << "\" score=" << std::setprecision(17) << got[i].second);
        REQUIRE(ref == expectedTop[i].reference);
        REQUIRE(wires == expectedTop[i].wirePerWinding);
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, kRelTol));
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS — captured 2026-05-19. Leave empty to auto-harvest.
// ----------------------------------------------------------------------------
// Phase 1 fix landed: previously the top-N had Heavy Build wires
// preferred over Single Build for the same conductor diameter due to
// non-deterministic stable-sort tie-break order. WireAdviser now breaks
// score ties by ascending outer diameter (thinner insulation wins), so
// every winding here now resolves to Single Build. Top-2 also swapped
// rankings because the new wire selection changed the magnetic-level
// scoring.
// Refreshed 2026-06-16 (ABT #10) after the 2026-06 saturation rework
// (5000909f inductor reclassification + 60fe7c79 gap-aware isat gate). The
// former winner — 96 E 8.3/4 2 stacks at only 10 turns — now FAILS the
// saturation gate (under-turned → B_peak too high), so the adviser promotes
// the 79 E 10/3 at 20 turns, which has the headroom to stay valid (top-1
// score exactly 2.0). NOTE: the stricter gate thinned the valid 3-winding
// pool for this 507 kHz spec to two designs, so slot 2 is now an
// INVALID-penalised candidate (score 0.85) — that is the magnetic adviser's
// documented thin-pool fallback, not a silenced failure; top-1/top-2 are
// both valid. If the valid pool should be deeper here, that is a separate
// investigation (the saturation gate is correct; the query is demanding).
const std::vector<MagneticEntry> kTopThreeWinding = {
    {"79 E 10/3, Turns: 20, Order: 012, Non-Interleaved, Margin Taped 00",
     "Round 33.0 - Single Build || Round 41.0 - Single Build || Round 41.0 - Single Build",
     2.0},
    {"79 E 10/3, Turns: 20, Order: 012, Non-Interleaved, Margin Taped 01",
     "Round 33.0 - Single Build || Round 41.0 - Heavy Build || Round 41.0 - Single Build",
     1.8855006474886866},
    {"INVALID (failed validity filters): 95 E 8.3/4 2 stacks, Turns: 34, Order: 102, Interleaved, Margin Taped 01",
     "Round 39.0 - Single Build || Round 44.0 - Single Build || Round 43.5 - Single Build",
     0.84999999999999998},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOT
// =============================================================================

TEST_CASE("MagneticAdviser 3-winding end-to-end top-3 snapshot",
          "[adviser][magnetic-adviser][characterisation][end-to-end]") {
    settings.reset();
    settings.set_coil_allow_margin_tape(true);
    settings.set_coil_allow_insulated_wire(false);
    settings.set_coil_fill_sections_with_margin_tape(true);

    auto inputs = make_three_winding_inputs();

    MagneticAdviser adviser;
    auto results = adviser.get_advised_magnetic(inputs, 3);
    check_top_n("3W_24_78_76_507kHz", results, kTopThreeWinding);

    settings.reset();
}

// =============================================================================
// BENCHMARK  (opt-in via [!benchmark])
// =============================================================================
//
// IMPORTANT: this benchmark exercises the entire MagneticAdviser pipeline
// (core selection + coil generation + scoring). Even a single sample takes
// many seconds. ALWAYS invoke with:
//
//   ./MKF_tests "[benchmark-magnetic-adviser]" --benchmark-samples 3 --benchmark-warmup-time 0
//

TEST_CASE("Benchmark MagneticAdviser 3-winding end-to-end (top-3)",
          "[adviser][magnetic-adviser][!benchmark][benchmark-magnetic-adviser]") {
    settings.reset();
    settings.set_coil_allow_margin_tape(true);
    settings.set_coil_allow_insulated_wire(false);
    settings.set_coil_fill_sections_with_margin_tape(true);

    auto inputs = make_three_winding_inputs();

    BENCHMARK("get_advised_magnetic 3W top-3") {
        MagneticAdviser adviser;
        return adviser.get_advised_magnetic(inputs, 3);
    };

    settings.reset();
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                          | mean (3 samples) | notes
// -----------+-----------------------------------+------------------+----------
// 2026-05-19 | 3-winding 24:78:76 / 507 kHz top-3| TBD              | initial
// 2026-05-25 | 3-winding 24:78:76 / 507 kHz top-3| 44.52 s ± 590 ms | first recorded baseline; uncontested host;
//            |                                   |                  | after Phase-7 dispatch unification +
//            |                                   |                  | sat-margin 1.0→1.2 default flip
//
// =============================================================================
