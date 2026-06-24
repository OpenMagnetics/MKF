// =============================================================================
// TestWireAdviserCharacterisation.cpp
// =============================================================================
// Characterisation tests for WireAdviser::get_advised_wire (standard wires)
// and ::get_advised_planar_wire.
//
// PURPOSE
//   Lock the CURRENT top-N (wire-name, score) output for a representative
//   winding fixture, so any refactor of the 675-line WireAdviser preserves
//   ranking AND scores to within 1e-6 relative tolerance.
//
// FIXTURE
//   Section 5 mm x 15 mm, 42 turns, triangular 2 A_rms @ 13.456 kHz, 22 °C.
//   This matches the historical TestWireAdviser.cpp top-of-file defaults.
//
// SCENARIOS
//   1. ALL_TYPES top-5            — round + litz + foil + rectangular enabled.
//   2. LITZ_ONLY_HIGH_FREQ top-5  — same section, 200 kHz, litz only.
//   3. PLANAR top-5               — get_advised_planar_wire path.
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true (or leave a snapshot empty), run, copy
//   BASELINE lines from stderr into the kSnapshots tables.
//
// BENCHMARKS  (tag: [!benchmark])
//   Time get_advised_wire end-to-end for the ALL_TYPES scenario (full DB).
// =============================================================================

#include <iomanip>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "advisers/WireAdviser.h"
#include "constructive_models/Coil.h"
#include "support/Settings.h"
#include "support/Utils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

constexpr bool kRegenerateBaselines = false;
constexpr double kRelTol = 1e-6;

struct TopEntry {
    std::string name;
    double score;
};

// --- Fixture builders -------------------------------------------------------

Winding make_winding(int64_t numberTurns) {
    Winding w;
    w.set_isolation_side(IsolationSide::PRIMARY);
    w.set_name("primary");
    w.set_number_parallels(1);
    w.set_number_turns(numberTurns);
    w.set_wire("Dummy");
    return w;
}

Section make_section(double width, double height) {
    Section s;
    s.set_dimensions({width, height});
    s.set_coordinate_system(CoordinateSystem::CARTESIAN);
    return s;
}

SignalDescriptor make_current(double rms, double effectiveFreq) {
    SignalDescriptor c;
    ProcessedWaveform p;
    p.set_peak_to_peak(rms * 2 * 1.4142);
    p.set_rms(rms);
    p.set_effective_frequency(effectiveFreq);
    p.set_offset(0);
    p.set_label(WaveformLabel::TRIANGULAR);
    c.set_processed(p);
    return c;
}

void enable_all_wire_types() {
    settings.set_wire_adviser_include_round(true);
    settings.set_wire_adviser_include_litz(true);
    settings.set_wire_adviser_include_foil(true);
    settings.set_wire_adviser_include_rectangular(true);
    settings.set_wire_adviser_include_planar(false);
}

void enable_only_litz() {
    settings.set_wire_adviser_include_round(false);
    settings.set_wire_adviser_include_litz(true);
    settings.set_wire_adviser_include_foil(false);
    settings.set_wire_adviser_include_rectangular(false);
    settings.set_wire_adviser_include_planar(false);
}

// --- Snapshot helper --------------------------------------------------------

void check_top_n(const std::string& label,
                 const std::vector<std::pair<Winding, double>>& got,
                 const std::vector<TopEntry>& expectedTop) {
    if (kRegenerateBaselines || expectedTop.empty()) {
        std::cerr << "\nBASELINE TOPN " << label << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            auto wire = OpenMagnetics::Coil::resolve_wire(got[i].first);
            std::string name = wire.get_name().value_or("<unnamed>");
            std::cerr << "  [" << i << "] name=\"" << name << "\" score="
                      << std::setprecision(17) << got[i].second << "\n";
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
        auto wire = OpenMagnetics::Coil::resolve_wire(got[i].first);
        std::string name = wire.get_name().value_or("<unnamed>");
        INFO("top[" << i << "] want=\"" << expectedTop[i].name
                    << "\" got=\"" << name << "\" score="
                    << std::setprecision(17) << got[i].second);
        REQUIRE(name == expectedTop[i].name);
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, kRelTol));
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS — captured 2026-05-19. Leave a vector empty to auto-harvest.
// ----------------------------------------------------------------------------
// Phase 1 fix landed: previously slots 2 & 3 were "Round 16.5 - Heavy Build"
// and "Round 16.5 - Single Build" tied at 3.9969312641903958 with no
// deterministic tie-break (both have identical conductor geometry so the
// resistance/skin/proximity filters score them identically). The
// WireAdviser now breaks score ties by ascending outer diameter, which
// promotes the thinner-insulation variant. Same fix promoted Round 17.0
// Single Build over Heavy Build in slot 4.
const std::vector<TopEntry> kTopAllTypes = {
    {"Round 16.0 - Single Build", 3.9999898124395115},
    {"Round 1.25 - Grade 1",      3.9983446101392519},
    {"Round 16.5 - Single Build", 3.9969312641903958},
    {"Round 16.5 - Heavy Build",  3.9969312641903958},
    {"Round 17.0 - Single Build", 3.9933077449224341},
};

// Phase 1 fix landed: previously top-2 tied at 3.9996749837491876 with
// Single-Served before Unserved with no deterministic tie-break. The
// outer-diameter tie-break now promotes Unserved (no serving = smaller
// outer diameter).
const std::vector<TopEntry> kTopLitzHighFreq = {
    {"Litz 10x0.3 - Grade 1 - Unserved",      3.9996749837491876},
    {"Litz 10x0.3 - Grade 1 - Single Served", 3.9996749837491876},
    {"Litz 6x0.4 - Grade 1 - Unserved",       3.9993320286548131},
    {"Litz 45x0.14 - Grade 1 - Unserved",     3.9953893064831778},
    {"Litz 4x0.5 - Grade 1 - Unserved",       3.995096563768834},
};

// Phase 1 fix landed: the planar dataset generator (create_planar_dataset)
// emitted the same catalogue wire at multiple parallel counts with the
// SAME display name (e.g. two "Planar 243.59 µm" rows with different
// scores at ranks 0 and 3). The parallels count is now encoded in the
// display name so users can disambiguate.
const std::vector<TopEntry> kTopPlanar = {
    {"Planar 243.59 µm x2 parallels", 1.3196636443022853},
    {"Planar 278.38 µm x2 parallels", 1.2926461763463761},
    {"Planar 208.79 µm x2 parallels", 1.2613406508784197},
    {"Planar 243.59 µm x3 parallels", 1.1647936870013549},
    {"Planar 313.18 µm x2 parallels", 1.1599709209384959},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("WireAdviser ALL_TYPES top-5 snapshot",
          "[adviser][wire-adviser][characterisation][all-types]") {
    settings.reset();
    clear_databases();
    enable_all_wire_types();

    auto winding = make_winding(42);
    auto section = make_section(0.005, 0.015);
    auto current = make_current(/*rms*/ 2, /*effFreq*/ 13456);

    WireAdviser adviser;
    auto results = adviser.get_advised_wire(
        winding, section, current, /*temp*/ 22, /*numberSections*/ 1, /*maxN*/ 5);
    check_top_n("ALL_TYPES_42T_2A_13k", results, kTopAllTypes);
}

TEST_CASE("WireAdviser LITZ_ONLY high-frequency top-5 snapshot",
          "[adviser][wire-adviser][characterisation][litz-high-freq]") {
    settings.reset();
    clear_databases();
    enable_only_litz();

    auto winding = make_winding(42);
    auto section = make_section(0.005, 0.015);
    auto current = make_current(/*rms*/ 2, /*effFreq*/ 200000);

    WireAdviser adviser;
    auto results = adviser.get_advised_wire(
        winding, section, current, /*temp*/ 22, /*numberSections*/ 1, /*maxN*/ 5);
    check_top_n("LITZ_ONLY_42T_2A_200k", results, kTopLitzHighFreq);
}

TEST_CASE("WireAdviser PLANAR top-5 snapshot",
          "[adviser][wire-adviser][characterisation][planar]") {
    settings.reset();
    clear_databases();
    settings.set_wire_adviser_include_planar(true);

    auto winding = make_winding(8);
    auto section = make_section(0.010, 0.0015);
    auto current = make_current(/*rms*/ 5, /*effFreq*/ 200000);

    WireAdviser adviser;
    auto results = adviser.get_advised_planar_wire(
        winding, section, current, /*temp*/ 22, /*numberSections*/ 1, /*maxN*/ 5);
    check_top_n("PLANAR_8T_5A_200k", results, kTopPlanar);
}

// =============================================================================
// BENCHMARKS  (opt-in via [!benchmark])
// =============================================================================
//
// Use:  ./MKF_tests "[benchmark-wire-adviser]" --benchmark-samples 5 --benchmark-warmup-time 0
//

TEST_CASE("Benchmark WireAdviser ALL_TYPES (top-20)",
          "[adviser][wire-adviser][!benchmark][benchmark-wire-adviser]") {
    settings.reset();
    clear_databases();
    enable_all_wire_types();

    auto winding = make_winding(42);
    auto section = make_section(0.005, 0.015);
    auto current = make_current(2, 13456);

    BENCHMARK("get_advised_wire ALL_TYPES top-20") {
        WireAdviser adviser;
        return adviser.get_advised_wire(winding, section, current, 22, 1, 20);
    };
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                | mean (5 samples) | notes
// -----------+-------------------------+------------------+----------
// 2026-05-19 | ALL_TYPES top-20        | TBD              | initial
//
// =============================================================================
