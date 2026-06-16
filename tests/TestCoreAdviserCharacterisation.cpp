// =============================================================================
// TestCoreAdviserCharacterisation.cpp
// =============================================================================
// Characterisation tests for the CoreAdviser end-to-end pipeline.
//
// PURPOSE
//   Lock the CURRENT top-N output of CoreAdviser::get_advised_core across the
//   matrix (mode × application), so that any refactor of the 3 196-line
//   CoreAdviser.cpp (Phase 5 in the audit) preserves the exact ranking and
//   scores to within 1e-6 relative tolerance.
//
// MATRIX
//   AVAILABLE_CORES  × POWER                       (typical inductor query)
//   STANDARD_CORES   × POWER                       (catalogue-only query)
//   AVAILABLE_CORES  × INTERFERENCE_SUPPRESSION    (toroidal CMC-style query)
//
// SHAPE OF EACH SNAPSHOT TEST
//   - Fixed inputs (frequency, voltage, current, inductance, etc.).
//   - Fixed cores fixture (tests/testData/test_cores.ndjson, ~4800 cores after
//     parse-skips) — already used by TestCoreAdviser.cpp.
//   - Strict equality on the ordered top-N core names.
//   - WithinRel(1e-6) on each score.
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true, run, copy the printed BASELINE lines
//   into the kSnapshots tables below, flip back.
//
// BENCHMARKS  (tag: [!benchmark])
//   Time get_advised_core on the full fixture for the three scenarios above.
//   Baselines recorded at the bottom of the file.
// =============================================================================

#include <source_location>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "advisers/CoreAdviser.h"
#include "processors/Inputs.h"
#include "physical_models/Impedance.h"
#include "support/Settings.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

constexpr bool kRegenerateBaselines = false;
constexpr double kRelTol = 1e-6;

// --- Fixture loader ---------------------------------------------------------
std::vector<Core> load_test_cores_fixture() {
    auto path = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "test_cores.ndjson");
    std::ifstream in(path);
    std::string line;
    std::vector<Core> cores;
    while (std::getline(in, line)) {
        try {
            json jf = json::parse(line);
            cores.emplace_back(jf, false, true, false);
        } catch (const CoreShapeNotFoundException&) {
            continue;
        } catch (const std::exception&) {
            continue;
        }
    }
    return cores;
}

// --- Inputs builders --------------------------------------------------------
OpenMagnetics::Inputs make_power_inputs() {
    // Same parameter set as TestCoreAdviser.cpp:Test_CoreAdviserAvailableCores_All_Cores.
    return OpenMagnetics::Inputs::create_quick_operating_point(
        /*frequency*/ 100000,
        /*magnetizingInductance*/ 100e-6,
        /*temperature*/ 25,
        /*label*/ WaveformLabel::SINUSOIDAL,
        /*peakToPeak*/ 600,
        /*dutyCycle*/ 0.5,
        /*dcCurrent*/ 0,
        /*turnsRatios*/ {});
}

OpenMagnetics::Inputs make_interference_suppression_inputs() {
    // Mirrors TestCoreAdviser.cpp:Test_CoreAdviserAvailableCores_Toroidal_Cores_With_Impedance.
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(
        /*frequency*/ 100000,
        /*magnetizingInductance*/ 1e-6,
        /*temperature*/ 25,
        /*label*/ WaveformLabel::SINUSOIDAL,
        /*peakToPeak*/ 60,
        /*dutyCycle*/ 0.5,
        /*dcCurrent*/ 0,
        /*turnsRatios*/ {});
    std::vector<ImpedanceAtFrequency> minimumImpedance;
    ImpedancePoint ip;
    ip.set_magnitude(1000);
    ImpedanceAtFrequency iaf;
    iaf.set_frequency(100000);
    iaf.set_impedance(ip);
    minimumImpedance.push_back(iaf);
    inputs.get_mutable_design_requirements().set_minimum_impedance(minimumImpedance);
    inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_minimum(1e-6);
    inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_nominal(std::nullopt);
    inputs.get_mutable_design_requirements().get_mutable_magnetizing_inductance().set_maximum(std::nullopt);
    return inputs;
}

// --- Snapshot type ----------------------------------------------------------
struct TopEntry {
    std::string name;
    double score;
};

void check_top_n(const std::string& label,
                 const std::vector<std::pair<OpenMagnetics::Mas, double>>& got,
                 const std::vector<TopEntry>& expectedTop) {
    // Print baseline lines whenever regenerating OR when no snapshot has
    // been recorded yet (expectedTop empty). Empty-snapshot still enforces
    // the monotonic-ordering structural contract.
    if (kRegenerateBaselines || expectedTop.empty()) {
        UNSCOPED_INFO("BASELINE TOPN " << label << " count=" << got.size());
        for (size_t i = 0; i < got.size(); ++i) {
            auto name = got[i].first.get_magnetic().get_core().get_name().value_or("<unnamed>");
            UNSCOPED_INFO("  [" << i << "] name=\"" << name << "\" score="
                          << std::setprecision(17) << got[i].second);
        }
        std::cerr << "\nBASELINE TOPN " << label << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            auto name = got[i].first.get_magnetic().get_core().get_name().value_or("<unnamed>");
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
    // Verify scores are monotonically non-increasing (the contract).
    for (size_t i = 1; i < got.size(); ++i) {
        REQUIRE(got[i].second <= got[i - 1].second);
    }
    // Strict equality on names; WithinRel on scores.
    for (size_t i = 0; i < expectedTop.size(); ++i) {
        auto name = got[i].first.get_magnetic().get_core().get_name().value_or("<unnamed>");
        INFO("top[" << i << "] expected=\"" << expectedTop[i].name << "\" got=\"" << name
                    << "\" score=" << std::setprecision(17) << got[i].second);
        REQUIRE(name == expectedTop[i].name);
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, kRelTol));
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS
// Baselines captured 2026-05-19. Top entries are ordering-sensitive; scores
// are within 1e-6 relative tolerance.
// ----------------------------------------------------------------------------
// Refreshed 2026-06-16 (ABT #10) after the proximity-aware re-rank landed
// (5bd45414 "feat(core-adviser): proximity-aware re-rank for single-winding
// inductors"). This POWER query is a single-winding inductor, so the adviser
// now winds+simulates the top-K and re-orders them by proximity-inclusive
// loss. That promotes EP 20 3C91 0.605 mm from slot 2 → slot 1 (its score
// jumps 3.830 → 3.992 once its lower real winding-proximity loss is counted)
// while PQ 20/20 3C97 is essentially unchanged (3.929 → 3.930) and slips to
// slot 2. The Kool Mµ Hƒ 40 toroid also enters the top-5. EP 20 winning by
// ~1.6 % over PQ 20/20 is the more-accurate ranking, not a regression.
// (Previous refresh was the saturation-margin 1.0 → 1.2 flip.)
// Refreshed 2026-06-16 (ABT #13) after the saturation-derating rework: the
// inductor saturation gate/sizing now evaluate at the 100 C hot junction corner
// with RAW B_sat (derating = hot temperature × margin, no 0.7 flux-proportion
// stacked). Raw B_sat at 100 C is ~1.14x looser than the former 0.7·B_sat@25 C
// gate, so the saturation-headroom scoring reshuffles the top-5 toward cores
// with real headroom (EP 20 3C96, EFD 25/13/9, E 25/13/7). All five are valid
// power ferrites; scores nudged < 1.5 %.
const std::vector<TopEntry> kTopAvailablePower = {
    {"EP 20 - 3C96 - Gapped 0.375 mm",             4.0394759597245846},
    {"EFD 25/13/9 - 97 - Gapped 0.5 mm",           4.0371918250916554},
    {"E 25/13/7 - 3C90 - Gapped 1.0 mm",           4.0259647054204928},
    {"PQ 20/20 - 3C94 - Gapped 0.472 mm",          3.9958256614081851},
    {"PQ 20/20 - 3C96 - Gapped 0.46900000000000003 mm", 3.9508508473488773},
};

// STANDARD_CORES x POWER: top-5 unique standard-shape ferrite candidates.
// Phase 1 fix landed: coreShapeDatabase had been keyed by canonical-name AND
// each alias, so iterating it for the STANDARD_CORES path yielded the same
// CoreShape (1 + #aliases) times and produced duplicate top-N entries
// (previously slots 0-2 were three identical copies of
// "95 PQ 27/15 gapped 0.36 mm"). See CoreAdviser.cpp get_advised_core
// dispatch path for the canonical-name dedupe.
// Refreshed 2026-06-16 (ABT #10): same proximity-aware re-rank as the
// AVAILABLE_CORES_POWER table above promotes 95 EP 20 (gapped 0.32 mm) to
// slot 1 over 95 PQ 27/15, and the post-rerank gap selection on PQ 27/15
// settles at 0.18 mm. Single-winding-inductor re-rank, intended.
// Refreshed 2026-06-16 (ABT #13): same saturation-derating rework as
// kTopAvailablePower above (RAW B_sat at the 100 C hot corner). The
// saturation-headroom rescoring promotes 95 RM 10/13 to slot 1 and reshuffles
// the standard-shape top-5; all valid, scores nudged < 1.5 %.
const std::vector<TopEntry> kTopStandardPower = {
    {"95 RM 10/13 gapped 0.27 mm",                 3.8963165812464418},
    {"98 E 19/8/9 2 stacks gapped 0.12 mm",        3.7753243780925869},
    {"98 EP 20 gapped 0.32 mm",                    3.7471654879366123},
    {"95 E 21/9/5 2 stacks gapped 0.21 mm",        3.715779819799347},
    {"98 E 21/9/5 2 stacks gapped 0.21 mm",        3.701426173763287},
};

// Refreshed 2026-06-16 (ABT #10) after landing the suppression returns-0 fix
// (3-layer: a905fd2 pipeline pass-through + 957c516 unscorable-material +
// 3be26c6 cross-ref loss-skip, cherry-picked from
// feat/wideband-impedance-resonances). Before the fix this scenario returned
// ZERO candidates on main — the losses filter rejected/threw on every complex-
// permeability suppression ferrite and the adviser produced an empty set. With
// the fix the SAME three cores as the 2026-05-19 baseline flow again (XFlux 26
// toroids, identical identities + order); only the scores nudged < 0.6 % from
// the 2026-06 loss/saturation recompute.
const std::vector<TopEntry> kTopAvailableInterference = {
    {"T 134/77/155 - XFlux 26 - Ungapped",         1.9451312382492831},
    {"T 134/77/78 - XFlux 26 - Ungapped",          1.8563034258168774},
    {"T 167/87/27 - XFlux 26 - Ungapped",          1.5057090764821521},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("CoreAdviser AVAILABLE_CORES x POWER top-5 snapshot",
          "[adviser][core-adviser][characterisation][available-cores][power]") {
    clear_databases();
    settings.reset();
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_power_inputs();
    auto cores = load_test_cores_fixture();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    CoreAdviser adviser;
    adviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    auto results = adviser.get_advised_core(inputs, weights, &cores, 5);

    check_top_n("AVAILABLE_CORES_POWER", results, kTopAvailablePower);
}

TEST_CASE("CoreAdviser STANDARD_CORES x POWER top-5 snapshot",
          "[adviser][core-adviser][characterisation][standard-cores][power]") {
    clear_databases();
    settings.reset();
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_power_inputs();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    CoreAdviser adviser;
    adviser.set_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto results = adviser.get_advised_core(inputs, weights, 5);

    check_top_n("STANDARD_CORES_POWER", results, kTopStandardPower);
}

TEST_CASE("CoreAdviser AVAILABLE_CORES x INTERFERENCE_SUPPRESSION top-3 snapshot",
          "[adviser][core-adviser][characterisation][available-cores][interference-suppression]") {
    clear_databases();
    settings.reset();
    settings.set_use_concentric_cores(false);
    settings.set_use_toroidal_cores(true);
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_interference_suppression_inputs();
    auto cores = load_test_cores_fixture();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 0;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0;

    CoreAdviser adviser;
    adviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    adviser.set_application(MAS::Application::INTERFERENCE_SUPPRESSION);
    auto results = adviser.get_advised_core(inputs, weights, &cores, 3);

    check_top_n("AVAILABLE_CORES_INTERFERENCE_SUPPRESSION", results, kTopAvailableInterference);
}

// =============================================================================
// BENCHMARKS  (opt-in via [!benchmark])
// =============================================================================
// The CoreAdviser pipeline is large; we benchmark end-to-end get_advised_core
// rather than internal stages — that's what users actually care about.
//
// IMPORTANT: Catch2 defaults to 100 samples per BENCHMARK, which means a
// single full-DB run can take several minutes. Always invoke with a small
// sample count for these, e.g.:
//
//   ./MKF_tests "[benchmark-available-power]" --benchmark-samples 3
//

TEST_CASE("Benchmark CoreAdviser AVAILABLE_CORES x POWER (full DB, top-50)",
          "[adviser][core-adviser][!benchmark][benchmark-available-power]") {
    clear_databases();
    settings.reset();
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_power_inputs();
    auto cores = load_test_cores_fixture();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    BENCHMARK("get_advised_core AVAILABLE x POWER full DB top-50") {
        CoreAdviser adviser;
        adviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
        return adviser.get_advised_core(inputs, weights, &cores, 50);
    };
}

TEST_CASE("Benchmark CoreAdviser STANDARD_CORES x POWER (top-50)",
          "[adviser][core-adviser][!benchmark][benchmark-standard-power]") {
    clear_databases();
    settings.reset();
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_power_inputs();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 1;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

    BENCHMARK("get_advised_core STANDARD x POWER top-50") {
        CoreAdviser adviser;
        adviser.set_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
        return adviser.get_advised_core(inputs, weights, 50);
    };
}

TEST_CASE("Benchmark CoreAdviser AVAILABLE_CORES x INTERFERENCE_SUPPRESSION (full DB, top-10)",
          "[adviser][core-adviser][!benchmark][benchmark-available-interference]") {
    clear_databases();
    settings.reset();
    settings.set_use_concentric_cores(false);
    settings.set_use_toroidal_cores(true);
    settings.set_use_only_cores_in_stock(false);

    auto inputs = make_interference_suppression_inputs();
    auto cores = load_test_cores_fixture();

    std::map<CoreAdviser::CoreAdviserFilters, double> weights;
    weights[CoreAdviser::CoreAdviserFilters::COST] = 0;
    weights[CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
    weights[CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0;

    BENCHMARK("get_advised_core AVAILABLE x IS full DB top-10") {
        CoreAdviser adviser;
        adviser.set_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
        adviser.set_application(MAS::Application::INTERFERENCE_SUPPRESSION);
        return adviser.get_advised_core(inputs, weights, &cores, 10);
    };
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                                | mean (1 sample) | notes
// -----------+-----------------------------------------+-----------------+----------
// 2026-05-19 | AVAILABLE_CORES x POWER (full, top-50)  | TBD             | initial
// 2026-05-19 | STANDARD_CORES x POWER (top-50)         | TBD             | initial
// 2026-05-19 | AVAILABLE x INTERFERENCE (full, top-10) | TBD             | initial
//
// Each row is "time per get_advised_core call". The full-DB scenarios are
// known to be expensive (tens of seconds each); when running for the first
// time, use `--benchmark-samples 3 --benchmark-warmup-time 0`.
// =============================================================================
