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
// Refreshed 2026-06-17 (ABT #14): the ABT #10 proximity-aware re-rank was
// removed (the post-gap fringing ceiling tightened to 1.3 prunes the
// catastrophic gap-fringing cores up front, making the per-advise re-rank
// winding redundant). Without the re-rank the ordering follows the raw
// cost/loss/dimensions score again: EP 20 3C96 keeps slot 1 (score 4.039 →
// 3.954 once the re-rank's proximity bonus is gone) and the tail fills with
// PQ 20/20 variants instead of EFD/E. All valid small-and-sane power ferrites
// (no tiny winding-killer cores).
// Refreshed 2026-07-07 (ABT #126 investigation): drift root-caused bidirectionally
// to two verified-correct changes — ABT #70 (dummy-coil parallels sized to current
// density, shifting the loss-filter renormalization) and the health-pass
// MagneticEnergy min/max interval fix (energy target at L_max -> 1.5x gap energy).
// Same top core; 3C96/3C94 swap in slots 1-2; scores -1.4%. Regenerated with
// kRegenerateBaselines on main 17e1f850.
const std::vector<TopEntry> kTopAvailablePower = {
    {"EP 20 - 3C96 - Gapped 0.375 mm",             3.8966756636009796},
    {"PQ 20/20 - 3C96 - Gapped 0.46900000000000003 mm", 3.854716197272996},
    {"PQ 20/20 - 3C94 - Gapped 0.472 mm",          3.8494928517386233},
    {"PQ 20/20 - 3C90 - Gapped 0.472 mm",          3.8341659260385046},
    {"PQ 20/20 - 3C97 - Gapped 0.477 mm",          3.8233252319348403},
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
// Refreshed 2026-06-17 (ABT #14): re-rank removal + 1.3 fringing ceiling (see
// kTopAvailablePower note). 95 EQ 25/6 and 95 RM 10/13 were tied on score
// (3.8963…); without the re-rank the lexicographic tiebreak now seats EQ 25/6
// at slot 1, and the tail fills with EP 20 / RM 10 variants in place of the
// E 19/E 21 2-stacks. All valid small-and-sane power ferrites.
// Refreshed 2026-07-07 (ABT #126 investigation): same two root causes as
// kTopAvailablePower. The 1.5x gap-energy requirement re-solves every standard
// gap larger (RM 10/ILP 0.24 -> 0.32 mm, exactly 1.5x on the E-stacks), which
// displaces EQ 25/6 and seats RM 10/ILP at slot 1; the paralleled dummy coil
// admits the 3/4-stack E 16/E 19 variants (physically sane for 100 uH /
// 100 kHz / 600 Vpp). Regenerated on main 17e1f850, user-approved.
const std::vector<TopEntry> kTopStandardPower = {
    {"98 RM 10/ILP gapped 0.32 mm",                3.8624574793397786},
    {"98 PQ 27/15 gapped 0.25 mm",                 3.7835392946594726},
    {"98 E 16/7/5 4 stacks gapped 0.09 mm",        3.7500041208232733},
    {"95 RM 10/ILP gapped 0.32 mm",                3.6867104806869007},
    {"98 E 19/8/9 3 stacks gapped 0.08 mm",        3.6749555384431662},
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
    adviser.set_application(MAS::MagneticApplication::INTERFERENCE_SUPPRESSION);
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
        adviser.set_application(MAS::MagneticApplication::INTERFERENCE_SUPPRESSION);
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
