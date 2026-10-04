// =============================================================================
// TestMagneticFilter.cpp
// =============================================================================
// Characterisation tests for the MagneticFilter family (src/advisers/MagneticFilter.{h,cpp}).
//
// PURPOSE
//   Lock the CURRENT behaviour of every factory-creatable MagneticFilter
//   subclass — including the known bugs documented in the audit — so that the
//   upcoming "no silent fallbacks" cleanup (Phase 1) can flip them to throws
//   loudly instead of silently changing scores. Also provides micro-benchmarks
//   so that any refactor preserves (or improves) performance.
//
// SHAPE OF EACH TEST
//   - Build a fixed reference Magnetic (inductor-style: E 35 / 3C97, 40+20 turns,
//     Round 1.00 - Grade 1 wire, ungapped — a design that FITS its window, ABT #785).
//   - Build a fixed reference Inputs (100 kHz triangular, 100 µH, 25 °C, ±√3 A).
//   - factory() the filter (where allowed).
//   - evaluate_magnetic(...) once.
//   - REQUIRE valid == snapshot && score within 1e-6 relative tolerance of snapshot.
//
// HOW TO REGENERATE SNAPSHOTS
//   When a deliberate behaviour change is made:
//     1. Set kRegenerateBaselines = true (below).
//     2. Build & run `./MKF_tests "[magnetic-filter][characterisation]"` —
//        the tests will print BASELINE lines instead of asserting.
//     3. Copy the printed numbers into the SNAPSHOTS table below.
//     4. Set kRegenerateBaselines = false and rerun to verify.
//
// BUG-LOCKING
//   Filters with KNOWN bugs (audit-tracked) have an explicit comment
//   "LOCKS BUG <id>" next to their snapshot. When Phase 1 fixes the bug, the
//   snapshot will (correctly) fail — that's the signal to update it.
//
// BENCHMARKS  (tag: [!benchmark], opt-in via Catch2)
//   For each filter, evaluate_magnetic is timed against tiered core sets
//   (10 / 100 / full test_cores.ndjson). Baseline numbers recorded at the
//   bottom of this file as comments — record new ones whenever a refactor
//   touches the hot path.
// =============================================================================

#include <source_location>
#include <cmath>
#include <numbers>
#include <filesystem>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>
#include <limits>
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "advisers/MagneticFilter.h"
#include "advisers/MagneticAdviser.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Core.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "support/Exceptions.h"
#include "support/Logger.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

// --- Knobs ------------------------------------------------------------------
constexpr bool kRegenerateBaselines = false;
constexpr double kRelTol = 1e-4;

// --- Fixtures ---------------------------------------------------------------
// A single deterministic Magnetic + Inputs pair used by all behaviour tests.
// Keeping it tiny + ungapped keeps filter evaluation fast and avoids dragging
// random external state into the snapshot.
//
// THE WIRE MUST KEEP THE DESIGN WINDABLE (ABT #785). The original fixture used
// Round 2.00 - Grade 1 (outer 2.074 mm) in E 35's 7.5 x 25 mm window, which holds
// floor(25/2.074) = 12 turns per layer and floor(7.5/2.074) = 3 layers — about 36
// turns against the 60 this fixture asks for. are_sections_and_layers_fitting()
// returned FALSE and wind() emitted a layout that cannot exist: winding 0 as a
// SINGLE layer 2.074 x 82.96 mm, 3.3x taller than the window containing it
// (section filling factor 1.84). Every layout-sensitive snapshot was therefore
// pinned to whatever an over-stuffed winder happened to emit, which is why
// DIMENSIONS / COST / AREA_NO_PARALLELS drifted with no code change anyone would
// call a regression, and why nine more sat inert at {false, 0.0}.
//
// Round 1.00 - Grade 1 (outer 1.062 mm) gives 23 turns per layer x 7 layers = 161
// turns of capacity for 60, so the coil genuinely fits and the filters measure a
// real winding. If you change the wire or the turn counts again, CHECK
// are_sections_and_layers_fitting() before re-harvesting anything.
OpenMagnetics::Magnetic make_reference_magnetic() {
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {40, 20};
    cfg.numberParallels = {1, 1};
    cfg.coreShapeName = "E 35";
    cfg.coreMaterialName = "3C97";
    cfg.wireNames = {"Round 1.00 - Grade 1", "Round 1.00 - Grade 1"};
    cfg.numberStacks = 1;
    auto m = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    // Wind the coil so filters that iterate turns/layers/sections (VOLUME, AREA,
    // HEIGHT, *_TIMES_VOLUME*, LOSSES family) don't trip COIL_NOT_PROCESSED.
    m.get_mutable_coil().wind();
    return m;
}

// The same magnetic with its coil explicitly UNWOUND, for the tests that pin the
// loud-failure contract: a filter needing turns/layers must THROW, never score
// silently (the 2026-05-20 "eliminate silent fallbacks" pass).
//
// Those tests used to lean on make_reference_magnetic() failing to wind — wind()
// above could not lay this design out when they were written, so the filters threw
// by accident. The winder lays it out now, so the throw stopped happening: six
// characterisation tests went red AND the contract they guard stopped being
// exercised at all. Unwinding on purpose restores it and makes the tests
// independent of whether any given design happens to fit.
//
// NOT for the tests whose contract is "wound but not delimited/compacted"
// (AREA_WITH_PARALLELS) or "wound, but the inputs lack a spec"
// (CORE_MINIMUM_IMPEDANCE) — those need the WOUND reference and still pass.
OpenMagnetics::Magnetic make_unwound_reference_magnetic() {
    auto m = make_reference_magnetic();
    m.get_mutable_coil().unwind();
    return m;
}

OpenMagnetics::Inputs make_reference_inputs() {
    OpenMagneticsTesting::QuickInputsConfig cfg;
    cfg.frequency = 100000;
    cfg.magnetizingInductance = 100e-6;
    cfg.temperature = 25;
    cfg.label = WaveformLabel::TRIANGULAR;
    cfg.peakToPeak = 2 * 1.73205;
    cfg.dutyCycle = 0.5;
    cfg.offset = 0;
    return OpenMagneticsTesting::create_quick_test_inputs(cfg);
}

// Print baseline if regenerating, otherwise assert against snapshot.
void check_or_print(const std::string& name, bool gotValid, double gotScore,
                    bool expectedValid, double expectedScore, double relTol = kRelTol) {
    if (kRegenerateBaselines) {
        // Emit a line that can be grepped + pasted back as a snapshot.
        UNSCOPED_INFO("BASELINE " << name
                      << " valid=" << (gotValid ? "true" : "false")
                      << " score=" << std::setprecision(17) << gotScore);
        CHECK(std::isfinite(gotScore));  // Always require finite even when regenerating.
        return;
    }
    INFO("filter=" << name << " score=" << std::setprecision(17) << gotScore);
    REQUIRE(gotValid == expectedValid);
    if (std::isnan(expectedScore)) {
        REQUIRE(std::isnan(gotScore));
    } else if (expectedScore == 0.0) {
        REQUIRE(gotScore == 0.0);
    } else {
        REQUIRE_THAT(gotScore, WithinRel(expectedScore, relTol));
    }
}

// Test fixture: pair of (valid, score) snapshots per filter against the
// reference magnetic + inputs. Snapshots are HARVESTED from a first run with
// kRegenerateBaselines=true; the numbers below are the captured truth, NOT
// hand-derived. They represent CURRENT behaviour (warts and all).
struct Snapshot {
    bool valid;
    double score;
};

// IMPORTANT: These snapshots are placeholders set to NaN. On the very first
// run of this file, kRegenerateBaselines must be flipped to true, the test
// run, and the printed BASELINE lines pasted in here.
const std::map<std::string, Snapshot> kSnapshots = {
    // Baselines captured 2026-05-19 against E 35 / 3C97 / 40+20 turns /
    // Round 2.00 - Grade 1 wire, 100 kHz triangular, ±√3 A, 25 °C.
    //
    // RE-HARVESTED 2026-08-17 (ABT #785) after the fixture wire went
    // Round 2.00 → Round 1.00 - Grade 1, because the old design did not fit its
    // winding window and every layout-sensitive entry was characterising an
    // impossible layout (see the note on make_reference_magnetic()). Only four
    // entries moved, and the direction is reassuring: the coil now lays up as
    // 3 conduction layers + 2 insulation layers, which is exactly the stack the
    // 2026-05-19 harvest saw — COST comes back to its ORIGINAL 7.0 on its own.
    //   DIMENSIONS         2.7616400000000005e-05 → 2.0178200000000003e-05
    //   COST               6.0 (drifted) → 7.0 (original value restored)
    //   AREA_NO_PARALLELS  valid true → false  (a real filter bug, see below)
    //   TEMPERATURE        27.158311798960909 → 28.050433964481037
    // TEMPERATURE rises because the thinner wire has more DC resistance; the
    // core-only entries (AREA_PRODUCT, ENERGY_STORED, ESTIMATED_COST,
    // MAGNETIZING_INDUCTANCE, SATURATION) are unmoved, as they should be.
    //
    // The two that moved to a NON-ZERO value are hand-derivable, so they are not
    // merely "whatever the code emitted" (the June 2026 review's trap). The layup
    // is 40 turns -> 2 layers, 20 turns -> 1 layer, plus 2 x 25 um insulation:
    //   COST       = numberLayers + sum(wire relative cost) = 5 + (1 + 1)      = 7
    //   DIMENSIONS = W x H x max(coreDepth, column0Depth + 2 x sum(layerWidths))
    //              = 0.035 x 0.035 x (0.010 + 2 x (3 x 1.062 + 2 x 0.025) mm)
    //              = 0.001225 x 0.016472                        = 2.017820e-05
    //
    // NOT fixed by the fixture change: the nine entries pinned at {false, 0.0}
    // below are dead for their own reasons (no DcResistance on the wound coil,
    // empty vectors, etc.), not because the reference could not be wound. Making
    // the design windable did not revive any of them.
    // ---------------------------------------------------------------
    {"AREA_PRODUCT",                                    {true,  1.8374733304713626e-08}},
    // ENERGY_STORED / TEMPERATURE / SATURATION refreshed 2026-06-16 (ABT #10).
    // The reference 40+20-turn fixture has no isolation sides set, so the
    // 2026-06 saturation rework now classifies it as an INDUCTOR
    // (5000909f "classify all-one-isolation-side magnetics as inductors")
    // and gates it by gap-aware saturation current (60fe7c79). An ungapped
    // E35 ferrite at ~1.73 A peak has isat << margin·ipeak, so SATURATION
    // flips valid true→false (was {true, 0.1185}). ENERGY_STORED (-0.4 %) and
    // TEMPERATURE (-1.4 %) drift from the same reclassified flux/B recompute.
    {"ENERGY_STORED",                                   {true,  0.00018781148528340766}},
    {"ESTIMATED_COST",                                  {true,  1.7774692926005085}},
    // COST re-derived 2026-09-30: the filter now scores the unit cost in US$ at 1000 pieces (CostBasis) instead
    // of numberLayers + wire relative cost (7.0). By hand, for this fixture:
    //   core      E 35 Ve 8.07 cm3 x 3C97 4800 kg/m3 = 38.74 g; MnZn law 21.03 x 0.03874^0.6407 = US$2.620
    //   conductor 60 turns, 3.004 m (50 mm/turn round the 10 x 10 mm centre leg) of 1.00 mm copper
    //             = 21.09 g x US$14.737/kg                                          = US$0.311
    //   labour    60 x 0.105 min / 0.5 = 12.6 min at US$7.27/h                      = US$1.527
    //   total                                                                       = US$4.457
    {"COST",                                            {true,  4.45732988417687}},
    // CORE_AND_DC_LOSSES returns (false, 0) here because the reference
    // magnetic has no DcResistance computed on the wound coil → ohmic
    // losses model rejects it. This locks the silent invalid-fallthrough
    // behaviour — Phase 1 will likely flip this to throw.
    {"CORE_AND_DC_LOSSES",                              {false, 0.0}},
    {"CORE_DC_AND_SKIN_LOSSES",                         {false, 0.0}},
    {"LOSSES",                                          {false, 0.0}},
    // LOCKS BUG: LossesNoProximity calls calculate_skin_effect_losses twice
    // (MagneticFilter.cpp:1070–1071). Score is 0 here because earlier
    // model rejects → the double-call cannot be observed at this fixture.
    {"LOSSES_NO_PROXIMITY",                             {false, 0.0}},
    {"DIMENSIONS",                                      {true,  2.0178200000000003e-05}},
    // LOCKS BUG (ABT #787): MagneticFilterAreaNoParallels compares
    // `wire.get_maximum_outer_width() < section.get_dimensions()[0]` — STRICTLY.
    // A section holding exactly one layer is compacted to exactly one wire width,
    // so the comparison is an exact float equality and the winding is rejected.
    // Here winding 1 (20 turns, one 1.062 mm layer) trips it, so the filter says
    // false for a coil that fits with room to spare. Not fixture-dependent: any
    // single-layer section hits it. Locked as {false, 0.0}; when #787 lands this
    // flips back to {true, 0.0} and that is the signal, not a regression.
    {"AREA_NO_PARALLELS",                               {false, 0.0}},
    {"EFFECTIVE_RESISTANCE",                            {false, 0.0}},
    {"PROXIMITY_FACTOR",                                {false, 0.0}},
    {"TURNS_RATIOS",                                    {true,  0.0}},
    {"MAXIMUM_DIMENSIONS",                              {true,  0.0}},
    {"SATURATION",                                      {false, 0.0}},  // 2026-06-16 ABT #10: now isat-gated inductor, rejected (see ENERGY_STORED note)
    {"DC_CURRENT_DENSITY",                              {false, 0.0}},
    {"EFFECTIVE_CURRENT_DENSITY",                       {false, 0.0}},
    {"IMPEDANCE",                                       {true,  0.0}},
    {"MAGNETIZING_INDUCTANCE",                          {false, 0.0057629428244944511}},
    {"SKIN_LOSSES_DENSITY",                             {false, 0.0}},
    {"TEMPERATURE_RISE",                                {true,  25.0}},
    // LOCKS BUG: MagnetomotiveForce only push_backs in INSULATION branch
    // (MagneticFilter.cpp:2310); CONDUCTION branch keeps vector empty →
    // the (false, 0) below is the bug speaking.
    {"MAGNETOMOTIVE_FORCE",                             {false, 0.0}},
    // LEAKAGE_INDUCTANCE no longer returns DBL_MAX sentinel — the error path
    // is now a throws-contract test below (see TEST_CASE "LEAKAGE_INDUCTANCE
    // throws on missing turns description"). No snapshot entry needed.
    {"TEMPERATURE",                                     {true,  26.425983545124080}},  // 2026-08-17 ABT #785: +3.3 % — Round 1.00 has more DC resistance than Round 2.00. 2026-09-28 ABT #1454/#1459: 28.0504 → 26.4260 (−5.8 %), the winding now conducts to bobbin/core and the core has its real exterior surface
    {"TURN_COUNT",                                      {true,  2.1000000000000001}},
    // FRINGING_FACTOR returns score=1.0 on every non-crashing path
    // (MagneticFilter.cpp:2079, 2082, 2089, 2092). After fix A the factory
    // no longer SEGV's; the score is still trivially 1.0 — locked here so
    // any future fix that makes the score meaningful trips the snapshot.
    {"FRINGING_FACTOR",                                 {true,  1.0}},
};

Snapshot lookup_snapshot(const std::string& key) {
    auto it = kSnapshots.find(key);
    if (it == kSnapshots.end()) {
        return {false, std::numeric_limits<double>::quiet_NaN()};
    }
    return it->second;
}

// Evaluate a filter and CHECK against snapshot. Used by every behaviour test.
void evaluate_and_check(const std::string& key, MagneticFilters which) {
    settings.reset();
    auto magnetic = make_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(which, inputs);
    REQUIRE(filter != nullptr);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);

    auto snap = lookup_snapshot(key);
    if (kRegenerateBaselines || std::isnan(snap.score)) {
        // No snapshot yet — print and require finite so we can collect.
        UNSCOPED_INFO("BASELINE " << key
                      << " valid=" << (valid ? "true" : "false")
                      << " score=" << std::setprecision(17) << score);
        REQUIRE(std::isfinite(score));
    } else {
        check_or_print(key, valid, score, snap.valid, snap.score);
    }
}

// --- Benchmark helper -------------------------------------------------------
std::vector<Core> load_test_cores(size_t limit = std::numeric_limits<size_t>::max()) {
    auto path = OpenMagneticsTesting::get_test_data_path(
        std::source_location::current(), "test_cores.ndjson");
    std::ifstream in(path);
    std::string line;
    std::vector<Core> cores;
    while (std::getline(in, line) && cores.size() < limit) {
        try {
            json jf = json::parse(line);
            cores.emplace_back(jf, false, true, false);
        } catch (const CoreShapeNotFoundException&) {
            // skip — same convention as TestCoreAdviser
        } catch (const std::exception&) {
            // skip malformed lines (some fixtures fail parsing on this branch)
        }
    }
    return cores;
}

} // namespace

// =============================================================================
// BEHAVIOUR SNAPSHOTS — one TEST_CASE per filter
// =============================================================================
//
// Tag layout:
//   [magnetic-filter]                always
//   [characterisation][heavy]        always (heavy: excluded from smoke runs)
//   [<filter-name-lowercased>]       so a specific filter can be re-run alone
//

// ABT #785: the snapshots below are only meaningful if the reference coil is a winding that can
// physically exist. It was not for three years (Round 2.00 in E 35 needed ~60 turns of capacity
// and had ~36), and the whole table silently characterised an over-stuffed layout. Guard it, so a
// future fixture edit fails HERE with an obvious reason rather than as a drifting snapshot.
TEST_CASE("MagneticFilter reference fixture is windable (ABT #785)",
          "[magnetic-filter][characterisation][fixture]") {
    settings.reset();
    auto magnetic = make_reference_magnetic();
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_sections_description());
    REQUIRE(coil.get_layers_description());
    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.get_turns_description().value().size() == 60);
    INFO("the reference design must fit its winding window — see the wire note on "
         "make_reference_magnetic()");
    REQUIRE(coil.are_sections_and_layers_fitting());
}

TEST_CASE("MagneticFilter AREA_PRODUCT snapshot",
          "[magnetic-filter][characterisation][heavy][area-product]") {
    evaluate_and_check("AREA_PRODUCT", MagneticFilters::AREA_PRODUCT);
}

TEST_CASE("MagneticFilter ENERGY_STORED snapshot",
          "[magnetic-filter][characterisation][heavy][energy-stored]") {
    evaluate_and_check("ENERGY_STORED", MagneticFilters::ENERGY_STORED);
}

TEST_CASE("MagneticFilter ESTIMATED_COST snapshot",
          "[magnetic-filter][characterisation][heavy][estimated-cost]") {
    evaluate_and_check("ESTIMATED_COST", MagneticFilters::ESTIMATED_COST);
}

TEST_CASE("MagneticFilter COST snapshot",
          "[magnetic-filter][characterisation][heavy][cost]") {
    evaluate_and_check("COST", MagneticFilters::COST);
}

TEST_CASE("MagneticFilter CORE_AND_DC_LOSSES snapshot",
          "[magnetic-filter][characterisation][heavy][core-and-dc-losses]") {
    evaluate_and_check("CORE_AND_DC_LOSSES", MagneticFilters::CORE_AND_DC_LOSSES);
}

TEST_CASE("MagneticFilter CORE_DC_AND_SKIN_LOSSES snapshot",
          "[magnetic-filter][characterisation][heavy][core-dc-and-skin-losses]") {
    evaluate_and_check("CORE_DC_AND_SKIN_LOSSES", MagneticFilters::CORE_DC_AND_SKIN_LOSSES);
}

TEST_CASE("MagneticFilter LOSSES snapshot",
          "[magnetic-filter][characterisation][heavy][losses]") {
    evaluate_and_check("LOSSES", MagneticFilters::LOSSES);
}

TEST_CASE("MagneticFilter LOSSES_NO_PROXIMITY snapshot",
          "[magnetic-filter][characterisation][heavy][losses-no-proximity]") {
    // LOCKS BUG: MagneticFilterLossesNoProximity calls
    // calculate_skin_effect_losses() TWICE (MagneticFilter.cpp:1070–1071).
    // The captured snapshot reflects this double-count; once fixed, snapshot
    // must be regenerated.
    evaluate_and_check("LOSSES_NO_PROXIMITY", MagneticFilters::LOSSES_NO_PROXIMITY);
}

TEST_CASE("MagneticFilter DIMENSIONS snapshot",
          "[magnetic-filter][characterisation][heavy][dimensions]") {
    evaluate_and_check("DIMENSIONS", MagneticFilters::DIMENSIONS);
}

TEST_CASE("MagneticFilter CORE_MINIMUM_IMPEDANCE snapshot",
          "[magnetic-filter][characterisation][heavy][core-minimum-impedance]") {
    // CONTRACT: CORE_MINIMUM_IMPEDANCE requires the design requirements to
    // carry a minimum-impedance spec. With our reference Inputs (a regular
    // power-converter setup, no impedance requirement) the filter throws —
    // this is the correct loud-failure behaviour.
    settings.reset();
    auto magnetic = make_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::CORE_MINIMUM_IMPEDANCE, inputs);
    REQUIRE_THROWS(filter->evaluate_magnetic(&magnetic, &inputs));
}

TEST_CASE("MagneticFilter AREA_NO_PARALLELS snapshot",
          "[magnetic-filter][characterisation][heavy][area-no-parallels]") {
    evaluate_and_check("AREA_NO_PARALLELS", MagneticFilters::AREA_NO_PARALLELS);
}

TEST_CASE("MagneticFilter AREA_WITH_PARALLELS snapshot",
          "[magnetic-filter][characterisation][heavy][area-with-parallels]") {
    // CONTRACT: AREA_WITH_PARALLELS needs the coil's conducting area, which
    // requires a fully delimited+compacted winding. Our reference magnetic
    // is only wound (not delimit_and_compact'd) so the filter throws
    // CoilNotProcessedException. This locks the precondition.
    settings.reset();
    auto magnetic = make_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::AREA_WITH_PARALLELS, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter EFFECTIVE_RESISTANCE snapshot",
          "[magnetic-filter][characterisation][heavy][effective-resistance]") {
    evaluate_and_check("EFFECTIVE_RESISTANCE", MagneticFilters::EFFECTIVE_RESISTANCE);
}

TEST_CASE("MagneticFilter PROXIMITY_FACTOR snapshot",
          "[magnetic-filter][characterisation][heavy][proximity-factor]") {
    evaluate_and_check("PROXIMITY_FACTOR", MagneticFilters::PROXIMITY_FACTOR);
}

TEST_CASE("MagneticFilter TURNS_RATIOS snapshot",
          "[magnetic-filter][characterisation][heavy][turns-ratios]") {
    evaluate_and_check("TURNS_RATIOS", MagneticFilters::TURNS_RATIOS);
}

TEST_CASE("MagneticFilter MAXIMUM_DIMENSIONS snapshot",
          "[magnetic-filter][characterisation][heavy][maximum-dimensions]") {
    evaluate_and_check("MAXIMUM_DIMENSIONS", MagneticFilters::MAXIMUM_DIMENSIONS);
}

TEST_CASE("MagneticFilter SATURATION snapshot",
          "[magnetic-filter][characterisation][heavy][saturation]") {
    evaluate_and_check("SATURATION", MagneticFilters::SATURATION);
}

TEST_CASE("MagneticFilter DC_CURRENT_DENSITY snapshot",
          "[magnetic-filter][characterisation][heavy][dc-current-density]") {
    evaluate_and_check("DC_CURRENT_DENSITY", MagneticFilters::DC_CURRENT_DENSITY);
}

TEST_CASE("MagneticFilter EFFECTIVE_CURRENT_DENSITY snapshot",
          "[magnetic-filter][characterisation][heavy][effective-current-density]") {
    evaluate_and_check("EFFECTIVE_CURRENT_DENSITY", MagneticFilters::EFFECTIVE_CURRENT_DENSITY);
}

TEST_CASE("MagneticFilter IMPEDANCE snapshot",
          "[magnetic-filter][characterisation][heavy][impedance]") {
    evaluate_and_check("IMPEDANCE", MagneticFilters::IMPEDANCE);
}

TEST_CASE("MagneticFilter MAGNETIZING_INDUCTANCE snapshot",
          "[magnetic-filter][characterisation][heavy][magnetizing-inductance]") {
    evaluate_and_check("MAGNETIZING_INDUCTANCE", MagneticFilters::MAGNETIZING_INDUCTANCE);
}

TEST_CASE("MagneticFilter SKIN_LOSSES_DENSITY snapshot",
          "[magnetic-filter][characterisation][heavy][skin-losses-density]") {
    evaluate_and_check("SKIN_LOSSES_DENSITY", MagneticFilters::SKIN_LOSSES_DENSITY);
}

TEST_CASE("MagneticFilter FRINGING_FACTOR snapshot",
          "[magnetic-filter][characterisation][heavy][fringing-factor]") {
    // Phase 1 fix applied: MagneticFilter::factory(FRINGING_FACTOR) now
    // routes through the Inputs-aware ctor at MagneticFilter.cpp:128 so
    // _reluctanceModel is initialised. Previously the factory used the
    // default ctor, left the model null, and segfaulted on first call.
    // NOTE: evaluate_magnetic still returns 1.0 on most code paths
    // (cpp:2079, 2082, 2089, 2092). This snapshot locks the current
    // (possibly trivial) score so a future fix that makes the score
    // meaningful will surface as a snapshot delta.
    evaluate_and_check("FRINGING_FACTOR", MagneticFilters::FRINGING_FACTOR);
}

TEST_CASE("MagneticFilter VOLUME snapshot",
          "[magnetic-filter][characterisation][heavy][volume]") {
    // CONTRACT: VOLUME iterates the wound coil's turns to compute the
    // bounding volume. Our reference magnetic doesn't have a fully processed
    // turns_description so it throws — Phase 1 should keep this loud.
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::VOLUME, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter AREA snapshot",
          "[magnetic-filter][characterisation][heavy][area]") {
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::AREA, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter HEIGHT snapshot",
          "[magnetic-filter][characterisation][heavy][height]") {
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::HEIGHT, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter TEMPERATURE_RISE snapshot",
          "[magnetic-filter][characterisation][heavy][temperature-rise]") {
    evaluate_and_check("TEMPERATURE_RISE", MagneticFilters::TEMPERATURE_RISE);
}

TEST_CASE("MagneticFilter LOSSES_TIMES_VOLUME snapshot",
          "[magnetic-filter][characterisation][heavy][losses-times-volume]") {
    // Same precondition as VOLUME — wound + delimited coil required.
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::LOSSES_TIMES_VOLUME, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter VOLUME_TIMES_TEMPERATURE_RISE snapshot",
          "[magnetic-filter][characterisation][heavy][volume-times-temperature-rise]") {
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE snapshot",
          "[magnetic-filter][characterisation][heavy][losses-times-volume-times-temperature-rise]") {
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(
        MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter MAGNETOMOTIVE_FORCE snapshot",
          "[magnetic-filter][characterisation][heavy][magnetomotive-force]") {
    // LOCKS BUG: MagnetomotiveForce only push_backs in the INSULATION branch
    // (MagneticFilter.cpp:2310); CONDUCTION branch updates a dead local. The
    // E35 / 100 µH reference is a conduction case; snapshot captures whatever
    // current behaviour produces (likely score = 0 due to empty vector → mean).
    evaluate_and_check("MAGNETOMOTIVE_FORCE", MagneticFilters::MAGNETOMOTIVE_FORCE);
}

TEST_CASE("MagneticFilter LEAKAGE_INDUCTANCE throws on missing turns description",
          "[magnetic-filter][characterisation][heavy][leakage-inductance]") {
    // Phase 1: removed the catch-all that returned (false, DBL_MAX) sentinel
    // in MagneticFilterLeakageInductance::evaluate_magnetic. The reference
    // 2-winding fixture has no fully processed turns/sections so the leakage
    // model throws CoilNotProcessedException — that exception must now
    // propagate instead of being swallowed and rewritten as DBL_MAX.
    settings.reset();
    auto magnetic = make_unwound_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto filter = MagneticFilter::factory(MagneticFilters::LEAKAGE_INDUCTANCE, inputs);
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs),
                      CoilNotProcessedException);
}

TEST_CASE("MagneticFilter TEMPERATURE snapshot",
          "[magnetic-filter][characterisation][heavy][temperature][smoke-test]") {
    // LOCKS BUG: MagneticFilter::factory(TEMPERATURE) goes through the
    // default ctor (MagneticFilter.cpp:158) leaving _coreLossesModel null,
    // which dereferences at cpp:2420. The default ctor path is therefore a
    // CRASH. We construct via the explicit 2-arg ctor to keep this test alive,
    // and assert that the factory variant THROWS or crashes — proving the bug.
    settings.reset();
    // Fixture correction: the temperature filter now solves the copper of every winding the
    // inputs excite, and a coil with a winding the inputs do not excite has no current to give
    // it (MagneticFilter::with_every_winding throws). The reference inputs excite one winding,
    // so the filter gets the reference magnetic's first winding alone (40 turns), wound.
    auto magnetic = make_reference_magnetic();
    {
        auto coil = magnetic.get_coil();
        auto windings = coil.get_functional_description();
        windings.resize(1);
        coil.set_functional_description(windings);
        coil.unwind();
        coil.wind();
        magnetic.set_coil(coil);
    }
    auto inputs = make_reference_inputs();
    MagneticFilterTemperature filter(inputs, 130.0);
    auto [valid, score] = filter.evaluate_magnetic(&magnetic, &inputs);

    auto snap = lookup_snapshot("TEMPERATURE");
    if (kRegenerateBaselines || std::isnan(snap.score)) {
        UNSCOPED_INFO("BASELINE TEMPERATURE valid=" << (valid ? "true" : "false")
                      << " score=" << std::setprecision(17) << score);
        REQUIRE(std::isfinite(score));
    } else {
        // 5% relative bracket: the coreOnly TEMPERATURE score legitimately shifts with
        // core-model refinements (column footprint area, conductivity resolution, radiation).
        check_or_print("TEMPERATURE", valid, score, snap.valid, snap.score, 0.05);
    }
}

TEST_CASE("MagneticFilter TURN_COUNT snapshot",
          "[magnetic-filter][characterisation][heavy][turn-count]") {
    evaluate_and_check("TURN_COUNT", MagneticFilters::TURN_COUNT);
}

// =============================================================================
// DATASHEET_LIMITS (ABT #19)
// =============================================================================
// Gate catalogue parts by their OWN datasheet electrical limits. The filter
// does not apply (applies_to) to designed/custom magnetics or to parts that
// publish no limit it checks; it scores UTILISATION (operating/limit, lower is
// better), not headroom. These tests build catalogue-style fixtures
// by attaching a DatasheetInfo with one MagneticDatasheetElectrical entry, and
// drive operating values straight into the excitations' processed signals (no
// physics path). The worked-example numbering matches the ABT #19 handoff.

namespace {

struct WindingExcitationSpec {
    std::optional<double> currentRms;
    std::optional<double> currentPeak;
    std::optional<double> voltageRms;
    std::optional<double> voltageDc;
};

// Build a single-operating-point Inputs whose per-winding excitations carry the
// requested processed RMS/peak/offset values directly.
OpenMagnetics::Inputs make_datasheet_inputs(const std::vector<WindingExcitationSpec>& windings) {
    OpenMagnetics::Inputs inputs;
    OperatingPoint operatingPoint;
    for (const auto& w : windings) {
        OperatingPointExcitation excitation;
        excitation.set_frequency(100000);
        if (w.currentRms || w.currentPeak) {
            SignalDescriptor current;
            ProcessedWaveform processed;
            if (w.currentRms) processed.set_rms(w.currentRms.value());
            if (w.currentPeak) processed.set_peak(w.currentPeak.value());
            current.set_processed(processed);
            excitation.set_current(current);
        }
        if (w.voltageRms || w.voltageDc) {
            SignalDescriptor voltage;
            ProcessedWaveform processed;
            if (w.voltageRms) processed.set_rms(w.voltageRms.value());
            if (w.voltageDc) processed.set_offset(w.voltageDc.value());
            voltage.set_processed(processed);
            excitation.set_voltage(voltage);
        }
        operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
    }
    inputs.get_mutable_operating_points().push_back(operatingPoint);
    return inputs;
}

// Attach a catalogue datasheet (one electrical entry) to the reference magnetic.
OpenMagnetics::Magnetic make_datasheet_magnetic(const MagneticDatasheetElectrical& electrical) {
    auto magnetic = make_reference_magnetic();
    DatasheetInfo datasheetInfo;
    datasheetInfo.set_electrical(std::vector<MagneticDatasheetElectrical>{electrical});
    MagneticManufacturerInfo manufacturerInfo;
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

}  // namespace

TEST_CASE("MagneticFilter DATASHEET_LIMITS rejects over-rated current (the bug)",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #1: ratedCurrents=[1.0], operating RMS 1.5 A ⇒ reject.
    settings.reset();
    MagneticDatasheetElectrical electrical;
    electrical.set_rated_currents(std::vector<double>{1.0});
    auto magnetic = make_datasheet_magnetic(electrical);
    auto inputs = make_datasheet_inputs({{1.5, std::nullopt, std::nullopt, std::nullopt}});
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE(valid == false);
    REQUIRE_THAT(score, WithinRel(1.5, kRelTol));  // utilisation above 1: over the limit
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS passes comfortable current and scores utilisation",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #2: ratedCurrents=[2.0], operating RMS 1.0 A ⇒ valid, utilisation 0.5.
    settings.reset();
    MagneticDatasheetElectrical electrical;
    electrical.set_rated_currents(std::vector<double>{2.0});
    auto magnetic = make_datasheet_magnetic(electrical);
    auto inputs = make_datasheet_inputs({{1.0, std::nullopt, std::nullopt, std::nullopt}});
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE(valid == true);
    REQUIRE_THAT(score, WithinRel(0.5, kRelTol));
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS is a no-op for custom magnetics",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #3: no manufacturerInfo/datasheetInfo ⇒ the filter does not apply, so the
    // adviser records no score for it: zero effect on designed parts even when
    // the operating current is enormous. Evaluating it anyway is a caller error.
    settings.reset();
    auto magnetic = make_reference_magnetic();  // no manufacturer info
    auto inputs = make_datasheet_inputs({{999.0, 999.0, std::nullopt, std::nullopt}});
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    REQUIRE_FALSE(filter->applies_to(&magnetic));
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs), InvalidInputException);
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS skips limits the datasheet omits",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #4: datasheet present but only publishes inductance (no current/voltage
    // limit) ⇒ nothing to gate on ⇒ the filter does not apply.
    settings.reset();
    MagneticDatasheetElectrical electrical;
    DimensionWithTolerance inductance;
    inductance.set_nominal(100e-6);
    electrical.set_inductance(inductance);
    auto magnetic = make_datasheet_magnetic(electrical);
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    REQUIRE_FALSE(filter->applies_to(&magnetic));
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS enforces rated AC voltage",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #5: ratedVoltageAc=250, operating 400 V RMS ⇒ reject.
    settings.reset();
    MagneticDatasheetElectrical electrical;
    electrical.set_rated_voltage_ac(250.0);
    auto magnetic = make_datasheet_magnetic(electrical);
    auto inputs = make_datasheet_inputs({{std::nullopt, std::nullopt, 400.0, std::nullopt}});
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE(valid == false);
    REQUIRE_THAT(score, WithinRel(400.0 / 250.0, kRelTol));
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS enforces saturation current peak",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #6: saturationCurrentPeak=1.5, operating peak 2.0 A ⇒ reject.
    settings.reset();
    MagneticDatasheetElectrical electrical;
    electrical.set_saturation_current_peak(1.5);
    auto magnetic = make_datasheet_magnetic(electrical);
    auto inputs = make_datasheet_inputs({{std::nullopt, 2.0, std::nullopt, std::nullopt}});
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);
    auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
    REQUIRE(valid == false);
    REQUIRE_THAT(score, WithinRel(2.0 / 1.5, kRelTol));
}

TEST_CASE("MagneticFilter DATASHEET_LIMITS gates each winding against its own rated current",
          "[magnetic-filter][datasheet-limits][smoke-test]") {
    // #7: ratedCurrents=[1.0, 3.0] (multi-entry ⇒ per winding).
    settings.reset();
    MagneticDatasheetElectrical electrical;
    electrical.set_rated_currents(std::vector<double>{1.0, 3.0});
    auto magnetic = make_datasheet_magnetic(electrical);
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);

    SECTION("both windings within their own rating ⇒ valid") {
        auto inputs = make_datasheet_inputs({
            {0.5, std::nullopt, std::nullopt, std::nullopt},
            {2.0, std::nullopt, std::nullopt, std::nullopt},
        });
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        REQUIRE(valid == true);
    }

    SECTION("second winding exceeds its own rating ⇒ invalid") {
        auto inputs = make_datasheet_inputs({
            {0.5, std::nullopt, std::nullopt, std::nullopt},
            {4.0, std::nullopt, std::nullopt, std::nullopt},
        });
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        REQUIRE(valid == false);
    }
}

// =============================================================================
// FACTORY CONTRACT TESTS
// =============================================================================
// Filters that require Inputs at construction must throw if none is provided.

// ABT #357 phase 1: the WE-MAPI catalogue corpus — 183 parts pulled from RedExpert's public
// JSON API (scripts/pull_we_mapi.py, schema-validated at generation, provenance carried per
// record). Pins: (a) the corpus stays complete and parseable; (b) DATASHEET_LIMITS gates on
// REAL vendor records end-to-end — 20% over the part's own rated current rejects, 50% under
// passes with headroom.
TEST_CASE("MagneticFilter DATASHEET_LIMITS gates the WE-MAPI corpus (ABT #357 phase 1)",
          "[magnetic-filter][datasheet-limits][we-mapi]") {
    settings.reset();
    auto path = std::filesystem::path{std::source_location::current().file_name()}
                    .parent_path().append("testData").append("we_mapi_datasheet_stubs.ndjson");
    std::ifstream in(path);
    REQUIRE(in.good());
    std::string line;
    size_t partCount = 0;
    std::optional<nlohmann::json> firstStub;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        nlohmann::json stub = nlohmann::json::parse(line);
        auto& electricalJson = stub.at("manufacturerInfo").at("datasheetInfo").at("electrical").at(0);
        REQUIRE(electricalJson.contains("inductance"));
        REQUIRE(electricalJson.contains("saturationCurrents"));
        REQUIRE(electricalJson.contains("ratedCurrents"));
        if (!firstStub) firstStub = stub;
        ++partCount;
    }
    CHECK(partCount == 183);

    DatasheetInfo datasheetInfo;
    from_json((*firstStub)["manufacturerInfo"]["datasheetInfo"], datasheetInfo);
    REQUIRE(datasheetInfo.get_electrical());
    auto electrical = datasheetInfo.get_electrical().value()[0];
    REQUIRE(electrical.get_rated_currents());
    double ratedCurrent = electrical.get_rated_currents().value()[0];
    auto magnetic = make_datasheet_magnetic(electrical);
    auto filter = MagneticFilter::factory(MagneticFilters::DATASHEET_LIMITS);

    auto overInputs = make_datasheet_inputs({{ratedCurrent * 1.2, std::nullopt, std::nullopt, std::nullopt}});
    auto [overValid, overScore] = filter->evaluate_magnetic(&magnetic, &overInputs);
    CHECK(overValid == false);

    auto underInputs = make_datasheet_inputs({{ratedCurrent * 0.5, std::nullopt, std::nullopt, std::nullopt}});
    auto [underValid, underScore] = filter->evaluate_magnetic(&magnetic, &underInputs);
    CHECK(underValid == true);
    CHECK(underScore > 0);
}

TEST_CASE("MagneticFilter factory requires Inputs for AREA_PRODUCT",
          "[magnetic-filter][characterisation][factory-contract][smoke-test]") {
    REQUIRE_THROWS_AS(MagneticFilter::factory(MagneticFilters::AREA_PRODUCT, std::nullopt),
                      InvalidInputException);
}

TEST_CASE("MagneticFilter factory requires Inputs for ENERGY_STORED",
          "[magnetic-filter][characterisation][factory-contract][smoke-test]") {
    REQUIRE_THROWS_AS(MagneticFilter::factory(MagneticFilters::ENERGY_STORED, std::nullopt),
                      InvalidInputException);
}

TEST_CASE("MagneticFilter factory requires Inputs for ESTIMATED_COST",
          "[magnetic-filter][characterisation][factory-contract][smoke-test]") {
    REQUIRE_THROWS_AS(MagneticFilter::factory(MagneticFilters::ESTIMATED_COST, std::nullopt),
                      InvalidInputException);
}

TEST_CASE("MagneticFilter factory requires Inputs for CORE_AND_DC_LOSSES",
          "[magnetic-filter][characterisation][factory-contract][smoke-test]") {
    REQUIRE_THROWS_AS(MagneticFilter::factory(MagneticFilters::CORE_AND_DC_LOSSES, std::nullopt),
                      InvalidInputException);
}

TEST_CASE("MagneticFilter factory requires Inputs for CORE_DC_AND_SKIN_LOSSES",
          "[magnetic-filter][characterisation][factory-contract][smoke-test]") {
    REQUIRE_THROWS_AS(MagneticFilter::factory(MagneticFilters::CORE_DC_AND_SKIN_LOSSES, std::nullopt),
                      InvalidInputException);
}

// =============================================================================
// BENCHMARKS — opt-in via [!benchmark]
// =============================================================================
// Tiered datasets: 10 / 100 / full (~4800 cores after parse-skips).

TEST_CASE("Benchmark MagneticFilter SATURATION (10 cores)",
          "[!benchmark][benchmark-saturation]") {
    settings.reset();
    auto inputs = make_reference_inputs();
    auto cores = load_test_cores(10);
    auto filter = MagneticFilter::factory(MagneticFilters::SATURATION, inputs);

    BENCHMARK("evaluate SATURATION on 10 cores") {
        double acc = 0.0;
        for (auto& core : cores) {
            OpenMagneticsTesting::QuickMagneticConfig cfg;
            cfg.numberTurns = {40};
            cfg.numberParallels = {1};
            cfg.coreShapeName = core.get_shape_name();
            cfg.wireNames = {"Round 2.00 - Grade 1"};
            OpenMagnetics::Magnetic m;
            m.set_core(core);
            auto coil = OpenMagnetics::Coil::create_quick_coil(
                core.get_shape_name(), {40}, {1},
                {find_wire_by_name("Round 2.00 - Grade 1")});
            m.set_coil(coil);
            auto [v, s] = filter->evaluate_magnetic(&m, &inputs);
            acc += s;
        }
        return acc;
    };
}

TEST_CASE("Benchmark MagneticFilter CORE_AND_DC_LOSSES (100 cores)",
          "[!benchmark][benchmark-core-and-dc-losses]") {
    settings.reset();
    auto inputs = make_reference_inputs();
    auto cores = load_test_cores(100);
    auto filter = MagneticFilter::factory(MagneticFilters::CORE_AND_DC_LOSSES, inputs);

    BENCHMARK("evaluate CORE_AND_DC_LOSSES on 100 cores") {
        double acc = 0.0;
        for (auto& core : cores) {
            OpenMagnetics::Magnetic m;
            m.set_core(core);
            auto coil = OpenMagnetics::Coil::create_quick_coil(
                core.get_shape_name(), {40}, {1},
                {find_wire_by_name("Round 2.00 - Grade 1")});
            m.set_coil(coil);
            auto [v, s] = filter->evaluate_magnetic(&m, &inputs);
            acc += s;
        }
        return acc;
    };
}

TEST_CASE("Benchmark MagneticFilter AREA_PRODUCT (full DB)",
          "[!benchmark][benchmark-area-product]") {
    settings.reset();
    auto inputs = make_reference_inputs();
    auto cores = load_test_cores();
    auto filter = MagneticFilter::factory(MagneticFilters::AREA_PRODUCT, inputs);

    BENCHMARK("evaluate AREA_PRODUCT on full test_cores.ndjson") {
        double acc = 0.0;
        for (auto& core : cores) {
            OpenMagnetics::Magnetic m;
            m.set_core(core);
            auto coil = OpenMagnetics::Coil::create_quick_coil(
                core.get_shape_name(), {40}, {1},
                {find_wire_by_name("Round 2.00 - Grade 1")});
            m.set_coil(coil);
            auto [v, s] = filter->evaluate_magnetic(&m, &inputs);
            acc += s;
        }
        return acc;
    };
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Filter                | n     | mean         | notes
// -----------+-----------------------+-------+--------------+--------------------
// 2026-05-19 | SATURATION            | 10    | 43.5 ms      | initial commit
// 2026-05-19 | CORE_AND_DC_LOSSES    | 100   | TBD          | initial commit
// 2026-05-19 | AREA_PRODUCT          | full  | TBD          | initial commit
//
// Note: SATURATION/10-cores at ~43 ms/iter is dominated by the per-iteration
// `Coil::create_quick_coil` rebuild (winding, not the filter itself). When
// optimising the filter, also extract a cached-coil benchmark to isolate the
// hot path. See "Performance" section of the Phase 6 plan.
// =============================================================================

// ABT #19 wiring: the DATASHEET_LIMITS filter existed and was factory-registered but was
// never part of the DEFAULT catalogue flow — so catalogue parts were only gated by
// MKF-simulated quantities and their own published limits were silently ignored. Pin the
// flow membership so it cannot fall out again.
TEST_CASE("Catalogue filter flow includes DATASHEET_LIMITS (ABT #19 wiring)",
          "[adviser][magnetic-filter][datasheet]") {
    OpenMagnetics::MagneticAdviser adviser;
    bool datasheetLimitsPresent = false;
    for (auto& operation : adviser._defaultCatalogueMagneticFilterFlow) {
        if (operation.get_filter() == MagneticFilters::DATASHEET_LIMITS) {
            datasheetLimitsPresent = true;
        }
    }
    CHECK(datasheetLimitsPresent);
}

TEST_CASE("Processed-only excitation advises without a bad optional access (ABT #825)",
          "[adviser][magnetic-filter]") {
    // ABT #825: an excitation given as {label, dutyCycle, offset, peakToPeak} — the ordinary way
    // to describe one without a waveform — short-circuits calculate_basic_processed_data, which is
    // the only place processed.peak is derived. So the processed block arrives with an rms, a thd
    // and harmonics but NO peak, and MagneticEnergy::calculate_required_magnetic_energy read it as
    // get_peak().value(). Via MagneticFilterEnergyStored that sits on the core adviser's critical
    // path, so advising ANY design described this way died with "bad optional access", naming
    // neither the field nor the excitation. A plain inductor with no turns ratio reproduces it.
    //
    // Note what is NOT asserted here: that processed.peak gets populated. Filling it in at the
    // source was the first fix and it was wrong — peak is read behind `if (get_peak())` in half a
    // dozen filters, and giving it a value where they had learned to expect none changed which
    // branch they take, emptying the candidate list for a flyback design
    // (Test_CoreAdviser_Flyback_From_Frontend_Inputs, deterministic across two full-suite runs).
    // The consumer that needs the number derives it instead, so this pins the BEHAVIOUR — the
    // energy filter builds — not the representation.
    json inputsJson = json();
    inputsJson["designRequirements"]["magnetizingInductance"]["nominal"] = 100e-6;
    inputsJson["designRequirements"]["turnsRatios"] = json::array();

    json excitation = json();
    excitation["frequency"] = 100000;
    excitation["current"]["processed"]["dutyCycle"] = 0.5;
    excitation["current"]["processed"]["label"] = "Triangular";
    excitation["current"]["processed"]["offset"] = 0;
    excitation["current"]["processed"]["peakToPeak"] = 10;

    json operatingPoint = json();
    operatingPoint["name"] = "Nominal";
    operatingPoint["conditions"]["ambientTemperature"] = 25;
    operatingPoint["excitationsPerWinding"] = json::array({excitation});
    inputsJson["operatingPoints"] = json::array({operatingPoint});

    OpenMagnetics::Inputs inputs(inputsJson);

    // The processed block still has no peak — that is the input shape this ticket is about.
    auto processedCurrent = inputs.get_operating_point(0).get_excitations_per_winding()[0]
                                  .get_current().value().get_processed().value();
    CHECK_FALSE(processedCurrent.get_peak());
    CHECK_THAT(processedCurrent.get_peak_to_peak().value(), Catch::Matchers::WithinRel(10.0, 1e-9));

    // ...and the energy filter builds anyway, deriving the peak it needs from the waveform.
    REQUIRE_NOTHROW(OpenMagnetics::MagneticFilterEnergyStored(inputs, {}));

    // The derived required energy is E = L*Ipk^2/2 with Ipk ~ 5 A for a 10 A pk-pk triangular
    // at no offset: 100e-6 * 25 / 2 = 1.25 mJ.
    OpenMagnetics::MagneticEnergy magneticEnergy;
    double requiredEnergy = OpenMagnetics::resolve_dimensional_values(
        magneticEnergy.calculate_required_magnetic_energy(inputs));
    CHECK_THAT(requiredEnergy, Catch::Matchers::WithinRel(1.25e-3, 0.05));
}

// ABT #1410: design mode holds maximumDimensions as the part can be mounted. A bobbin-wound
// part keeps its height (laying it down needs a horizontal bobbin the design doesn't have) but
// may be turned on the board (width and depth swap). A toroid mounts flat or on edge on a
// vertical header, so it fits if either orientation does.
TEST_CASE("Design-mode envelope fit: height strict, in-plane rotation, toroid flat or on edge",
          "[adviser][magnetic-filter][abt-1410]") {
    auto inputs_with = [](std::optional<double> width, std::optional<double> height, std::optional<double> depth) {
        MaximumDimensions maximumDimensions;
        if (width) maximumDimensions.set_width(width.value());
        if (height) maximumDimensions.set_height(height.value());
        if (depth) maximumDimensions.set_depth(depth.value());
        DesignRequirements designRequirements;
        designRequirements.set_maximum_dimensions(maximumDimensions);
        OpenMagnetics::Inputs inputs;
        inputs.set_design_requirements(designRequirements);
        return inputs;
    };

    auto eCore = Core::create_quick_core("E 65/32/27", "N87");       // 65.15 x 65.0 x 27.0 mm
    auto toroid = Core::create_quick_core("T 58/35/15", "N87");      // 58.04 x 58.04 x 14.9 mm

    SECTION("a 30 mm height does not admit an E 65 lying on its 27 mm side") {
        CHECK_FALSE(MagneticFilterMaximumDimensions::core_fits(eCore, inputs_with(std::nullopt, 0.030, std::nullopt)));
        CHECK(MagneticFilterMaximumDimensions::core_fits(eCore, inputs_with(std::nullopt, 0.066, std::nullopt)));
    }
    SECTION("turning the part on the board is allowed: width and depth may swap") {
        CHECK(MagneticFilterMaximumDimensions::core_fits(eCore, inputs_with(0.030, 0.066, 0.070)));
        CHECK_FALSE(MagneticFilterMaximumDimensions::core_fits(eCore, inputs_with(0.030, 0.066, 0.060)));
    }
    SECTION("a toroid fits flat (height = axial length) or on edge (height = OD)") {
        CHECK(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(std::nullopt, 0.020, std::nullopt)));
        CHECK(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(0.060, 0.020, 0.060)));
        CHECK_FALSE(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(0.050, 0.020, 0.060)));
        CHECK_FALSE(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(std::nullopt, 0.014, std::nullopt)));
        // Only on edge fits: a 60 x 16 mm footprint is too narrow flat, 60 mm is tall enough standing.
        CHECK(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(0.060, 0.060, 0.016)));
        CHECK_FALSE(MagneticFilterMaximumDimensions::core_fits(toroid, inputs_with(0.060, 0.050, 0.016)));
    }
}

// ABT #1456: a catalogue built on a MHz grade (TP5H: Steinmetz fitted over 1-5 MHz only) advised at 100 kHz.
// Its loss filters now throw rather than extrapolate, and the catalogue path aborts after 8 identical throws
// in a row, which a catalogue holding many parts of one material produces. The always-on
// LOSS_MODEL_FREQUENCY_SPAN gate must drop those parts with a verdict instead: the run completes, returns
// the N87 part, and returns no TP5H part.
TEST_CASE("MagneticAdviser catalogue path drops parts whose material has no loss fit at the operating frequency",
          "[adviser][magnetic-adviser][abt-1456]") {
    settings.reset();
    clear_databases();

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 2, 0.5, 1);
    auto part = [&](const std::string& material, double gap, int64_t turns, const std::string& reference) {
        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15", OpenMagneticsTesting::get_ground_gap(gap),
                                                                 std::vector<int64_t>{turns}, 1, material);
        MAS::MagneticManufacturerInfo manufacturerInfo;
        manufacturerInfo.set_name("ABT 1456 test");
        manufacturerInfo.set_reference(reference);
        magnetic.set_manufacturer_info(manufacturerInfo);
        return OpenMagnetics::magnetic_autocomplete(magnetic);
    };

    auto tp5hPart = part("TP5H", 0.001, 20, "TP5H-0");
    auto n87Part = part("N87", 0.001, 20, "N87-0");
    auto spanFilter = MagneticFilter::factory(MagneticFilters::LOSS_MODEL_FREQUENCY_SPAN);
    CHECK(spanFilter->evaluate_magnetic(&tp5hPart, &inputs) == std::pair<bool, double>{false, 0.0});
    CHECK(spanFilter->evaluate_magnetic(&n87Part, &inputs) == std::pair<bool, double>{true, 1.0});

    std::vector<OpenMagnetics::Magnetic> catalogue;
    for (int i = 0; i < 10; ++i) {
        catalogue.push_back(part("TP5H", 0.0005 + 0.0001 * i, 15 + i, "TP5H-" + std::to_string(i)));
    }
    catalogue.push_back(n87Part);

    MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    // A strictly-required loss filter: without the gate, each TP5H part makes it throw the same
    // out-of-span message, and the 10 consecutive identical throws abort the run.
    std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::CORE_AND_DC_LOSSES, true, false, true, 1.0)};
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, catalogue, flow, 5, false));
    REQUIRE_FALSE(results.empty());
    for (auto& [mas, scoring] : results) {
        INFO(mas.get_magnetic().get_reference());
        CHECK(mas.get_magnetic().get_core().get_material_name() != "TP5H");
    }
    settings.reset();
}

// ABT #1652 (Alf): advisers NEVER extrapolate, even with Settings::allowMaterialDataExtrapolation on (the web
// turns it on for hand edits). The same TP5H-at-100 kHz catalogue as above, flag ON: the span gate must still
// drop every TP5H part, no extrapolated loss may be computed (no WARNING from module MaterialDataExtrapolation),
// and once the adviser is gone the flag works again, so the setup really had it on.
TEST_CASE("MagneticAdviser catalogue path keeps barring out-of-span parts with material extrapolation allowed",
          "[adviser][magnetic-adviser][abt-1652]") {
    settings.reset();
    clear_databases();
    Logger::getInstance().disableCollector();
    Logger::getInstance().enableCollector(LogLevel::WARNING);

    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 2, 0.5, 1);
    auto part = [&](const std::string& material, double gap, int64_t turns, const std::string& reference) {
        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15", OpenMagneticsTesting::get_ground_gap(gap),
                                                                 std::vector<int64_t>{turns}, 1, material);
        MAS::MagneticManufacturerInfo manufacturerInfo;
        manufacturerInfo.set_name("ABT 1652 test");
        manufacturerInfo.set_reference(reference);
        magnetic.set_manufacturer_info(manufacturerInfo);
        return OpenMagnetics::magnetic_autocomplete(magnetic);
    };
    std::vector<OpenMagnetics::Magnetic> catalogue;
    for (int i = 0; i < 10; ++i) {
        catalogue.push_back(part("TP5H", 0.0005 + 0.0001 * i, 15 + i, "TP5H-" + std::to_string(i)));
    }
    catalogue.push_back(part("N87", 0.001, 20, "N87-0"));

    settings.set_allow_material_data_extrapolation(true);
    Logger::getInstance().drainCollected();
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    {
        MagneticAdviser adviser;
        std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::CORE_AND_DC_LOSSES, true, false, true, 1.0)};
        REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, catalogue, flow, 11, false));
    }
    auto collected = Logger::getInstance().drainCollected();
    size_t extrapolationRecords = 0;
    for (auto& record : collected) {
        if (record.moduleOfOrigin == kMaterialDataExtrapolationModule) {
            UNSCOPED_INFO("extrapolated under the adviser: " << record.message);
            extrapolationRecords++;
        }
    }
    CHECK(extrapolationRecords == 0);
    REQUIRE_FALSE(results.empty());
    for (auto& [mas, scoring] : results) {
        INFO(mas.get_magnetic().get_reference());
        CHECK(mas.get_magnetic().get_core().get_material_name() != "TP5H");
    }

    // The adviser is gone: the caller's flag is back, and the same gate now lets TP5H through.
    auto spanFilter = MagneticFilter::factory(MagneticFilters::LOSS_MODEL_FREQUENCY_SPAN);
    CHECK(settings.get_allow_material_data_extrapolation());
    CHECK(spanFilter->evaluate_magnetic(&catalogue[0], &inputs).first);

    Logger::getInstance().disableCollector();
    settings.reset();
}

// ---------------------------------------------------------------------------
// COST: the unit cost in US$ at 1000 pieces = core + conductor + labour (see CostBasis).
// Every expected value below is derived by hand from the reference magnetic, not read back
// from the filter.
// ---------------------------------------------------------------------------
namespace {
double sum_of_turn_lengths(OpenMagnetics::Magnetic& magnetic) {
    double length = 0;
    auto turns = magnetic.get_coil().get_turns_description().value();
    for (const auto& turn : turns) {
        length += turn.get_length();
    }
    return length;
}
}  // namespace

TEST_CASE("MagneticFilter COST breakdown of the reference magnetic",
          "[magnetic-filter][cost]") {
    auto magnetic = make_reference_magnetic();
    auto cost = MagneticFilterCost().calculate_cost(magnetic);

    // Labour: 60 turns x 0.105 min / 0.5 winding share = 12.6 min at US$7.27/h.
    CHECK_THAT(cost.labour, WithinRel(12.6 / 60 * 7.27, 1e-12));

    // Conductor: Round 1.00 mm, conducting area pi x (0.5 mm)^2, copper 8940 kg/m3 at US$14.737/kg,
    // no litz premium, over the wound length of all 60 turns.
    double area = std::numbers::pi * 0.5e-3 * 0.5e-3;
    CHECK_THAT(cost.conductor, WithinRel(sum_of_turn_lengths(magnetic) * area * 8940 * 14.737, 1e-9));

    // Core: 3C97 is a MnZn ferrite, price = 21.03 x mass^0.6407, mass = Ve(E 35) x 4800 kg/m3.
    double mass = magnetic.get_core().get_effective_volume() * 4800;
    REQUIRE(cost.priced());
    CHECK_THAT(*cost.core, WithinRel(21.03 * std::pow(mass, 0.6407), 1e-9));

    CHECK_THAT(cost.total(), WithinRel(*cost.core + cost.conductor + cost.labour, 1e-12));
}

TEST_CASE("MagneticFilter COST scores a priced core as valid",
          "[magnetic-filter][cost]") {
    auto magnetic = make_reference_magnetic();
    auto inputs = make_reference_inputs();
    auto [valid, score] = MagneticFilterCost().evaluate_magnetic(&magnetic, &inputs);
    CHECK(valid);
    CHECK_THAT(score, WithinRel(MagneticFilterCost().calculate_cost(magnetic).total(), 1e-12));
}

TEST_CASE("MagneticFilter COST charges Sullivan's strand premium on litz",
          "[magnetic-filter][cost]") {
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {40, 20};
    cfg.numberParallels = {1, 1};
    cfg.coreShapeName = "E 35";
    cfg.coreMaterialName = "3C97";
    cfg.wireNames = {"Litz 4x0.1 - Grade 1 - Single Served", "Litz 4x0.1 - Grade 1 - Single Served"};
    auto magnetic = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    magnetic.get_mutable_coil().wind();
    auto cost = MagneticFilterCost().calculate_cost(magnetic);
    // 4 strands of 0.1 mm copper; C_m(0.1 mm) = 1 + 6e-26/(1e-4)^6 + 2.7e-9/(1e-4)^2 = 1 + 0.06 + 0.27.
    double area = 4 * std::numbers::pi * 0.05e-3 * 0.05e-3;
    CHECK_THAT(cost.conductor, WithinRel(sum_of_turn_lengths(magnetic) * area * 8940 * 14.737 * 1.33, 1e-9));
}

TEST_CASE("MagneticFilter COST flags a core with no price law instead of estimating or dropping it",
          "[magnetic-filter][cost]") {
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {40, 20};
    cfg.numberParallels = {1, 1};
    cfg.coreShapeName = "E 35";
    cfg.coreMaterialName = "46";  // Fair-Rite MgZn ferrite: no distributor prices fitted
    cfg.wireNames = {"Round 1.00 - Grade 1", "Round 1.00 - Grade 1"};
    auto magnetic = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    magnetic.get_mutable_coil().wind();
    auto cost = MagneticFilterCost().calculate_cost(magnetic);
    CHECK_FALSE(cost.priced());
    CHECK_THAT(cost.coreUnpricedReason, Catch::Matchers::ContainsSubstring("no price law for ferrite/MgZn"));
    // Copper and labour are still known; only the total is not.
    CHECK_THAT(cost.labour, WithinRel(12.6 / 60 * 7.27, 1e-12));
    CHECK_THROWS_WITH(cost.total(), Catch::Matchers::ContainsSubstring("unit cost unknown"));
    // The filter reports it as a failed scoring, which a non-strict adviser flow ranks worst for cost.
    auto inputs = make_reference_inputs();
    CHECK(MagneticFilterCost().evaluate_magnetic(&magnetic, &inputs) == std::pair<bool, double>{false, 0.0});

    // The same MnZn reference is unpriced when the basis carries no law for its class.
    auto reference = make_reference_magnetic();
    auto basis = MagneticFilterCost::default_basis();
    basis.corePrice.erase("ferrite/MnZn");
    auto unpriced = MagneticFilterCost(basis).calculate_cost(reference);
    CHECK_FALSE(unpriced.priced());
    CHECK_THAT(unpriced.coreUnpricedReason, Catch::Matchers::ContainsSubstring("no price law for ferrite/MnZn"));
}

TEST_CASE("MagneticAdviser keeps a core COST cannot price and ranks it last for cost",
          "[magnetic-filter][cost][magnetic-adviser]") {
    settings.reset();
    clear_databases();
    // Two ~100 uH parts on an E 42/21/15 (the same inputs as the ABT #1456 catalogue test), both gapped
    // 1 mm at 21 turns so the gap sets the inductance: N87 (priced, ferrite/MnZn) and Fair-Rite 46 (MgZn
    // ferrite, which has no price law).
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        100000, 100e-6, 25, WaveformLabel::TRIANGULAR, 2, 0.5, 1);
    auto part = [&](const std::string& material, double gap, int64_t turns, const std::string& reference) {
        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 42/21/15", OpenMagneticsTesting::get_ground_gap(gap),
                                                                 std::vector<int64_t>{turns}, 1, material);
        MAS::MagneticManufacturerInfo manufacturerInfo;
        manufacturerInfo.set_name("COST test");
        manufacturerInfo.set_reference(reference);
        magnetic.set_manufacturer_info(manufacturerInfo);
        return OpenMagnetics::magnetic_autocomplete(magnetic);
    };
    std::vector<OpenMagnetics::Magnetic> catalogue{part("46", 0.001, 21, "unpriced MgZn"), part("N87", 0.001, 21, "priced MnZn")};
    std::vector<MagneticFilterOperation> flow{MagneticFilterOperation(MagneticFilters::COST, true, false, 1.0)};
    MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, catalogue, flow, 5, false));
    REQUIRE(results.size() == 2);
    CHECK(results[0].first.get_magnetic().get_reference() == "priced MnZn");
    CHECK(results[1].first.get_magnetic().get_reference() == "unpriced MgZn");
    CHECK(results[1].second < results[0].second);
    settings.reset();
}

TEST_CASE("MagneticFilter COST prices a proprietary grade only with its own manufacturer's law",
          "[magnetic-filter][cost]") {
    auto make = [](const std::string& material) {
        OpenMagneticsTesting::QuickMagneticConfig cfg;
        cfg.numberTurns = {40, 20};
        cfg.numberParallels = {1, 1};
        cfg.coreShapeName = "E 35";
        cfg.coreMaterialName = material;
        cfg.wireNames = {"Round 1.00 - Grade 1", "Round 1.00 - Grade 1"};
        auto magnetic = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
        magnetic.get_mutable_coil().wind();
        return magnetic;
    };
    // Micrometals OC 60: powder/proprietary, priced by the law fitted on Micrometals' own grades.
    auto micrometals = make("OC 60");
    CHECK(MagneticFilterCost().calculate_cost(micrometals).priced());
    // Poco GPC 60 is also powder/proprietary, but no Poco part was priced: the Micrometals law must not be
    // borrowed for it.
    auto poco = make("GPC 60");
    auto cost = MagneticFilterCost().calculate_cost(poco);
    CHECK_FALSE(cost.priced());
    CHECK_THAT(cost.coreUnpricedReason, Catch::Matchers::ContainsSubstring("no price law for powder/proprietary/Poco"));
}

// =============================================================================
// Impedance display |Z| outside the core material's mu(f) span; adviser failure report
// =============================================================================
namespace {
OpenMagnetics::Magnetic load_wound_cmc_catalogue_part(const std::string& partNumber) {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc_catalogue_" + partNumber + ".json");
    std::ifstream file(path);
    REQUIRE(file.good());
    return magnetic_autocomplete(OpenMagnetics::Magnetic(json::parse(file)));
}

// The same part without its datasheet's measured impedance curve: IMPEDANCE then judges it on the
// core material's mu(f) range alone (the measured band judges points outside it on that curve).
OpenMagnetics::Magnetic without_measured_impedance(OpenMagnetics::Magnetic magnetic) {
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
    auto electrical = datasheetInfo.get_electrical().value();
    for (auto& entry : electrical) {
        entry.set_impedance_points(std::nullopt);
    }
    datasheetInfo.set_electrical(electrical);
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);
    REQUIRE_FALSE(MagneticFilterImpedance::get_measured_impedance_curve(magnetic));
    return magnetic;
}

OpenMagnetics::Inputs make_cmc_inputs(double frequency) {
    // Two windings (turns ratio 1), sinusoidal line current: a common-mode choke's operating point.
    return OpenMagnetics::Inputs::create_quick_operating_point_only_current(
        frequency, 0.001, 25, WaveformLabel::SINUSOIDAL, {1.0, 1.0}, 0.5, 0, {1.0});
}

void set_minimum_impedance(OpenMagnetics::Inputs& inputs, double frequency, double magnitude) {
    ImpedancePoint impedancePoint;
    impedancePoint.set_magnitude(magnitude);
    ImpedanceAtFrequency impedanceAtFrequency;
    impedanceAtFrequency.set_frequency(frequency);
    impedanceAtFrequency.set_impedance(impedancePoint);
    inputs.get_mutable_design_requirements().set_minimum_impedance(std::vector<ImpedanceAtFrequency>{impedanceAtFrequency});
}
}

// El Choker hands the adviser a CMC's noise point (150 kHz) and its 50 Hz mains point (line
// current, thermal). The IMPEDANCE filter's display |Z| at every operating point asked the core
// material for mu(50 Hz); every CMC ferrite is tabulated from above 50 Hz, so that threw and the
// adviser dropped the part: 4 of 299 chokes survived. The display |Z| is now written only where
// the material has mu(f) data and is absent elsewhere; the requirement and the score are unchanged.
TEST_CASE("MagneticFilter_Impedance_Display_Output_Only_Inside_Material_Mu_Span", "[magnetic-filter][impedance][mu-span]") {
    settings.reset();
    auto magnetic = load_wound_cmc_catalogue_part("7448052502");
    auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(magnetic.get_core().resolve_material());
    REQUIRE(minimumMaterialFrequency > 50);
    REQUIRE(minimumMaterialFrequency < 150000);
    REQUIRE(maximumMaterialFrequency > 150000);

    auto noiseOnly = make_cmc_inputs(150000);
    set_minimum_impedance(noiseOnly, 150000, 100);
    auto noiseAndMains = noiseOnly;
    noiseAndMains.get_mutable_operating_points().push_back(make_cmc_inputs(50).get_operating_points()[0]);
    REQUIRE(noiseAndMains.get_operating_points().size() == 2);
    REQUIRE(noiseAndMains.get_operating_points()[1].get_excitations_per_winding()[0].get_frequency() == 50);

    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, noiseOnly);

    std::vector<OpenMagnetics::Outputs> outputsNoiseOnly;
    auto [validNoiseOnly, scoringNoiseOnly] = filter->evaluate_magnetic(&magnetic, &noiseOnly, &outputsNoiseOnly);

    std::vector<OpenMagnetics::Outputs> outputs;
    auto [valid, scoring] = filter->evaluate_magnetic(&magnetic, &noiseAndMains, &outputs);

    CHECK(valid);
    CHECK(valid == validNoiseOnly);
    CHECK(scoring == scoringNoiseOnly);
    REQUIRE(outputs.size() == 2);
    REQUIRE(outputs[0].get_impedance());
    REQUIRE(outputs[0].get_impedance()->get_impedance_matrix());
    auto impedanceMatrix = outputs[0].get_impedance()->get_impedance_matrix().value();
    REQUIRE(impedanceMatrix.size() == 1);
    CHECK(impedanceMatrix[0].get_frequency() == 150000);
    auto expectedMatrix = outputsNoiseOnly[0].get_impedance()->get_impedance_matrix().value()[0];
    auto windingName = magnetic.get_coil().get_functional_description()[0].get_name();
    CHECK(impedanceMatrix[0].get_magnitude().at(windingName).at(windingName).get_nominal().value() ==
          expectedMatrix.get_magnitude().at(windingName).at(windingName).get_nominal().value());
    CHECK_FALSE(outputs[1].get_impedance());

    // A minimum-impedance REQUIREMENT with no point inside the material's data cannot be judged at
    // all, and is an error (points outside are skipped only when some point is judged: covered-band).
    auto requirementOutsideSpan = noiseAndMains;
    set_minimum_impedance(requirementOutsideSpan, 50, 1);
    std::vector<OpenMagnetics::Outputs> outputsRequirementOutsideSpan;
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &requirementOutsideSpan, &outputsRequirementOutsideSpan),
                      RequirementOutsideMaterialDataException);
}

// A candidate whose evaluation raises is still kept out of the ranking, but the adviser reports it
// with its error instead of returning fewer parts with no word why. 7448229004's A07 is tabulated
// from 9.9 kHz, so a 1 kHz impedance requirement cannot be evaluated on it; 7448052502's SC-1K107
// starts at 418 Hz and can.
TEST_CASE("MagneticAdviser_Reports_Candidates_Whose_Evaluation_Throws", "[magnetic-adviser][failed-candidates][mu-span]") {
    settings.reset();
    auto evaluable = load_wound_cmc_catalogue_part("7448052502");
    // Without its measured curve (which starts at 1 kHz and would judge the point: measured-band).
    auto notEvaluable = without_measured_impedance(load_wound_cmc_catalogue_part("7448229004"));
    const double requirementFrequency = 1000;
    REQUIRE(ComplexPermeability().get_frequency_range(evaluable.get_core().resolve_material()).first < requirementFrequency);
    REQUIRE(ComplexPermeability().get_frequency_range(notEvaluable.get_core().resolve_material()).first > requirementFrequency);

    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedance(inputs, requirementFrequency, 1e-3);

    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{evaluable.get_reference(), evaluable},
                                                             {notEvaluable.get_reference(), notEvaluable}};
    bool strictlyRequired = GENERATE(true, false);
    INFO("strictlyRequired " << strictlyRequired);
    std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, true, strictlyRequired, 1.0)};

    MagneticAdviser adviser(false);
    auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
    REQUIRE(results.size() == 1);
    CHECK(results[0].first.get_magnetic().get_reference() == evaluable.get_reference());

    const auto& failedCandidates = adviser.get_failed_candidates();
    REQUIRE(failedCandidates.size() == 1);
    CHECK(failedCandidates[0].first == notEvaluable.get_reference());
    CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("IMPEDANCE"));
    CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("A07"));
}

// =============================================================================
// Rank on the covered part: IMPEDANCE judges a minimumImpedance point only inside the core
// material's tabulated mu(f) range (El Choker CISPR-band request up to 30 MHz; 142 of 299 chokes
// were excluded because their ferrite's data stops at 1-16 MHz).
// =============================================================================
namespace {
void set_minimum_impedances(OpenMagnetics::Inputs& inputs, const std::vector<std::pair<double, double>>& frequencyAndMagnitude) {
    std::vector<ImpedanceAtFrequency> requirement;
    for (auto [frequency, magnitude] : frequencyAndMagnitude) {
        ImpedancePoint impedancePoint;
        impedancePoint.set_magnitude(magnitude);
        ImpedanceAtFrequency impedanceAtFrequency;
        impedanceAtFrequency.set_frequency(frequency);
        impedanceAtFrequency.set_impedance(impedancePoint);
        requirement.push_back(impedanceAtFrequency);
    }
    inputs.get_mutable_design_requirements().set_minimum_impedance(requirement);
}
}

// 7448229004's A07 is tabulated to ~12.5 MHz. A requirement at 150 kHz, 1 MHz and 30 MHz ranks it on
// the two points it covers, scores it as the mean over those two exactly as a requirement made of
// only those two would, and reports them as its judged frequencies. The 30 MHz point is not judged.
TEST_CASE("MagneticAdviser_Impedance_Ranks_On_The_Covered_Requirement_Points", "[magnetic-adviser][impedance][covered-band]") {
    settings.reset();
    // Without its measured curve: a part that carries one is judged above the mu(f) range on it (measured-band).
    auto part = without_measured_impedance(load_wound_cmc_catalogue_part("7448229004"));
    auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(part.get_core().resolve_material());
    REQUIRE(minimumMaterialFrequency < 150000);
    REQUIRE(maximumMaterialFrequency > 1e6);
    REQUIRE(maximumMaterialFrequency < 30e6);

    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedances(inputs, {{150000, 1e-3}, {1e6, 1e-3}, {30e6, 1e-3}});
    auto coveredOnly = make_cmc_inputs(150000);
    set_minimum_impedances(coveredOnly, {{150000, 1e-3}, {1e6, 1e-3}});

    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, inputs);
    auto [valid, scoring] = filter->evaluate_magnetic(&part, &inputs);
    auto [validCoveredOnly, scoringCoveredOnly] = filter->evaluate_magnetic(&part, &coveredOnly);
    CHECK(valid);
    CHECK(validCoveredOnly);
    CHECK(scoring > 0);
    CHECK(scoring == scoringCoveredOnly);
    REQUIRE(filter->get_judged_frequencies(&part, &inputs));
    CHECK(filter->get_judged_frequencies(&part, &inputs).value() == std::vector<double>{150000, 1e6});

    // A point inside the range that the part fails still invalidates it; one outside is not judged.
    auto failsInside = make_cmc_inputs(150000);
    set_minimum_impedances(failsInside, {{150000, 1e9}, {30e6, 1e-3}});
    CHECK_FALSE(filter->evaluate_magnetic(&part, &failsInside).first);
    auto hugeOutside = make_cmc_inputs(150000);
    set_minimum_impedances(hugeOutside, {{150000, 1e-3}, {30e6, 1e9}});
    CHECK(filter->evaluate_magnetic(&part, &hugeOutside).first);

    bool strictlyRequired = GENERATE(true, false);
    INFO("strictlyRequired " << strictlyRequired);
    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{part.get_reference(), part}};
    std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, true, strictlyRequired, 1.0)};
    MagneticAdviser adviser(false);
    auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
    REQUIRE(results.size() == 1);
    CHECK(results[0].first.get_magnetic().get_reference() == part.get_reference());
    CHECK(adviser.get_failed_candidates().empty());
    const auto& judgedFrequencies = adviser.get_judged_frequencies();
    REQUIRE(judgedFrequencies.contains(part.get_reference()));
    REQUIRE(judgedFrequencies.at(part.get_reference()).contains(MagneticFilters::IMPEDANCE));
    CHECK(judgedFrequencies.at(part.get_reference()).at(MagneticFilters::IMPEDANCE) == std::vector<double>{150000, 1e6});
    CHECK(adviser.get_measured_frequencies().at(part.get_reference()).at(MagneticFilters::IMPEDANCE).empty());
}

// A part whose material has no data at any requirement frequency cannot be judged on it: it is not
// ranked, and the adviser reports it, naming its range and the requirement's.
TEST_CASE("MagneticAdviser_Impedance_Reports_A_Part_With_No_Covered_Requirement_Point", "[magnetic-adviser][impedance][covered-band]") {
    settings.reset();
    // Without its measured curve (measured-band covers 20 and 30 MHz on a part that has one).
    auto part = without_measured_impedance(load_wound_cmc_catalogue_part("7448229004"));
    REQUIRE(ComplexPermeability().get_frequency_range(part.get_core().resolve_material()).second < 20e6);

    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedances(inputs, {{20e6, 1e-3}, {30e6, 1e-3}});

    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, inputs);
    CHECK(filter->get_judged_frequencies(&part, &inputs).value().empty());
    REQUIRE_THROWS_AS(filter->evaluate_magnetic(&part, &inputs), RequirementOutsideMaterialDataException);

    bool strictlyRequired = GENERATE(true, false);
    INFO("strictlyRequired " << strictlyRequired);
    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{part.get_reference(), part}};
    std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, true, strictlyRequired, 1.0)};
    MagneticAdviser adviser(false);
    auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
    CHECK(results.empty());
    const auto& failedCandidates = adviser.get_failed_candidates();
    REQUIRE(failedCandidates.size() == 1);
    CHECK(failedCandidates[0].first == part.get_reference());
    CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("IMPEDANCE"));
    CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("A07"));
    CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("no minimumImpedance frequency"));
}

// =============================================================================
// Gate only when losses count (ABT #1679): the loss-model frequency-span gate excludes a part only
// when the flow or the final simulation computes core losses, and reports every part it excludes.
// =============================================================================
namespace {
// A common-mode choke on P47, whose Steinmetz fit spans 100 kHz to 1 MHz: its losses at the 50 Hz
// mains point would be an extrapolation of the fit (WE 7448640406-0418 and 7448680200 are of this kind).
OpenMagnetics::Magnetic make_p47_cmc() {
    auto part = load_wound_cmc_catalogue_part("7448229004");
    part.get_mutable_core().get_mutable_functional_description().set_material("P47");
    auto manufacturerInfo = part.get_manufacturer_info().value();
    manufacturerInfo.set_reference("P47 choke");
    part.set_manufacturer_info(manufacturerInfo);
    return part;
}

OpenMagnetics::Inputs make_cmc_inputs_with_mains() {
    auto inputs = make_cmc_inputs(150000);
    inputs.get_mutable_operating_points().push_back(make_cmc_inputs(50).get_operating_points()[0]);
    return inputs;
}
}

TEST_CASE("MagneticAdviser_Loss_Span_Gate_Keeps_A_Part_When_No_Losses_Are_Computed", "[magnetic-adviser][loss-span-gate][covered-band]") {
    settings.reset();
    auto part = make_p47_cmc();
    auto inputs = make_cmc_inputs_with_mains();
    REQUIRE_FALSE(MagneticFilterLossModelFrequencySpan::is_material_evaluable(part.get_core().resolve_material(), inputs));

    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{part.get_reference(), part}};
    std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::VOLUME, true, true, false, 1.0),
                                                    MagneticFilterOperation(MagneticFilters::TURNS_RATIOS, true, false, true, 1.0)};
    MagneticAdviser adviser(false);
    auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
    REQUIRE(results.size() == 1);
    CHECK(results[0].first.get_magnetic().get_reference() == "P47 choke");
    CHECK(adviser.get_failed_candidates().empty());
}

TEST_CASE("MagneticAdviser_Loss_Span_Gate_Reports_A_Part_When_Losses_Are_Computed", "[magnetic-adviser][loss-span-gate][covered-band]") {
    settings.reset();
    auto part = make_p47_cmc();
    auto inputs = make_cmc_inputs_with_mains();
    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{part.get_reference(), part}};

    SECTION("a loss-based filter in the flow") {
        std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::VOLUME, true, true, false, 1.0),
                                                        MagneticFilterOperation(MagneticFilters::LOSSES_NO_PROXIMITY, true, true, false, 1.0)};
        MagneticAdviser adviser(false);
        auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
        CHECK(results.empty());
        const auto& failedCandidates = adviser.get_failed_candidates();
        REQUIRE(failedCandidates.size() == 1);
        CHECK(failedCandidates[0].first == "P47 choke");
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("LOSS_MODEL_FREQUENCY_SPAN"));
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("P47"));
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("50.000000 Hz"));
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("a filter of the flow computes core losses"));
    }
    SECTION("the final simulation, which computes core losses at every operating point") {
        std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::VOLUME, true, true, false, 1.0)};
        MagneticAdviser adviser(true);
        auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
        CHECK(results.empty());
        const auto& failedCandidates = adviser.get_failed_candidates();
        REQUIRE(failedCandidates.size() == 1);
        CHECK(failedCandidates[0].first == "P47 choke");
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("LOSS_MODEL_FREQUENCY_SPAN"));
        CHECK_THAT(failedCandidates[0].second, Catch::Matchers::ContainsSubstring("final simulation"));
    }
}

// =============================================================================
// A solid wire stored without numberConductors (optional in MAS basicWire) is one conductor.
// WE-CMB / WE-CMBNC chokes (e.g. 7448012002, 744842565) store their round wire that way; the Albach
// skin-effect model threw "Missing number of conductors" on it, so El Choker's final simulation
// excluded them, while Wire::calculate_conducting_area already read such a wire as one conductor.
// =============================================================================
TEST_CASE("Wire_Without_Number_Of_Conductors_Is_One_Conductor_Unless_Litz", "[wire-conductors][covered-band]") {
    settings.reset();
    auto counted = OpenMagnetics::find_wire_by_name("Round 0.475 - Grade 1");
    counted.set_number_conductors(1);
    auto uncounted = counted;
    uncounted.set_number_conductors(std::nullopt);
    REQUIRE_FALSE(uncounted.get_number_conductors());
    CHECK(uncounted.resolve_number_conductors() == 1);

    WindingSkinEffectLossesAlbachModel albach;
    for (double frequency : {1e5, 1e6, 1e7}) {
        INFO("frequency " << frequency);
        CHECK(albach.calculate_skin_factor(uncounted, frequency, 25) == albach.calculate_skin_factor(counted, frequency, 25));
    }

    auto litz = OpenMagnetics::find_wire_by_name("Litz 225x0.04 - Grade 1 - Double Served");
    litz.set_number_conductors(std::nullopt);
    CHECK_THROWS_AS(litz.resolve_number_conductors(), InvalidInputException);
    CHECK_THROWS_AS(albach.calculate_skin_factor(litz, 1e6, 25), InvalidInputException);

    // End to end: the winding losses of a catalogue choke whose round wires carry no count.
    auto part = load_wound_cmc_catalogue_part("7448229004");
    for (const auto& winding : part.get_coil().get_functional_description()) {
        REQUIRE_FALSE(OpenMagnetics::Coil::resolve_wire(winding).get_number_conductors());
    }
    auto inputs = make_cmc_inputs(150000);
    auto operatingPoint = inputs.get_operating_points()[0];
    WindingLossesOutput losses;
    REQUIRE_NOTHROW(losses = WindingLosses().calculate_losses(part, operatingPoint, 25));
    CHECK(losses.get_winding_losses() > 0);
}

// =============================================================================
// Measured band: a minimumImpedance point outside the core material's mu(f) range is judged on the
// part's own measured common-mode |Z| (datasheet impedancePoints, zero DC bias), interpolated
// log|Z| vs log f between the two bracketing measured points, never extrapolated. ACME's A07 mu(f)
// ends at ~12.5 MHz (the plot stops at mu = 10) while 7448229004 is measured from 1 kHz to 1 GHz.
// =============================================================================
namespace {
// 7448229004's measured points bracketing 22 MHz (tests/testData/cmc_catalogue_7448229004.json).
constexpr double kBracketLowFrequency = 21256900.0;
constexpr double kBracketLowMagnitude = 741.031;
constexpr double kBracketHighFrequency = 22387200.0;
constexpr double kBracketHighMagnitude = 728.922;
// Hand-computed: exp(ln 741.031 + (ln 22e6 - ln 21256900) / (ln 22387200 - ln 21256900) * (ln 728.922 - ln 741.031)).
constexpr double kMeasuredMagnitudeAt22MHz = 732.9775932440663;

OpenMagnetics::Magnetic with_impedance_points(OpenMagnetics::Magnetic magnetic, std::vector<DatasheetImpedancePoint> points) {
    auto manufacturerInfo = magnetic.get_manufacturer_info().value();
    auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
    auto electrical = datasheetInfo.get_electrical().value();
    REQUIRE(electrical.size() == 1);
    electrical[0].set_impedance_points(points);
    datasheetInfo.set_electrical(electrical);
    manufacturerInfo.set_datasheet_info(datasheetInfo);
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

std::vector<DatasheetImpedancePoint> impedance_points_of(const OpenMagnetics::Magnetic& magnetic) {
    return magnetic.get_manufacturer_info()->get_datasheet_info()->get_electrical()->at(0).get_impedance_points().value();
}

DatasheetImpedancePoint make_impedance_point(double frequency, double magnitude) {
    ImpedancePoint impedance;
    impedance.set_magnitude(magnitude);
    DatasheetImpedancePoint point;
    point.set_frequency(frequency);
    point.set_impedance(impedance);
    return point;
}
}

TEST_CASE("MagneticFilter_Impedance_Judges_Above_The_Mu_Span_On_The_Measured_Curve", "[magnetic-filter][impedance][measured-band]") {
    settings.reset();
    auto part = load_wound_cmc_catalogue_part("7448229004");
    auto [minimumMaterialFrequency, maximumMaterialFrequency] = ComplexPermeability().get_frequency_range(part.get_core().resolve_material());
    REQUIRE(maximumMaterialFrequency > 1e6);
    REQUIRE(maximumMaterialFrequency < 13e6);

    auto curve = MagneticFilterImpedance::get_measured_impedance_curve(part);
    REQUIRE(curve);
    CHECK(curve->front().first == 1000);
    CHECK(curve->back().first == 1e9);
    auto upper = std::find_if(curve->begin(), curve->end(), [](const auto& point) { return point.first >= 22e6; });
    REQUIRE(upper != curve->begin());
    REQUIRE(upper != curve->end());
    CHECK(std::prev(upper)->first == kBracketLowFrequency);
    CHECK(std::prev(upper)->second == kBracketLowMagnitude);
    CHECK(upper->first == kBracketHighFrequency);
    CHECK(upper->second == kBracketHighMagnitude);

    double measuredAt22MHz = MagneticFilterImpedance::interpolate_measured_impedance(curve.value(), 22e6);
    CHECK_THAT(measuredAt22MHz, WithinRel(kMeasuredMagnitudeAt22MHz, 1e-12));
    // A measured frequency returns its own point; outside the measured range nothing is extrapolated.
    CHECK(MagneticFilterImpedance::interpolate_measured_impedance(curve.value(), kBracketHighFrequency) == kBracketHighMagnitude);
    CHECK_THROWS_AS(MagneticFilterImpedance::interpolate_measured_impedance(curve.value(), 2e9), InvalidInputException);
    CHECK_THROWS_AS(MagneticFilterImpedance::interpolate_measured_impedance(curve.value(), 500), InvalidInputException);

    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, make_cmc_inputs(150000));

    // Same pass/fail and the same per-point score as a model point: valid just under the measured |Z|,
    // invalid just over it, score log10(Zmeasured / Zreq).
    auto under = make_cmc_inputs(150000);
    set_minimum_impedances(under, {{22e6, 100}});
    auto [validUnder, scoringUnder] = filter->evaluate_magnetic(&part, &under);
    CHECK(validUnder);
    CHECK_THAT(scoringUnder, WithinRel(std::log10(kMeasuredMagnitudeAt22MHz / 100), 1e-9));
    auto justUnder = make_cmc_inputs(150000);
    set_minimum_impedances(justUnder, {{22e6, kMeasuredMagnitudeAt22MHz * 0.999}});
    CHECK(filter->evaluate_magnetic(&part, &justUnder).first);
    auto justOver = make_cmc_inputs(150000);
    set_minimum_impedances(justOver, {{22e6, kMeasuredMagnitudeAt22MHz * 1.001}});
    CHECK_FALSE(filter->evaluate_magnetic(&part, &justOver).first);

    // A mixed requirement: model points inside the mu(f) range, measured points above it, and a point
    // above the measured range that is not judged. The score is the mean over the judged points.
    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedances(inputs, {{150000, 1e-3}, {1e6, 1e-3}, {22e6, 100}, {30e6, 1e-3}, {2e9, 1e-3}});
    auto modelOnly = make_cmc_inputs(150000);
    set_minimum_impedances(modelOnly, {{150000, 1e-3}, {1e6, 1e-3}});
    auto thirtyOnly = make_cmc_inputs(150000);
    set_minimum_impedances(thirtyOnly, {{30e6, 1e-3}});
    auto [valid, scoring] = filter->evaluate_magnetic(&part, &inputs);
    auto [validModelOnly, scoringModelOnly] = filter->evaluate_magnetic(&part, &modelOnly);
    auto [validThirtyOnly, scoringThirtyOnly] = filter->evaluate_magnetic(&part, &thirtyOnly);
    CHECK(valid);
    CHECK(validModelOnly);
    CHECK(validThirtyOnly);
    CHECK_THAT(scoring, WithinRel((2 * scoringModelOnly + scoringUnder + scoringThirtyOnly) / 4, 1e-9));
    CHECK(filter->get_judged_frequencies(&part, &inputs).value() == std::vector<double>{150000, 1e6, 22e6, 30e6});
    CHECK(filter->get_measured_frequencies(&part, &inputs).value() == std::vector<double>{22e6, 30e6});

    // The adviser reports both lists.
    bool strictlyRequired = GENERATE(true, false);
    INFO("strictlyRequired " << strictlyRequired);
    std::map<std::string, OpenMagnetics::Magnetic> catalogue{{part.get_reference(), part}};
    std::vector<MagneticFilterOperation> filterFlow{MagneticFilterOperation(MagneticFilters::IMPEDANCE, true, true, strictlyRequired, 1.0)};
    MagneticAdviser adviser(false);
    auto results = adviser.get_advised_magnetic(inputs, catalogue, filterFlow, 10);
    REQUIRE(results.size() == 1);
    CHECK(adviser.get_failed_candidates().empty());
    CHECK(adviser.get_judged_frequencies().at(part.get_reference()).at(MagneticFilters::IMPEDANCE) == std::vector<double>{150000, 1e6, 22e6, 30e6});
    CHECK(adviser.get_measured_frequencies().at(part.get_reference()).at(MagneticFilters::IMPEDANCE) == std::vector<double>{22e6, 30e6});
}

// A point measured under DC bias is a different curve: only the zero-bias points are read.
TEST_CASE("MagneticFilter_Impedance_Measured_Curve_Ignores_DC_Biased_Points", "[magnetic-filter][impedance][measured-band]") {
    settings.reset();
    auto part = load_wound_cmc_catalogue_part("7448229004");
    auto points = impedance_points_of(part);
    auto biased = make_impedance_point(22e6, 1e6);
    biased.set_current(2.0);
    points.push_back(biased);
    auto partWithBiased = with_impedance_points(part, points);
    auto curve = MagneticFilterImpedance::get_measured_impedance_curve(partWithBiased);
    REQUIRE(curve);
    CHECK(curve->size() == impedance_points_of(part).size());
    CHECK_THAT(MagneticFilterImpedance::interpolate_measured_impedance(curve.value(), 22e6), WithinRel(kMeasuredMagnitudeAt22MHz, 1e-12));

    // Only biased points: no zero-bias curve, so the part is judged on the mu(f) range alone.
    auto onlyBiased = with_impedance_points(part, {biased});
    CHECK_FALSE(MagneticFilterImpedance::get_measured_impedance_curve(onlyBiased));
}

// A part without a measured curve keeps the covered-band behaviour: points above the mu(f) range
// are not judged, and the score is the one of its covered points.
TEST_CASE("MagneticFilter_Impedance_Without_Measured_Curve_Keeps_The_Covered_Band", "[magnetic-filter][impedance][measured-band]") {
    settings.reset();
    auto part = without_measured_impedance(load_wound_cmc_catalogue_part("7448229004"));
    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedances(inputs, {{150000, 1e-3}, {1e6, 1e-3}, {22e6, 1e9}});
    auto coveredOnly = make_cmc_inputs(150000);
    set_minimum_impedances(coveredOnly, {{150000, 1e-3}, {1e6, 1e-3}});
    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, inputs);
    auto [valid, scoring] = filter->evaluate_magnetic(&part, &inputs);
    auto [validCoveredOnly, scoringCoveredOnly] = filter->evaluate_magnetic(&part, &coveredOnly);
    CHECK(valid);
    CHECK(scoring == scoringCoveredOnly);
    CHECK(filter->get_judged_frequencies(&part, &inputs).value() == std::vector<double>{150000, 1e6});
    CHECK(filter->get_measured_frequencies(&part, &inputs).value().empty());

    // With its curve the 22 MHz point is judged, and the 1 GOhm request fails it.
    auto measuredPart = load_wound_cmc_catalogue_part("7448229004");
    CHECK_FALSE(filter->evaluate_magnetic(&measuredPart, &inputs).first);
}

// Malformed measured data throws a specific error; nothing is skipped or guessed.
TEST_CASE("MagneticFilter_Impedance_Malformed_Measured_Curve_Throws", "[magnetic-filter][impedance][measured-band]") {
    settings.reset();
    auto part = load_wound_cmc_catalogue_part("7448229004");
    auto inputs = make_cmc_inputs(150000);
    set_minimum_impedances(inputs, {{150000, 1e-3}, {22e6, 1e-3}});
    auto filter = MagneticFilter::factory(MagneticFilters::IMPEDANCE, inputs);
    auto points = impedance_points_of(part);

    SECTION("a non-positive magnitude") {
        auto bad = points;
        bad[100].get_mutable_impedance().set_magnitude(0);
        auto badPart = with_impedance_points(part, bad);
        CHECK_THROWS_AS(filter->evaluate_magnetic(&badPart, &inputs), InvalidDatasheetImpedanceException);
        CHECK_THROWS_AS(MagneticFilterImpedance::get_measured_impedance_curve(badPart), InvalidDatasheetImpedanceException);
    }
    SECTION("a non-finite frequency") {
        auto bad = points;
        bad[100].set_frequency(std::numeric_limits<double>::quiet_NaN());
        auto badPart = with_impedance_points(part, bad);
        CHECK_THROWS_AS(filter->evaluate_magnetic(&badPart, &inputs), InvalidDatasheetImpedanceException);
    }
    SECTION("two values at one frequency") {
        auto bad = points;
        bad.push_back(make_impedance_point(kBracketLowFrequency, 2 * kBracketLowMagnitude));
        auto badPart = with_impedance_points(part, bad);
        CHECK_THROWS_AS(filter->evaluate_magnetic(&badPart, &inputs), InvalidDatasheetImpedanceException);
    }
    SECTION("two commonModeChoke entries with a zero-bias curve") {
        auto manufacturerInfo = part.get_manufacturer_info().value();
        auto datasheetInfo = manufacturerInfo.get_datasheet_info().value();
        auto electrical = datasheetInfo.get_electrical().value();
        electrical.push_back(electrical[0]);
        datasheetInfo.set_electrical(electrical);
        manufacturerInfo.set_datasheet_info(datasheetInfo);
        auto badPart = part;
        badPart.set_manufacturer_info(manufacturerInfo);
        CHECK_THROWS_AS(filter->evaluate_magnetic(&badPart, &inputs), InvalidDatasheetImpedanceException);
    }
    SECTION("a point without magnitude does not even parse") {
        json magneticJson;
        to_json(magneticJson, part);
        magneticJson["manufacturerInfo"]["datasheetInfo"]["electrical"][0]["impedancePoints"][100]["impedance"].erase("magnitude");
        CHECK_THROWS(OpenMagnetics::Magnetic(magneticJson));
    }
}
