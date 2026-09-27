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

#include "advisers/MagneticFilter.h"
#include <filesystem>
#include <fstream>
#include <source_location>
#include <chrono>
#include <cstdint>
#include <iostream>
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
#include "processors/MagneticSimulator.h"
#include "support/Settings.h"
#include "support/Utils.h"

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
// the 79 E 10/3 at 20 turns, which has the headroom to stay valid.
// Refreshed 2026-06-16 (ABT #13) after the saturation-derating rework (RAW
// B_sat at the 100 C hot corner, derating = hot temperature × margin). Top-1
// and top-2 are the SAME 79 E 10/3 / 20-turn designs as before — only their
// composite scores shifted with the recomputed saturation-headroom term
// (top-1 2.0 → 2.997). Slot 3 IMPROVED: the raw-B_sat gate admits a deeper
// valid pool for this 507 kHz spec, so the former INVALID-penalised fallback
// (95 E 8.3/4, 0.85) is replaced by a genuinely VALID third design
// (79 E 8/2 3 stacks, 0.75) — the thin-pool note below no longer applies.
// Refreshed 2026-06-27 (ABT #10): the core names now annotate the solved gap
// ("gapped 0.03 mm" on the E 10/3, "gapped 0.00 mm" on the E 8/2 — i.e. still
// physically ungapped, just rendered with the gap field). This is a pure
// naming/representation change from the gap-labelling rework: the shapes, turn
// counts, winding order, margins, wires AND scores are byte-identical to the
// previous snapshot, so the advised designs are unchanged — only the reference
// strings gained the gap suffix.
const std::vector<MagneticEntry> kTopThreeWinding = {
    // Refreshed 2026-07-07 (ABT #126 investigation): the health-pass MagneticEnergy
    // min/max interval fix (energy target at L_max) re-solves the E 10/3 to a
    // 0.04 mm gap needing 21 turns (conservative sizing for L at +tolerance) and
    // re-ranks slot 3 to the 95-material 9-turn stack. Regenerated on main
    // 17e1f850, user-approved.
    // Slot 0 re-pinned 2026-07-30 (ABT #377), user-approved. The value moved +0.061%
    // (2.30960636439206812 -> 2.31102462762095229) while EVERYTHING ELSE held: same winning
    // reference, same wire triple, and slots 1 and 2 unchanged (1.88518538772831645 exactly,
    // 0.70004398 vs 0.70004337 — inside kRelTol). A global normalisation shift would have moved
    // all three, so this is one candidate's score, not a re-ranking.
    // NOT caused by the drum/molded family work: reverting the TAK/SDE/B45 material additions
    // reproduces the identical new value, so does reverting the ABT #358 per-shape-family DC-bias
    // fix, and the windingOrder propagation cannot reach it (no catalogue record carries
    // windingOrder). ~50 physics/constructive commits landed between the previous pin and this
    // one — among them the ABT #339 permeability-curve fixes and the #358 DC-bias work — and the
    // exact one responsible was not isolated; bisect recipe if it ever matters: check out each
    // candidate, rebuild, run this test, watch for 2.30960636439206812.
    // Refreshed 2026-08-20 (ABT #834, user-approved). A REAL re-ranking, and the one
    // snapshot in that ticket whose head actually changed identity: the 95-material
    // E 8/2 2-stack at 13 turns (gapped 0.02 mm) now wins over the 79-material
    // E 10/3 at 21 turns, which drops to slot 2. NOT the winding-loss physics: with
    // the #832/#835/#837 loss fixes reverted in isolation this scenario produced the
    // IDENTICAL new head ("want 79 E 10/3 ... got 95 E 8/2" both ways), so the driver
    // is the same catalogue evolution as the sibling tables (135 MAS data/ commits
    // since 2026-06-16). The family was already competitive — the old slot 3 was a
    // 95-material E 8/2 3-stack at 9 turns — and the wire triple is unchanged across
    // all slots. Slots 0 and 1 are the same design differing only in winding order
    // (021 vs 012) with scores 0.0005% apart: a documented near-tie, not a meaningful
    // preference. Score scale moved (2.31 -> 1.70) with the normalization pool.
    {"95 E 8/2 2 stacks gapped 0.02 mm, Turns: 13, Order: 021, Non-Interleaved, Margin Taped 01",
     "Round 33.0 - Single Build || Round 41.0 - Single Build || Round 41.0 - Single Build",
     1.7000769237774718},
    {"95 E 8/2 2 stacks gapped 0.02 mm, Turns: 13, Order: 012, Non-Interleaved, Margin Taped 01",
     "Round 33.0 - Single Build || Round 41.0 - Single Build || Round 41.0 - Single Build",
     1.7},
    {"79 E 10/3 gapped 0.04 mm, Turns: 21, Order: 021, Non-Interleaved, Margin Taped 00",
     "Round 33.0 - Single Build || Round 41.0 - Single Build || Round 41.0 - Single Build",
     1},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOT
// =============================================================================

// ABT #1328: El Magnetic's /v1/design/custom with coreMode "available cores" on an ordinary buck
// inductor (11 V -> 5 V, 10 W, 100 kHz, ripple 0.4; inputs from its Kirchhoff buck wizard) grew
// past 16 GB in PyOpenMagnetics, whose build left cores_stock.ndjson out and so advised over all
// 18943 catalogue cores. Under the default stock setting the adviser works on the 1573-core stock
// catalogue: a few hundred MB and a few seconds here.
TEST_CASE("MagneticAdviser available cores advises a buck inductor from the stock catalogue",
          "[adviser][magnetic-adviser][available-cores][abt-1328]") {
    settings.reset();
    clear_databases();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "abt1328_buck_inductor_inputs.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    OpenMagnetics::Inputs inputs(nlohmann::json::parse(file));

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    auto results = adviser.get_advised_magnetic(inputs, 3);
    CHECK(results.size() == 3);

    auto stockPath = std::filesystem::path(__FILE__).parent_path().parent_path() / "MAS" / "data" / "cores_stock.ndjson";
    std::ifstream stockFile(stockPath);
    REQUIRE(stockFile.is_open());
    size_t stockRecords = 0;
    std::string line;
    while (std::getline(stockFile, line)) {
        if (!line.empty()) {
            stockRecords++;
        }
    }
    CHECK(coreDatabase.size() == stockRecords);
    settings.reset();
}

// Field report (ABT #1410/#1411): the secondary resonant inductor of a 30 kW CLLC, one
// winding, 33.5 uH +/-15 %, 41 x 43 x 44 mm envelope, two measured operating points of the
// resonant tank current (118-143 kHz, no DC).
static OpenMagnetics::Inputs load_cllc_resonant_inductor_inputs() {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "abt1410_cllc_resonant_inductor_inputs.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    return OpenMagnetics::Inputs(nlohmann::json::parse(file));
}

// ABT #1410: standard-cores design mode returned an ER 51/10/38 (51 mm wide) inside a
// 41 x 43 x 44 mm envelope, because only the height was ever checked.
TEST_CASE("MagneticAdviser standard cores keeps the wound magnetic inside the maximum dimensions",
          "[adviser][magnetic-adviser][standard-cores][abt-1410]") {
    settings.reset();
    clear_databases();
    auto inputs = load_cllc_resonant_inductor_inputs();
    REQUIRE(inputs.get_design_requirements().get_maximum_dimensions());

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto results = adviser.get_advised_magnetic(inputs, 10);

    REQUIRE(!results.empty());
    for (auto& [mas, scoring] : results) {
        auto& magnetic = mas.get_mutable_magnetic();
        auto dimensions = magnetic.get_maximum_dimensions();
        INFO(magnetic.get_reference() << ": " << dimensions[0] * 1000 << " x " << dimensions[1] * 1000 << " x " << dimensions[2] * 1000 << " mm");
        CHECK(MagneticFilterMaximumDimensions::magnetic_fits(magnetic, inputs));
    }
    settings.reset();
}

// ABT #1411: available-cores mode returned powder toroids at 4.6-10.4 uH against the
// 28.475-38.525 uH requirement. Two defects: the single-winding resonant inductor was
// classified as a transformer by its converter topology (turns seeded from volt-seconds),
// and the available-cores power path had no inductance gate to catch the result.
TEST_CASE("MagneticAdviser available cores returns only designs inside the magnetizing inductance band",
          "[adviser][magnetic-adviser][available-cores][abt-1411]") {
    settings.reset();
    clear_databases();
    auto inputs = load_cllc_resonant_inductor_inputs();
    const auto& requirement = inputs.get_design_requirements().get_magnetizing_inductance();
    double minimum = requirement.get_minimum().value();
    double maximum = requirement.get_maximum().value();

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    auto results = adviser.get_advised_magnetic(inputs, 10);

    REQUIRE(!results.empty());
    for (auto& [mas, scoring] : results) {
        REQUIRE(mas.get_outputs().size() == inputs.get_operating_points().size());
        for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
            auto inductance = mas.get_outputs()[operatingPointIndex].get_inductance();
            REQUIRE(inductance);
            double magnetizingInductance = resolve_dimensional_values(inductance->get_magnetizing_inductance().get_magnetizing_inductance());
            INFO(mas.get_mutable_magnetic().get_reference() << " at operating point " << operatingPointIndex << ": " << magnetizingInductance * 1e6 << " uH");
            CHECK(magnetizingInductance >= minimum);
            CHECK(magnetizingInductance <= maximum);
        }
    }
    settings.reset();
}

// ABT #1426: the resonant inductor of the same CLLC (ABT #1410/#1411) carries pure AC, so it
// swings the full flux every cycle and CORE LOSS, not saturation, decides its turns. The
// adviser took the fewest turns that hold L and clear the saturation margin and stopped
// there: 8-9 turns, 0.34 T, 59-62 W of core loss against 8-13 W of copper, 590-620 C.
//
// Criterion. At a fixed inductance B ~ 1/N, so core loss falls as N^-beta (Steinmetz
// beta ~ 2.5 for these power ferrites) while window-filled copper rises as N^2; their sum
// is least where Pcore = (2 / beta) Pcu, below Pcu. A design whose core dissipates more than
// three times its copper is therefore far on the few-turns side of that minimum (one more
// turn would lower the total), unless it stopped at a buildability limit; the 3x bound
// leaves that room plus the proximity loss the full simulation adds and the adviser's
// estimate leaves out. Each returned design is simulated in full here (core and winding
// losses, every operating point) and held to it.
TEST_CASE("MagneticAdviser standard cores sizes an AC-dominated inductor by losses, not at the saturation floor",
          "[adviser][magnetic-adviser][standard-cores][abt-1426]") {
    settings.reset();
    clear_databases();
    auto inputs = load_cllc_resonant_inductor_inputs();

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto results = adviser.get_advised_magnetic(inputs, 10);

    REQUIRE(!results.empty());
    MagneticSimulator simulator;
    for (auto& [mas, scoring] : results) {
        auto magnetic = mas.get_magnetic();
        auto simulated = simulator.simulate(inputs, magnetic);
        double coreLosses = 0;
        double windingLosses = 0;
        for (auto& output : simulated.get_outputs()) {
            REQUIRE(output.get_core_losses());
            REQUIRE(output.get_winding_losses());
            coreLosses += output.get_core_losses()->get_core_losses();
            windingLosses += output.get_winding_losses()->get_winding_losses();
        }
        INFO(magnetic.get_reference() << ": N = " << magnetic.get_coil().get_functional_description()[0].get_number_turns()
             << ", core " << coreLosses << " W, winding " << windingLosses << " W");
        CHECK(coreLosses <= 3 * windingLosses);
    }
    settings.reset();
}

// Hidden ABT #1426 driver (no suite tag, so no suite run pulls it in): times the adviser on
// the CLLC resonant inductor (standard cores with and without the envelope, available cores)
// and prints every returned design fully simulated, per operating point.
static void report_abt1426_designs(const std::string& label, OpenMagnetics::Inputs inputs,
                                   CoreAdviser::CoreAdviserModes mode) {
    settings.reset();
    clear_databases();
    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(mode);
    auto start = std::chrono::steady_clock::now();
    auto results = adviser.get_advised_magnetic(inputs, 10);
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cerr << "ABT1426 " << label << ": " << results.size() << " results in " << seconds << " s" << std::endl;
    MagneticSimulator simulator;
    for (auto& [mas, scoring] : results) {
        auto magnetic = mas.get_magnetic();
        auto core = magnetic.get_core();
        std::cerr << "ABT1426   " << core.get_name().value_or("?") << " | N="
                  << magnetic.get_coil().get_functional_description()[0].get_number_turns() << " | gaps mm [";
        for (auto& gap : core.get_gapping()) {
            if (gap.get_type() != GapType::RESIDUAL) {
                std::cerr << gap.get_length() * 1000 << " ";
            }
        }
        std::cerr << "]" << std::endl;
        auto simulated = simulator.simulate(inputs, magnetic);
        for (size_t operatingPointIndex = 0; operatingPointIndex < simulated.get_outputs().size(); ++operatingPointIndex) {
            auto& output = simulated.get_outputs()[operatingPointIndex];
            double bPeak = output.get_core_losses()->get_magnetic_flux_density()->get_processed()->get_peak().value();
            std::cerr << "ABT1426     OP" << operatingPointIndex << ": Bpk=" << bPeak
                      << " T core=" << output.get_core_losses()->get_core_losses()
                      << " W winding=" << output.get_winding_losses()->get_winding_losses()
                      << " W Tmax=" << output.get_temperature()->get_maximum_temperature() << " C";
            auto perWinding = output.get_winding_losses()->get_winding_losses_per_winding().value()[0];
            double skin = 0, proximity = 0;
            auto skinLosses = perWinding.get_skin_effect_losses().value();
            auto proximityLosses = perWinding.get_proximity_effect_losses().value();
            for (auto v : skinLosses.get_losses_per_harmonic()) skin += v;
            for (auto v : proximityLosses.get_losses_per_harmonic()) proximity += v;
            std::cerr << " (ohmic " << perWinding.get_ohmic_losses()->get_losses() << " skin " << skin << " prox " << proximity << ") wire "
                      << magnetic.get_mutable_coil().get_wires()[0].get_name().value_or("?")
                      << " x" << magnetic.get_coil().get_functional_description()[0].get_number_parallels() << std::endl;
        }
    }
    settings.reset();
}

TEST_CASE("ABT 1426 driver: CLLC resonant inductor designs and adviser timing", "[.][abt-1426-driver]") {
    auto inputs = load_cllc_resonant_inductor_inputs();
    report_abt1426_designs("standard cores, envelope", inputs, CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    auto withoutEnvelope = inputs;
    auto requirements = withoutEnvelope.get_design_requirements();
    requirements.set_maximum_dimensions(std::nullopt);
    withoutEnvelope.set_design_requirements(requirements);
    report_abt1426_designs("standard cores, no envelope", withoutEnvelope, CoreAdviser::CoreAdviserModes::STANDARD_CORES);
    report_abt1426_designs("available cores, envelope", inputs, CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
}

TEST_CASE("MagneticAdviser 3-winding end-to-end top-3 snapshot",
          "[adviser][magnetic-adviser][characterisation][heavy][end-to-end]") {
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
          "[!benchmark][benchmark-magnetic-adviser]") {
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

// ABT #1369: the web CMC wizard's design (230 V, 10 A line current, 2 windings, 500 ohm at
// 150 kHz — tests/testData/cmc/cmc_default_web_inputs.json) came back with ZERO magnetics from
// calculate_advised_magnetics while the CoreAdviser alone returned toroids. The final saturation
// gate in process_wound_candidate classified the CMC as an inductor (every winding on the line
// side) and compared I_sat, derived from the mH-range common-mode inductance, against the full
// 10 A line current — which cancels in a CMC's core. It must gate on the common-mode magnetizing
// current, as MagnetizingInductance and the saturation filter already do.
TEST_CASE("MagneticAdviser advises a toroidal CMC for the web wizard's default design",
          "[adviser][magnetic-adviser][available-cores][cmc][suppression][abt-1369]") {
    settings.reset();
    clear_databases();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "cmc/cmc_default_web_inputs.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    OpenMagnetics::Inputs inputs(nlohmann::json::parse(file));
    REQUIRE(OpenMagnetics::Inputs::can_be_common_mode_choke(inputs.get_operating_points()[0]));

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(CoreAdviser::CoreAdviserModes::AVAILABLE_CORES);
    auto results = adviser.get_advised_magnetic(inputs, 3);

    REQUIRE(results.size() > 0);
    for (auto& [mas, scoring] : results) {
        CHECK(mas.get_mutable_magnetic().get_mutable_core().get_type() == CoreType::TOROIDAL);
        CHECK(mas.get_mutable_magnetic().get_mutable_coil().get_functional_description().size() == 2);
    }
    settings.reset();
    clear_databases();
}
