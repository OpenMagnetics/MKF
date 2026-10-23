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
#include <catch2/matchers/catch_matchers_string.hpp>

#include "advisers/MagneticAdviser.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "physical_models/InitialPermeability.h"
#include "support/Exceptions.h"
#include "processors/MagneticSimulator.h"
#include "processors/Sweeper.h"
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
    // Refreshed 2026-09-30 (user-approved): the COST filter now scores the unit cost in US$ at 1000
    // pieces (core price law + copper + winding labour) instead of numberLayers + wire relative cost, so
    // the composite re-ranks. The 79-material E 10/3 at 21 turns returns to slot 0; the 95-material
    // E 8/2 2-stack leaves the top 3 and a 79 E 8.3/4 2-stack at 12 turns takes slots 1-2 (the same
    // design in two winding orders). Both 79 and 95 are priced MnZn ferrites, so no candidate was
    // dropped by a COST throw; the old winner was out-scored. The per-design US$ were not broken down.
    {"79 E 10/3 gapped 0.04 mm, Turns: 21, Order: 021, Non-Interleaved, Margin Taped 01",
     "Round 35.5 - Single Build || Round 41.0 - Single Build || Round 41.0 - Single Build",
     2},
    {"79 E 8.3/4 2 stacks gapped 0.01 mm, Turns: 12, Order: 021, Non-Interleaved, Margin Taped 00",
     "Round 33.0 - Single Build || Round 38.0 - Single Build || Round 38.0 - Single Build",
     0.85007868429843947},
    {"79 E 8.3/4 2 stacks gapped 0.01 mm, Turns: 12, Order: 012, Non-Interleaved, Margin Taped 00",
     "Round 33.0 - Single Build || Round 38.0 - Single Build || Round 38.0 - Single Build",
     0.84999999999999998},
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

// The same inductor at a fraction of its measured tank current (voltage scaled with it, the
// inductor's own V = L di/dt). At full current no standard core inside the 41 x 43 x 44 mm
// envelope stays under the 130 C temperature gate once the copper is scored as the coil stage
// would wind it (coolest: 95 E 41/16.5/12.5 x3, 136.7 C), so the tests below that check what
// the adviser does with a design it can build (envelope, inductance band, loss balance) run it
// at half current (2026-10-04, ABT #1412).
static OpenMagnetics::Inputs load_cllc_resonant_inductor_inputs_at_current_fraction(double currentFraction) {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "abt1410_cllc_resonant_inductor_inputs.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    auto inputsJson = nlohmann::json::parse(file);
    for (auto& operatingPoint : inputsJson["operatingPoints"]) {
        for (auto& excitation : operatingPoint["excitationsPerWinding"]) {
            for (const char* signal : {"current", "voltage"}) {
                REQUIRE(excitation[signal]["waveform"].contains("data"));
                REQUIRE_FALSE(excitation[signal].contains("processed"));
                for (auto& value : excitation[signal]["waveform"]["data"]) {
                    value = value.get<double>() * currentFraction;
                }
            }
        }
    }
    return OpenMagnetics::Inputs(inputsJson);
}

// ABT #1410: standard-cores design mode returned an ER 51/10/38 (51 mm wide) inside a
// 41 x 43 x 44 mm envelope, because only the height was ever checked.
TEST_CASE("MagneticAdviser standard cores keeps the wound magnetic inside the maximum dimensions",
          "[adviser][magnetic-adviser][standard-cores][abt-1410]") {
    settings.reset();
    clear_databases();
    auto inputs = load_cllc_resonant_inductor_inputs_at_current_fraction(0.5);
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
    auto inputs = load_cllc_resonant_inductor_inputs_at_current_fraction(0.5);
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
    auto inputs = load_cllc_resonant_inductor_inputs_at_current_fraction(0.5);

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

// Hidden ABT #1449 driver: one adviser case of the CLLC resonant inductor, timed, with a
// fingerprint of every returned design (core, gaps, turns, parallels, wires, score and a hash
// of the full Mas JSON) so two engine builds can be diffed for identical results.
// ABT1449_CASE = envelope | noenvelope (default) | available; ABT1449_N = number of results (10).
TEST_CASE("ABT 1449 driver: CLLC resonant inductor, one adviser case, fingerprinted", "[.][abt-1449-driver]") {
    auto inputs = load_cllc_resonant_inductor_inputs();
    std::string adviserCase = std::getenv("ABT1449_CASE") ? std::getenv("ABT1449_CASE") : "noenvelope";
    size_t numberResults = std::getenv("ABT1449_N") ? std::stoul(std::getenv("ABT1449_N")) : 10;
    auto mode = CoreAdviser::CoreAdviserModes::STANDARD_CORES;
    if (adviserCase == "noenvelope") {
        auto requirements = inputs.get_design_requirements();
        requirements.set_maximum_dimensions(std::nullopt);
        inputs.set_design_requirements(requirements);
    }
    else if (adviserCase == "available") {
        mode = CoreAdviser::CoreAdviserModes::AVAILABLE_CORES;
    }
    else {
        REQUIRE(adviserCase == "envelope");
    }
    settings.reset();
    clear_databases();
    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(mode);
    auto start = std::chrono::steady_clock::now();
    auto results = adviser.get_advised_magnetic(inputs, numberResults);
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cerr << "ABT1449 " << adviserCase << " n=" << numberResults << ": " << results.size() << " results in " << seconds << " s" << std::endl;
    for (auto& [mas, scoring] : results) {
        auto magnetic = mas.get_magnetic();
        std::ostringstream line;
        line << std::setprecision(17) << "ABT1449   " << magnetic.get_core().get_name().value_or("?") << " | gaps [";
        for (auto& gap : magnetic.get_core().get_gapping()) {
            line << gap.get_length() << " ";
        }
        line << "] |";
        auto wires = magnetic.get_mutable_coil().get_wires();
        const auto& functionalDescription = magnetic.get_coil().get_functional_description();
        for (size_t windingIndex = 0; windingIndex < functionalDescription.size(); ++windingIndex) {
            line << " N=" << functionalDescription[windingIndex].get_number_turns()
                 << " x" << functionalDescription[windingIndex].get_number_parallels()
                 << " " << wires[windingIndex].get_name().value_or("?");
        }
        json masJson;
        to_json(masJson, mas);
        line << " | score " << scoring << " | mas hash " << std::hash<std::string>{}(masJson.dump());
        std::cerr << line.str() << std::endl;
    }
    settings.reset();
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

namespace {

// Replays WebLibMKF calculate_advised_magnetics_impl for an interferenceSuppression
// input: weights normalised to sum 1, then the explicit suppression filter flow the
// binding builds (CORE_MINIMUM_IMPEDANCE, COST, LOSSES, DIMENSIONS, LEAKAGE_INDUCTANCE,
// TURN_COUNT) instead of the weights map.
std::vector<std::pair<OpenMagnetics::Mas, double>> replay_web_suppression_magnetic_adviser(const std::string& fixtureName) {
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/" + fixtureName);
    std::ifstream file(path);
    REQUIRE(file.is_open());
    auto fixture = nlohmann::json::parse(file);
    OpenMagneticsTesting::apply_web_engine_settings(fixture.at("engineSettingsAtCall"));
    Settings::GetInstance().set_coil_delimit_and_compact(true);

    OpenMagnetics::Inputs inputs(fixture.at("inputs"));
    REQUIRE(inputs.get_design_requirements().get_application().value() == "interferenceSuppression");
    CoreAdviser::CoreAdviserModes coreMode;
    from_json(fixture.at("coreMode"), coreMode);

    std::map<std::string, double> weightsByName = fixture.at("weights");
    double weightSum = 0;
    for (auto const& [name, weight] : weightsByName) {
        weightSum += weight;
    }
    std::map<MagneticFilters, double> weights;
    for (auto const& [name, weight] : weightsByName) {
        MagneticFilters filter;
        from_json(name, filter);
        weights[filter] = weight / weightSum;
    }
    double wCost = weights.at(MagneticFilters::COST);
    double wLosses = weights.at(MagneticFilters::LOSSES);
    double wDims = weights.at(MagneticFilters::DIMENSIONS);
    std::vector<MagneticFilterOperation> filterFlow{
        MagneticFilterOperation(MagneticFilters::CORE_MINIMUM_IMPEDANCE, true, true, true, std::max(1.0, wDims * 2.0)),
        MagneticFilterOperation(MagneticFilters::COST, true, true, false, wCost),
        MagneticFilterOperation(MagneticFilters::LOSSES, true, true, false, wLosses),
        MagneticFilterOperation(MagneticFilters::DIMENSIONS, true, true, false, wDims),
        MagneticFilterOperation(MagneticFilters::LEAKAGE_INDUCTANCE, true, true, false, wDims),
        MagneticFilterOperation(MagneticFilters::TURN_COUNT, true, false, false, std::max(wLosses, wDims)),
    };

    OpenMagnetics::MagneticAdviser adviser;
    adviser.set_core_mode(coreMode);
    return adviser.get_advised_magnetic(inputs, filterFlow, fixture.at("maximumNumberResults").get<size_t>());
}

}  // namespace

// ABT #1699: the web's CMC wizard -> Magnetic Adviser call, captured at the engine proxy.
// It threw "[COIL_NOT_PROCESSED] ... Missing turns description to evaluate cost filter":
// a candidate reached the COST filter without a wound coil.
// The fixture's inputs are Kirchhoff design_cmc's output for the wizard spec (230 V, 10 A, 50 Hz, 2
// windings, 500 ohm at 150 kHz, C_par 10 pF, dV/dt 50 V/ns). Its "Nominal" operating point was captured
// when Kirchhoff drove each winding with I_cm = C*dV/dt = 0.5 A as a continuous 150 kHz sinusoid, so each
// winding saw v = L*w*(sum I_cm) = 0.5305 mH * 2pi*150 kHz * 1 A = 500 V peak (~0.7 T; no toroid stayed
// under the adviser's 130 C gate). C*dV/dt only lasts the few-ns switching edge. The CM current at 150 kHz
// is the switch node's 150 kHz component through C_par: V_bus = sqrt(2)*230 = 325.3 V, t_r = V_bus/(dV/dt)
// = 6.5 ns, V_1 = (2/pi)*V_bus*sinc(pi*f*t_r) = 207.1 V, I_cm = 2pi*f*C_par*V_1 = 1.95 mA in total, i.e.
// 0.98 mA per winding on its +-10 A line (DM) level, and the voltage on every winding is the CM drop on
// the choke's own CM impedance, |Z|*I_cm = 500 ohm * 1.95 mA = 0.98 V peak. That "Nominal" operating point
// was regenerated with the corrected Kirchhoff (cmc_common_mode_current_peak); the line-frequency operating
// point and the design requirements are unchanged.
TEST_CASE("MagneticAdviser web CMC suppression flow returns wound toroidal CMCs",
          "[adviser][magnetic-adviser][available-cores][cmc][suppression][abt-1699]") {
    settings.reset();
    clear_databases();
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = replay_web_suppression_magnetic_adviser("cmc_calculate_advised_magnetics.json"));
    REQUIRE(results.size() > 0);
    for (auto& [mas, scoring] : results) {
        CHECK(mas.get_mutable_magnetic().get_mutable_core().get_type() == CoreType::TOROIDAL);
        CHECK(mas.get_magnetic().get_coil().get_turns_description().has_value());
    }
    settings.reset();
    clear_databases();
}

// ABT #1699: the web's DMC wizard -> Magnetic Adviser call (standard cores, 6 results).
TEST_CASE("MagneticAdviser web DMC suppression flow returns wound chokes",
          "[adviser][magnetic-adviser][standard-cores][dmc][suppression][abt-1699]") {
    settings.reset();
    clear_databases();
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = replay_web_suppression_magnetic_adviser("dmc_calculate_advised_magnetics.json"));
    REQUIRE(results.size() > 0);
    for (auto& [mas, scoring] : results) {
        CHECK(mas.get_magnetic().get_coil().get_turns_description().has_value());
    }
    settings.reset();
    clear_databases();
}

// ABT #1794: the candidates of the DMC switching-ripple call below include Nanoperm 2000 toroids.
// To reach 450 ohm at 150 kHz on T 25/16/10 the filter raises the turns until the line current biases
// the core to 1554 A/m, past the knee of that material, where its permeability model gives mu_r < 1
// (0.963) and the impedance model refuses the bias with CALCULATION_CORE_SATURATES (ABT #1740). The
// filter caught only ModelNotAvailableException, so one saturating candidate threw out of the whole
// adviser. It is rejected instead and the adviser goes on with the others (Alf, 2026-10-12).
TEST_CASE("CoreMinimumImpedance filter rejects a candidate whose core saturates at its bias",
          "[adviser][magnetic-filter][dmc][suppression][abt-1794]") {
    settings.reset();
    clear_databases();
    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json");
    std::ifstream file(path);
    REQUIRE(file.is_open());
    auto fixture = nlohmann::json::parse(file);
    OpenMagneticsTesting::apply_web_engine_settings(fixture.at("engineSettingsAtCall"));
    OpenMagnetics::Inputs inputs(fixture.at("inputs"));
    REQUIRE(inputs.get_design_requirements().get_minimum_impedance().has_value());

    // The premise, outside the filter: at that bias the material is saturated, and at a small one it
    // is a material. Until ABT #1855 the permeability refused that bias (SaturationException); since then
    // (Alf approved this rewrite 2026-10-14, ABT #1855) the saturated material takes the vacuum's
    // permeability, 1, and the filter rejects the candidate through the shared saturation check
    // (ResultsAnalyzer::check_saturation on B at the bias).
    auto coreMaterial = find_core_material_by_name("Nanoperm 2000");
    CHECK(InitialPermeability::get_initial_permeability(coreMaterial, 25, 1554.28, 150000) == 1);
    CHECK(InitialPermeability::get_initial_permeability(coreMaterial, 25, 10, 150000) > 1);

    // The adviser's standard-cores candidates: each toroid shape with the material, ungapped.
    for (std::string shapeName : {"T 25/16/10", "T 30/20/10", "T 40/32/15"}) {
        INFO("core " << shapeName << " - Nanoperm 2000 - Ungapped");
        json coreJson = json::parse(R"({"functionalDescription": {"gapping": [], "material": "Nanoperm 2000", "numberStacks": 1, "shape": "", "type": "toroidal"}, "name": "candidate"})");
        coreJson["functionalDescription"]["shape"] = shapeName;
        OpenMagnetics::Core core(coreJson);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(OpenMagnetics::Coil(json::parse(R"({"bobbin": "Dummy", "functionalDescription": [{"isolationSide": "primary", "name": "Primary", "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}, {"isolationSide": "primary", "name": "Secondary", "numberParallels": 1, "numberTurns": 1, "wire": "Dummy"}], "layersDescription": null, "sectionsDescription": null, "turnsDescription": null})")));

        MagneticFilterCoreMinimumImpedance filter;
        std::pair<bool, double> verdict{true, -1};
        REQUIRE_NOTHROW(verdict = filter.evaluate_magnetic(&magnetic, &inputs));
        CHECK_FALSE(verdict.first);
        CHECK(verdict.second == 0);
    }
    settings.reset();
    clear_databases();
}

TEST_CASE("MagneticAdviser web DMC whose switching ripple is resolved returns wound chokes",
          "[adviser][magnetic-adviser][standard-cores][dmc][suppression][bug]") {
    // The web's DMC wizard -> Magnetic Adviser call of 2026-10-11: two windings, each carrying the
    // 50 Hz line current (4.8 A peak) with a 0.44 A, 100 kHz switching ripple, sampled at 16384
    // points so the ripple is resolved (the 1024-point fixture above aliases it to 2.4 kHz).
    // The magnetizing-inductance model took the CMC branch for it (matching winding currents on a
    // toroid) and drove the core with a triangle of the line current's whole swing at 100 kHz:
    // ~400 W of core losses, every candidate over the 130 C limit, no result. The DMC topology
    // now selects the DMC magnetizing current, the sum of the winding currents as sampled.
    settings.reset();
    clear_databases();
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = replay_web_suppression_magnetic_adviser("dmc_calculate_advised_magnetics_switching_ripple.json"));
    REQUIRE(results.size() > 0);
    for (auto& [mas, scoring] : results) {
        CHECK(mas.get_magnetic().get_coil().get_turns_description().has_value());
    }

    auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "web_adviser_flows/dmc_calculate_advised_magnetics_switching_ripple.json");
    std::ifstream file(path);
    OpenMagnetics::Inputs inputs(nlohmann::json::parse(file).at("inputs"));
    auto operatingPoint = inputs.get_operating_point(0);
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    REQUIRE(excitations.size() == 2);
    auto lineCurrent = excitations[0].get_current()->get_waveform()->get_data();
    auto neutralCurrent = excitations[1].get_current()->get_waveform()->get_data();
    REQUIRE(lineCurrent.size() == 16384);

    MagnetizingInductance magnetizingInductance;
    magnetizingInductance.set_topology(inputs.get_design_requirements().get_topology());
    magnetizingInductance.calculate_inductance_and_magnetic_flux_density(results[0].first.get_magnetic(), &operatingPoint);
    auto magnetizingCurrent = operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform()->get_data();
    REQUIRE(magnetizingCurrent.size() == lineCurrent.size());
    for (size_t index = 0; index < lineCurrent.size(); index += 1024) {
        CHECK_THAT(magnetizingCurrent[index], Catch::Matchers::WithinAbs(lineCurrent[index] + neutralCurrent[index], 1e-9));
    }

    // The core-losses-over-frequency plot drives the core the same way (ABT #1754): at the
    // operating frequency it must report what the simulator reports for this DMC. Without the
    // topology the model cannot tell the DMC from a CMC and takes the CMC triangle, which is what
    // the sweep's own copy of the branches always did: ~550 W against ~2 W.
    auto magnetic = results[0].first.get_magnetic();
    auto topology = inputs.get_design_requirements().get_topology();
    double operatingFrequency = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    MagneticSimulator simulator;
    simulator.set_core_losses_model_name(CoreLossesModels::STEINMETZ);
    simulator.set_topology(topology);
    auto simulatedOperatingPoint = inputs.get_operating_point(0);
    auto simulatedCoreLosses = simulator.calculate_core_losses(simulatedOperatingPoint, magnetic);
    // The simulator solves the losses with the temperature they heat the core to (ABT #838: 76.4 C
    // here, from 25 C ambient), so the sweep is run at that solved temperature and both sides price
    // the same core (at ambient the sweep reads 6.836 W against the simulator's 6.932 W). Approved
    // by Alf 2026-10-13.
    double solvedTemperature = simulatedCoreLosses.get_temperature();
    CHECK(solvedTemperature > inputs.get_operating_point(0).get_conditions().get_ambient_temperature());

    auto dmcSweep = Sweeper::sweep_core_losses_over_frequency(magnetic, inputs.get_operating_point(0), topology,
                                                              operatingFrequency, 100 * operatingFrequency, 3, solvedTemperature);
    auto cmcBranchSweep = Sweeper::sweep_core_losses_over_frequency(magnetic, inputs.get_operating_point(0), std::nullopt,
                                                                    operatingFrequency, 100 * operatingFrequency, 3, solvedTemperature);
    REQUIRE_THAT(dmcSweep.get_x_points()[0], Catch::Matchers::WithinRel(operatingFrequency, 1e-9));
    INFO("solved core temperature " << solvedTemperature << " C");
    INFO("f = " << operatingFrequency << " Hz: simulator " << simulatedCoreLosses.get_core_losses() << " W, DMC sweep "
         << dmcSweep.get_y_points()[0] << " W, CMC branch " << cmcBranchSweep.get_y_points()[0] << " W");
    CHECK_THAT(dmcSweep.get_y_points()[0], Catch::Matchers::WithinRel(simulatedCoreLosses.get_core_losses(), 0.01));
    CHECK_FALSE(cmcBranchSweep.get_y_points()[0] < 10 * simulatedCoreLosses.get_core_losses());

    // Each switching cycle swings at its own instantaneous line bias (ABT #1764): with the permeability of the
    // line swing's peak (~104) for every cycle the ripple near the zero crossings was priced ~4x low, 1.74 W.
    // Reference: the record split at its current's valleys, each cycle's ripple B = mu0 * mu(H of the cycle's
    // mean current) * N * di / l_e on the material's DC-bias curve, Steinmetz per cycle at its own rate: 6.84 W
    // (first-principles script on this magnetic, the ungapped T 48/23/19 MPP 300 with 30 + 30 turns).
    REQUIRE(magnetic.get_core().get_shape_name() == "T 48/23/19");
    REQUIRE(magnetic.get_core().get_material_name() == "MPP 300");
    REQUIRE(magnetic.get_core().get_functional_description().get_gapping().empty());
    REQUIRE(magnetic.get_coil().get_functional_description()[0].get_number_turns() == 30);
    CHECK_THAT(simulatedCoreLosses.get_core_losses(), Catch::Matchers::WithinRel(6.84, 0.05));

    // The operating point as stored (waveforms, no harmonics) is the same excitation as its processed
    // copy, so the sweep must price it the same (ABT #1764): it read 16.3 W against 1.74 W, because
    // the switching frequency came from the cached harmonics and fell back to 50 Hz without them,
    // which dropped the line swing's DC-bias roll-off of the MPP 300 permeability.
    std::ifstream storedFile(path);
    OperatingPoint storedOperatingPoint(nlohmann::json::parse(storedFile)["inputs"]["operatingPoints"][0]);
    REQUIRE_FALSE(storedOperatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics());
    auto storedSweep = Sweeper::sweep_core_losses_over_frequency(magnetic, storedOperatingPoint, topology,
                                                                 operatingFrequency, 100 * operatingFrequency, 3, solvedTemperature);
    INFO("stored operating point: " << storedSweep.get_y_points()[0] << " W, processed: " << dmcSweep.get_y_points()[0] << " W");
    CHECK_THAT(storedSweep.get_y_points()[0], Catch::Matchers::WithinRel(dmcSweep.get_y_points()[0], 0.01));

    // Swept to 500 Hz and 5 kHz the record is the same line current with its ripple, faster: the line
    // swing still biases the core (ABT #1764). A fixed f < 400 Hz gate on the switching-frequency
    // detection gave those points no bias: mu 300 and 1.02 T on an MPP 300 core that saturates at
    // 0.8 T, and 685 W / 28754 W in the plot.
    auto fluxDensityPeakAt = [&](double lineFrequency) {
        auto swept = inputs.get_operating_point(0);
        OpenMagnetics::Inputs::scale_time_to_frequency(swept, lineFrequency, false, true);
        for (auto& excitation : swept.get_mutable_excitations_per_winding()) {
            excitation.set_magnetizing_current(std::nullopt);
            excitation.set_magnetic_flux_density(std::nullopt);
        }
        MagnetizingInductance sweptModel;
        sweptModel.set_topology(topology);
        return sweptModel.calculate_inductance_and_magnetic_flux_density(magnetic, &swept).second.get_processed()->get_peak().value();
    };
    double lineFrequencyFluxDensityPeak = fluxDensityPeakAt(operatingFrequency);
    auto saturation = magnetic.get_core().resolve_material().get_saturation();
    REQUIRE_FALSE(saturation.empty());
    double saturationFluxDensity = saturation[0].get_magnetic_flux_density();
    CHECK(lineFrequencyFluxDensityPeak < saturationFluxDensity);
    for (double lineFrequency : {10 * operatingFrequency, 100 * operatingFrequency}) {
        double fluxDensityPeak = fluxDensityPeakAt(lineFrequency);
        INFO("line frequency " << lineFrequency << " Hz: B peak " << fluxDensityPeak << " T, at " << operatingFrequency << " Hz "
             << lineFrequencyFluxDensityPeak << " T, saturation " << saturationFluxDensity << " T");
        CHECK(fluxDensityPeak < saturationFluxDensity);
        CHECK_THAT(fluxDensityPeak, Catch::Matchers::WithinRel(lineFrequencyFluxDensityPeak, 0.02));
    }
    settings.reset();
    clear_databases();
}

static void strip_nulls_abt1942_adviser(nlohmann::json& value) {
    if (value.is_array()) {
        for (auto& element : value) {
            strip_nulls_abt1942_adviser(element);
        }
        return;
    }
    if (!value.is_object()) {
        return;
    }
    for (auto it = value.begin(); it != value.end();) {
        if (it->is_null()) {
            it = value.erase(it);
        }
        else {
            strip_nulls_abt1942_adviser(*it);
            ++it;
        }
    }
}

// ABT #1942 (Alf, 2026-10-22): a wound candidate one of whose operating points has no thermal operating point (here a
// 95/97 ferrite core whose losses heat it towards its Curie point, where the magnetisation curve can no longer carry the
// voltage's flux) gets no temperature from simulate (ABT #1876). The MagneticAdviser's temperature filter threw
// "[MISSING_DATA] simulate() returned no temperature" and took the whole call with it; such a candidate is rejected,
// recorded in get_failed_candidates with the solve's reason, and the search goes on.
TEST_CASE("MagneticAdviser rejects a wound candidate with no thermal operating point instead of failing the call",
          "[adviser][magnetic-adviser][abt-1942]") {
    settings.reset();
    clear_databases();
    std::ifstream file(OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "catalogue_missing_turns.json"));
    REQUIRE(file.good());
    auto masJson = nlohmann::json::parse(file);
    strip_nulls_abt1942_adviser(masJson);
    OpenMagnetics::Inputs inputs(masJson.at("inputs"));

    OpenMagnetics::MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inputs, 3));
    CHECK(!results.empty());

    const auto& failedCandidates = adviser.get_failed_candidates();
    bool noThermalOperatingPointRecorded = false;
    for (const auto& [reference, reason] : failedCandidates) {
        INFO(reference << ": " << reason);
        if (reason.find("no thermal operating point") != std::string::npos) {
            noThermalOperatingPointRecorded = true;
            CHECK_THAT(reason, Catch::Matchers::ContainsSubstring("No thermal operating point"));
        }
    }
    CHECK(noThermalOperatingPointRecorded);
    for (const auto& [mas, scoring] : results) {
        for (const auto& output : mas.get_outputs()) {
            CHECK(output.get_temperature());
        }
    }
    settings.reset();
    clear_databases();
}
