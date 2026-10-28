// =============================================================================
// TestLossesNotEvaluableCatalogue.cpp
// =============================================================================
// A catalogue part whose core material carries no core-loss data any model can
// run (volumetricLosses: {}) must not abort the catalogue MagneticAdviser, and
// must not vanish from it either. Like a datasheet-only part:
//   - the loss-based filters do not apply to it (applies_to) and record no
//     score -- never a stand-in loss;
//   - it is ranked on the filters that do apply;
//   - the adviser flags it, with the reason, in get_losses_not_evaluable().
// =============================================================================

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "advisers/MagneticAdviser.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Magnetic.h"
#include "physical_models/CoreLosses.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include "TestingUtils.h"

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::ContainsSubstring;

namespace {

const std::string losslessMaterialName = "Test Ferrite Without Loss Data";

OpenMagnetics::Magnetic catalogue_part(const std::string& reference, std::optional<CoreMaterial> material = std::nullopt) {
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {10};
    cfg.numberParallels = {1};
    cfg.coreShapeName = "E 35";
    cfg.coreMaterialName = "3C97";
    cfg.wireNames = {"Round 1.00 - Grade 1"};
    cfg.numberStacks = 1;
    auto magnetic = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    if (material) {
        magnetic.get_mutable_core().get_mutable_functional_description().set_material(material.value());
    }
    magnetic.get_mutable_coil().wind();
    MagneticManufacturerInfo manufacturerInfo;
    manufacturerInfo.set_name("Test Maker");
    manufacturerInfo.set_reference(reference);
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

// A catalogue part on another core with another number of turns, wound the same way.
OpenMagnetics::Magnetic catalogue_part_with(const std::string& reference, const std::string& shapeName, int64_t numberTurns) {
    OpenMagneticsTesting::QuickMagneticConfig cfg;
    cfg.numberTurns = {numberTurns};
    cfg.numberParallels = {1};
    cfg.coreShapeName = shapeName;
    cfg.coreMaterialName = "3C97";
    cfg.wireNames = {"Round 1.00 - Grade 1"};
    cfg.numberStacks = 1;
    auto magnetic = OpenMagneticsTesting::create_quick_test_magnetic(cfg);
    magnetic.get_mutable_coil().wind();
    MagneticManufacturerInfo manufacturerInfo;
    manufacturerInfo.set_name("Test Maker");
    manufacturerInfo.set_reference(reference);
    magnetic.set_manufacturer_info(manufacturerInfo);
    return magnetic;
}

// 3C97 with its core-loss data removed, as MAS stores a material nobody has measured losses for.
CoreMaterial lossless_material() {
    auto material = find_core_material_by_name("3C97");
    material.set_name(losslessMaterialName);
    material.set_volumetric_losses({});
    material.set_mass_losses(std::nullopt);
    return material;
}

OpenMagnetics::Inputs inductor_inputs() {
    return OpenMagnetics::Inputs::create_quick_operating_point_only_current(100e3, 100e-6, 25, WaveformLabel::TRIANGULAR, 0.4, 0.5, 1.0);
}

}  // namespace

TEST_CASE("A material with no core-loss model makes the loss filters not applicable, and says why",
          "[losses-not-evaluable][adviser][smoke-test]") {
    settings.reset();
    auto normal = catalogue_part("NORMAL");
    auto lossless = catalogue_part("LOSSLESS", lossless_material());
    REQUIRE(CoreLossesModel::get_methods(lossless.get_core().resolve_material()).empty());

    REQUIRE_FALSE(MagneticFilter::core_losses_not_evaluable_reason(&normal));
    auto reason = MagneticFilter::core_losses_not_evaluable_reason(&lossless);
    REQUIRE(reason);
    REQUIRE_THAT(reason.value(), ContainsSubstring(losslessMaterialName));

    for (auto filter : {MagneticFilters::LOSSES, MagneticFilters::LOSSES_NO_PROXIMITY, MagneticFilters::CORE_AND_DC_LOSSES,
                        MagneticFilters::CORE_DC_AND_SKIN_LOSSES, MagneticFilters::TEMPERATURE_RISE, MagneticFilters::LOSSES_TIMES_VOLUME}) {
        auto magneticFilter = MagneticFilter::factory(filter, inductor_inputs());
        CHECK(magneticFilter->applies_to(&normal));
        CHECK_FALSE(magneticFilter->applies_to(&lossless));
    }
    // Filters that do not use the core losses still judge it.
    CHECK(MagneticFilter::factory(MagneticFilters::VOLUME, inductor_inputs())->applies_to(&lossless));
    // The loss-span gate has no span to judge on it, so it does not apply either (it would drop the part).
    MagneticFilterLossModelFrequencySpan spanFilter;
    CHECK(spanFilter.applies_to(&normal));
    CHECK_FALSE(spanFilter.applies_to(&lossless));
    settings.reset();
}

TEST_CASE("The catalogue adviser ranks a part whose core losses are not evaluable, and flags it",
          "[losses-not-evaluable][adviser]") {
    settings.reset();
    std::vector<OpenMagnetics::Magnetic> catalogue{catalogue_part("NORMAL-A"), catalogue_part("LOSSLESS", lossless_material()),
                                                   catalogue_part("NORMAL-B")};
    std::vector<MagneticFilterOperation> flow{
        MagneticFilterOperation(MagneticFilters::LOSSES_NO_PROXIMITY, true, false, false, 1.0),
        MagneticFilterOperation(MagneticFilters::VOLUME, true, false, false, 1.0),
    };
    MagneticAdviser adviser;
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inductor_inputs(), catalogue, flow, 10, false));

    std::set<std::string> references;
    for (auto& [mas, score] : results) {
        references.insert(mas.get_magnetic().get_reference());
        CHECK(std::isfinite(score));
    }
    REQUIRE(references.count("NORMAL-A") == 1);
    REQUIRE(references.count("NORMAL-B") == 1);
    // Neither dropped nor fatal: the part is in the results ...
    REQUIRE(references.count("LOSSLESS") == 1);

    // ... flagged, with the reason ...
    auto flagged = adviser.get_losses_not_evaluable();
    REQUIRE(flagged.size() == 1);
    REQUIRE(flagged.count("LOSSLESS") == 1);
    REQUIRE_THAT(flagged.at("LOSSLESS"), ContainsSubstring(losslessMaterialName));
    REQUIRE_THAT(flagged.at("LOSSLESS"), ContainsSubstring("no core-loss data"));

    // ... and ranked only on what applies to it: no loss score, never a stand-in.
    auto scorings = adviser.get_scorings();
    REQUIRE(scorings["LOSSLESS"].count(MagneticFilters::LOSSES_NO_PROXIMITY) == 0);
    REQUIRE(scorings["LOSSLESS"].count(MagneticFilters::VOLUME) == 1);
    REQUIRE(scorings["NORMAL-A"].count(MagneticFilters::LOSSES_NO_PROXIMITY) == 1);

    // A later run over parts that all have loss models leaves no stale flag behind.
    std::vector<OpenMagnetics::Magnetic> normalOnly{catalogue_part("NORMAL-A"), catalogue_part("NORMAL-B")};
    REQUIRE_NOTHROW(adviser.get_advised_magnetic(inductor_inputs(), normalOnly, flow, 10, false));
    REQUIRE(adviser.get_losses_not_evaluable().empty());
    settings.reset();
}

TEST_CASE("The catalogue adviser evaluates losses in one pass, its filters and its hits alike",
          "[adviser][catalogue][loss-temperature-loop][cost]") {
    // ABT #2124 (Alf, 2026-10-28): screening a catalogue ran the loss <-> temperature loop (ABT #838, with its
    // lowest-stable-point search, ABT #1840) for every part, 20x slower than before it. The catalogue search runs one
    // pass now, its filters and the final simulation of the hits alike (this test first asserted that the hits kept
    // the loop; Alf then decided one pass for them too).
    settings.reset();
    std::vector<OpenMagnetics::Magnetic> catalogue{catalogue_part("PART-A"), catalogue_part("PART-B")};
    std::vector<MagneticFilterOperation> flow{
        MagneticFilterOperation(MagneticFilters::LOSSES_NO_PROXIMITY, true, false, false, 1.0),
    };
    MagneticAdviser adviser(true);
    const size_t onePassBefore = MagneticSimulator::get_loss_solves(false);
    const size_t loopBefore = MagneticSimulator::get_loss_solves(true);
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inductor_inputs(), catalogue, flow, 10, true));
    const size_t onePass = MagneticSimulator::get_loss_solves(false) - onePassBefore;
    const size_t loop = MagneticSimulator::get_loss_solves(true) - loopBefore;
    INFO(onePass << " one-pass loss solves, " << loop << " looped loss solves, " << results.size() << " hits");
    REQUIRE(results.size() == 2);
    for (auto& [mas, score] : results) {
        REQUIRE(mas.get_outputs().size() == 1);
        CHECK(mas.get_outputs()[0].get_core_losses());
    }
    // The filter stage: one pass per part; the final simulation: one pass per hit.
    CHECK(onePass == catalogue.size() + results.size());
    // Never the loop.
    CHECK(loop == 0);
    // Outside the catalogue search the advisers' forced loop is back.
    CHECK_FALSE(CatalogueSearchOnePassScope::active());
    CHECK(settings.get_simulator_loss_temperature_loop());
    settings.reset();
}

TEST_CASE("The catalogue adviser refills its hits when a candidate fails its final simulation",
          "[adviser][catalogue][refill]") {
    // ABT #2124 (Alf, 2026-10-28): the adviser cut its ranking to the requested number of hits BEFORE the final
    // simulation, so every hit that failed there shortened the list (a real 400-part slice returned 2 of 5). The
    // next-ranked candidate takes a failed hit's place now, and the failure is still reported with its reason.
    settings.reset();
    auto overflowing = catalogue_part_with("OVERFLOWING", "E 25/13/7", 120);
    std::vector<OpenMagnetics::Magnetic> catalogue{overflowing, catalogue_part("FITS-A"), catalogue_part("FITS-B"),
                                                   catalogue_part("FITS-C")};
    // Smallest volume first: the overflowing part on the smaller core ranks first.
    std::vector<MagneticFilterOperation> flow{
        MagneticFilterOperation(MagneticFilters::VOLUME, true, false, false, 1.0),
    };
    MagneticAdviser adviser(true);
    std::vector<std::pair<OpenMagnetics::Mas, double>> results;
    REQUIRE_NOTHROW(results = adviser.get_advised_magnetic(inductor_inputs(), catalogue, flow, 2, true));
    std::string failures;
    for (const auto& [reference, reason] : adviser.get_failed_candidates()) {
        failures += "\n  " + reference + ": " + reason;
    }
    INFO("failed candidates:" << failures);
    std::vector<std::string> references;
    for (auto& [mas, score] : results) {
        references.push_back(mas.get_magnetic().get_reference());
    }
    INFO("hits: " << references.size());
    // The premise: the overflowing part ranks first and fails its final simulation.
    const auto& failed = adviser.get_failed_candidates();
    auto overflowingFailure = std::find_if(failed.begin(), failed.end(), [](const auto& entry) { return entry.first == "OVERFLOWING"; });
    REQUIRE(overflowingFailure != failed.end());
    REQUIRE_THAT(overflowingFailure->second, ContainsSubstring("final simulation"));
    // Two hits were asked for and two fitting parts remain: two hits, the overflowing one not among them.
    REQUIRE(results.size() == 2);
    CHECK(std::find(references.begin(), references.end(), "OVERFLOWING") == references.end());
    settings.reset();
}
