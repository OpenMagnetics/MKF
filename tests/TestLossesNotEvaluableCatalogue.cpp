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
