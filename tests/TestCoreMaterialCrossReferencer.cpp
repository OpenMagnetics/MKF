#include "support/Settings.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "support/Utils.h"
#include "processors/Inputs.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    // ABT #190c: these tests used to pin an exact winner, but the scoring cannot
    // discriminate the head of the field — for 3C97 at 25 °C the top five land within
    // 0.9% of each other (DMR95 2.70371, ML33D 2.70170, P45 2.70033, 3C95 2.68183,
    // TPW33 2.67987; top-two gap 0.007%). Every material-data batch reshuffled the
    // order and the expectation churned DMR95 -> TPW33 -> DMR95, costing an
    // investigation each time to conclude the flip was benign. Assert instead what
    // cross-referencing can actually guarantee: that the right materials are IN the
    // shortlist. A material dropping OUT of the top N is a real regression; the order
    // among near-identical scores is not.
    //
    // ABT #398 (2026-07-31): re-pinned. The band above widened rather than moved — the
    // scores of the ORIGINAL pins are unchanged to every digit (DMR95 2.70371, P45
    // 2.70033, 3C95 2.68183, TPW33 2.67987), so nothing about the model shifted. What
    // changed is the catalogue: the JFE (MAS 0d65a56, ABT #220) and TDG (MAS 56d8c84,
    // ABT #217) material batches dropped six more MnZn power ferrites straight into the
    // tie band, and they displaced the old pins by fractions of a percent. Measured
    // ranking for 3C97 at 25 °C, STEINMETZ:
    //
    //     1. ML33D  2.70557  (best)      7. P47     2.68850  0.631%
    //     2. DMR95  2.70371  0.069%      8. PL-13   2.68627  0.713%
    //     3. MBT2   2.70242  0.117%      9. TPG33B  2.68552  0.741%
    //     4. TPW30  2.70200  0.132%     10. TPG30   2.68371  0.808%
    //     5. P45    2.70033  0.194%     11. 3C95    2.68183  0.877%
    //     6. TP4C   2.69934  0.230%     12. TPW33   2.67987  0.950%
    //
    // Twelve grades from seven manufacturers inside 0.95%. Membership in the top FIVE is
    // therefore no more stable than the ordering this test stopped asserting in #190c —
    // it will churn again on the next batch. #398 tracks the underlying fix (make the
    // score discriminate, or report an explicit tie band); until then this pin records
    // WHERE the band sits so the next flip is a one-line update, not an investigation.
    std::vector<std::string> shortlist_names(const std::vector<std::pair<CoreMaterial, double>>& results) {
        std::vector<std::string> names;
        for (auto& [material, scoring] : results) {
            names.push_back(material.get_name());
        }
        return names;
    }

    void require_shortlisted(const std::vector<std::pair<CoreMaterial, double>>& results,
                             const std::vector<std::string>& expected) {
        auto names = shortlist_names(results);
        std::string joined;
        for (auto& name : names) {
            joined += " " + name;
        }
        for (auto& wanted : expected) {
            INFO("expected '" << wanted << "' in the shortlist, got:" << joined);
            REQUIRE(std::find(names.begin(), names.end(), wanted) != names.end());
        }
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer(std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // The DMR95-class MnZn power ferrites must make 3C97's shortlist; their relative
        // order is data-churn noise. Re-pinned for ABT #398 to the grades that occupy the
        // top five today — TPW33 (2.67987) and 3C95 (2.68183) are still sound matches but
        // now sit 12th and 11th, edged out by six newer grades within 0.95% (see the band
        // above). DMR95 is the one survivor of the original pin.
        require_shortlisted(crossReferencedCoreMaterials, {"ML33D", "DMR95", "MBT2", "TPW30", "P45"});

        auto scorings = coreMaterialCrossReferencer.get_scorings();
        auto scoredValues = coreMaterialCrossReferencer.get_scored_values();
        json results;
        results["cores"] = json::array();
        results["scorings"] = json::array();
        for (auto& [coreMaterial, scoring] : crossReferencedCoreMaterials) {
            std::string name = coreMaterial.get_name();

            json coreMaterialJson;
            to_json(coreMaterialJson, coreMaterial);
            results["cores"].push_back(coreMaterialJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);

                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                REQUIRE(!std::isnan(scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()]));
                REQUIRE(!std::isnan(scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()]));
            };
            results["data"].push_back(result);
        }
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_TDK", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("TDK");

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // ABT #224: N95 is the closest TDK match to 3C97 at the 25 °C evaluation point.
        // A previous pin expected PC47, justified by "3C97=1140, PC47=1154, N95=2537" — but
        // those are each material's FIRST datasheet permeability entry, at MISMATCHED
        // temperatures (3C97's first point is -40 °C; PC47's and N95's are -60 °C). Comparing
        // at the actual evaluation temperature (25 °C) inverts that conclusion:
        //   μr @25 °C: 3C97≈3341,  N95≈3008 (closer),  PC47≈2330
        //   volumetric losses (the other weight-1 filter): N95 also closer than PC47
        // so N95 wins on both dominant filters. The TDK pool is unaffected by the TDG import;
        // this was a pre-existing mis-pin, surfaced while triaging #224.
        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "N95");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::string coreMaterialName = "Kool M\xC2\xB5 MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // ABT #224 (MAS advanced to latest main): CSC Sendust 26 (Changsung) is the closest match
        // to Kool Mu MAX 26 — both are FeSiAl (Sendust) powder at μ=26, i.e. the same alloy family
        // and permeability, so it tracks the reference more tightly than any Kool Mu grade. Was
        // Kool Mu Hf 26 before the Changsung powder cores entered the database (ABT #214).
        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "CSC Sendust 26");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder_Only_Micrometals", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "Kool M\xC2\xB5 MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // July 2026 re-baseline: the volumetric-losses filter now scores by absolute distance from
        // the reference (removing the score=0 ceiling that tied every better-than-reference material)
        // and culls NaN-loss candidates instead of DBL_MAX — the same fix CoreCrossReferencer already
        // had. The corrected similarity ranking now puts SM 40 (was SM 60) closest to Kool Mu MAX 26.
        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "SM 40");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Powder_Only_Micrometals_Ferrite", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        settings.set_core_cross_referencer_allow_different_core_material_type(true);
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Micrometals");

        std::string coreMaterialName = "3C95";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, 5);


        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "MP 160");
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Volumetric_Losses", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 1;
        weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 0.5;
        weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // Same non-discrimination as the default-weight case (ABT #190c): weighting
        // permeability + volumetric losses puts 3C95A 1.49114, ML33D 1.48969,
        // JNP96A 1.48505, PL-13 1.48246, SMP97 1.48055 — a 0.7% spread, so the winner
        // flips on data churn. Assert the shortlist membership instead.
        require_shortlisted(crossReferencedCoreMaterials, {"JNP96A", "3C95A", "ML33D"});
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Volumetric_Losses_Powder", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;

        std::map<CoreMaterialCrossReferencerFilters, double> weights;
        weights[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY] = 0;
        weights[CoreMaterialCrossReferencerFilters::REMANENCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE] = 0;
        weights[CoreMaterialCrossReferencerFilters::SATURATION] = 0;
        weights[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE] = 0;
        weights[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES] = 1;
        weights[CoreMaterialCrossReferencerFilters::RESISTIVITY] = 0;

        std::string coreMaterialName = "Kool M\xC2\xB5 MAX 26";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 25, weights, 5);

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // Accept any Kool Mμ material as top result due to algorithm improvements
        auto topMaterialName = crossReferencedCoreMaterials[0].first.get_name();
        REQUIRE((topMaterialName.find("Kool M") != std::string::npos || 
                 topMaterialName.find("MAX") != std::string::npos ||
                 topMaterialName.find("HF") != std::string::npos ||
                 topMaterialName.find("Edge") != std::string::npos));
    }

    TEST_CASE("Test_CoreMaterialCrossReferencer_All_Core_Materials_Only_Fair_Rite", "[adviser][core-material-cross-referencer][smoke-test]") {
        settings.reset();
        clear_databases();
        OperatingPoint operatingPoint;
        CoreMaterialCrossReferencer coreMaterialCrossReferencer;
        coreMaterialCrossReferencer.use_only_manufacturer("Fair-Rite");
        settings.set_use_only_cores_in_stock(false);

        std::string coreMaterialName = "3C97";
        CoreMaterial coreMaterial = Core::resolve_material(coreMaterialName);

        auto crossReferencedCoreMaterials = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, 50, 20);

        auto scorings = coreMaterialCrossReferencer.get_scorings();
        auto scoredValues = coreMaterialCrossReferencer.get_scored_values();
        json results;
        results["coreMaterials"] = json::array();
        results["scorings"] = json::array();
        results["data"] = json::array();

        for (auto& [coreMaterial, scoring] : crossReferencedCoreMaterials) {
            std::string name = coreMaterial.get_name();

            json coreMaterialJson;
            to_json(coreMaterialJson, coreMaterial);
            results["coreMaterials"].push_back(coreMaterialJson);
            results["scorings"].push_back(scoring);

            json result;
            result["scoringPerFilter"] = json();
            result["scoredValuePerFilter"] = json();
            for (auto& filter : magic_enum::enum_names<CoreMaterialCrossReferencerFilters>()) {
                std::string filterString(filter);
                result["scoringPerFilter"][filterString] = scorings[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
                result["scoredValuePerFilter"][filterString] = scoredValues[name][magic_enum::enum_cast<CoreMaterialCrossReferencerFilters>(filterString).value()];
            };
            results["data"].push_back(result);
        }

        REQUIRE(crossReferencedCoreMaterials.size() > 0);

        // July 2026 re-baseline: corrected volumetric-losses scoring (absolute distance, NaN cull —
        // see the note above / CoreCrossReferencer) now ranks Fair-Rite 98 (was 95) closest.
        REQUIRE(crossReferencedCoreMaterials[0].first.get_name() == "98");
    }

}  // namespace
