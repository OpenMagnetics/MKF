// =============================================================================
// TestCoreMaterialCrossReferencerCharacterisation.cpp
// =============================================================================
// Characterisation tests for CoreMaterialCrossReferencer.
//
// PURPOSE
//   Lock the CURRENT top-N (material, score) output across a small matrix of
//   reference materials and configurations, so any refactor of the 578-line
//   CoreMaterialCrossReferencer preserves both ranking and scores to within
//   1e-6 relative tolerance.
//
// MATRIX
//   ferrite-default    : 3C97 @ 25 °C (Steinmetz loss model)
//   ferrite-only-TDK   : 3C97 @ 25 °C restricted to TDK
//   powder-default     : Kool Mµ MAX 26 @ 25 °C
//   powder-only-MM     : Kool Mµ MAX 26 @ 25 °C restricted to Micrometals
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true (or leave a snapshot vector empty), run,
//   copy BASELINE lines from stderr into the kSnapshots tables, flip back.
//
// BENCHMARKS  (tag: [!benchmark])
//   Time get_cross_referenced_core_material end-to-end. The material set is
//   small (~hundreds), so Catch2 defaults are tolerable here, but for
//   consistency with the other characterisation files we still recommend
//      --benchmark-samples 5 --benchmark-warmup-time 0
// =============================================================================

#include <cmath>
#include <iomanip>
#include <map>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "advisers/CoreMaterialCrossReferencer.h"
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

void check_top_n(const std::string& label,
                 const std::vector<std::pair<CoreMaterial, double>>& got,
                 const std::vector<TopEntry>& expectedTop) {
    if (kRegenerateBaselines || expectedTop.empty()) {
        std::cerr << "\nBASELINE TOPN " << label << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            std::cerr << "  [" << i << "] name=\"" << got[i].first.get_name()
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
        INFO("top[" << i << "] want=\"" << expectedTop[i].name << "\" got=\""
                    << got[i].first.get_name() << "\" score="
                    << std::setprecision(17) << got[i].second);
        REQUIRE(got[i].first.get_name() == expectedTop[i].name);
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, kRelTol));
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS — captured 2026-05-19. Leave a vector empty to auto-harvest.
// ----------------------------------------------------------------------------
const std::vector<TopEntry> kTopFerriteDefault = {
    {"DMR95", 2.7229156168842046},
    {"P45",   2.7216465717490372},
    {"3C95",  2.7144645242903547},
    {"P492",  2.7032348169719453},
    {"P452",  2.6934204788417202},
};

const std::vector<TopEntry> kTopFerriteOnlyTdk = {
    {"PC47", 2.6839499384443108},
    {"N95",  2.651360761361282},
    {"N97",  2.5578230120884067},
    {"N27",  2.4991932937788723},
    {"N92",  2.4814549828497379},
};

const std::vector<TopEntry> kTopPowderDefault = {
    {"Kool Mµ Hƒ 26",    2.7299983338570506},
    {"Kool Mµ Ultra 26", 2.7299790326971172},
    {"Kool Mµ 26",       2.7211505563536287},
    {"Kool Mµ MAX 19",   2.7171791666466776},
    {"Kool Mµ MAX 40",   2.7166417605966702},
};

const std::vector<TopEntry> kTopPowderOnlyMicrometals = {
    {"SM 60", 2.6337004497830838},
    {"SM 40", 2.6308271217456958},
    {"SP 26", 2.6057980867714345},
    {"SM 26", 2.5918014945962127},
    {"OC 26", 2.5610286222647689},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("CoreMaterialCrossReferencer 3C97 default top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][ferrite-default]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref(
        std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});
    auto ref = Core::resolve_material("3C97");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("FERRITE_DEFAULT_3C97", results, kTopFerriteDefault);
}

TEST_CASE("CoreMaterialCrossReferencer 3C97 only-TDK top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][only-tdk]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref;
    xref.use_only_manufacturer("TDK");
    auto ref = Core::resolve_material("3C97");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("FERRITE_ONLY_TDK_3C97", results, kTopFerriteOnlyTdk);
}

TEST_CASE("CoreMaterialCrossReferencer Kool Mµ MAX 26 default top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][powder-default]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref;
    auto ref = Core::resolve_material("Kool M\xC2\xB5 MAX 26");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("POWDER_DEFAULT_KOOLMU_MAX_26", results, kTopPowderDefault);
}

TEST_CASE("CoreMaterialCrossReferencer Kool Mµ MAX 26 only-Micrometals top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][only-micrometals]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref;
    xref.use_only_manufacturer("Micrometals");
    auto ref = Core::resolve_material("Kool M\xC2\xB5 MAX 26");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("POWDER_ONLY_MM_KOOLMU_MAX_26", results, kTopPowderOnlyMicrometals);
}

// =============================================================================
// BENCHMARKS  (opt-in via [!benchmark])
// =============================================================================

TEST_CASE("Benchmark CoreMaterialCrossReferencer 3C97 (top-20)",
          "[adviser][core-material-cross-referencer][!benchmark][benchmark-mat-xref]") {
    settings.reset();
    clear_databases();

    auto ref = Core::resolve_material("3C97");

    BENCHMARK("get_cross_referenced_core_material 3C97 top-20") {
        CoreMaterialCrossReferencer xref(
            std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});
        return xref.get_cross_referenced_core_material(ref, 25, 20);
    };
}

TEST_CASE("Benchmark CoreMaterialCrossReferencer Kool Mµ MAX 26 (top-20)",
          "[adviser][core-material-cross-referencer][!benchmark][benchmark-mat-xref-powder]") {
    settings.reset();
    clear_databases();

    auto ref = Core::resolve_material("Kool M\xC2\xB5 MAX 26");

    BENCHMARK("get_cross_referenced_core_material Kool Mu MAX 26 top-20") {
        CoreMaterialCrossReferencer xref;
        return xref.get_cross_referenced_core_material(ref, 25, 20);
    };
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                    | mean (5 samples) | notes
// -----------+-----------------------------+------------------+----------
// 2026-05-19 | 3C97 top-20                 | TBD              | initial
// 2026-05-19 | Kool Mu MAX 26 top-20       | TBD              | initial
//
// =============================================================================
