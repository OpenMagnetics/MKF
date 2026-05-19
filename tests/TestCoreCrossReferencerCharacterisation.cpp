// =============================================================================
// TestCoreCrossReferencerCharacterisation.cpp
// =============================================================================
// Characterisation tests for CoreCrossReferencer::get_cross_referenced_core.
//
// PURPOSE
//   Lock the CURRENT top-N (core, score) output across a small matrix of
//   reference cores and filter configurations, so any refactor of the
//   740-line CoreCrossReferencer preserves both the ranking AND the scores
//   to within 1e-6 relative tolerance.
//
// MATRIX
//   default          : ferrite reference (EC 35/17/10 - 3C91 - Gapped 1 mm)
//   same-material    : restricted to reference material (3C91)
//   only-manufacturer: restricted to TDK
//   powder           : powder reference (E 25/9.5/6.3 - XFlux 60 - Ungapped)
//
//   Inputs are the OpenMagneticsTesting::create_quick_test_inputs() defaults
//   used throughout TestCoreCrossReferencer.cpp, 28 turns.
//
// REGENERATING SNAPSHOTS
//   Flip kRegenerateBaselines = true (or leave a snapshot vector empty), run,
//   copy BASELINE lines from stderr into the kSnapshots tables, flip back.
//
// BENCHMARKS  (tag: [!benchmark])
//   Time get_cross_referenced_core end-to-end. Catch2 defaults to 100 samples
//   per BENCHMARK; for the full-DB scenarios use
//      --benchmark-samples 3 --benchmark-warmup-time 0
// =============================================================================

#include <cmath>
#include <iomanip>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "advisers/CoreCrossReferencer.h"
#include "processors/Inputs.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include "TestingUtils.h"

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
                 const std::vector<std::pair<Core, double>>& got,
                 const std::vector<TopEntry>& expectedTop) {
    if (kRegenerateBaselines || expectedTop.empty()) {
        std::cerr << "\nBASELINE TOPN " << label << " count=" << got.size() << "\n";
        for (size_t i = 0; i < got.size(); ++i) {
            auto name = got[i].first.get_name().value_or("<unnamed>");
            std::cerr << "  [" << i << "] name=\"" << name << "\" score="
                      << std::setprecision(17) << got[i].second << "\n";
        }
        // Structural contract still enforced even with no snapshot.
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
        auto name = got[i].first.get_name().value_or("<unnamed>");
        INFO("top[" << i << "] want=\"" << expectedTop[i].name << "\" got=\""
                    << name << "\" score=" << std::setprecision(17)
                    << got[i].second);
        REQUIRE(name == expectedTop[i].name);
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, kRelTol));
    }
}

// ----------------------------------------------------------------------------
// SNAPSHOTS — captured 2026-05-19. Leave a vector empty to auto-harvest on
// next run (the helper prints BASELINE lines to stderr).
// ----------------------------------------------------------------------------

// Phase 1 fix landed: the SATURATION filter (weight 0.5) was completely
// inactive — the apply_filters switch had no SATURATION case and no
// MagneticCoreFilterSaturation existed, so candidates were ranked on
// geometry + losses only. Powder cores with very high Bsat competed
// equally with ferrites despite hugely different saturation headroom,
// biasing every ferrite reference toward powder cross-references.
// Saturation is now scored by |refBratio - candBratio| inside
// MagneticCoreFilterCoreLosses (where the values were already computed)
// and normalised with weight 0.5. The reference ferrite shape +
// adjacent material (3C94) now correctly tops the list; powders still
// appear in slots 2-5 because they remain geometrically and
// magnetically competitive (low Bpeak/Bsat in both cases).
const std::vector<TopEntry> kTopDefault = {
    {"EC 35/17/10 - 3C94 - Gapped 1.000 mm",                       3.0915469239085107},
    {"T 28/14/12 - epoxy coated - Edge 75 - Ungapped",             2.6023810299281194},
    {"T 28/14/12 - epoxy coated - Kool Mµ 75 - Ungapped",          2.5598799018962883},
    {"T 28/14/12 - epoxy coated - Kool Mµ MAX 75 - Ungapped",      2.4488472314630836},
    {"T 28/14/12 - epoxy coated - XFlux 75 - Ungapped",            2.2430466357109586},
};

const std::vector<TopEntry> kTopSameMaterial = {
    {"EP 20 - 3C91 - Gapped 0.605 mm",                             2.3260884808139326},
    {"EC 41/19/12 - 3C91 - Gapped 1.000 mm",                       2.2891873964853477},
    {"EP 17 - 3C91 - Gapped 0.414 mm",                             1.6374042631675698},
    {"RM 8/I - 3C91 - Gapped 0.480 mm",                            1.5642380256468966},
    {"EP 17 - 3C91 - Gapped 0.255 mm",                             1.414304619475435},
};

const std::vector<TopEntry> kTopOnlyTdk = {
    {"ETD 29/16/10 - N87 - Gapped 1.000 mm",                       2.8695155571782585},
    {"E 32/16/9 - N87 - Gapped 1.000 mm",                          2.7089256051581794},
    {"ETD 29/16/10 - N27 - Gapped 1.000 mm",                       2.6299427704000218},
    {"ETD 29/16/10 - N27 - Distributed gapped 0.500 mm",           2.5625439304740238},
    {"EC 35/17/10 - N27 - Distributed gapped 0.500 mm",            2.232510837412359},
};

// Phase 1 fix landed: with SATURATION now active, the powder ordering
// shifts. FS 60 (closest losses to reference) no longer monopolises the
// top; MS 60 wins on combined losses + saturation distance. The
// previous CORE_LOSSES ceiling-clamp removal also still applies — no
// candidates tie at the artificial 2.6 maximum.
const std::vector<TopEntry> kTopPowder = {
    {"E 25/9.5/6.3 - MS 60 - Ungapped",                            2.9338806063397938},
    {"E 25/9.5/6.3 - FS 60 - Ungapped",                            2.7454801535568123},
    {"E 25/9.5/6.3 - Kool Mµ Hƒ 60 - Ungapped",                    2.6147453693692406},
    {"E 25/9.5/6.3 - Kool Mµ 60 - Ungapped",                       2.4001098422524034},
    {"E 25/9.5/6.3 - Kool Mµ MAX 60 - Ungapped",                   2.2461078440549671},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("CoreCrossReferencer default ferrite top-5 snapshot",
          "[adviser][core-cross-referencer][characterisation][default]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    CoreCrossReferencer xref;
    Core ref = find_core_by_name("EC 35/17/10 - 3C91 - Gapped 1.000 mm");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    auto results = xref.get_cross_referenced_core(ref, 28, inputs, 5);
    check_top_n("DEFAULT_FERRITE", results, kTopDefault);
}

TEST_CASE("CoreCrossReferencer same-material top-5 snapshot",
          "[adviser][core-cross-referencer][characterisation][same-material]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    CoreCrossReferencer xref;
    xref.use_only_reference_material(true);
    Core ref = find_core_by_name("EC 35/17/10 - 3C91 - Gapped 1.000 mm");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    auto results = xref.get_cross_referenced_core(ref, 28, inputs, 5);
    check_top_n("SAME_MATERIAL_3C91", results, kTopSameMaterial);
}

TEST_CASE("CoreCrossReferencer only-TDK top-5 snapshot",
          "[adviser][core-cross-referencer][characterisation][only-tdk]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    CoreCrossReferencer xref;
    xref.use_only_manufacturer("TDK");
    Core ref = find_core_by_name("EC 35/17/10 - 3C91 - Gapped 1.000 mm");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    auto results = xref.get_cross_referenced_core(ref, 28, inputs, 5);
    check_top_n("ONLY_TDK", results, kTopOnlyTdk);
}

TEST_CASE("CoreCrossReferencer powder reference top-5 snapshot",
          "[adviser][core-cross-referencer][characterisation][powder]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    CoreCrossReferencer xref;
    Core ref = find_core_by_name("E 25/9.5/6.3 - XFlux 60 - Ungapped");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();
    auto results = xref.get_cross_referenced_core(ref, 28, inputs, 5);
    check_top_n("POWDER_XFLUX60", results, kTopPowder);
}

// =============================================================================
// BENCHMARKS  (opt-in via [!benchmark])
// =============================================================================
//
// Use:  ./MKF_tests "[benchmark-core-xref]" --benchmark-samples 3 --benchmark-warmup-time 0
//

TEST_CASE("Benchmark CoreCrossReferencer default ferrite (top-50)",
          "[adviser][core-cross-referencer][!benchmark][benchmark-core-xref]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    Core ref = find_core_by_name("EC 35/17/10 - 3C91 - Gapped 1.000 mm");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();

    BENCHMARK("get_cross_referenced_core default top-50") {
        CoreCrossReferencer xref;
        return xref.get_cross_referenced_core(ref, 28, inputs, 50);
    };
}

TEST_CASE("Benchmark CoreCrossReferencer powder (top-50)",
          "[adviser][core-cross-referencer][!benchmark][benchmark-core-xref-powder]") {
    settings.reset();
    settings.set_use_only_cores_in_stock(false);
    clear_databases();

    Core ref = find_core_by_name("E 25/9.5/6.3 - XFlux 60 - Ungapped");
    auto inputs = OpenMagneticsTesting::create_quick_test_inputs();

    BENCHMARK("get_cross_referenced_core powder top-50") {
        CoreCrossReferencer xref;
        return xref.get_cross_referenced_core(ref, 28, inputs, 50);
    };
}

// =============================================================================
// BASELINE BENCHMARKS (record after each refactor)
// =============================================================================
//
// Date       | Scenario                              | mean (3 samples) | notes
// -----------+---------------------------------------+------------------+--------
// 2026-05-19 | default ferrite top-50                | TBD              | initial
// 2026-05-19 | powder XFlux60 top-50                 | TBD              | initial
//
// =============================================================================
