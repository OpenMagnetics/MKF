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
// Refreshed 2026-05-24 after commit 3c851744 (fix(coil): two OOB reads in
// Coil::equalize_margins). That fix corrected margin sizing on multi-section
// bridge topologies, which shifts the WINDING_WINDOW_AREA filter's
// contribution downward by ~3 % for several powder toroids. The top-1 ferrite
// is unaffected; powder slots 1-4 keep the same identities and order but
// with lower scores.
// Refreshed 2026-06-16 (ABT #10): the 2026-06 saturation/loss recompute
// reshuffled the powder cross-references below the (unchanged, score-stable)
// EC 35/17/10 3C94 top-1 ferrite — Kool Mµ MAX 75 rises to slot 1, Edge 75
// drops to slot 2, and two T 27 parylene powders (OE 75 / Mix 66) displace
// the former Kool Mµ 75 / XFlux 75 in slots 3-4. The actual cross-reference
// winner (the ferrite) is unchanged; only the powder ordering moved.
const std::vector<TopEntry> kTopDefault = {
    {"EC 35/17/10 - 3C94 - Gapped 1.000 mm",                       3.0915469239085107},
    {"T 28/14/12 - epoxy coated - Kool Mµ MAX 75 - Ungapped",      2.5524218280737472},
    {"T 28/14/12 - epoxy coated - Edge 75 - Ungapped",             2.4306644789891148},
    {"T 27/14.7/14.0 - parylene coated - OE 75 - Ungapped",        2.3440327755042603},
    {"T 27/14.5/14.6 - parylene coated - Mix 66 - Ungapped",       2.052278908151667},
};

// Refreshed 2026-06-16 (ABT #10): identical ordering, scores nudged < 0.3 %
// by the 2026-06 loss/saturation recompute.
const std::vector<TopEntry> kTopSameMaterial = {
    {"EP 20 - 3C91 - Gapped 0.605 mm",                             2.326855381298027},
    {"EC 41/19/12 - 3C91 - Gapped 1.000 mm",                       2.2856989762412012},
    {"EP 17 - 3C91 - Gapped 0.414 mm",                             1.6353493768514515},
    {"RM 8/I - 3C91 - Gapped 0.480 mm",                            1.5594438851744803},
    {"EP 17 - 3C91 - Gapped 0.255 mm",                             1.4153519258586909},
};

// Refreshed 2026-05-24: 3c851744 reshuffled slots 1 and 4 — the
// EC 35/17/10 distributed-gap N27 moved from slot 4 up past
// E 32/16/9 N87, which dropped to slot 4. ETD 29/16/10 stays on top.
// Refreshed 2026-06-16 (ABT #10): identical ordering, scores nudged < 0.1 %
// by the 2026-06 loss/saturation recompute.
const std::vector<TopEntry> kTopOnlyTdk = {
    {"ETD 29/16/10 - N87 - Gapped 1.000 mm",                       2.8471783659710788},
    {"EC 35/17/10 - N27 - Distributed gapped 0.500 mm",            2.6358944598653742},
    {"ETD 29/16/10 - N27 - Gapped 1.000 mm",                       2.609290926957466},
    {"ETD 29/16/10 - N27 - Distributed gapped 0.500 mm",           2.5223183281041655},
    {"E 32/16/9 - N27 - Gapped 1.000 mm",                          1.9187849608198064},
};

// Phase 1 fix landed: with SATURATION now active, the powder ordering
// shifts. FS 60 (closest losses to reference) no longer monopolises the
// top; MS 60 wins on combined losses + saturation distance. The
// previous CORE_LOSSES ceiling-clamp removal also still applies — no
// candidates tie at the artificial 2.6 maximum.
// Refreshed 2026-05-24: 3c851744 moved Edge 60 into slot 2 (was absent
// from top-5) and pushed Kool Mµ MAX 60 out of top-5. MS / FS / Kool Mµ
// keep their identities; Kool Mµ Hƒ drops two slots.
const std::vector<TopEntry> kTopPowder = {
    {"E 25/9.5/6.3 - MS 60 - Ungapped",                            2.9338806063397938},
    {"E 25/9.5/6.3 - FS 60 - Ungapped",                            2.7454801535568123},
    {"E 25/9.5/6.3 - Edge 60 - Ungapped",                          2.5837166775568221},
    {"E 25/9.5/6.3 - Kool Mµ 60 - Ungapped",                       2.4001098422524034},
    {"E 25/9.5/6.3 - Kool Mµ Hƒ 60 - Ungapped",                    2.1310286918124186},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("CoreCrossReferencer default ferrite top-5 snapshot",
          "[adviser][core-cross-referencer][characterisation][heavy][default]") {
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
          "[adviser][core-cross-referencer][characterisation][heavy][same-material]") {
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
          "[adviser][core-cross-referencer][characterisation][heavy][only-tdk]") {
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
          "[adviser][core-cross-referencer][characterisation][heavy][powder]") {
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
