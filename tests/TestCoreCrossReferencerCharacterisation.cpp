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
// DEFAULT_FERRITE only: 1e-4, approved by Alf on 2026-11-04 ("Re-pin all 5, kRelTol 1e-4"). Its slot-0 score
// had moved by 1.3e-6 relative with the same core, which failed the 1e-6 pin on its own. The core NAMES are
// still compared exactly, so an identity flip fails as before.
constexpr double kRelTolDefault = 1e-4;

struct TopEntry {
    std::string name;
    double score;
};

void check_top_n(const std::string& label,
                 const std::vector<std::pair<Core, double>>& got,
                 const std::vector<TopEntry>& expectedTop,
                 double relTol = kRelTol) {
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
        REQUIRE_THAT(got[i].second, WithinRel(expectedTop[i].score, relTol));
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
// Re-pinned 2026-07-31 (user-approved). Same phenomenon as ABT #398: these scores are MIN-MAX
// NORMALIZED over the candidate population, so they move whenever the set of candidates changes,
// even when nothing about the model or the winning core does. The top-1 core NAME is unchanged in
// all three scenarios here — only the scores moved.
//   DEFAULT_FERRITE   3.09154692390851 -> 3.09155142828103  (0.00015%)
//   SAME_MATERIAL     2.32685538129803 -> 2.38560929615197  (2.5%)
//   ONLY_TDK          2.84717836597108 -> 2.84601006322993  (0.04%)
// DEFAULT_FERRITE was a hard "[GAP_INVALID_DIMENSIONS] Gap Area is not set" exception until the
// seven Magnetics gap lengths were corrected in MAS 9561a60 (ABT #407); those cores could not be
// built before, so they were absent from the population and are now in it. The other two carried
// these exact numbers both before and after that fix, so they are older drift, not its
// consequence. Tolerance here is 0.0001%, which no population-normalized score can hold across a
// growing catalogue — #398 tracks making the score independent of the population.
// Refreshed 2026-08-20 (ABT #834, user-approved): the same catalogue-growth
// normalization shift as the CoreAdviser tables (135 MAS data/ commits since the
// last refresh; physics-reverted byte-identity verified — slot-1 score
// 2.5615519709070664 with AND without the #832 loss fixes). Slots 0-1 keep their
// identities with scores +0.35%; the Sinomag E 35 - MS 75 (6a97c76) enters at
// slot 3, displacing the Edge 75 toroid from the five.
// Re-pinned 2026-11-04 (Alf approved: "Re-pin all 5, kRelTol 1e-4"; see kRelTolDefault).
// Moved by c0a840c0 (2026-10-11, ABT #1739: mu_i(T) read from the measured tables, no floor at 1; MAS -> 8e03cd0
// with Ferroxcube's measured mu_i(T) and new Curie temperatures). Bisected on clean builds over 15ce7e25..8e935f20:
// fac0d777 passes the old pin, c0a840c0 gives these values exactly. The slot-0 core is unchanged (score -1.3e-6);
// Kool Mu MAX 75 and OE 75 drop out because only the top int(1.1 x 5) = 5 geometric candidates reach the loss and
// saturation scoring (they are 6th and later now; the order depends on the requested count, see apply_filters).
// Old: EC 35/17/10 3C94 3.0915514 | T 28/14/12 Kool Mµ MAX 75 2.5615520 | T 27/14.7/14.0 OE 75 2.4040048 |
//      E 35 MS 75 2.3430703 | T 27/14.5/14.6 Mix 66 2.1095916.
// New: EC 35/17/10 3C94 3.0915474 | T 27/14.5/14.6 Mix 66 2.5660300 | E 35 MS 75 2.5342217 |
//      E 35 Mix 52 2.3019516 | E 35 Mix 26 2.0388255.
const std::vector<TopEntry> kTopDefault = {
    {"EC 35/17/10 - 3C94 - Gapped 1.000 mm",                     3.0915474073038323},
    {"T 27/14.5/14.6 - parylene coated - Mix 66 - Ungapped",     2.5660300148730584},
    {"E 35 - MS 75 - Ungapped",                                  2.5342217224600243},
    {"E 35 - Mix 52 - Ungapped",                                 2.3019515870195209},
    {"E 35 - Mix 26 - Ungapped",                                 2.038825529217986},
};

// Refreshed 2026-06-16 (ABT #10): identical ordering, scores nudged < 0.3 %
// by the 2026-06 loss/saturation recompute.
// Re-pinned 2026-10-31 (Alf approved): pool drift of the min-max normalised scores (#398). The same
// five cores, EP 20 top-1 (2.3856 -> 2.3646, -0.88 %); only the near-tied RM 8/I and EP 17 0.414 mm
// (1.640 vs 1.657) swap slots 2 and 3.
const std::vector<TopEntry> kTopSameMaterial = {
    {"EP 20 - 3C91 - Gapped 0.605 mm",         2.364573591123257},
    {"EC 41/19/12 - 3C91 - Gapped 1.000 mm",   2.2429106060342221},
    {"EP 17 - 3C91 - Gapped 0.414 mm",         1.6570886597615659},
    {"RM 8/I - 3C91 - Gapped 0.480 mm",        1.6397777808006964},
    {"EP 17 - 3C91 - Gapped 0.255 mm",         1.4313624007456922},
};

// Refreshed 2026-05-24: 3c851744 reshuffled slots 1 and 4 — the
// EC 35/17/10 distributed-gap N27 moved from slot 4 up past
// E 32/16/9 N87, which dropped to slot 4. ETD 29/16/10 stays on top.
// Refreshed 2026-06-16 (ABT #10): identical ordering, scores nudged < 0.1 %
// by the 2026-06 loss/saturation recompute.
// Re-pinned 2026-10-31 (Alf approved): pool drift of the min-max normalised scores (#398). The same
// five cores, ETD 29 N87 top-1 (2.8460 -> 2.8477, +0.06 %); only the near-tied EC 35 N27 distributed
// and ETD 29 N27 gapped (2.560 vs 2.568) swap slots 1 and 2.
const std::vector<TopEntry> kTopOnlyTdk = {
    {"ETD 29/16/10 - N87 - Gapped 1.000 mm",              2.8476771864923007},
    {"ETD 29/16/10 - N27 - Gapped 1.000 mm",              2.5682675365229759},
    {"EC 35/17/10 - N27 - Distributed gapped 0.500 mm",   2.5595460929649891},
    {"ETD 29/16/10 - N27 - Distributed gapped 0.500 mm",  2.5264140805586299},
    {"E 32/16/9 - N27 - Gapped 1.000 mm",                 1.9153230627127971},
};

// Phase 1 fix landed: with SATURATION now active, the powder ordering
// shifts. FS 60 (closest losses to reference) no longer monopolises the
// top; MS 60 wins on combined losses + saturation distance. The
// previous CORE_LOSSES ceiling-clamp removal also still applies — no
// candidates tie at the artificial 2.6 maximum.
// Refreshed 2026-05-24: 3c851744 moved Edge 60 into slot 2 (was absent
// from top-5) and pushed Kool Mµ MAX 60 out of top-5. MS / FS / Kool Mµ
// keep their identities; Kool Mµ Hƒ drops two slots.
// Re-pinned 2026-10-31 (Alf approved): pool drift of the min-max normalised scores (#398). The same
// five cores in the same order; MS 60 top-1 2.9339 -> 2.8271 (-3.6 %).
// Re-pinned 2026-11-03 (ABT #2234, Alf-approved): MAS c7a7cbc (Magnetics provenance pass) replaced the
// Kool Mµ Hƒ 60 loss coefficient set and the XFlux 60 / Kool Mµ 60 densities. Kool Mµ Hƒ 60 climbs to slot 2,
// Edge 60 drops to slot 4; the same five cores, MS 60 top-1 2.82713 -> 2.82817.
const std::vector<TopEntry> kTopPowder = {
    {"E 25/9.5/6.3 - MS 60 - Ungapped",                            2.8281707080848451},
    {"E 25/9.5/6.3 - FS 60 - Ungapped",                            2.7197819506350833},
    {"E 25/9.5/6.3 - Kool Mµ Hƒ 60 - Ungapped",                    2.539263790197857},
    {"E 25/9.5/6.3 - Kool Mµ 60 - Ungapped",                       2.3378351612339765},
    {"E 25/9.5/6.3 - Edge 60 - Ungapped",                          2.1019014034044576},
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
    check_top_n("DEFAULT_FERRITE", results, kTopDefault, kRelTolDefault);
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
          "[!benchmark][benchmark-core-xref]") {
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
          "[!benchmark][benchmark-core-xref-powder]") {
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
