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
// Refreshed 2026-07-07 (ABT #118 loss-factor band clamp, user-approved re-pin):
// the interpolation no longer extrapolates a material's loss factor beyond its
// measured frequency band, so materials whose scores rode on extrapolated
// (lower) loss values lose that unearned advantage. FerriteDefault: same
// ranking, ~1e-6 drifts. FerriteOnlyTdk: slots 4-5 change (N27/N92 -> N51/N41),
// scores ~2%. PowderDefault: Kool Mu Ultra 26 drops from slot 2 to slot 5.
// PowderOnlyMicrometals: SM 40 and SM 60 swap slots 1-2. Regenerated with
// kRegenerateBaselines on main 17e1f850.
// Re-baselined 2026-07-14 (ABT #224 + advance MAS to latest main 489ddd9, user-approved re-pin).
// The MAS bump is 232c2f0 (TDG import) -> latest main: it also adds Changsung (CSC) powder cores,
// Ferroxcube 3E/3F37, TDK PC40/PC44/PC90, VAC nanocrystalline (ABT #213-216, #189). Effects:
// FerriteDefault: TPW33 (new TDG DMR95-class MnZn) is slot 0 — at 25 °C its μr (3325) tracks 3C97
//   (3341) even closer than DMR95 (3480) with identical saturation & losses. Slots 1-3 unchanged;
//   slot 4 is now "T" (new Magnetics MnZn power ferrite), a marginally closer match than P492.
// FerriteOnlyTdk: N95 overtakes PC47 at slot 0. The previous PC47 pin's justification compared
//   permeabilities at MISMATCHED temperatures (each material's first datasheet point: 3C97 at
//   -40 °C, PC47/N95 at -60 °C). At the actual 25 °C evaluation 3C97≈3341, N95≈3008 (closer),
//   PC47≈2330, and N95 is also closer on volumetric losses (both weight-1 filters) — so N95 is
//   the correct 25 °C match. Slots 1-3 unchanged; slot 4 is now PC44 (new TDK power ferrite).
// PowderDefault: CSC Sendust 26 (Changsung Sendust, same FeSiAl alloy and μ=26 as the Kool Mµ MAX
//   26 reference) takes slot 0 as a near-identical same-family match. NPH-L 26 (Poco High-Flux)
//   also participates now that the curie-NaN crash is fixed (6 Poco NPN/NPU materials gained
//   curieTemperature in MAS). Kool Mµ MAX 40 falls off the top-5.
// Re-baselined 2026-07-27 (ABT #190c, user-approved re-pin) after advancing MAS past
// 489ddd9. Only the two FERRITE tables move; both POWDER tables reproduce to the last
// digit, which is the evidence that this is a material-set effect and not a scoring
// change. What changed, and why it is benign:
//  * FerriteDefault reshuffles inside a statistical tie — DMR95 2.70371, ML33D 2.70170,
//    P45 2.70033, 3C95 2.68183, TPW33 2.67987 span 0.9%, with 0.007% between the top
//    two. DMR95 leading is if anything the more defensible result: it is the classic
//    3C97 equivalent (mu_i ~3480 vs 3341, same 0.53 T saturation, comparable losses).
//    Because the head of this list is a coin flip, the smoke-level tests now assert
//    shortlist MEMBERSHIP rather than an exact winner; this file keeps the full pin on
//    purpose, as a refactor tripwire.
//  * FerriteOnlyTdk keeps N95 at slot 0 with the same 25 °C justification as before;
//    slots 1-4 gain PEM95 and PCL47, TDK materials that entered MAS after the last
//    baseline, pushing N51/PC44 out of the top five.
//  * Every score drifts (1e-5 to 0.9%) because scorings are normalised across the whole
//    candidate set, so any material added to — or corrected in — MAS moves the
//    denominators for everyone, including the A10/A102 permeability extension (ABT #178)
//    and the POCO toroid fix (ABT #306) that landed the same day.
// Re-pinned 2026-07-31 (ABT #398, user-approved). The JFE (MAS 0d65a56) and TDG
// (MAS 56d8c84) material batches dropped six more MnZn power ferrites into a tie
// band that is only 0.95% wide, so the head of this list reshuffled: MBT2 and
// TPW30 enter at slots 3-4, pushing 3C95 and TPW33 out of the top five (they now
// sit 11th and 12th, at 2.68183 and 2.67987). DMR95 and P45 keep their scores to
// the last digit; ML33D's rose 2.70170 -> 2.70557 without its own data changing,
// because these scores are MIN-MAX NORMALIZED over the whole candidate population
// (CrossReferencerCommon.h::compute_normalized_scorings) — every material added to
// the catalogue moves the normalization, and therefore everyone else's score. That
// is why an absolute-score snapshot cannot hold still while the catalogue grows;
// #398 tracks making the score discriminate (or reporting an explicit tie band).
// Refreshed 2026-08-20 (ABT #834, user-approved): same catalogue-growth normalization
// shift (135 MAS data/ commits; physics-reverted byte-identity verified — slot-0 score
// 2.7030235692090656 with AND without the #832 loss fixes). The head coin flip flipped
// again: DMR95 over ML33D by 0.003% — both materials verified to carry real Steinmetz
// loss models and initial permeability, so neither wins by missing data. Same five
// names; TPW30/MBT2 swap. The tie band remains ~0.1% wide; #398 still tracks making
// the score discriminate instead of pinning a coin flip.
const std::vector<TopEntry> kTopFerriteDefault = {
    {"DMR95", 2.7030235692090656},
    {"ML33D", 2.7029322665532569},
    {"TPW30", 2.7012756176446233},
    {"MBT2",  2.6999765288775177},
    {"P45",   2.6996579563546494},
};

// Refreshed 2026-08-20 (ABT #834): identical five names and order; scores +0.27%
// from the same normalization shift.
const std::vector<TopEntry> kTopFerriteOnlyTdk = {
    {"N95",   2.5632063030546171},
    {"PEM95", 2.5509949258323354},
    {"PCL47", 2.5055538730926847},
    {"PC47",  2.4923030895288774},
    {"N97",   2.4920803486041487},
};

// Refreshed 2026-08-20 (ABT #834, user-approved): same normalization shift. Kool Mµ
// MAX 40 enters at slot 0 over CSC Sendust 26 on a 0.06% margin (another tie-band
// head; both carry real 'magnetics' loss models), Kool Mµ MAX 19 climbs to slot 2,
// and NPH-L 26 leaves the five.
const std::vector<TopEntry> kTopPowderDefault = {
    {"Kool Mµ MAX 40", 2.7166417605966702},
    {"CSC Sendust 26", 2.7150314971044129},
    {"Kool Mµ MAX 19", 2.7046770516050276},
    {"Kool Mµ Hƒ 26",  2.7034810625383088},
    {"Kool Mµ 26",     2.7031674272762469},
};

const std::vector<TopEntry> kTopPowderOnlyMicrometals = {
    {"SM 40", 2.6346308391155557},
    {"SM 60", 2.6285960369736614},
    {"SP 26", 2.6095998852559847},
    {"SM 26", 2.5954523252716264},
    {"OC 26", 2.5646142299287216},
};

} // namespace

// =============================================================================
// CHARACTERISATION SNAPSHOTS
// =============================================================================

TEST_CASE("CoreMaterialCrossReferencer 3C97 default top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][heavy][ferrite-default]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref(
        std::map<std::string, std::string>{{"coreLosses", "STEINMETZ"}});
    auto ref = Core::resolve_material("3C97");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("FERRITE_DEFAULT_3C97", results, kTopFerriteDefault);
}

TEST_CASE("CoreMaterialCrossReferencer 3C97 only-TDK top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][heavy][only-tdk]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref;
    xref.use_only_manufacturer("TDK");
    auto ref = Core::resolve_material("3C97");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("FERRITE_ONLY_TDK_3C97", results, kTopFerriteOnlyTdk);
}

TEST_CASE("CoreMaterialCrossReferencer Kool Mµ MAX 26 default top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][heavy][powder-default]") {
    settings.reset();
    clear_databases();

    CoreMaterialCrossReferencer xref;
    auto ref = Core::resolve_material("Kool M\xC2\xB5 MAX 26");
    auto results = xref.get_cross_referenced_core_material(ref, 25, 5);
    check_top_n("POWDER_DEFAULT_KOOLMU_MAX_26", results, kTopPowderDefault);
}

TEST_CASE("CoreMaterialCrossReferencer Kool Mµ MAX 26 only-Micrometals top-5 snapshot",
          "[adviser][core-material-cross-referencer][characterisation][heavy][only-micrometals]") {
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
          "[!benchmark][benchmark-mat-xref]") {
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
          "[!benchmark][benchmark-mat-xref-powder]") {
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
