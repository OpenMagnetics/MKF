// ABT #933 — WE-HCC ("High Current Cube"): a rod core inside a shell, modelled with the
// EXISTING `molded` family. No new shape family, no schema change.
//
// THE PART. WE-HCC is an air-wound solenoid slipped over a cylindrical rod core; coil and rod
// are then potted into a rectangular block that forms the outer flux return. Two construction
// variants ship under the series: 7443 31xxxx uses an iron-powder rod in an iron-powder block
// (ACE), while 7443 32/33/34xxxx use a NiZn rod in a MnZn block (TAK). RedExpert labels these
// `Core_Material` "Iron Powder" and "Ferrite" respectively; they are the two MATERIAL CLASSES
// below. Three body sizes exist: 8070, 1090 and 1210 (the last in both material classes).
//
// WHY `molded` AND NOT A NEW FAMILY. `molded` (ABT #357) is already exactly this topology: a
// coil inside a magnetic body, post + two plates + rectangular return shell, one closed solid
// with no discrete gaps. Its letters map one-to-one onto the HCC construction:
//     A  body width                 -> package length  a
//     B  body height (COIL AXIS)    -> package width   b   (the coil axis lies in the PCB
//                                                           plane, along b — see below)
//     C  body depth                 -> package height  c
//     D  cavity height              -> coil axial length
//     E  cavity outer diameter      -> coil outer diameter
//     F  cavity inner diameter      -> THE ROD (the winding post)
// The glue lines between rod, coil and block are a distributed gap, which is precisely what
// `molded` represents ("the distributed gap lives in the MATERIAL, mu_eff ~15-40"). The pooled
// permeabilities fitted below (23.4 and 27.7) land inside that band, so the family's own stated
// domain of validity covers this part. A two-material split (rod grade != block grade) would be
// physically truer for the ferrite variant, but nothing published separates the two grades, so
// fitting them would mean two free permeabilities where the data supports one. Not done.
//
// WHERE THE DIMENSIONS COME FROM — and the trap. The list-of-parts sheets carry columns a,b,c
// (package L/W/H, which match RedExpert's published package dims on all 47 parts) plus d,e,f.
// It is TEMPTING to read d as the rod diameter and f as the coil length; that reading is WRONG.
// Measured off the WE mechanical drawings (which are exact: the 1090 drawing is 2:1 to within
// 0.07%), d/e/f are all TERMINAL/CLIP geometry:
//     d = width of the clip band, measured along b
//     e = clip length along a, with e + <ref> + e == a EXACTLY on all three sizes
//         (1.5+5.4+1.5=8.4, 1.6+7.7+1.6=10.9, 2.0+8.1+2.0=12.1)
//     f = (b - d)/2 EXACTLY on all three sizes
//         (7.9-2.3)/2=2.8, (10.0-3.0)/2=3.5, (11.4-3.5)/2=3.95
// f being an exact function of b and d proves it carries no independent information, so it
// cannot also be a coil length. Two further independent falsifications of "d = rod diameter":
// the coil is DRAWN on the drawings but never dimensioned (and the drawn coil does not even
// scale between sizes — it is a generic symbol); and the saturation data gives an identical
// Bsat*Ae for the 1090 and 1210 ferrite bodies (8.78 vs 8.73 uWb, 0.7% apart) although their d
// differs by 17%, i.e. 36% in area. The internal magnetic geometry is simply not published.
//
// THE RECONSTRUCTION (this is the recipe consumers must use):
//     A = a,  B = b,  C = c                        (published package dimensions)
//     F = ROD_DIAMETER_OVER_PACKAGE_HEIGHT * c     (one structural constant, see below)
//     E = F + 2 * wireDiameter                     (a single layer of wire over the rod)
//     D = numberTurns * wireDiameter               (single-layer solenoid axial length)
// N and the wire diameter are published per part in the list of parts, so D and E vary
// part-to-part exactly as the real coil does. Only F needed to be supplied.
//
// F IS NOT FITTED TO INDUCTANCE. It is anchored on the SATURATION data, which is independent of
// the inductance the test below checks. Bsat*Ae = L*Isat/N is measurable per part; dividing by
// pi/4*(k*c)^2 gives the implied Bsat as a function of k. At k = 0.60 the three ferrite bodies
// return 0.303 T (8070), 0.359 T (1090) and 0.342 T (1210) — mutually consistent and a textbook
// MnZn/NiZn saturation flux density; the iron-powder body returns 0.568 T, which is right for a
// 30%-roll-off point on a soft-saturating powder. Neighbouring k values do not do this: k=0.45
// implies 0.54-0.64 T for ferrite (too high) and k=0.68 implies 0.24-0.28 T (too low).
// The choice is in any case not load-bearing: sweeping k from 0.35 to 0.68 (a 2x range in rod
// diameter) keeps every group median inside 0.91..1.07 and all 47 parts inside the per-part
// gate. The cross-body scaling is carried by the package dimensions, not by k.
//
// WHY THE COIL AXIS IS b. The WE drawings show the coil through the body in the top view: the
// turns stack along the package WIDTH b, so the rod lies in the PCB plane, not standing up. The
// choice barely matters numerically (swapping the a/b roles moves every group median by <0.03),
// but B is the coil-axis letter in `molded` and this is what the drawings show.
//
// WHAT THIS MODEL CANNOT DO, stated plainly. Within one body and one material class, parts with
// IDENTICAL turns, IDENTICAL wire and IDENTICAL published material differ in inductance by up
// to 1.5x (7443320022 vs 7443320033: both 2.5 turns of 1.5 mm wire, 0.22 vs 0.33 uH). WE
// evidently grades the rod permeability part by part, and that grade is not published anywhere.
// No single-permeability model can reproduce per-part inductance for this family, and one that
// appeared to would be fitting noise. The model targets the body-level median and accepts the
// residual grade scatter, which is what the per-part 0.5..2.0 band below is sized for.
//
// The permeabilities are PINNED constants, not refitted at run time, so this is a genuine
// regression test on CorePieceMolded's geometry: if the sectioning changes, the ratios move and
// this test fails.
#include "constructive_models/Core.h"
#include "physical_models/Reluctance.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <string>
#include <vector>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

// The one structural constant: rod diameter as a fraction of the package height, anchored on
// the saturation data (see header).
constexpr double ROD_DIAMETER_OVER_PACKAGE_HEIGHT = 0.60;

// Pooled effective permeability per material class — the geometric mean over every part of the
// class of the permeability each part would need. ONE value per class, never per part.
constexpr double MU_EFFECTIVE_FERRITE = 27.7;
constexpr double MU_EFFECTIVE_IRON_POWDER = 23.4;

// Acceptance gates (the heimdall gates).
constexpr double PER_PART_RATIO_MINIMUM = 0.5;
constexpr double PER_PART_RATIO_MAXIMUM = 2.0;
constexpr double GROUP_MEDIAN_RATIO_MINIMUM = 0.7;
constexpr double GROUP_MEDIAN_RATIO_MAXIMUM = 1.4;

enum class HccMaterialClass { FERRITE, IRON_POWDER };

struct HccPart {
    std::string orderCode;
    double packageLength;      // a, mm
    double packageWidth;       // b, mm  (the coil axis)
    double packageHeight;      // c, mm
    double numberTurns;
    double wireDiameter;       // mm
    double vendorInductance;   // uH
    HccMaterialClass materialClass;
};

constexpr auto FERRITE = HccMaterialClass::FERRITE;
constexpr auto IRON_POWDER = HccMaterialClass::IRON_POWDER;

// All 47 catalogue parts. Package dimensions from RedExpert (identical to the list-of-parts
// a/b/c on 47/47 parts); turns, wire diameter and inductance from the list of parts.
const std::vector<HccPart> WE_HCC_PARTS = {
    {"7443310022", 12.1, 11.4, 9.5,  2.5, 1.5,  0.22, IRON_POWDER},
    {"7443310033", 12.1, 11.4, 9.5,  2.5, 1.5,  0.33, IRON_POWDER},
    {"7443310047", 12.1, 11.4, 9.5,  3.5, 1.5,  0.47, IRON_POWDER},
    {"7443310068", 12.1, 11.4, 9.5,  3.5, 1.5,  0.68, IRON_POWDER},
    {"7443310082", 12.1, 11.4, 9.5,  4.5, 1.3,  0.82, IRON_POWDER},
    {"7443310100", 12.1, 11.4, 9.5,  4.5, 1.3,  1.00, IRON_POWDER},
    {"7443310150", 12.1, 11.4, 9.5,  5.5, 1.1,  1.50, IRON_POWDER},
    {"7443310220", 12.1, 11.4, 9.5,  6.5, 1.0,  2.20, IRON_POWDER},
    {"7443310330", 12.1, 11.4, 9.5,  8.5, 0.8,  3.30, IRON_POWDER},
    {"7443310390", 12.1, 11.4, 9.5,  9.5, 0.7,  3.90, IRON_POWDER},
    {"7443310470", 12.1, 11.4, 9.5,  9.5, 0.7,  4.70, IRON_POWDER},
    {"7443320022", 12.1, 11.4, 9.5,  2.5, 1.5,  0.22, FERRITE},
    {"7443320033", 12.1, 11.4, 9.5,  2.5, 1.5,  0.33, FERRITE},
    {"7443320047", 12.1, 11.4, 9.5,  3.5, 1.4,  0.47, FERRITE},
    {"7443320068", 12.1, 11.4, 9.5,  3.5, 1.4,  0.68, FERRITE},
    {"7443320082", 12.1, 11.4, 9.5,  4.5, 1.3,  0.82, FERRITE},
    {"7443320100", 12.1, 11.4, 9.5,  4.5, 1.3,  1.00, FERRITE},
    {"7443320150", 12.1, 11.4, 9.5,  5.5, 1.1,  1.50, FERRITE},
    {"7443320220", 12.1, 11.4, 9.5,  6.5, 1.0,  2.20, FERRITE},
    {"7443320330", 12.1, 11.4, 9.5,  7.5, 0.9,  3.30, FERRITE},
    {"7443320470", 12.1, 11.4, 9.5,  8.5, 0.8,  4.70, FERRITE},
    {"7443320680", 12.1, 11.4, 9.5,  9.5, 0.7,  6.80, FERRITE},
    {"7443320820", 12.1, 11.4, 9.5,  9.5, 0.7,  8.20, FERRITE},
    {"7443321000", 12.1, 11.4, 9.5, 10.5, 0.6, 10.00, FERRITE},
    {"7443330022", 10.9, 10.0, 9.3,  2.5, 1.4,  0.22, FERRITE},
    {"7443330033", 10.9, 10.0, 9.3,  2.5, 1.4,  0.33, FERRITE},
    {"7443330047", 10.9, 10.0, 9.3,  3.5, 1.4,  0.47, FERRITE},
    {"7443330068", 10.9, 10.0, 9.3,  4.5, 1.2,  0.68, FERRITE},
    {"7443330082", 10.9, 10.0, 9.3,  4.5, 1.2,  0.82, FERRITE},
    {"7443330100", 10.9, 10.0, 9.3,  4.5, 1.2,  1.00, FERRITE},
    {"7443330150", 10.9, 10.0, 9.3,  5.5, 1.0,  1.50, FERRITE},
    {"7443330220", 10.9, 10.0, 9.3,  6.5, 0.9,  2.20, FERRITE},
    {"7443330330", 10.9, 10.0, 9.3,  7.5, 0.8,  3.30, FERRITE},
    {"7443330470", 10.9, 10.0, 9.3,  8.5, 0.7,  4.70, FERRITE},
    {"7443330680", 10.9, 10.0, 9.3,  9.5, 0.6,  6.80, FERRITE},
    {"7443330820", 10.9, 10.0, 9.3,  9.5, 0.6,  8.20, FERRITE},
    {"7443331000", 10.9, 10.0, 9.3, 10.5, 0.5, 10.00, FERRITE},
    {"7443340030",  8.4,  7.9, 7.2,  3.5, 0.90,  0.30, FERRITE},
    {"7443340047",  8.4,  7.9, 7.2,  4.5, 0.90,  0.47, FERRITE},
    {"7443340068",  8.4,  7.9, 7.2,  4.5, 0.90,  0.68, FERRITE},
    {"7443340100",  8.4,  7.9, 7.2,  5.5, 0.80,  1.00, FERRITE},
    {"7443340150",  8.4,  7.9, 7.2,  6.5, 0.70,  1.50, FERRITE},
    {"7443340220",  8.4,  7.9, 7.2,  6.5, 0.70,  2.20, FERRITE},
    {"7443340330",  8.4,  7.9, 7.2,  7.5, 0.60,  3.30, FERRITE},
    {"7443340470",  8.4,  7.9, 7.2,  8.5, 0.50,  4.70, FERRITE},
    {"7443340680",  8.4,  7.9, 7.2, 10.5, 0.40,  6.80, FERRITE},
    {"7443341000",  8.4,  7.9, 7.2, 12.5, 0.35, 10.00, FERRITE},
};

// Build the inline custom `molded` shape for one HCC part. Everything in metres, as everywhere
// in MAS. This function IS the recipe consumers must reproduce.
json build_we_hcc_molded_shape(const HccPart& part) {
    double rodDiameter = ROD_DIAMETER_OVER_PACKAGE_HEIGHT * part.packageHeight;
    double cavityOuterDiameter = rodDiameter + 2 * part.wireDiameter;
    double cavityHeight = part.numberTurns * part.wireDiameter;
    return json{
        {"magneticCircuit", "closed"},
        {"type", "custom"},
        {"family", "molded"},
        {"aliases", json::array()},
        {"name", "WE-HCC " + part.orderCode},
        {"dimensions", {
            {"A", {{"nominal", part.packageLength * 1e-3}}},
            {"B", {{"nominal", part.packageWidth * 1e-3}}},
            {"C", {{"nominal", part.packageHeight * 1e-3}}},
            {"D", {{"nominal", cavityHeight * 1e-3}}},
            {"E", {{"nominal", cavityOuterDiameter * 1e-3}}},
            {"F", {{"nominal", rodDiameter * 1e-3}}}}}};
}

json build_we_hcc_core(const HccPart& part) {
    // The material NAME is irrelevant to this test: the reluctance is evaluated with the pinned
    // pooled permeability below, not with the named material's own permeability model. A real
    // consumer supplies a material whose initial permeability is the pooled value.
    return json{
        {"name", "WE-HCC " + part.orderCode},
        {"functionalDescription", {
            {"type", "closedShape"},
            {"material", "Kool Mµ 26"},
            {"shape", build_we_hcc_molded_shape(part)},
            {"gapping", json::array()},
            {"numberStacks", 1}}}};
}

double pooled_permeability(HccMaterialClass materialClass) {
    switch (materialClass) {
        case HccMaterialClass::FERRITE: return MU_EFFECTIVE_FERRITE;
        case HccMaterialClass::IRON_POWDER: return MU_EFFECTIVE_IRON_POWDER;
    }
    throw std::runtime_error("WE-HCC: unhandled material class");
}

double median_of(std::vector<double> values) {
    if (values.empty()) {
        throw std::runtime_error("WE-HCC: median of an empty sample");
    }
    std::sort(values.begin(), values.end());
    size_t middle = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[middle];
    }
    return (values[middle - 1] + values[middle]) / 2.0;
}

// L = N^2 / R, with R the ungapped-core reluctance MKF derives from the molded sectioning.
double modelled_inductance_microhenry(const HccPart& part) {
    Core core(build_we_hcc_core(part));
    auto reluctanceModel = ReluctanceModel::factory();
    double reluctance = reluctanceModel->get_ungapped_core_reluctance(
        core, pooled_permeability(part.materialClass));
    if (!(reluctance > 0)) {
        throw std::runtime_error("WE-HCC: non-positive reluctance for " + part.orderCode);
    }
    return part.numberTurns * part.numberTurns / reluctance * 1e6;
}

std::string group_key(const HccPart& part) {
    std::string size = part.orderCode.substr(0, 6);
    std::string material = part.materialClass == HccMaterialClass::FERRITE ? "Ferrite" : "Iron Powder";
    return size + " (" + material + ")";
}

}  // namespace

// The molded reconstruction must reproduce every published WE-HCC inductance to within the
// heimdall gates, using ONE pooled permeability per material class across ALL THREE bodies.
TEST_CASE("Test_We_Hcc_Molded_Reconstruction_Reproduces_Vendor_Inductance",
          "[core][molded][we-hcc][abt933]") {
    settings.reset();

    std::map<std::string, std::vector<double>> ratiosByGroup;
    for (const auto& part : WE_HCC_PARTS) {
        double modelled = modelled_inductance_microhenry(part);
        double ratio = modelled / part.vendorInductance;
        INFO("part=" << part.orderCode << " vendor=" << part.vendorInductance
             << " uH modelled=" << modelled << " uH ratio=" << ratio);
        CHECK(ratio >= PER_PART_RATIO_MINIMUM);
        CHECK(ratio <= PER_PART_RATIO_MAXIMUM);
        ratiosByGroup[group_key(part)].push_back(ratio);
    }

    REQUIRE(ratiosByGroup.size() == 4);
    for (const auto& [group, ratios] : ratiosByGroup) {
        double median = median_of(ratios);
        INFO("group=" << group << " n=" << ratios.size() << " median ratio=" << median);
        CHECK(median >= GROUP_MEDIAN_RATIO_MINIMUM);
        CHECK(median <= GROUP_MEDIAN_RATIO_MAXIMUM);
    }
    settings.reset();
}

// The pooled permeabilities must stay inside the band `molded` claims for a distributed-gap
// body (~15-40), and inside the 10-100 band that makes sense for an iron-powder rod-and-block
// assembly. A refit that walked outside these would mean the geometry recipe had drifted.
TEST_CASE("Test_We_Hcc_Pooled_Permeabilities_Are_Physical", "[core][molded][we-hcc][abt933]") {
    for (double permeability : {MU_EFFECTIVE_FERRITE, MU_EFFECTIVE_IRON_POWDER}) {
        CHECK(permeability > 10.0);
        CHECK(permeability < 100.0);
    }
}

// The reconstruction must stay inside `molded`'s own geometric constraints for every part, with
// real margin — not perched on the E < min(A, C) edge where the answer would be an artefact of
// the constraint rather than of the physics.
TEST_CASE("Test_We_Hcc_Reconstruction_Respects_Molded_Constraints", "[core][molded][we-hcc][abt933]") {
    settings.reset();
    for (const auto& part : WE_HCC_PARTS) {
        double rodDiameter = ROD_DIAMETER_OVER_PACKAGE_HEIGHT * part.packageHeight;
        double cavityOuterDiameter = rodDiameter + 2 * part.wireDiameter;
        double cavityHeight = part.numberTurns * part.wireDiameter;
        INFO("part=" << part.orderCode << " rod=" << rodDiameter << " cavityOD="
             << cavityOuterDiameter << " cavityHeight=" << cavityHeight);
        CHECK(cavityOuterDiameter > rodDiameter);
        CHECK(cavityHeight < part.packageWidth);
        CHECK(cavityOuterDiameter < std::min(part.packageLength, part.packageHeight));
        // The return path must not be the bottleneck: the shell cross-section `molded` works
        // with (footprint minus cavity circle) has to comfortably exceed the rod's. Note this
        // is an AREA criterion, not a wall-thickness one — the shell wraps the cavity on all
        // sides, so the thinnest wall (0.54 mm, on the 8070 with 0.9 mm wire) is not the
        // quantity that sets the reluctance.
        double shellArea = part.packageLength * part.packageHeight
                           - std::numbers::pi / 4 * std::pow(cavityOuterDiameter, 2);
        double rodArea = std::numbers::pi / 4 * std::pow(rodDiameter, 2);
        // The criterion is shellArea > rodArea; the 1.5x floor leaves the assertion able to
        // catch real drift. The tightest part in the catalogue (7443330033/47 — the 1090 with
        // 1.4 mm wire, the fattest wire relative to its body) achieves 1.89x.
        CHECK(shellArea > 1.5 * rodArea);
        // And the core must actually build.
        Core core(build_we_hcc_core(part));
        auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
        CHECK(effectiveParameters.get_effective_area() > 0);
        CHECK(effectiveParameters.get_effective_length() > 0);
    }
    settings.reset();
}
