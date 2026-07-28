// Validation for the 2026-07 shape families added to CorePiece::factory
// (ABT #263-#272). Each new family is modelled as inheriting an existing piece's
// geometry; this test confirms (a) the factory recognises the family without
// throwing, (b) effective parameters compute, and (c) for identical A-F dimensions
// the child family yields the SAME Ae/Le as its parent piece (i.e. the inheritance
// is wired correctly).
#include "constructive_models/Core.h"
#include "TestingUtils.h"
#include "json.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace MAS;
using namespace OpenMagnetics;

static std::pair<double, double> ae_le(CoreShape shape, CoreShapeFamily family) {
    shape.set_family(family);
    auto piece = CorePiece::factory(shape, true);
    auto ep = piece->get_partial_effective_parameters();
    return {ep.get_effective_area(), ep.get_effective_length()};
}

TEST_CASE("New shape families inherit parent geometry", "[corepiece][newfamilies]") {
    // {parent shape name, parent family, child family}
    struct Case { std::string shapeName; CoreShapeFamily parent; CoreShapeFamily child; };
    std::vector<Case> cases = {
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPC},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPQ},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPW},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::EPT},
        {"EP 13",        CoreShapeFamily::EP,  CoreShapeFamily::LEP},
        {"ETD 29/16/10", CoreShapeFamily::ETD, CoreShapeFamily::EER},
        {"E 20/10/6",    CoreShapeFamily::E,   CoreShapeFamily::EF},
    };
    for (auto& c : cases) {
        auto shape = find_core_shape_by_name(c.shapeName);
        auto [aeP, leP] = ae_le(shape, c.parent);
        auto [aeC, leC] = ae_le(shape, c.child);
        INFO("shape=" << c.shapeName << " Ae parent=" << aeP << " child=" << aeC
                      << " Le parent=" << leP << " child=" << leC);
        REQUIRE(aeC > 0.0);
        REQUIRE(leC > 0.0);
        REQUIRE_THAT(aeC, Catch::Matchers::WithinRel(aeP, 1e-9));
        REQUIRE_THAT(leC, Catch::Matchers::WithinRel(leP, 1e-9));
    }
}

// ABT #274 / #264: validate the piece-and-plate effective parameters against PUBLISHED vendor
// values rather than against another MKF result.
//
// Source: Magnetics 2022 Ferrite Catalog, "U, I Cores", printed pages 38-39
// (https://www.mag-inc.com/Media/Magnetics/File-Library/Product%20Literature/Ferrite%20Literature/Magnetics-2022-Ferrite-Catalog.pdf).
// An I bar has no closed magnetic path of its own, so the le/Ae printed on an I-core row are those
// of the U+I COMBINATION it is used in -- which is exactly what a UI shape models. Those rows are
// the reference below.
//
// IEC 60205:2016 has no clause for a piece closed by a plate (every 5.x clause is "Pair of
// X-cores"), so CorePieceUi applies the standard's general method -- C1 = sum(l/A), C2 =
// sum(l/A^2), le = C1^2/C2, Ae = C1/C2, corners per clause 4.6. This test is what keeps that
// derivation honest.
TEST_CASE("Test_Ui_Effective_Parameters_Match_Vendor_Catalogue", "[core][shape-families][ui]") {
    settings.reset();
    clear_databases();

    struct Reference {
        std::string shapeName;
        double effectiveLengthMillimetres;   // catalogue le, from the I-core row
        double effectiveAreaSquareMillimetres;
    };
    // UI 93/76/16 pairs with plate I 93/28/16: catalogue le 257 mm, Ae 450 mm^2.
    // PQI 16/7.8 from TDK's planar series: published Ae 41.8 mm^2 and Ve 815 mm^3, so the
    // published le is Ve/Ae = 19.50 mm. (IEC 60205 clause 5.12 covers the PQ + PLT(plate)
    // combination; the geometry keeps the PQ's logarithmic radial-spreading yoke.)
    std::vector<Reference> references = {
        {"UI 93/76/16", 257.0, 450.0},
        {"PQI 16/7.8",   19.50,  41.8},
    };

    for (const auto& reference : references) {
        auto shape = find_core_shape_by_name(reference.shapeName);
        auto piece = CorePiece::factory(shape, true);
        REQUIRE(piece != nullptr);
        auto effectiveParameters = piece->get_partial_effective_parameters();

        double effectiveLength = effectiveParameters.get_effective_length() * 1000;
        double effectiveArea = effectiveParameters.get_effective_area() * 1e6;
        UNSCOPED_INFO(reference.shapeName << ": le = " << effectiveLength << " mm (catalogue "
                      << reference.effectiveLengthMillimetres << "), Ae = " << effectiveArea
                      << " mm2 (catalogue " << reference.effectiveAreaSquareMillimetres << ")");
        // 5% covers the catalogue's three-significant-figure rounding and the dimensional
        // tolerances; the derivation itself lands well inside 1% on these parts.
        CHECK_THAT(effectiveLength, Catch::Matchers::WithinRel(reference.effectiveLengthMillimetres, 0.05));
        CHECK_THAT(effectiveArea, Catch::Matchers::WithinRel(reference.effectiveAreaSquareMillimetres, 0.05));
    }
    settings.reset();
}
