// ABT #958 / windability: the minimum bend radius of a wound conductor, and what the geometry
// does when the former's corner is tighter than the wire can bend.
//
// Every expected value here is transcribed from the standard named beside it, NOT captured from
// this implementation's output. IEC 60317-0-1:2013 Tables 6 and 7 (round copper), IEC 60317-0-2
// Table 6 and clause 9 (rectangular copper), and IEC 60851-3 ed.3.1 §5.1.1 for the fact that the
// wire is wound ON the mandrel, which is what puts its centreline one wire radius outside it.
#include "physical_models/WireBend.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Magnetic.h"
#include "advisers/MagneticFilter.h"
#include "support/Utils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <vector>

using json = nlohmann::json;
using namespace MAS;
using namespace OpenMagnetics;

namespace {

// Grade 1 enamel, because Wire::calculate_outer_diameter DERIVES the outer diameter from the
// conductor and the coating (it does not read a stored outerDiameter), and the bend model must
// use the very same outer dimension that turn placement uses for the standoff.
OpenMagnetics::Wire round_wire(double conductingDiameter) {
    OpenMagnetics::Wire wire;
    wire.set_type(WireType::ROUND);
    wire.set_nominal_value_conducting_diameter(conductingDiameter);
    InsulationWireCoating coating;
    coating.set_type(InsulationWireCoatingType::ENAMELLED);
    coating.set_grade(1);
    wire.set_coating(coating);
    return wire;
}

OpenMagnetics::Wire rectangular_wire(double conductingWidth, double conductingHeight) {
    OpenMagnetics::Wire wire;
    wire.set_type(WireType::RECTANGULAR);
    wire.set_nominal_value_conducting_width(conductingWidth);
    wire.set_nominal_value_conducting_height(conductingHeight);
    InsulationWireCoating coating;
    coating.set_type(InsulationWireCoatingType::ENAMELLED);
    coating.set_grade(1);
    wire.set_coating(coating);
    return wire;
}

constexpr double kTight = 1e-9;

} // namespace

TEST_CASE("Test_WireBend_Round_Flexibility_Mandrel_Is_The_Conductor", "[wire][wire-bend]") {
    // IEC 60317-0-1 Table 6, last row: 0,140 mm < d <= 1,600 mm, no pre-stretch, mandrel = d.
    auto wire = round_wire(0.400e-3);
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.400e-3, kTight));

    // IEC 60851-3 §5.1.1: wound ON the mandrel, so centreline = mandrel/2 + outer/2.
    // ...and Round 0.4 grade 1 measures 0,430 mm over the enamel, so the centreline may not come
    // inside 0,200 + 0,215 mm. The tolerance is the enamel table's own resolution, not slack.
    REQUIRE_THAT(WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.200e-3 + 0.215e-3, 1e-6));
}

TEST_CASE("Test_WireBend_Round_Fine_Wire_Uses_The_Fixed_Mandrel", "[wire][wire-bend]") {
    // Table 6: every row up to and including 0,140 mm winds on a fixed 0,150 mm mandrel.
    auto fine = round_wire(0.100e-3);
    REQUIRE_THAT(WireBend::get_mandrel_diameter(fine, BendCriterion::FLEXIBILITY, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.150e-3, kTight));
    // ...and Table 7's note sends the heat-shock case to Table 6 for those sizes too.
    REQUIRE_THAT(WireBend::get_mandrel_diameter(fine, BendCriterion::HEAT_SHOCK, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.150e-3, kTight));
}

TEST_CASE("Test_WireBend_Round_Heat_Shock_Table_Rows", "[wire][wire-bend]") {
    // IEC 60317-0-1 Table 7, spot rows across all three of its plateaus.
    struct Row { double conductor; double mandrel; };
    const std::vector<Row> rows = {
        {0.200e-3, 0.315e-3}, {0.250e-3, 0.400e-3}, {0.280e-3, 0.630e-3},
        {0.400e-3, 0.900e-3}, {1.000e-3, 2.240e-3}, {1.120e-3, 3.550e-3},
        {1.600e-3, 5.000e-3},
    };
    for (const auto& row : rows) {
        auto wire = round_wire(row.conductor);
        REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::HEAT_SHOCK, BendAxis::ROUND),
                     Catch::Matchers::WithinAbs(row.mandrel, kTight));
    }
}

TEST_CASE("Test_WireBend_Round_Intermediate_Size_Takes_The_Next_Larger_Row", "[wire][wire-bend]") {
    // The note under Table 7: an intermediate conductor takes the mandrel of the NEXT LARGER
    // tabulated diameter. 0,375 mm is between the 0,355 and 0,400 rows, so it gets 0,900 mm.
    auto wire = round_wire(0.375e-3);
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::HEAT_SHOCK, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.900e-3, kTight));
}

TEST_CASE("Test_WireBend_Round_Above_The_Table_Throws", "[wire][wire-bend]") {
    // Above 1,600 mm the standard replaces the mandrel test with a stretching test (8.2), so
    // there is no bend radius to give. It must say so rather than extrapolate.
    auto wire = round_wire(2.000e-3);
    REQUIRE_THROWS(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND));
}

TEST_CASE("Test_WireBend_Rectangular_Mandrels", "[wire][wire-bend]") {
    // IEC 60317-0-2 Table 6: on the width, 4 x width up to 10 mm; on the thickness, 4 x thickness.
    auto wire = rectangular_wire(3.0e-3, 0.5e-3);
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::EDGEWISE),
                 Catch::Matchers::WithinAbs(4.0 * 3.0e-3, kTight));
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::FLATWISE),
                 Catch::Matchers::WithinAbs(4.0 * 0.5e-3, kTight));

    // Clause 9: heat shock is bent flatwise on six times the thickness...
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::HEAT_SHOCK, BendAxis::FLATWISE),
                 Catch::Matchers::WithinAbs(6.0 * 0.5e-3, kTight));
    // ...and edgewise is simply not specified, which must be an explicit refusal.
    REQUIRE_THROWS(WireBend::get_mandrel_diameter(wire, BendCriterion::HEAT_SHOCK, BendAxis::EDGEWISE));
}

TEST_CASE("Test_WireBend_Rectangular_Over_Ten_Millimetres_Is_Five_Times", "[wire][wire-bend]") {
    // Table 6's second row: sizes over 10 mm go to 5 x width.
    auto wire = rectangular_wire(12.0e-3, 1.0e-3);
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::EDGEWISE),
                 Catch::Matchers::WithinAbs(5.0 * 12.0e-3, kTight));
}

TEST_CASE("Test_WireBend_Axis_Follows_The_Dimension_In_The_Bend_Plane", "[wire][wire-bend]") {
    // The standard's "width" is the larger dimension whichever MAS axis holds it.
    auto flat = rectangular_wire(3.0e-3, 0.5e-3);
    REQUIRE(WireBend::axis_from_bend_plane_dimension(flat, /*bendPlaneIsWidth=*/true) == BendAxis::EDGEWISE);
    REQUIRE(WireBend::axis_from_bend_plane_dimension(flat, /*bendPlaneIsWidth=*/false) == BendAxis::FLATWISE);

    auto onEdge = rectangular_wire(0.5e-3, 3.0e-3);
    REQUIRE(WireBend::axis_from_bend_plane_dimension(onEdge, /*bendPlaneIsWidth=*/true) == BendAxis::FLATWISE);
    REQUIRE(WireBend::axis_from_bend_plane_dimension(onEdge, /*bendPlaneIsWidth=*/false) == BendAxis::EDGEWISE);
}

TEST_CASE("Test_WireBend_Conforming_Leaves_The_Turn_Where_It_Was", "[wire][wire-bend]") {
    // Round 0.4 - Grade 1 on a former whose corner is at the flexibility floor (mandrel radius).
    auto wire = round_wire(0.400e-3);
    const double standoff = 0.215e-3;
    auto solution = WireBend::solve(0.200e-3, 0.25 * std::numbers::pi, wire, standoff, BendAxis::ROUND);

    REQUIRE(solution.regime == BendRegime::CONFORMING);
    REQUIRE_THAT(solution.faceStandoff, Catch::Matchers::WithinAbs(standoff, kTight));
    REQUIRE_THAT(solution.bendRadius, Catch::Matchers::WithinAbs(0.415e-3, 1e-6));
    REQUIRE_THAT(solution.cornerCentreInset, Catch::Matchers::WithinAbs(0.0, kTight));
    // At the flexibility floor but below the 0,900 mm heat-shock mandrel, so it must not read OK.
    REQUIRE(solution.verdict == BendVerdict::BELOW_HEAT_SHOCK);
}

TEST_CASE("Test_WireBend_Regime_Boundary_Is_Continuous", "[wire][wire-bend]") {
    // Exactly at Rmin the two branches must agree, or the placement jumps at the crossover.
    auto wire = round_wire(0.400e-3);
    const double standoff = 0.215e-3;
    const double flexibilityRadius =
        WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND);

    auto atBoundary = WireBend::solve(flexibilityRadius - standoff, 0.25 * std::numbers::pi, wire,
                                      standoff, BendAxis::ROUND);
    auto justInside = WireBend::solve(flexibilityRadius - standoff - 1e-9, 0.25 * std::numbers::pi,
                                      wire, standoff, BendAxis::ROUND);

    REQUIRE(atBoundary.regime == BendRegime::CONFORMING);
    REQUIRE(justInside.regime == BendRegime::LIFTED);
    REQUIRE_THAT(justInside.faceStandoff,
                 Catch::Matchers::WithinAbs(atBoundary.faceStandoff, 1e-9));
}

TEST_CASE("Test_WireBend_Sharp_Former_Lifts_The_Winding_Off", "[wire][wire-bend]") {
    // The 02_flyback case, measured: EFD25 rectangular column modelled with a sharp corner,
    // Round 0.4 - Grade 1 standing off by its own coated radius. The wire cannot make that
    // corner, so it rides the corners and stands off the faces by more than its radius.
    auto wire = round_wire(0.400e-3);
    const double standoff = 0.215e-3;
    auto solution = WireBend::solve(0.0, 0.25 * std::numbers::pi, wire, standoff, BendAxis::ROUND);

    REQUIRE(solution.regime == BendRegime::LIFTED);
    REQUIRE_THAT(solution.bendRadius, Catch::Matchers::WithinAbs(0.415e-3, 1e-6));
    REQUIRE_THAT(solution.cornerCentreInset, Catch::Matchers::WithinAbs(0.200e-3, 1e-6));
    // s + (1 - sin 45deg) * 0.200 mm = 0.215 + 0.0585786 = 0.2735786 mm
    REQUIRE_THAT(solution.faceStandoff, Catch::Matchers::WithinAbs(0.2735786e-3, 1e-6));
    REQUIRE(solution.verdict == BendVerdict::BELOW_HEAT_SHOCK);
}

TEST_CASE("Test_WireBend_No_Corner_Means_No_Lift_Off", "[wire][wire-bend]") {
    // As the corner opens out to a straight junction the lift-off term must vanish, whatever
    // the wire: (1 - sin(pi/2)) = 0.
    auto wire = round_wire(0.400e-3);
    auto solution = WireBend::solve(0.0, 0.5 * std::numbers::pi, wire, 0.215e-3, BendAxis::ROUND);
    REQUIRE(solution.regime == BendRegime::LIFTED);
    REQUIRE_THAT(solution.faceStandoff, Catch::Matchers::WithinAbs(0.215e-3, kTight));
}

TEST_CASE("Test_WireBend_Evaluate_Judges_As_Built_Geometry", "[wire][wire-bend]") {
    // What MVB++ draws today on 02_flyback: the turn wraps a sharp corner at its own coated
    // radius. That is tighter than the wire is qualified to bend at all.
    auto wire = round_wire(0.400e-3);
    auto asBuilt = WireBend::evaluate(0.2147e-3, wire, BendAxis::ROUND);
    REQUIRE(asBuilt.verdict == BendVerdict::BELOW_FLEXIBILITY);
    // Outer-fibre strain r/R, the number that says how far past the qualification bend it is.
    REQUIRE_THAT(asBuilt.outerFibreStrain, Catch::Matchers::WithinRel(0.200e-3 / 0.2147e-3, 1e-9));

    // The heat-shock radius is the comfortable one, and it must be well clear.
    auto comfortable = WireBend::evaluate(0.665e-3, wire, BendAxis::ROUND);
    REQUIRE(comfortable.verdict == BendVerdict::OK);
}

TEST_CASE("Test_WireBend_Refuses_Wire_Types_The_Standards_Do_Not_Cover", "[wire][wire-bend]") {
    // Litz bundle bending is not in IEC 60317; it must refuse rather than borrow the round rule.
    OpenMagnetics::Wire litz;
    litz.set_type(WireType::LITZ);
    REQUIRE_THROWS(WireBend::get_minimum_bend_radius(litz, BendCriterion::FLEXIBILITY, BendAxis::ROUND));
}

TEST_CASE("Test_WireBend_Catalogue_Wire_Matches_The_Hand_Numbers", "[wire][wire-bend]") {
    // End to end on the real catalogue entry rather than a hand-built wire: Round 0.4 - Grade 1
    // has no nominal outer diameter, only 0,421/0,439 mm, so this also pins that the resolver is
    // what collapses it (0,430 mm) rather than any hand-read bound.
    auto wire = find_wire_by_name("Round 0.4 - Grade 1");
    REQUIRE_THAT(WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.400e-3, 1e-9));
    REQUIRE_THAT(WireBend::get_minimum_bend_radius(wire, BendCriterion::FLEXIBILITY, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.200e-3 + 0.215e-3, 1e-6));
    REQUIRE_THAT(WireBend::get_minimum_bend_radius(wire, BendCriterion::HEAT_SHOCK, BendAxis::ROUND),
                 Catch::Matchers::WithinAbs(0.450e-3 + 0.215e-3, 1e-6));
}

// CENSUS. Read-only: it changes no geometry and asserts almost nothing. It walks the MAS example
// corpus and reports, per design, whether each winding's wire can actually follow the corner of
// the column it is wound on -- and, where it cannot, how far the winding really stands off.
// Run it with:  ./MKF_tests "[wire-bend-census]" -s
TEST_CASE("Test_WireBend_Corpus_Census", "[wire][wire-bend][wire-bend-census]") {
    const std::filesystem::path examples = "/home/alf/OpenMagnetics/MKF/MAS/examples";
    if (!std::filesystem::exists(examples)) {
        SKIP("MAS examples not available in this checkout");
    }

    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(examples)) {
        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());

    size_t lifted = 0;
    size_t belowFlexibility = 0;
    size_t examined = 0;
    std::cout << "\n=== WIRE BEND CENSUS (corner radius derived, no MAS datum yet) ===\n";
    for (const auto& file : files) {
        std::ifstream stream(file);
        if (!stream.good()) {
            continue;
        }
        json masJson;
        try {
            masJson = json::parse(stream);
        }
        catch (const std::exception&) {
            continue;
        }
        json magneticJson = masJson.contains("magnetic") ? masJson["magnetic"] : masJson;
        try {
            OpenMagnetics::Magnetic magnetic(magneticJson);
            auto coil = magnetic.get_mutable_coil();
            auto bobbin = coil.resolve_bobbin();
            if (!bobbin.get_processed_description()) {
                // The examples carry a bobbin by NAME ("Basic"/"Dummy"), which is what the rest
                // of MKF resolves into a synthesised one -- same call Utils makes when winding.
                bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(magnetic.get_mutable_core(), false);
            }
            if (!bobbin.get_processed_description()) {
                std::cout << file.stem().string() << "  SKIPPED (bobbin not processed)\n";
                continue;
            }
            const auto shape = bobbin.get_processed_description()->get_column_shape();
            if (shape != ColumnShape::RECTANGULAR && shape != ColumnShape::IRREGULAR) {
                // Nothing to lift off: a round or oblong column is all corner.
                std::cout << file.stem().string() << "  column is "
                          << (shape == ColumnShape::ROUND ? "ROUND" : "OBLONG")
                          << ", every turn conforms by construction\n";
                continue;
            }
            const double cornerRadius = bobbin.get_column_corner_radius();
            const double halfAngle = bobbin.get_column_corner_half_angle();

            for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size();
                 ++windingIndex) {
                auto wire = coil.resolve_wire(windingIndex);
                if (wire.get_type() != WireType::ROUND && wire.get_type() != WireType::RECTANGULAR) {
                    std::cout << file.stem().string() << "  winding " << windingIndex
                              << "  SKIPPED (no standardised bend data for this wire type)\n";
                    continue;
                }
                // A first-layer turn stands off by its own outer radius, which is how MKF
                // places it against the column today.
                const double standoff = 0.5 * (wire.get_type() == WireType::ROUND
                                                   ? OpenMagnetics::Wire::calculate_outer_diameter(wire)
                                                   : OpenMagnetics::Wire::calculate_outer_width(wire));
                const auto axis = wire.get_type() == WireType::ROUND
                                      ? BendAxis::ROUND
                                      : WireBend::axis_from_bend_plane_dimension(wire, true);
                const auto asBuilt = WireBend::evaluate(standoff, wire, axis);
                const auto solved = WireBend::solve(cornerRadius, halfAngle, wire, standoff, axis);
                ++examined;
                if (solved.regime == BendRegime::LIFTED) {
                    ++lifted;
                }
                if (asBuilt.verdict == BendVerdict::BELOW_FLEXIBILITY) {
                    ++belowFlexibility;
                }
                std::cout << file.stem().string() << "  winding " << windingIndex
                          << "  asBuilt=" << WireBend::to_string(asBuilt.verdict)
                          << " strain=" << asBuilt.outerFibreStrain
                          << "  formerCorner=" << cornerRadius * 1e3 << "mm"
                          << "  ->  " << WireBend::to_string(solved.regime)
                          << " bend=" << solved.bendRadius * 1e3 << "mm"
                          << " standoff=" << standoff * 1e3 << "->" << solved.faceStandoff * 1e3
                          << "mm  " << WireBend::to_string(solved.verdict) << "\n";
            }
        }
        catch (const std::exception& e) {
            std::cout << file.stem().string() << "  SKIPPED (" << e.what() << ")\n";
        }
    }
    std::cout << "=== " << examined << " windings examined, " << lifted << " lift off the former, "
              << belowFlexibility << " are drawn today tighter than the wire may be bent ===\n";
    REQUIRE(examined > 0);
}

TEST_CASE("Test_WireBend_Filter_Judges_The_Wire_Against_Its_Former", "[wire][wire-bend][filter]") {
    // The filter reads the MAS datum, so this also exercises the new columnCornerRadius field.
    const std::filesystem::path example =
        "/home/alf/OpenMagnetics/MKF/MAS/examples/02_flyback_efd25_3c95.json";
    std::ifstream stream(example);
    if (!stream.good()) {
        SKIP("MAS example not available in this checkout");
    }
    json masJson = json::parse(stream);
    json magneticJson = masJson.contains("magnetic") ? masJson["magnetic"] : masJson;

    auto withCornerRadius = [&](double cornerRadius) {
        OpenMagnetics::Magnetic magnetic(magneticJson);
        auto coil = magnetic.get_mutable_coil();
        auto bobbin = coil.resolve_bobbin();
        if (!bobbin.get_processed_description()) {
            bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(magnetic.get_mutable_core(), false);
        }
        auto processedDescription = bobbin.get_processed_description().value();
        processedDescription.set_column_corner_radius(cornerRadius);
        bobbin.set_processed_description(processedDescription);
        coil.set_bobbin(bobbin);
        magnetic.set_coil(coil);
        return magnetic;
    };

    OpenMagnetics::Inputs inputs;
    auto filter = OpenMagnetics::MagneticFilter::factory(MagneticFilters::WINDABILITY);

    // A sharp former -- what MKF has always modelled -- cannot be wound with this wire.
    auto sharp = withCornerRadius(0.0);
    auto sharpVerdict = filter->evaluate_magnetic(&sharp, &inputs);
    REQUIRE_FALSE(sharpVerdict.first);

    // The floor is set by the THICKEST wire on the former: this design also carries a 0,500 mm
    // conductor, whose Table 6 mandrel is 0,500 mm and so whose corner floor is 0,250 mm. Legal
    // there, but the insulation is not qualified to be heated afterwards, so it must still score
    // worse than...
    auto atFloor = withCornerRadius(0.250e-3);
    auto atFloorVerdict = filter->evaluate_magnetic(&atFloor, &inputs);
    REQUIRE(atFloorVerdict.first);

    // ...a former with the corner the heat-shock mandrel asks for.
    auto comfortable = withCornerRadius(0.700e-3);
    auto comfortableVerdict = filter->evaluate_magnetic(&comfortable, &inputs);
    REQUIRE(comfortableVerdict.first);
    REQUIRE(comfortableVerdict.second < atFloorVerdict.second);
}

TEST_CASE("Test_WireBend_Toroid_Edge_Comes_From_The_Coating", "[wire][wire-bend][filter]") {
    // A toroid has no moulded former: the turn is a closed loop around the ring cross-section and
    // bends at that section's edges. MAS says so directly -- cornerRadius is "the radius of the
    // core cross-section edges (or of its coating, when coated), which is what the wire is pulled
    // over" -- and create_quick_bobbin now resolves it instead of leaving the zero that the
    // rectangular-column moulding fallback produced from a toroid's zero wall thickness.
    //
    // The coating is the only quantity on that surface anyone dimensions. IEC 62317-12 gives a
    // ring core as A/B/C and admits the chamfer only through the effective height; the MMPA/IMA
    // toroid specification says only that the corners "shall not be sharp or rough"; Fair-Rite
    // burnishes to break the edge; TDK gives tumbling or a chamfer by size tier. None publishes a
    // number. The jacket does: Ferroxcube draws PA11 at ~0,3 mm, TDK specifies epoxy below 0,4 mm
    // and parylene at 12,7 or 25 um. A conformal jacket cannot form an edge sharper than itself,
    // so its thickness is a sourced lower bound on what the wire rides over.
    for (const std::string name : {"05_pfc_inductor_t4020_hf60", "07_cmc_t2515_w800",
                                   "12_boost_inductor_t5026_26"}) {
        const std::filesystem::path example =
            std::filesystem::path(std::string(__FILE__).substr(0, std::string(__FILE__).rfind('/')))
            / "../MAS/examples" / (name + ".json");
        std::ifstream stream(example);
        if (!stream.good()) {
            SKIP("MAS example not available in this checkout");
        }
        json masJson = json::parse(stream);
        json magneticJson = masJson.contains("magnetic") ? masJson["magnetic"] : masJson;
        OpenMagnetics::Magnetic magnetic(magneticJson);

        auto core = magnetic.get_mutable_core();
        REQUIRE(core.get_shape_family() == CoreShapeFamily::T);
        // create_quick_bobbin processes its own copy, so process this one to query it directly.
        core.process_data();
        auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, false);
        const double edgeRadius = bobbin.get_column_corner_radius();

        // The edge is the jacket's, not the zero a wall-less part used to resolve.
        REQUIRE(edgeRadius > 0);
        REQUIRE_THAT(edgeRadius, Catch::Matchers::WithinRel(core.get_toroid_edge_radius(), 1e-12));

        // It is its OWN datum, not the coating thickness: that one is the dielectric path normal
        // to the flat faces (what StrayCapacitance integrates over), this one is the curvature at
        // the corner. They are resolved separately so moving one cannot silently move the other.
        REQUIRE(edgeRadius >= core.get_coating_thickness());

        // And it never exceeds what the ring section can carry: a corner radius past half the
        // smaller dimension of the section would have eaten the face between the two edges.
        // Named local: get_processed_description() returns the optional BY VALUE, so binding a
        // reference through it would point into a temporary that is already gone.
        const auto coreProcessed = core.get_processed_description().value();
        const auto& ringSection = coreProcessed.get_columns()[0];
        REQUIRE(edgeRadius <=
                0.5 * std::min(ringSection.get_width(), ringSection.get_depth()) + 1e-12);

        // And it is a real bend the filter can judge. IEC 60851-3 winds the wire ON the mandrel,
        // so the wire's own outer radius appears on both sides of the comparison and cancels:
        // the edge is legal for a given wire exactly when it is at least the mandrel radius.
        for (const std::string wireName : {"Round 0.2 - Grade 1", "Round 0.4 - Grade 1"}) {
            auto wire = find_wire_by_name(wireName);
            const double standoff = 0.5 * OpenMagnetics::Wire::calculate_outer_diameter(wire);
            const double mandrelRadius =
                0.5 * WireBend::get_mandrel_diameter(wire, BendCriterion::FLEXIBILITY,
                                                     BendAxis::ROUND);
            const auto asLaid = WireBend::evaluate(edgeRadius + standoff, wire, BendAxis::ROUND);
            const bool legal = asLaid.verdict != BendVerdict::BELOW_FLEXIBILITY;
            REQUIRE(legal == (edgeRadius >= mandrelRadius - 1e-9));
            std::cout << name << "  " << wireName << "  edge=" << edgeRadius * 1e3
                      << " mm  mandrelRadius=" << mandrelRadius * 1e3 << " mm  -> "
                      << (legal ? "windable" : "BELOW_FLEXIBILITY") << "\n";
        }
    }
}

TEST_CASE("Test_WireBend_Windability_Filter_Survives_A_Json_Round_Trip", "[wire][wire-bend][filter]") {
    // MagneticFilters has two hand-maintained serializers, and a new enumerator is only reachable
    // over JSON once it is registered in BOTH. Every adviser pipeline crosses that boundary -- a
    // MagneticFilterOperation is how a flow is configured and how it is echoed back -- so an
    // unregistered value throws on the way in ("does not conform to MagneticFilters schema") and
    // again on the way out ("Unexpected value in enumeration"). The C++ factory never notices,
    // which is exactly why this needs its own test.
    const json asJson = "Windability";
    REQUIRE(asJson.get<MagneticFilters>() == MagneticFilters::WINDABILITY);

    json roundTripped = MagneticFilters::WINDABILITY;
    REQUIRE(roundTripped.get<std::string>() == "Windability");

    // And through the operation that actually carries it into an adviser.
    OpenMagnetics::MagneticFilterOperation operation;
    operation.set_filter(MagneticFilters::WINDABILITY);
    operation.set_invert(false);
    operation.set_log(false);
    operation.set_weight(1.0);
    json operationJson = operation;
    REQUIRE(operationJson["filter"].get<std::string>() == "Windability");
    REQUIRE(operationJson.get<OpenMagnetics::MagneticFilterOperation>().get_filter() ==
            MagneticFilters::WINDABILITY);
}


TEST_CASE("Test_WireBend_Toroid_Edge_Follows_The_Size_Tier", "[wire][wire-bend][filter]") {
    // Ring cores are not all finished the same way, and TDK's ring-core general information says
    // so by size: small cores get "edges rounded by tumbling", medium and large ones get a
    // chamfer. A tumbled edge is far tighter than a drawn chamfer, so the two tiers carry
    // different numbers and the tier is read off the same size boundary the coating crosses at.
    auto toroid = [](const std::string& shape) {
        json coreJson;
        coreJson["functionalDescription"]["type"] = "toroidal";
        coreJson["functionalDescription"]["material"] = "N97";
        coreJson["functionalDescription"]["shape"] = shape;
        coreJson["functionalDescription"]["gapping"] = json::array();
        coreJson["functionalDescription"]["numberStacks"] = 1;
        return OpenMagnetics::Core(coreJson, true);
    };
    auto defaults = Defaults();

    // Chamfered tier: the drawn ring-core jacket.
    auto large = toroid("T 25/15/10");
    REQUIRE_FALSE(large.get_default_toroid_coating_is_parylene());
    REQUIRE_THAT(large.get_toroid_edge_radius(),
                 Catch::Matchers::WithinAbs(defaults.defaultToroidRingEdgeRadius, 1e-12));

    // Tumbled tier: an order of magnitude tighter, and the reason a small toroid is judged much
    // more harshly than a large one. See Defaults for how weakly sourced that number is.
    auto small = toroid("T 2.5/1.5/1");
    REQUIRE(small.get_default_toroid_coating_is_parylene());
    REQUIRE_THAT(small.get_toroid_edge_radius(),
                 Catch::Matchers::WithinAbs(defaults.defaultTumbledToroidEdgeRadius, 1e-12));
    REQUIRE(small.get_toroid_edge_radius() < large.get_toroid_edge_radius());

    // Whatever the tier says, the edge can never exceed what the ring section can carry: a corner
    // radius past half the smaller dimension of a section has eaten the face between its two
    // edges. No catalogue toroid is thin enough for that clamp to bite today (nothing has a wall
    // or height under 0,6 mm), so this asserts it as a property across the whole family rather
    // than picking a part that happens to trip it.
    size_t examined = 0;
    for (const auto& shape : get_shapes(true)) {
        if (shape.get_family() != CoreShapeFamily::T || !shape.get_name()) {
            continue;
        }
        OpenMagnetics::Core core = [&] {
            json coreJson;
            coreJson["functionalDescription"]["type"] = "toroidal";
            coreJson["functionalDescription"]["material"] = "N97";
            coreJson["functionalDescription"]["shape"] = shape.get_name().value();
            coreJson["functionalDescription"]["gapping"] = json::array();
            coreJson["functionalDescription"]["numberStacks"] = 1;
            return OpenMagnetics::Core(coreJson, true);
        }();
        // Named local: get_processed_description() returns the optional BY VALUE, so binding a
        // reference through it would point into a temporary that is already gone.
        const auto coreProcessed = core.get_processed_description().value();
        const auto& ringSection = coreProcessed.get_columns()[0];
        const double edge = core.get_toroid_edge_radius();
        REQUIRE(edge > 0);
        REQUIRE(edge <= 0.5 * std::min(ringSection.get_width(), ringSection.get_depth()) + 1e-12);
        ++examined;
    }
    REQUIRE(examined > 0);
    std::cout << "toroid edge radius checked against the ring section on " << examined
              << " catalogue T shapes\n";

    // A bobbin-wound core has no ring edge to ask about, and says so rather than inventing one.
    json eCoreJson;
    eCoreJson["functionalDescription"]["type"] = "two-piece set";
    eCoreJson["functionalDescription"]["material"] = "N97";
    eCoreJson["functionalDescription"]["shape"] = "E 55/28/21";
    eCoreJson["functionalDescription"]["gapping"] = json::array();
    eCoreJson["functionalDescription"]["numberStacks"] = 1;
    OpenMagnetics::Core eCore(eCoreJson, true);
    REQUIRE_THROWS(eCore.get_toroid_edge_radius());
}
