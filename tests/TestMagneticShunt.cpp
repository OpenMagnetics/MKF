// Magnetic shunts (MAS-RFC 0015, ABT #1176).
//
// Literature cases, encoded from the papers' own tables; every value not published is named below
// as an assumption and is never presented as the paper's:
//
//   [Li]    M. Li, Z. Ouyang, M. A. E. Andersen, "High Frequency LLC Resonant Converter with Magnetic
//           Shunt Integrated Planar Transformer", IEEE Trans. Power Electron. 34(3), 2019
//           (https://backend.orbit.dtu.dk/ws/files/149117514/08369125.pdf), Tables IV and V:
//           E32/6/20 3F46, kp 2, Np 4, hp 0.07 mm, h_dp 0.25 mm, ks 1, Ns 4, hs 0.07 mm, h_ds 0.25 mm,
//           shunt IFL04 mus 45, tsh 0.1 mm, air gap 2la 0.18 mm, xp = xs = 2 mm; Lm 26 uH, Lr 3.5 uH
//           (resonant inductance integrated as the transformer leakage), measured at 1 MHz.
//   [Zhang] J. Zhang, Z. Ouyang, M. C. Duffy, M. A. E. Andersen, W. G. Hurley, "Leakage Inductance
//           Calculation for Planar Transformers with a Magnetic Shunt", IEEE Trans. Ind. Appl. 50(6),
//           2014 (https://backend.orbit.dtu.dk/ws/files/97832800/06810816.pdf), Table II and Fig. 12:
//           ELP 43/10/28, kp = ks = 1, hp = hs = 0.15 mm; PCB1 4 layers h_d 0.4 mm, PCB2 2 layers
//           h_d 1.5 mm; Trans. 1 PCB1/PCB1, Trans. 2 PCB2/PCB1, Trans. 3 PCB2/PCB2; shunt mus 30,
//           0.5 mm and 1 mm; measured at 100 kHz: 1402/2508, 360/636, 335/663 nH. Fig. 9 (FEA and
//           calculation, Table I = Trans. 1, N87, mus 10..200, h 0.1..2 mm) spans 0..35 uH.
//
// Assumptions (not published):
//   - PCB track width: the turns of a layer fill the window width at MKF's minimum border and
//     wire-to-wire distances (Defaults). Neither paper gives a track width.
//   - Both papers insert the sheet between the two core halves ("inserted in between the two magnetic
//     cores", [Li] Sec. III), so the halves are separated by the sheet plus the air: the core gap per
//     column is 2la + tsh = 0.28 mm for [Li], and tsh alone for [Zhang], which reports no air gap.
//   - [Zhang] gives no distance between the boards and the sheet; the outer copper layers are taken as
//     touching it. [Li] gives xp and xs, taken as copper-to-sheet distances.
//   - Sheet width and depth = the core's nominal width and depth (the sheet spans the core).
//   - Material names in the fixtures ("IFL04" for [Li]) describe the papers' parts; the permeability
//     used is the paper's mus, passed explicitly, because neither sheet grade is in the database
//     (owner decision 2026-09-13: validate materials on TDK FPC C350/C351 only).

#include "physical_models/LeakageInductance.h"
#include "physical_models/MagneticShunt.h"
#include "physical_models/MagnetizingInductance.h"
#include "advisers/MagneticFilter.h"
#include "support/Painter.h"
#include "support/Utils.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <source_location>

using namespace MAS;
using namespace OpenMagnetics;
using Catch::Matchers::WithinRel;

namespace {

// CMC Lk/Lm score of the LEAKAGE_INDUCTANCE filter on the T 36/23/15 3C90 12:12 fixture below, measured
// with the filter files as they were before ABT #1176 (base 137a2485, MagneticFilter.h/.cpp and
// MagneticFilterAdvanced.cpp restored, same binary otherwise): 0.0010388323822073957.
const double kCmcLeakageRatioBeforeAbt1176 = 0.0010388323822073957;

// Relative error against a published measurement, referred to the measurement (Catch's WithinRel
// refers to the larger of the two values, which is looser when the model overshoots).
double relative_error(double value, double reference) {
    return std::fabs(value - reference) / reference;
}

auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");

OpenMagnetics::Wire make_planar_track(double width, double thickness) {
    OpenMagnetics::Wire wire;
    wire.set_nominal_value_conducting_width(width);
    wire.set_nominal_value_conducting_height(thickness);
    wire.set_nominal_value_outer_width(width);
    wire.set_nominal_value_outer_height(thickness);
    wire.set_number_conductors(1);
    wire.set_material("copper");
    wire.set_type(WireType::PLANAR);
    return wire;
}

struct PlanarShuntCase {
    std::string shapeName;
    std::string coreMaterial;
    double coreGap;             // separation of the two halves, per column
    int64_t primaryTurnsPerLayer;
    size_t primaryLayers;
    double primaryInsulation;   // between primary layers
    int64_t secondaryTurnsPerLayer;
    size_t secondaryLayers;
    double secondaryInsulation; // between secondary layers
    double copperThickness;
    double interfaceDistance;   // last primary copper to first secondary copper, sheet included
    double sheetThickness;
    std::string sheetMaterialName;
};

// Planar transformer: primary layers above, secondary below, the sheet at the core mating plane.
OpenMagnetics::Magnetic make_planar_shunt_transformer(const PlanarShuntCase& spec, bool withShunt = true) {
    auto core = OpenMagneticsTesting::get_quick_core(spec.shapeName, OpenMagneticsTesting::get_spacer_gap(spec.coreGap), 1, spec.coreMaterial);
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);
    auto windowDimensions = bobbin.get_winding_window_dimensions();
    double windowWidth = windowDimensions[0];
    double windowHeight = windowDimensions[1];

    double border = Defaults().minimumBorderToWireDistance;
    double wireToWire = Defaults().minimumWireToWireDistance;
    auto trackWidth = [&](int64_t turnsPerLayer) {
        return (windowWidth - 2 * border - static_cast<double>(turnsPerLayer - 1) * wireToWire) / static_cast<double>(turnsPerLayer);
    };

    OpenMagnetics::Coil coil;
    OpenMagnetics::Winding primary;
    primary.set_name("Primary");
    primary.set_number_turns(spec.primaryTurnsPerLayer * static_cast<int64_t>(spec.primaryLayers));
    primary.set_number_parallels(1);
    primary.set_isolation_side(IsolationSide::PRIMARY);
    primary.set_wire(make_planar_track(trackWidth(spec.primaryTurnsPerLayer), spec.copperThickness));
    OpenMagnetics::Winding secondary;
    secondary.set_name("Secondary");
    secondary.set_number_turns(spec.secondaryTurnsPerLayer * static_cast<int64_t>(spec.secondaryLayers));
    secondary.set_number_parallels(1);
    secondary.set_isolation_side(IsolationSide::SECONDARY);
    secondary.set_wire(make_planar_track(trackWidth(spec.secondaryTurnsPerLayer), spec.copperThickness));
    coil.get_mutable_functional_description().push_back(primary);
    coil.get_mutable_functional_description().push_back(secondary);
    coil.set_bobbin(bobbin);

    std::vector<size_t> stackUp;
    for (size_t layer = 0; layer < spec.primaryLayers; ++layer) {
        stackUp.push_back(0);
    }
    for (size_t layer = 0; layer < spec.secondaryLayers; ++layer) {
        stackUp.push_back(1);
    }
    std::map<std::pair<size_t, size_t>, double> insulation;
    insulation[{0, 0}] = spec.primaryInsulation;
    insulation[{1, 1}] = spec.secondaryInsulation;
    insulation[{0, 1}] = spec.interfaceDistance;
    double stackHeight = static_cast<double>(spec.primaryLayers + spec.secondaryLayers) * spec.copperThickness +
                         static_cast<double>(spec.primaryLayers - 1) * spec.primaryInsulation +
                         static_cast<double>(spec.secondaryLayers - 1) * spec.secondaryInsulation + spec.interfaceDistance;
    double coreToLayerDistance = (windowHeight - stackHeight) / 2;
    REQUIRE(coreToLayerDistance > 0);
    settings.set_coil_wind_even_if_not_fit(true);
    coil.set_section_alignment(CoilAlignment::CENTERED);
    coil.wind_planar(stackUp, border, {{0, wireToWire}, {1, wireToWire}}, insulation, coreToLayerDistance);
    REQUIRE(coil.get_turns_description());

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    if (withShunt) {
        MagneticShunt shunt;
        shunt.set_name("Shunt");
        shunt.set_placement(MagneticShuntPlacement::BETWEEN_SECTIONS);
        shunt.set_coordinates({0, 0, 0});
        shunt.set_dimensions({core.get_width(), spec.sheetThickness, core.get_depth()});
        shunt.set_material(spec.sheetMaterialName);
        magnetic.set_shunts(std::vector<MagneticShunt>{shunt});
    }
    return magnetic;
}

// [Li] Table V.
PlanarShuntCase li_2018_case() {
    PlanarShuntCase spec;
    spec.shapeName = "E 32/6/20";
    spec.coreMaterial = "3F46";
    spec.sheetThickness = 0.1e-3;
    spec.coreGap = 0.18e-3 + spec.sheetThickness;
    spec.primaryTurnsPerLayer = 2;
    spec.primaryLayers = 4;
    spec.primaryInsulation = 0.25e-3;
    spec.secondaryTurnsPerLayer = 1;
    spec.secondaryLayers = 4;
    spec.secondaryInsulation = 0.25e-3;
    spec.copperThickness = 0.07e-3;
    spec.interfaceDistance = 2e-3 + spec.sheetThickness + 2e-3;
    spec.sheetMaterialName = "IFL04";
    return spec;
}

// [Zhang] Table II; boards: PCB1 4 layers h_d 0.4 mm, PCB2 2 layers h_d 1.5 mm.
PlanarShuntCase zhang_2014_case(bool primaryIsPcb2, bool secondaryIsPcb2, double sheetThickness) {
    PlanarShuntCase spec;
    spec.shapeName = "E 43/10/28";
    spec.coreMaterial = "N87";
    spec.sheetThickness = sheetThickness;
    spec.coreGap = sheetThickness;
    spec.primaryTurnsPerLayer = 1;
    spec.primaryLayers = primaryIsPcb2 ? 2 : 4;
    spec.primaryInsulation = primaryIsPcb2 ? 1.5e-3 : 0.4e-3;
    spec.secondaryTurnsPerLayer = 1;
    spec.secondaryLayers = secondaryIsPcb2 ? 2 : 4;
    spec.secondaryInsulation = secondaryIsPcb2 ? 1.5e-3 : 0.4e-3;
    spec.copperThickness = 0.15e-3;
    spec.interfaceDistance = sheetThickness;
    spec.sheetMaterialName = "Zhang 2014 shunt (mus 30)";
    return spec;
}

double copper_pitch_of_winding(OpenMagnetics::Magnetic& magnetic, const std::string& windingName) {
    std::vector<double> heights;
    auto turns = magnetic.get_coil().get_turns_description().value();
    for (auto& turn : turns) {
        if (turn.get_winding() == windingName && std::find_if(heights.begin(), heights.end(), [&](double y) { return std::fabs(y - turn.get_coordinates()[1]) < 1e-9; }) == heights.end()) {
            heights.push_back(turn.get_coordinates()[1]);
        }
    }
    std::sort(heights.begin(), heights.end());
    REQUIRE(heights.size() >= 2);
    return heights[1] - heights[0];
}

} // namespace

TEST_CASE("Magnetic shunt: Li 2018 E32/6/20 integrated leakage inductance within 15 % of 3.5 uH", "[magnetic-shunt][leakage-inductance][physical-model]") {
    settings.reset();
    auto spec = li_2018_case();
    auto magnetic = make_planar_shunt_transformer(spec);

    // The fixture is the paper's stack: 0.07 mm copper on 0.25 mm insulation.
    CHECK_THAT(copper_pitch_of_winding(magnetic, "Primary"), WithinRel(0.32e-3, 1e-6));
    CHECK_THAT(copper_pitch_of_winding(magnetic, "Secondary"), WithinRel(0.32e-3, 1e-6));

    double frequency = 1e6;
    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, frequency, 0, 1, 1, std::vector<double>{45.0});
    INFO("winding (Energy) Lk = " << result.windingLeakageInductance << " H");
    INFO("shunt network Lk = " << result.shunts[0].networkLeakageInductance << " H, displaced air " << result.shunts[0].displacedAirLeakageInductance << " H");
    INFO("total Lk = " << result.leakageInductance << " H");
    REQUIRE(result.shunts.size() == 1);
    REQUIRE(result.shunts[0].crossings.size() == 2);
    CHECK(result.shunts[0].enclosedMagnetomotiveForcePerAmpere == 8.0);
    CHECK(relative_error(result.leakageInductance, 3.5e-6) <= 0.15);

    // B_shunt is reported, and is the network's own flux over the sheet section.
    auto& crossing = result.shunts[0].crossings[0];
    double expectedFluxDensityPerAmpere = 8.0 / crossing.totalReluctance / (spec.sheetThickness * crossing.depth);
    CHECK(result.shunts[0].magneticFluxDensityPerAmpere > 0);
    CHECK_THAT(result.shunts[0].magneticFluxDensityPerAmpere, WithinRel(expectedFluxDensityPerAmpere, 1e-9));
    // No material record was read: no losses are claimed.
    CHECK(!result.shunts[0].lossesPerAmpereSquared);
}

TEST_CASE("Magnetic shunt: the sheet network reproduces Li 2018 eq. (17) within 5 %", "[magnetic-shunt][leakage-inductance][physical-model]") {
    settings.reset();
    auto spec = li_2018_case();
    auto magnetic = make_planar_shunt_transformer(spec);
    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, 1e6, 0, 1, 1, std::vector<double>{45.0});

    // Li eqs. (13), (16), (17) with the E 32/6/20 nominal dimensions of the MAS shape record
    // (A 31.75, C 20.325, E 25.5, F 6.35 mm): b_d = (A - E)/2, 2 b_d = F for the centre leg as Li
    // draws it, b_w = (E - F)/2, l_w = C, A_c = F C, l_a = 0.09 mm, t_sh = 0.1 mm, mus 45, n_p = 8.
    // The core sections R_c1, R_c2, R_cc (mur of order 10^3, about 1e5 A/Wb against 1.7e8 A/Wb in the
    // loop) are left out of this reference; they move it by about 0.1 %.
    double mu0 = Constants().vacuumPermeability;
    double A = 31.75e-3, C = 20.325e-3, E = 25.5e-3, F = 6.35e-3;
    double bd = (A - E) / 2, bw = (E - F) / 2, lw = C, Ac = F * C, la = 0.09e-3, tsh = 0.1e-3, mus = 45, np = 8;
    double Rs1 = tsh / (2 * mu0 * mus * bd * lw);
    double Rss = tsh / (2 * mu0 * mus * Ac);
    double Rs2 = bw / (mu0 * mus * tsh * lw);
    double Rg1 = la / (mu0 * (bd + la) * (lw + la));
    double Rg2 = la / (mu0 * (2 * bd + la) * (lw + la));
    double Rm = Rs1 + 2 * Rss;
    double publishedShuntLeakage = 4 * np * np / (Rm + Rg1 + 2 * Rg2 + 2 * Rs2);
    INFO("MKF network " << result.shunts[0].networkLeakageInductance << " H, Li eq. (17) " << publishedShuntLeakage << " H");
    CHECK(relative_error(result.shunts[0].networkLeakageInductance, publishedShuntLeakage) <= 0.05);
}

TEST_CASE("Magnetic shunt: Li 2018 magnetizing inductance within 10 % of 26 uH with the sheet in the gap", "[magnetic-shunt][magnetizing-inductance][physical-model]") {
    settings.reset();
    auto spec = li_2018_case();
    auto magnetic = make_planar_shunt_transformer(spec);
    double frequency = 1e6;
    auto core = MagneticShuntModel::apply_shunts_to_gapping(magnetic, Defaults().ambientTemperature, frequency, std::vector<double>{45.0});
    for (auto& gap : core.get_functional_description().get_gapping()) {
        // 0.28 mm of separation of which 0.1 mm is sheet at mus 45.
        CHECK_THAT(gap.get_length(), WithinRel(0.18e-3 + 0.1e-3 / 45.0, 1e-9));
    }
    double magnetizingInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(core, magnetic.get_coil()).get_magnetizing_inductance().get_nominal().value();
    INFO("Lm = " << magnetizingInductance << " H");
    CHECK(relative_error(magnetizingInductance, 26e-6) <= 0.10);
}

TEST_CASE("Magnetic shunt: Zhang 2014 ELP 43/10/28 prototypes within 15 % of the measured leakage", "[magnetic-shunt][leakage-inductance][physical-model]") {
    settings.reset();
    struct Prototype {
        std::string name;
        bool primaryIsPcb2;
        bool secondaryIsPcb2;
        double thickness;
        double measured;
    };
    // Fig. 12 (a) 0.5 mm and (b) 1 mm sheets, experimental bars.
    std::vector<Prototype> prototypes = {
        {"Trans. 1, 0.5 mm", false, false, 0.5e-3, 1402e-9},
        {"Trans. 2, 0.5 mm", true, false, 0.5e-3, 360e-9},
        {"Trans. 3, 0.5 mm", true, true, 0.5e-3, 335e-9},
        {"Trans. 1, 1 mm", false, false, 1e-3, 2508e-9},
        {"Trans. 2, 1 mm", true, false, 1e-3, 636e-9},
        {"Trans. 3, 1 mm", true, true, 1e-3, 663e-9},
    };
    for (auto& prototype : prototypes) {
        auto magnetic = make_planar_shunt_transformer(zhang_2014_case(prototype.primaryIsPcb2, prototype.secondaryIsPcb2, prototype.thickness));
        auto result = LeakageInductance().calculate_shunt_leakage(magnetic, 100e3, 0, 1, 1, std::vector<double>{30.0});
        INFO(prototype.name << ": winding " << result.windingLeakageInductance << " H, total " << result.leakageInductance << " H, measured " << prototype.measured << " H");
        CHECK(relative_error(result.leakageInductance, prototype.measured) <= 0.15);
    }
}

TEST_CASE("Magnetic shunt: Zhang 2014 sweep grows with mus and thickness and stays inside the Fig. 9 range", "[magnetic-shunt][leakage-inductance][physical-model]") {
    settings.reset();
    std::vector<double> permeabilities = {10, 40, 100, 150, 200};
    std::vector<double> thicknesses = {0.1e-3, 0.2e-3, 0.4e-3, 1e-3, 2e-3};
    std::map<std::pair<size_t, size_t>, double> leakage;
    for (size_t thicknessIndex = 0; thicknessIndex < thicknesses.size(); ++thicknessIndex) {
        auto magnetic = make_planar_shunt_transformer(zhang_2014_case(false, false, thicknesses[thicknessIndex]));
        for (size_t permeabilityIndex = 0; permeabilityIndex < permeabilities.size(); ++permeabilityIndex) {
            auto result = LeakageInductance().calculate_shunt_leakage(magnetic, 100e3, 0, 1, 1, std::vector<double>{permeabilities[permeabilityIndex]});
            leakage[{thicknessIndex, permeabilityIndex}] = result.leakageInductance;
            INFO("mus " << permeabilities[permeabilityIndex] << " h " << thicknesses[thicknessIndex] << " -> " << result.leakageInductance);
            CHECK(result.leakageInductance > 0);
            CHECK(result.leakageInductance <= 35e-6);
        }
    }
    for (size_t thicknessIndex = 0; thicknessIndex < thicknesses.size(); ++thicknessIndex) {
        for (size_t permeabilityIndex = 0; permeabilityIndex < permeabilities.size(); ++permeabilityIndex) {
            if (permeabilityIndex > 0) {
                CHECK(leakage[{thicknessIndex, permeabilityIndex}] > leakage[{thicknessIndex, permeabilityIndex - 1}]);
            }
            if (thicknessIndex > 0) {
                CHECK(leakage[{thicknessIndex, permeabilityIndex}] > leakage[{thicknessIndex - 1, permeabilityIndex}]);
            }
        }
    }
    // Ordered by mus * t across the whole sweep (products differing by at least 10 %).
    for (auto& [keyA, valueA] : leakage) {
        for (auto& [keyB, valueB] : leakage) {
            double productA = permeabilities[keyA.second] * thicknesses[keyA.first];
            double productB = permeabilities[keyB.second] * thicknesses[keyB.first];
            if (productA > 1.1 * productB) {
                CHECK(valueA > valueB);
            }
        }
    }
}

TEST_CASE("Magnetic shunt: TDK C350 by name selects the Shunt method and reports B and losses inside the data span", "[magnetic-shunt][leakage-inductance][physical-model]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);

    auto output = LeakageInductance().calculate_leakage_inductance(magnetic, 1e6, 0, 1);
    CHECK(output.get_method_used() == "Shunt");

    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, 1e6, 0, 1);
    CHECK_THAT(output.get_leakage_inductance_per_winding()[0].get_nominal().value(), WithinRel(result.leakageInductance, 1e-12));
    // C350 initial permeability 9 (TDK FPC datasheet, in core_materials.ndjson).
    CHECK_THAT(result.shunts[0].relativePermeability, WithinRel(9.0, 1e-9));
    CHECK(result.shunts[0].magneticFluxDensityPerAmpere > 0);
    REQUIRE(result.shunts[0].lossesPerAmpereSquared);
    CHECK(result.shunts[0].lossesPerAmpereSquared.value() > 0);

    // Below the measured complex-permeability span (1 MHz up) no losses are extrapolated.
    auto lowFrequency = LeakageInductance().calculate_shunt_leakage(magnetic, 100e3, 0, 1);
    CHECK(!lowFrequency.shunts[0].lossesPerAmpereSquared);

    // A thicker sheet of the same material adds leakage; removing it returns the Energy value.
    auto withoutShunt = magnetic;
    withoutShunt.set_shunts(std::nullopt);
    auto energy = LeakageInductance().calculate_leakage_inductance(withoutShunt, 1e6, 0, 1);
    CHECK(energy.get_method_used() == "Energy");
    CHECK_THAT(energy.get_leakage_inductance_per_winding()[0].get_nominal().value(), WithinRel(result.windingLeakageInductance, 1e-12));
    CHECK(result.leakageInductance > result.windingLeakageInductance);
}

TEST_CASE("Magnetic shunt: a sheet material without initial permeability throws", "[magnetic-shunt][leakage-inductance]") {
    settings.reset();
    auto spec = li_2018_case();
    auto magnetic = make_planar_shunt_transformer(spec);
    // Deliberately incomplete in-memory record to exercise the guard (never persisted).
    auto material = find_core_material_by_name("C350");
    auto permeability = material.get_permeability();
    permeability.set_initial(std::vector<PermeabilityPoint>{});
    material.set_permeability(permeability);
    auto shunts = magnetic.get_shunts().value();
    shunts[0].set_material(material);
    magnetic.set_shunts(shunts);
    REQUIRE_THROWS_AS(LeakageInductance().calculate_leakage_inductance(magnetic, 1e6, 0, 1), MaterialDataMissingException);
}

TEST_CASE("Magnetic shunt: a sheet material name not in the database throws", "[magnetic-shunt][leakage-inductance]") {
    settings.reset();
    auto magnetic = make_planar_shunt_transformer(li_2018_case());  // "IFL04" is not filed
    REQUIRE_THROWS(LeakageInductance().calculate_leakage_inductance(magnetic, 1e6, 0, 1));
}

TEST_CASE("Magnetic shunt: geometry that contradicts the declared gaps throws", "[magnetic-shunt][leakage-inductance]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);
    auto core = magnetic.get_core();
    auto shunts = magnetic.get_shunts().value();

    SECTION("a sheet running into a column where the column has no gap") {
        shunts[0].set_coordinates({0, 1e-3, 0});
        magnetic.set_shunts(shunts);
        REQUIRE_THROWS_AS(LeakageInductance().calculate_shunt_leakage(magnetic, 1e6, 0, 1), InvalidInputException);
    }
    SECTION("a sheet short of the columns without gapToColumns") {
        double innerFace = core.get_columns()[core.get_main_column_index()].get_width() / 2;
        shunts[0].set_coordinates({innerFace + 2e-3, 0, 0});
        shunts[0].set_dimensions({2e-3, spec.sheetThickness, core.get_depth()});
        magnetic.set_shunts(shunts);
        REQUIRE_THROWS_AS(LeakageInductance().calculate_shunt_leakage(magnetic, 1e6, 0, 1), InvalidInputException);
    }
    SECTION("an outsideWindow sheet is not modelled") {
        shunts[0].set_placement(MagneticShuntPlacement::OUTSIDE_WINDOW);
        magnetic.set_shunts(shunts);
        REQUIRE_THROWS_AS(LeakageInductance().calculate_leakage_inductance(magnetic, 1e6, 0, 1), NotImplementedException);
    }
}

TEST_CASE("Magnetic shunt: sheets with gaps to the columns add the gap reluctance in series", "[magnetic-shunt][leakage-inductance]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);
    auto core = magnetic.get_core();
    auto columns = core.get_columns();
    auto main = columns[core.get_main_column_index()];
    double innerFace = main.get_width() / 2;
    double outerFace = 0;
    for (auto& column : columns) {
        if (column.get_coordinates()[0] > 0) {
            outerFace = column.get_coordinates()[0] - column.get_width() / 2;
        }
    }
    double gap = 0.1e-3;
    std::vector<MagneticShunt> sheets;
    for (int side : {1, -1}) {
        MagneticShunt sheet;
        sheet.set_placement(MagneticShuntPlacement::IN_WINDOW);
        sheet.set_coordinates({side * (innerFace + outerFace) / 2, 0, 0});
        sheet.set_dimensions({outerFace - innerFace - 2 * gap, spec.sheetThickness, core.get_depth()});
        MagneticShuntGapToColumns gaps;
        gaps.set_inner(gap);
        gaps.set_outer(gap);
        sheet.set_gap_to_columns(gaps);
        sheet.set_material("C350");
        sheets.push_back(sheet);
    }
    magnetic.set_shunts(sheets);
    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, 1e6, 0, 1);
    REQUIRE(result.shunts.size() == 2);
    for (auto& contribution : result.shunts) {
        REQUIRE(contribution.crossings.size() == 1);
        auto& crossing = contribution.crossings[0];
        double mu0 = Constants().vacuumPermeability;
        double sheetReluctance = (outerFace - innerFace - 2 * gap) / (mu0 * 9.0 * spec.sheetThickness * core.get_depth());
        double gapReluctance = gap / (mu0 * (spec.sheetThickness + gap) * (core.get_depth() + gap));
        CHECK_THAT(crossing.shuntBranchReluctance, WithinRel(sheetReluctance + 2 * gapReluctance, 1e-9));
    }
    // Two separate sheets with gaps carry less flux than one continuous sheet through the gap.
    auto continuous = make_planar_shunt_transformer(spec);
    auto continuousResult = LeakageInductance().calculate_shunt_leakage(continuous, 1e6, 0, 1);
    CHECK(result.leakageInductance < continuousResult.leakageInductance);
}

TEST_CASE("Magnetic shunt: shunts survive magnetic autocomplete and the JSON round trip", "[magnetic-shunt][masautocomplete][support]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);
    auto expected = magnetic.get_shunts().value();

    auto check = [&](const OpenMagnetics::Magnetic& candidate) {
        REQUIRE(candidate.get_shunts());
        REQUIRE(candidate.get_shunts()->size() == expected.size());
        auto shunt = candidate.get_shunts().value()[0];
        CHECK(shunt.get_name() == expected[0].get_name());
        CHECK(shunt.get_placement() == expected[0].get_placement());
        CHECK(shunt.get_dimensions() == expected[0].get_dimensions());
        CHECK(shunt.get_coordinates() == expected[0].get_coordinates());
        CHECK(std::get<std::string>(shunt.get_material()) == "C350");
    };

    auto autocompleted = magnetic_autocomplete(magnetic);
    check(autocompleted);

    json magneticJson;
    OpenMagnetics::to_json(magneticJson, magnetic);
    REQUIRE(magneticJson.contains("shunts"));
    OpenMagnetics::Magnetic fromJson;
    OpenMagnetics::from_json(magneticJson, fromJson);
    check(fromJson);

    MAS::Magnetic masMagnetic(magneticJson);
    check(OpenMagnetics::Magnetic(masMagnetic));

    OpenMagnetics::Mas mas;
    mas.set_magnetic(magnetic);
    mas.set_inputs(OpenMagnetics::Inputs::create_quick_operating_point(1e6, 26e-6, 25, WaveformLabel::SINUSOIDAL, 10, 0.5, 0, {2.0}));
    auto autocompletedMas = mas_autocomplete(mas, false);
    check(autocompletedMas.get_magnetic());
}

TEST_CASE("Magnetic shunt: LEAKAGE_INDUCTANCE filter keeps the CMC Lk/Lm score", "[magnetic-shunt][adviser][magnetic-filter]") {
    settings.reset();
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("T 36/23/15", json::array(), {12, 12}, 1, "3C90");
    magnetic = magnetic_autocomplete(magnetic);
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(150e3, 1e-3, 25, WaveformLabel::SINUSOIDAL, 1, 0.5, 0, {1.0});

    auto factoryFilter = MagneticFilter::factory(MagneticFilters::LEAKAGE_INDUCTANCE, inputs);
    MagneticFilterLeakageInductance explicitFilter(LeakageInductanceFilterMode::MINIMIZE_LEAKAGE_RATIO);
    MagneticFilterLeakageInductance defaultFilter;
    auto [factoryValid, factoryScore] = factoryFilter->evaluate_magnetic(&magnetic, &inputs);
    auto [explicitValid, explicitScore] = explicitFilter.evaluate_magnetic(&magnetic, &inputs);
    auto [defaultValid, defaultScore] = defaultFilter.evaluate_magnetic(&magnetic, &inputs);
    CHECK(defaultFilter.get_mode() == LeakageInductanceFilterMode::MINIMIZE_LEAKAGE_RATIO);
    CHECK(factoryScore == explicitScore);
    CHECK(defaultScore == explicitScore);
    CHECK(factoryValid);
    UNSCOPED_INFO("CMC Lk/Lm score " << std::setprecision(17) << factoryScore);
    // Pinned from the pre-change filter (ABT #1176 base 137a2485), same fixture.
    CHECK_THAT(factoryScore, WithinRel(kCmcLeakageRatioBeforeAbt1176, 1e-12));
}

TEST_CASE("Magnetic shunt: LEAKAGE_INDUCTANCE_TARGET scores the distance to the leakage band", "[magnetic-shunt][adviser][magnetic-filter]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);
    double leakage = LeakageInductance().calculate_leakage_inductance(magnetic, 1e6, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    auto inputs = OpenMagnetics::Inputs::create_quick_operating_point(1e6, 26e-6, 25, WaveformLabel::SINUSOIDAL, 1, 0.5, 0, {2.0});
    auto filter = MagneticFilter::factory(MagneticFilters::LEAKAGE_INDUCTANCE_TARGET, inputs);

    SECTION("inside the band") {
        DimensionWithTolerance band;
        band.set_minimum(0.9 * leakage);
        band.set_maximum(1.1 * leakage);
        inputs.get_mutable_design_requirements().set_leakage_inductance(std::vector<DimensionWithTolerance>{band});
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        CHECK(valid);
        CHECK(score == 0.0);
    }
    SECTION("above the band") {
        DimensionWithTolerance band;
        band.set_minimum(0.4 * leakage);
        band.set_maximum(0.5 * leakage);
        inputs.get_mutable_design_requirements().set_leakage_inductance(std::vector<DimensionWithTolerance>{band});
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        CHECK(!valid);
        CHECK_THAT(score, WithinRel((leakage - 0.5 * leakage) / (0.45 * leakage), 1e-9));
    }
    SECTION("nominal only") {
        DimensionWithTolerance band;
        band.set_nominal(2 * leakage);
        inputs.get_mutable_design_requirements().set_leakage_inductance(std::vector<DimensionWithTolerance>{band});
        auto [valid, score] = filter->evaluate_magnetic(&magnetic, &inputs);
        CHECK_THAT(score, WithinRel(0.5, 1e-9));
        (void) valid;
    }
    SECTION("no requirement throws") {
        REQUIRE_THROWS_AS(filter->evaluate_magnetic(&magnetic, &inputs), InvalidInputException);
    }
}

TEST_CASE("Magnetic shunt: size_shunt_for_leakage reaches the target below 0.7 Bs", "[magnetic-shunt][leakage-inductance][adviser]") {
    settings.reset();
    CHECK(!settings.get_coil_adviser_size_magnetic_shunts());
    auto spec = li_2018_case();
    auto magnetic = make_planar_shunt_transformer(spec, false);
    auto material = find_core_material_by_name("C350");
    double frequency = 1e6;
    double windingLeakage = LeakageInductance().calculate_leakage_inductance(magnetic, frequency, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    double target = 1.3 * windingLeakage;
    double gap = 0.1e-3;

    auto sheets = MagneticShuntModel::size_shunt_for_leakage(magnetic, target, material, gap, gap, 1.0, frequency, 0, 1);
    REQUIRE(sheets.size() == 2);
    magnetic.set_shunts(sheets);
    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, frequency, 0, 1);
    CHECK_THAT(result.leakageInductance, WithinRel(target, 1e-4));
    double thickness = sheets[0].get_dimensions()[1];
    CHECK(thickness > 0);
    CHECK(thickness <= 4.1e-3);
    CHECK(result.shunts[0].magneticFluxDensityPerAmpere * 1.0 < 0.7 * Core::get_magnetic_flux_density_saturation(material, 25, false));

    auto unshunted = make_planar_shunt_transformer(spec, false);
    SECTION("target out of reach") {
        REQUIRE_THROWS_AS(MagneticShuntModel::size_shunt_for_leakage(unshunted, 1000 * windingLeakage, material, gap, gap, 1.0, frequency, 0, 1), InvalidInputException);
    }
    SECTION("saturation") {
        REQUIRE_THROWS_AS(MagneticShuntModel::size_shunt_for_leakage(unshunted, target, material, gap, gap, 1e4, frequency, 0, 1), InvalidInputException);
    }
}

TEST_CASE("Magnetic shunt: the 2D painter draws the sheet", "[magnetic-shunt][painter]") {
    settings.reset();
    auto spec = li_2018_case();
    spec.sheetMaterialName = "C350";
    auto magnetic = make_planar_shunt_transformer(spec);
    auto outFile = outputFilePath;
    outFile.append("Test_Magnetic_Shunt_Li_2018.svg");
    std::filesystem::remove(outFile);
    Painter painter(outFile);
    painter.paint_magnetic(magnetic);
    auto svg = painter.export_svg();
    CHECK_THAT(svg, Catch::Matchers::ContainsSubstring("class=\"shunt\""));
    CHECK_THAT(svg, Catch::Matchers::ContainsSubstring("<title>Shunt</title>"));

    auto withoutShunt = magnetic;
    withoutShunt.set_shunts(std::nullopt);
    Painter plainPainter(outFile);
    plainPainter.paint_magnetic(withoutShunt);
    CHECK_THAT(plainPainter.export_svg(), !Catch::Matchers::ContainsSubstring("class=\"shunt\""));
}
