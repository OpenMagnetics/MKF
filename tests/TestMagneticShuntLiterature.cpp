// Magnetic shunts (MAS-RFC 0015): more published prototypes with MEASURED leakage inductance for the
// "Shunt" leakage method (ABT #1252; model ABT #1176). TestMagneticShunt.cpp holds Li 2018 and Zhang 2014.
//
// Every number below is the paper's, with its table or figure. What a paper does not publish is listed
// per case as an assumption. The reference is always the measurement, never the paper's own model.
//
//   [Ansari] S. A. Ansari, J. N. Davidson, M. P. Foster, "Fully-Integrated Planar Transformer With a
//            Segmental Shunt for LLC Resonant Converters", IEEE Trans. Ind. Electron., 2021
//            (IEEE Xplore 9559762; accepted version White Rose eprint 178974,
//            https://eprints.whiterose.ac.uk/id/eprint/178974/1/All-21-TIE-2060.pdf).
//            Fig. 2 (topologies), Table I (specification), Table II (mus 800 for both shunts),
//            Table III (measured at 200 kHz, Omicron Bode 100):
//              five-segment: Lm 29.6 uH, Llk 9.23 uH; two-segment: Lm 28.95 uH, Llk 8.95 uH.
//            Core E32/6/20/R-3F4 (planar E + E); shunt 3F4 ferrite, mus 800, t_sh 1.5 mm; N_P 10 (n_P 5
//            layers x k_P 2 turns), N_S 2 (n_S 2 layers x k_S 1 turn); copper h_P = h_S 35 um, insulation
//            h_dP = h_dS 30 um; primary track width 3.1 mm, secondary 5.7 mm.
//              five-segment: horizontal air gap (shunt to each core half) l_fg1 0.13 mm, vertical gaps
//              l_fg2 0.5 mm, four of them, at the column faces (eq. (12): R_fS2 = (b_w - 2 l_fg2)/...),
//              windings x_fp = x_fs = 1.5 mm from the shunt.
//              two-segment: horizontal gap between the halves l_tg1 0.28 mm, one shunt segment per
//              window with l_tg2 0.5 mm gaps to both columns, windings x_tp = x_ts = 0.76 mm from it.
//   [Tan]    W. Tan, X. Margueron, L. Taylor, N. Idir, "Leakage Inductance Analytical Calculation for
//            Planar Components With Leakage Layers", IEEE Trans. Power Electron. 31(6), 4462-4473, 2016
//            (HAL hal-01886543). Sec. IV.A and Fig. 13(c): planar CM choke, two 3C90 E43 cores, window
//            13.3 mm x 10.8 mm, 8 copper layers of 70 um with 4 turns each (4 layers per winding, 16
//            turns), track width 2.1 mm, 0.7 mm between tracks, first track 1.9 mm from the window wall,
//            0.7 mm between the two layers of a board, 1.75 mm between boards, 2.2 mm between the
//            windings around a 0.96 mm FPC C350 layer (mu_r 9) spanning the window; FPC inside the core
//            only (Sec. IV.D). Measured leakage: 14.6 uH, read from Table III: the model values
//            13.4 / 14.7 / 12.1 uH are given with errors of 8.2 % / < 1 % / 17.1 % against the
//            measurement, which all three fit only for a measured 14.6 uH (13.4/(1-0.082) = 14.60,
//            12.1/(1-0.171) = 14.60, and 14.7 within 1 %); Fig. 16 shows the measured curve above the
//            13.4 uH line, flat from 1 kHz to 1 MHz.
//
// Considered and not encoded (ABT #1252 comment has the details):
//   - Tan 2016 LLC transformer (E38, 0.2 mm FPC): the copper-to-FPC distance is not dimensioned, and the
//     winding term, about 80 % of the leakage there, scales with it.
//   - Ansari, Davidson, Foster, IEEE TIE 70(3) 2023 (eprint 186238) and IEEE OJPEL 3, 2022 (eprint
//     181798): measured, but the shunts lie across the front and back faces of the core (outsideWindow
//     in RFC 0015), which the Shunt method does not model (NotImplementedException).
//   - Ansari et al. IECON 2020 (eprint 186257) and PEMD 2022 (eprint 188638): FEA only, no measurement.
//
// Cases the model misses by more than the 15 % used for Zhang (10 % for Lm, as for Li) are findings, not
// bugs of the test: they are hidden ([.]) and tagged only [shunt-literature], so neither the default run
// nor "[magnetic-shunt]" includes them; run them with "[shunt-literature]". At MKF e31c73ac the three
// leakage cases fail (+93 % to +123 %) and the two Lm cases pass.
// See ABT #1250 (the field model ignores the sheet) and ABT #1240 (field model under-estimates).

#include "physical_models/LeakageInductance.h"
#include "physical_models/MagneticShunt.h"
#include "physical_models/MagnetizingInductance.h"
#include "constructive_models/Bobbin.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdio>
#include <numbers>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

double relative_error(double value, double reference) {
    return std::fabs(value - reference) / reference;
}

OpenMagnetics::Wire literature_planar_track(double width, double thickness) {
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

// One copper layer as drawn in the paper, in the +x window.
struct LiteratureLayer {
    size_t winding;          // 0 = primary (source), 1 = secondary
    double copperCentre;     // y of the copper centre (m), core mating plane at 0
};

struct LiteratureTransformer {
    std::string shapeName;
    std::string coreMaterial;
    json gapping;
    double copperThickness;
    std::vector<int64_t> turnsPerLayer;   // per winding
    std::vector<double> trackWidth;       // per winding
    std::vector<double> trackClearance;   // per winding, between adjacent tracks of a layer
    // Distance from the central column face to the first track, per winding; nullopt centres the
    // tracks of a layer in the window.
    std::vector<std::optional<double>> firstTrackFromCentralColumn;
    // Measure firstTrackFromCentralColumn from the lateral column face instead (sensitivity check).
    bool firstTrackFromLateralColumn = false;
    std::vector<LiteratureLayer> layers;  // top to bottom
};

// Builds the magnetic with MKF's planar winder and then puts every layer's copper where the paper draws
// it. The winder takes one insulation value per winding pair; the papers' stacks alternate distances
// (board core, board-to-board), so the drawn coordinates are written onto the turns and layers.
OpenMagnetics::Magnetic make_literature_transformer(const LiteratureTransformer& spec) {
    auto core = OpenMagneticsTesting::get_quick_core(spec.shapeName, spec.gapping, 1, spec.coreMaterial);
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);
    auto windowDimensions = bobbin.get_winding_window_dimensions();
    double windowWidth = windowDimensions[0];
    double windowHeight = windowDimensions[1];
    double windowInnerEdge = bobbin.get_processed_description()->get_winding_windows()[0].get_coordinates().value()[0] - windowWidth / 2;

    std::vector<size_t> layersPerWinding(2, 0);
    std::vector<size_t> stackUp;
    for (auto& layer : spec.layers) {
        stackUp.push_back(layer.winding);
        layersPerWinding[layer.winding]++;
    }

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < 2; ++windingIndex) {
        OpenMagnetics::Winding winding;
        winding.set_name(windingIndex == 0 ? "Primary" : "Secondary");
        winding.set_number_turns(spec.turnsPerLayer[windingIndex] * static_cast<int64_t>(layersPerWinding[windingIndex]));
        winding.set_number_parallels(1);
        winding.set_isolation_side(windingIndex == 0 ? IsolationSide::PRIMARY : IsolationSide::SECONDARY);
        winding.set_wire(literature_planar_track(spec.trackWidth[windingIndex], spec.copperThickness));
        coil.get_mutable_functional_description().push_back(winding);
    }
    coil.set_bobbin(bobbin);
    settings.set_coil_wind_even_if_not_fit(true);
    coil.set_section_alignment(CoilAlignment::CENTERED);
    // Any small uniform insulation: the vertical positions are overwritten below.
    std::map<std::pair<size_t, size_t>, double> insulation;
    insulation[{0, 0}] = 0.1e-3;
    insulation[{1, 1}] = 0.1e-3;
    insulation[{0, 1}] = 0.1e-3;
    coil.wind_planar(stackUp, 0, {{0, spec.trackClearance[0]}, {1, spec.trackClearance[1]}}, insulation, 0);
    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.get_layers_description());

    // Conduction layers, top to bottom, as the winder laid them.
    auto layers = coil.get_layers_description().value();
    std::vector<size_t> conductionLayerIndexes;
    for (size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
        if (layers[layerIndex].get_type() == ElectricalType::CONDUCTION) {
            conductionLayerIndexes.push_back(layerIndex);
        }
    }
    REQUIRE(conductionLayerIndexes.size() == spec.layers.size());
    std::sort(conductionLayerIndexes.begin(), conductionLayerIndexes.end(), [&](size_t a, size_t b) {
        return layers[a].get_coordinates()[1] > layers[b].get_coordinates()[1];
    });

    auto turns = coil.get_turns_description().value();
    for (size_t order = 0; order < conductionLayerIndexes.size(); ++order) {
        auto& layer = layers[conductionLayerIndexes[order]];
        auto& drawn = spec.layers[order];
        REQUIRE(coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding()) == drawn.winding);
        auto coordinates = layer.get_coordinates();
        coordinates[1] = drawn.copperCentre;
        layer.set_coordinates(coordinates);

        std::vector<size_t> turnIndexes;
        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
            if (turns[turnIndex].get_layer() == layer.get_name()) {
                turnIndexes.push_back(turnIndex);
            }
        }
        REQUIRE(static_cast<int64_t>(turnIndexes.size()) == spec.turnsPerLayer[drawn.winding]);
        std::sort(turnIndexes.begin(), turnIndexes.end(), [&](size_t a, size_t b) {
            return turns[a].get_coordinates()[0] < turns[b].get_coordinates()[0];
        });
        double width = spec.trackWidth[drawn.winding];
        double clearance = spec.trackClearance[drawn.winding];
        double blockWidth = static_cast<double>(turnIndexes.size()) * width + static_cast<double>(turnIndexes.size() - 1) * clearance;
        double firstTrackEdge = windowInnerEdge + (windowWidth - blockWidth) / 2;
        if (spec.firstTrackFromCentralColumn[drawn.winding]) {
            double offset = spec.firstTrackFromCentralColumn[drawn.winding].value();
            firstTrackEdge = spec.firstTrackFromLateralColumn ? windowInnerEdge + windowWidth - offset - blockWidth : windowInnerEdge + offset;
        }
        for (size_t position = 0; position < turnIndexes.size(); ++position) {
            auto& turn = turns[turnIndexes[position]];
            auto turnCoordinates = turn.get_coordinates();
            double x = firstTrackEdge + static_cast<double>(position) * (width + clearance) + width / 2;
            // Rectangular column: the turn length grows by 2 pi per metre of radial offset.
            turn.set_length(turn.get_length() + 2 * std::numbers::pi * (x - turnCoordinates[0]));
            turnCoordinates[0] = x;
            turnCoordinates[1] = drawn.copperCentre;
            turn.set_coordinates(turnCoordinates);
            REQUIRE(x - width / 2 >= windowInnerEdge - 1e-9);
            REQUIRE(x + width / 2 <= windowInnerEdge + windowWidth + 1e-9);
            REQUIRE(std::fabs(drawn.copperCentre) + spec.copperThickness / 2 <= windowHeight / 2 + 1e-9);
        }
    }
    coil.set_layers_description(layers);
    coil.set_turns_description(turns);

    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

// Copper centres of a stack of layers of thickness h, listed from the top, given the y of the first
// copper's upper face and the clear distances between consecutive coppers.
std::vector<double> stack_downwards(double topFace, double copperThickness, const std::vector<double>& clearances) {
    std::vector<double> centres;
    double face = topFace;
    centres.push_back(face - copperThickness / 2);
    for (double clearance : clearances) {
        face -= copperThickness + clearance;
        centres.push_back(face - copperThickness / 2);
    }
    return centres;
}

// One in-window sheet per winding window, centred on the mating plane.
std::vector<MagneticShunt> window_sheets(Core core, double thickness, double innerGap, double outerGap, const std::string& material) {
    auto columns = core.get_columns();
    double innerFace = columns[core.get_main_column_index()].get_width() / 2;
    double outerFace = 0;
    for (auto& column : columns) {
        if (column.get_coordinates()[0] > 0) {
            outerFace = column.get_coordinates()[0] - column.get_width() / 2;
        }
    }
    std::vector<MagneticShunt> sheets;
    for (int side : {1, -1}) {
        MagneticShunt sheet;
        sheet.set_name(side > 0 ? "Shunt +x" : "Shunt -x");
        sheet.set_placement(MagneticShuntPlacement::IN_WINDOW);
        // A sheet abutting a column (gap 0) is drawn 1 pm clear of it: centre +- width/2 otherwise rounds
        // a fraction of a femtometre into the column, which the model reads as running into it.
        const double abutmentRounding = 1e-12;
        double start = innerFace + std::max(innerGap, abutmentRounding);
        double end = outerFace - std::max(outerGap, abutmentRounding);
        sheet.set_coordinates({side * (start + end) / 2, 0, 0});
        sheet.set_dimensions({end - start, thickness, core.get_depth()});
        MagneticShuntGapToColumns gaps;
        gaps.set_inner(innerGap);
        gaps.set_outer(outerGap);
        sheet.set_gap_to_columns(gaps);
        sheet.set_material(material);
        sheets.push_back(sheet);
    }
    return sheets;
}

struct LiteratureLeakage {
    double winding = 0;
    double sheetNetwork = 0;
    double displacedAir = 0;
    double total = 0;
    double noShunt = 0;
};

LiteratureLeakage evaluate_leakage(const std::string& label, OpenMagnetics::Magnetic magnetic, double frequency, const std::vector<double>& permeabilityPerShunt, double measured) {
    auto result = LeakageInductance().calculate_shunt_leakage(magnetic, frequency, 0, 1, 1, permeabilityPerShunt);
    LiteratureLeakage leakage;
    leakage.winding = result.windingLeakageInductance;
    for (auto& contribution : result.shunts) {
        leakage.sheetNetwork += contribution.networkLeakageInductance;
        leakage.displacedAir += contribution.displacedAirLeakageInductance;
    }
    leakage.total = result.leakageInductance;
    auto withoutShunt = magnetic;
    withoutShunt.set_shunts(std::nullopt);
    auto noShunt = LeakageInductance().calculate_leakage_inductance(withoutShunt, frequency, 0, 1);
    REQUIRE(noShunt.get_method_used() == "Energy");
    leakage.noShunt = noShunt.get_leakage_inductance_per_winding()[0].get_nominal().value();
    std::printf("[shunt-literature] %s: winding %.4g uH, sheet network %.4g uH, displaced air %.4g uH, total %.4g uH, "
                "no-shunt %.4g uH, measured %.4g uH, error %+.1f %%\n",
                label.c_str(), leakage.winding * 1e6, leakage.sheetNetwork * 1e6, leakage.displacedAir * 1e6, leakage.total * 1e6,
                leakage.noShunt * 1e6, measured * 1e6, 100 * (leakage.total - measured) / measured);
    return leakage;
}

// [Ansari] Table I, both topologies: primary 5 layers above the shunt, secondary 2 layers below.
LiteratureTransformer ansari_2021_transformer(double coreSeparation, double windingToShunt, double shuntThickness) {
    LiteratureTransformer spec;
    spec.shapeName = "E 32/6/20";
    spec.coreMaterial = "3F4";
    spec.gapping = OpenMagneticsTesting::get_spacer_gap(coreSeparation);
    spec.copperThickness = 35e-6;
    spec.turnsPerLayer = {2, 1};
    spec.trackWidth = {3.1e-3, 5.7e-3};
    // Assumption: the clearance between the two primary tracks of a layer is not published; MKF's
    // minimum wire-to-wire distance is used and each layer's tracks are centred in the window.
    spec.trackClearance = {Defaults().minimumWireToWireDistance, Defaults().minimumWireToWireDistance};
    spec.firstTrackFromCentralColumn = {std::nullopt, std::nullopt};
    double insulation = 30e-6;
    auto primary = stack_downwards(shuntThickness / 2 + windingToShunt + 5 * spec.copperThickness + 4 * insulation, spec.copperThickness,
                                   {insulation, insulation, insulation, insulation});
    auto secondary = stack_downwards(-shuntThickness / 2 - windingToShunt, spec.copperThickness, {insulation});
    for (double y : primary) {
        spec.layers.push_back({0, y});
    }
    for (double y : secondary) {
        spec.layers.push_back({1, y});
    }
    return spec;
}

// [Tan] Fig. 13(c): upper winding (source) above the FPC layer, lower winding below, symmetric.
LiteratureTransformer tan_2016_cm_choke() {
    LiteratureTransformer spec;
    spec.shapeName = "E 43/10/28";
    spec.coreMaterial = "3C90";
    spec.gapping = OpenMagneticsTesting::get_residual_gap();
    spec.copperThickness = 70e-6;
    spec.turnsPerLayer = {4, 4};
    spec.trackWidth = {2.1e-3, 2.1e-3};
    spec.trackClearance = {0.7e-3, 0.7e-3};
    // Assumption: Fig. 13(c) does not say which window wall is the central column; the 1.9 mm is taken
    // from the central column face.
    spec.firstTrackFromCentralColumn = {1.9e-3, 1.9e-3};
    // Assumptions: 0.7 mm, 1.75 mm and 2.2 mm are clear copper-to-copper distances, and the FPC sits in
    // the middle of the 2.2 mm (drawn so), on the core mating plane.
    double h = spec.copperThickness;
    double topFace = 1.1e-3 + h + 0.7e-3 + h + 1.75e-3 + h + 0.7e-3 + h;
    auto upper = stack_downwards(topFace, h, {0.7e-3, 1.75e-3, 0.7e-3});
    auto lower = stack_downwards(-1.1e-3, h, {0.7e-3, 1.75e-3, 0.7e-3});
    for (double y : upper) {
        spec.layers.push_back({0, y});
    }
    for (double y : lower) {
        spec.layers.push_back({1, y});
    }
    return spec;
}

} // namespace

TEST_CASE("Magnetic shunt literature: Ansari 2021 five-segment shunt, E32/6/20 3F4, measured Llk 9.23 uH", "[.][shunt-literature]") {
    settings.reset();
    // The shunt lies between the halves with l_fg1 = 0.13 mm of air on each side: the halves are
    // 1.5 + 2 x 0.13 = 1.76 mm apart at every column.
    double shuntThickness = 1.5e-3;
    auto spec = ansari_2021_transformer(shuntThickness + 2 * 0.13e-3, 1.5e-3, shuntThickness);
    auto magnetic = make_literature_transformer(spec);
    // The two in-window segments, each with its l_fg2 = 0.5 mm gaps at the column faces. The three
    // segments inside the column gaps are not sheets the Shunt method takes (a segmented sheet running
    // into a column gap throws NotImplemented); on the leakage path they are the column face the gap
    // reluctance sees. Their series 0.13 mm air gaps to the halves (Ansari R_fg1, R_fgg) are not in
    // MKF's network.
    magnetic.set_shunts(window_sheets(magnetic.get_core(), shuntThickness, 0.5e-3, 0.5e-3, "3F4"));
    double measured = 9.23e-6;
    // Table II: mus 800 for the 3F4 shunt, passed as the paper's value.
    auto leakage = evaluate_leakage("Ansari 2021 five-segment", magnetic, 200e3, {800.0, 800.0}, measured);
    CHECK_THAT(leakage.noShunt, Catch::Matchers::WithinRel(leakage.winding, 1e-12));
    CHECK(relative_error(leakage.total, measured) <= 0.15);
}

TEST_CASE("Magnetic shunt literature: Ansari 2021 five-segment shunt, measured Lm 29.6 uH", "[shunt-literature][magnetic-shunt][magnetizing-inductance]") {
    settings.reset();
    double shuntThickness = 1.5e-3;
    auto spec = ansari_2021_transformer(shuntThickness + 2 * 0.13e-3, 1.5e-3, shuntThickness);
    auto magnetic = make_literature_transformer(spec);
    // For Lm the ferrite segments inside the three column gaps matter: a continuous 1.5 mm sheet on the
    // mating plane gives each gap MKF's covered length g - t + t/mus; its in-window part does not enter
    // the magnetizing path in MKF.
    MagneticShunt sheet;
    sheet.set_name("Column segments");
    sheet.set_placement(MagneticShuntPlacement::BETWEEN_SECTIONS);
    sheet.set_coordinates({0, 0, 0});
    sheet.set_dimensions({magnetic.get_core().get_width(), shuntThickness, magnetic.get_core().get_depth()});
    sheet.set_material("3F4");
    magnetic.set_shunts(std::vector<MagneticShunt>{sheet});
    auto core = MagneticShuntModel::apply_shunts_to_gapping(magnetic, Defaults().ambientTemperature, 200e3, std::vector<double>{800.0});
    double magnetizingInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(core, magnetic.get_coil()).get_magnetizing_inductance().get_nominal().value();
    double measured = 29.6e-6;
    std::printf("[shunt-literature] Ansari 2021 five-segment: Lm %.4g uH, measured %.4g uH, error %+.1f %%\n", magnetizingInductance * 1e6, measured * 1e6,
                100 * (magnetizingInductance - measured) / measured);
    CHECK(relative_error(magnetizingInductance, measured) <= 0.10);
}

TEST_CASE("Magnetic shunt literature: Ansari 2021 two-segment shunt, E32/6/20 3F4, measured Llk 8.95 uH", "[.][shunt-literature]") {
    settings.reset();
    double shuntThickness = 1.5e-3;
    auto spec = ansari_2021_transformer(0.28e-3, 0.76e-3, shuntThickness);
    auto magnetic = make_literature_transformer(spec);
    // One segment per window, l_tg2 = 0.5 mm to both columns.
    magnetic.set_shunts(window_sheets(magnetic.get_core(), shuntThickness, 0.5e-3, 0.5e-3, "3F4"));
    double measured = 8.95e-6;
    auto leakage = evaluate_leakage("Ansari 2021 two-segment", magnetic, 200e3, {800.0, 800.0}, measured);
    CHECK_THAT(leakage.noShunt, Catch::Matchers::WithinRel(leakage.winding, 1e-12));
    CHECK(relative_error(leakage.total, measured) <= 0.15);
}

TEST_CASE("Magnetic shunt literature: Ansari 2021 two-segment shunt, measured Lm 28.95 uH", "[shunt-literature][magnetic-shunt][magnetizing-inductance]") {
    settings.reset();
    auto spec = ansari_2021_transformer(0.28e-3, 0.76e-3, 1.5e-3);
    auto magnetic = make_literature_transformer(spec);
    // The segments sit in the window, not in the 0.28 mm column gaps: MKF's Lm is that of the gapped core.
    double magnetizingInductance = MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance().get_nominal().value();
    double measured = 28.95e-6;
    std::printf("[shunt-literature] Ansari 2021 two-segment: Lm %.4g uH, measured %.4g uH, error %+.1f %%\n", magnetizingInductance * 1e6, measured * 1e6,
                100 * (magnetizingInductance - measured) / measured);
    CHECK(relative_error(magnetizingInductance, measured) <= 0.10);
}

TEST_CASE("Magnetic shunt literature: Tan 2016 planar CM choke, E43 3C90 with 0.96 mm FPC C350, measured Llk 14.6 uH", "[.][shunt-literature]") {
    settings.reset();
    auto spec = tan_2016_cm_choke();
    auto magnetic = make_literature_transformer(spec);
    // The FPC layer spans the window wall to wall (Fig. 13(c)), inside the core only, depth = core depth.
    // C350 is in the core-material database (mu_i 9); the paper's mu_r 9 is passed explicitly.
    magnetic.set_shunts(window_sheets(magnetic.get_core(), 0.96e-3, 0.0, 0.0, "C350"));
    double measured = 14.6e-6;
    // Fig. 16: flat from 1 kHz to 1 MHz; evaluated at 100 kHz.
    auto leakage = evaluate_leakage("Tan 2016 CM choke", magnetic, 100e3, {9.0, 9.0}, measured);
    CHECK_THAT(leakage.noShunt, Catch::Matchers::WithinRel(leakage.winding, 1e-12));
    CHECK(relative_error(leakage.total, measured) <= 0.15);

    // Sensitivity to the one geometric assumption, the side the 1.9 mm is measured from (reported only).
    spec.firstTrackFromLateralColumn = true;
    auto mirrored = make_literature_transformer(spec);
    mirrored.set_shunts(window_sheets(mirrored.get_core(), 0.96e-3, 0.0, 0.0, "C350"));
    evaluate_leakage("Tan 2016 CM choke, 1.9 mm from the lateral column", mirrored, 100e3, {9.0, 9.0}, measured);
}

TEST_CASE("Magnetic shunt literature: rewriting a planar stack to the coordinates the winder produced leaves the Energy leakage unchanged", "[magnetic-shunt][leakage-inductance]") {
    settings.reset();
    // Guards the fixture builder above: the Li 2018 stack wound by MKF (as in TestMagneticShunt.cpp),
    // then rebuilt by make_literature_transformer from the winder's own layer heights.
    auto core = OpenMagneticsTesting::get_quick_core("E 32/6/20", OpenMagneticsTesting::get_spacer_gap(0.28e-3), 1, "3F46");
    auto bobbin = OpenMagnetics::Bobbin::create_quick_bobbin(core, true);
    double windowWidth = bobbin.get_winding_window_dimensions()[0];
    double border = Defaults().minimumBorderToWireDistance;
    double wireToWire = Defaults().minimumWireToWireDistance;
    auto trackWidth = [&](int64_t turnsPerLayer) {
        return (windowWidth - 2 * border - static_cast<double>(turnsPerLayer - 1) * wireToWire) / static_cast<double>(turnsPerLayer);
    };

    OpenMagnetics::Coil coil;
    for (size_t windingIndex = 0; windingIndex < 2; ++windingIndex) {
        OpenMagnetics::Winding winding;
        winding.set_name(windingIndex == 0 ? "Primary" : "Secondary");
        winding.set_number_turns(windingIndex == 0 ? 8 : 4);
        winding.set_number_parallels(1);
        winding.set_isolation_side(windingIndex == 0 ? IsolationSide::PRIMARY : IsolationSide::SECONDARY);
        winding.set_wire(literature_planar_track(trackWidth(windingIndex == 0 ? 2 : 1), 0.07e-3));
        coil.get_mutable_functional_description().push_back(winding);
    }
    coil.set_bobbin(bobbin);
    settings.set_coil_wind_even_if_not_fit(true);
    coil.set_section_alignment(CoilAlignment::CENTERED);
    std::map<std::pair<size_t, size_t>, double> insulation;
    insulation[{0, 0}] = 0.25e-3;
    insulation[{1, 1}] = 0.25e-3;
    insulation[{0, 1}] = 4.1e-3;
    coil.wind_planar({0, 0, 0, 0, 1, 1, 1, 1}, border, {{0, wireToWire}, {1, wireToWire}}, insulation, 0);
    REQUIRE(coil.get_turns_description());
    OpenMagnetics::Magnetic wound;
    wound.set_core(core);
    wound.set_coil(coil);

    LiteratureTransformer spec;
    spec.shapeName = "E 32/6/20";
    spec.coreMaterial = "3F46";
    spec.gapping = OpenMagneticsTesting::get_spacer_gap(0.28e-3);
    spec.copperThickness = 0.07e-3;
    spec.turnsPerLayer = {2, 1};
    spec.trackWidth = {trackWidth(2), trackWidth(1)};
    spec.trackClearance = {wireToWire, wireToWire};
    spec.firstTrackFromCentralColumn = {std::nullopt, std::nullopt};
    std::vector<std::pair<double, size_t>> heights;
    auto woundCoilTurns = coil.get_turns_description().value();
    for (auto& turn : woundCoilTurns) {
        double y = turn.get_coordinates()[1];
        size_t windingIndex = turn.get_winding() == "Primary" ? 0 : 1;
        if (std::none_of(heights.begin(), heights.end(), [&](auto& h) { return std::fabs(h.first - y) < 1e-12; })) {
            heights.push_back({y, windingIndex});
        }
    }
    std::sort(heights.begin(), heights.end(), [](auto& a, auto& b) { return a.first > b.first; });
    for (auto& [y, windingIndex] : heights) {
        spec.layers.push_back({windingIndex, y});
    }
    auto rebuilt = make_literature_transformer(spec);

    auto woundTurns = wound.get_coil().get_turns_description().value();
    auto rebuiltTurns = rebuilt.get_coil().get_turns_description().value();
    REQUIRE(woundTurns.size() == rebuiltTurns.size());
    for (size_t turnIndex = 0; turnIndex < woundTurns.size(); ++turnIndex) {
        CHECK_THAT(rebuiltTurns[turnIndex].get_coordinates()[0], Catch::Matchers::WithinAbs(woundTurns[turnIndex].get_coordinates()[0], 1e-9));
        CHECK_THAT(rebuiltTurns[turnIndex].get_coordinates()[1], Catch::Matchers::WithinAbs(woundTurns[turnIndex].get_coordinates()[1], 1e-9));
        CHECK_THAT(rebuiltTurns[turnIndex].get_length(), Catch::Matchers::WithinAbs(woundTurns[turnIndex].get_length(), 1e-9));
    }
    double woundLeakage = LeakageInductance().calculate_leakage_inductance(wound, 1e6, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    double rebuiltLeakage = LeakageInductance().calculate_leakage_inductance(rebuilt, 1e6, 0, 1).get_leakage_inductance_per_winding()[0].get_nominal().value();
    std::printf("[shunt-literature] fixture guard, Li 2018 stack without the sheet: wound %.4g uH, rebuilt %.4g uH\n", woundLeakage * 1e6, rebuiltLeakage * 1e6);
    CHECK_THAT(rebuiltLeakage, Catch::Matchers::WithinRel(woundLeakage, 1e-6));
}
