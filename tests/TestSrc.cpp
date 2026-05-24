/**
 * @file TestSrc.cpp
 * @brief Phase 1+2 unit tests for the Series Resonant Converter (SRC).
 *
 * Coverage (SRC_PLAN.md §11 Phase 1+2):
 *   - JSON parsing and run_checks for valid / malformed inputs
 *   - process_design_requirements: turns ratio, Lr/Cr derivation
 *     (Steigerwald 1988 / FHA conventions)
 *   - process_operating_points: steady-state at resonance (Λ=1) gives
 *     M = 1 (HB) and Vout matches the expected n·Vbridge_dc
 *   - Above-resonance CCM: gain decreases monotonically with Λ (FHA)
 *   - Below-resonance: throws (deferred to Phase 3)
 *   - User overrides for Lr / Cr take precedence over Q-derived values
 *   - AdvancedSrc round-trip sets desired turns ratio and resonant components
 */
#include "converter_models/Src.h"
#include "converter_models/Topology.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <variant>

using namespace OpenMagnetics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Build a minimal half-bridge SRC JSON. fr is implicit via fmin/fmax bounds
// (geometric mean) unless resonantFrequency is set.
json make_src_json(double Vin, double Vout, double Iout,
                   double fsw, double fr,
                   double Q = 2.0,
                   const std::string& bridge = "halfBridge") {
    json j;
    j["inputVoltage"]          = {{"nominal", Vin}};
    j["minSwitchingFrequency"] = fr * 0.5;
    j["maxSwitchingFrequency"] = fr * 2.0;
    j["resonantFrequency"]     = fr;
    j["qualityFactor"]         = Q;
    j["bridgeType"]            = bridge;
    j["isolated"]              = true;

    json op;
    op["outputVoltages"]      = {Vout};
    op["outputCurrents"]      = {Iout};
    op["switchingFrequency"]  = fsw;
    op["ambientTemperature"]  = 25.0;
    j["operatingPoints"]      = json::array({op});
    return j;
}

} // namespace


TEST_CASE("SRC: JSON parses and run_checks accepts a valid HB design", "[src-topology]") {
    auto j = make_src_json(/*Vin*/400, /*Vout*/48, /*Iout*/10,
                           /*fsw*/100e3, /*fr*/100e3);
    Src src(j);
    REQUIRE(src.run_checks(/*assert*/false));
    REQUIRE(src.is_bridge_topology());
    REQUIRE(src.topology_kind() == MAS::Topologies::SERIES_RESONANT_CONVERTER);
    REQUIRE(src.get_bridge_voltage_factor() == 0.5);  // HB
}


TEST_CASE("SRC: full-bridge factor is 1.0", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0, "fullBridge");
    Src src(j);
    REQUIRE(src.get_bridge_voltage_factor() == 1.0);
}


TEST_CASE("SRC: process_design_requirements derives turns ratio and tank", "[src-topology]") {
    // HB, fr=100kHz, Q=2, Vin=400, Vout=48, Iout=10A → Rload=4.8Ω
    // n = k_bridge·Vin/Vout = 0.5·400/48 = 4.1667
    // Rac = (8·n²/π²)·Rload = 8·17.36·4.8/π² ≈ 67.5
    // Zr = Q·Rac = 135;  Lr = Zr/(2π·fr) = 215 µH;  Cr = 1/(2π·fr·Zr) = 11.8 nF
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, /*Q*/2.0);
    Src src(j);
    auto dr = src.process_design_requirements();
    REQUIRE(dr.get_topology() == Topologies::SERIES_RESONANT_CONVERTER);

    REQUIRE(dr.get_turns_ratios().size() == 1);
    auto n = dr.get_turns_ratios()[0].get_nominal().value();
    REQUIRE_THAT(n, WithinRel(0.5 * 400.0 / 48.0, 0.02));

    double Lr = src.get_computed_resonant_inductance();
    double Cr = src.get_computed_resonant_capacitance();
    REQUIRE(Lr > 0);
    REQUIRE(Cr > 0);
    double fr_computed = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));
    REQUIRE_THAT(fr_computed, WithinRel(100e3, 0.01));
}


TEST_CASE("SRC: at resonance (Λ=1) FHA gain M = 1 (HB)", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, 2.0);
    Src src(j);
    auto ops = src.process_operating_points(std::vector<double>{}, /*Lm*/0);
    REQUIRE(ops.size() >= 1);

    // FHA gain at Λ=1: M = 1 / sqrt(1 + Q²·(Λ-1/Λ)²) = 1
    REQUIRE_THAT(src.get_last_gain(), WithinRel(1.0, 1e-6));
    REQUIRE_THAT(src.get_last_normalized_fsw(), WithinRel(1.0, 1e-6));
    REQUIRE(src.get_last_is_above_resonance());
    REQUIRE(src.get_last_ir_peak()  > 0);
    REQUIRE(src.get_last_vcr_peak() > 0);
}


TEST_CASE("SRC: gain decreases monotonically above resonance", "[src-topology]") {
    // Same design, sweep fsw above resonance and verify M(Λ) is monotone ↓.
    double fr = 100e3;
    double prev_M = 2.0;
    for (double Lambda : {1.0, 1.1, 1.3, 1.6, 2.0}) {
        auto j = make_src_json(400, 48, 10, /*fsw*/Lambda * fr, fr, 2.0);
        Src src(j);
        src.process_operating_points(std::vector<double>{}, 0);
        double M = src.get_last_gain();
        REQUIRE(M <= prev_M + 1e-9);
        prev_M = M;
    }
    REQUIRE(prev_M < 1.0);  // far above resonance must give M < 1
}


TEST_CASE("SRC: below-resonance operation throws (Phase 3 deferred)", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, /*fsw*/80e3, /*fr*/100e3, 2.0);
    Src src(j);
    REQUIRE_THROWS_AS(
        src.process_operating_points(std::vector<double>{}, 0),
        std::runtime_error);
}


TEST_CASE("SRC: user override of Lr/Cr replaces Q-derived values", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.set_user_resonant_inductance(150e-6);
    src.set_user_resonant_capacitance(10e-9);
    src.process_design_requirements();
    REQUIRE_THAT(src.get_computed_resonant_inductance(),  WithinRel(150e-6, 1e-9));
    REQUIRE_THAT(src.get_computed_resonant_capacitance(), WithinRel(10e-9,  1e-9));
}


TEST_CASE("SRC: missing operating points fails run_checks", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["operatingPoints"] = json::array();
    Src src(j);
    REQUIRE(!src.run_checks(false));
    REQUIRE_THROWS_AS(src.run_checks(true), std::runtime_error);
}


TEST_CASE("SRC: AdvancedSrc round-trips desired turns ratio and Lr/Cr", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["desiredTurnsRatios"]         = {3.5};
    j["desiredResonantInductance"]  = 200e-6;
    j["desiredResonantCapacitance"] = 12e-9;
    AdvancedSrc adv(j);
    auto dr = adv.process_design_requirements();
    REQUIRE(dr.get_turns_ratios().size() == 1);
    REQUIRE_THAT(dr.get_turns_ratios()[0].get_nominal().value(), WithinRel(3.5, 1e-6));
    REQUIRE_THAT(adv.get_computed_resonant_inductance(),  WithinRel(200e-6, 1e-9));
    REQUIRE_THAT(adv.get_computed_resonant_capacitance(), WithinRel(12e-9,  1e-9));
}


// =====================================================================
// CURRENT_DOUBLER rectifier (analytical, Phase-2 extension)
// =====================================================================

TEST_CASE("SRC: currentDoubler emits 1 secondary per output (FB-like winding count)",
          "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["rectifierType"] = "currentDoubler";
    Src src(j);
    REQUIRE(src.get_effective_rectifier_type() == SrcRectifierType::CURRENT_DOUBLER);
    REQUIRE(src.is_center_tapped() == false);
    REQUIRE(src.windings_per_output() == 1);

    auto ops = src.process_operating_points(std::vector<double>{}, 0);
    REQUIRE(!ops.empty());
    // 1 primary + 1 secondary per output (the 2 output inductors Lo1/Lo2 are
    // emitted via get_extra_components_inputs, not as additional windings).
    REQUIRE(ops[0].get_excitations_per_winding().size() == 2);
}


TEST_CASE("SRC: currentDoubler stashes Lo1/Lo2 waveforms for extras", "[src-topology]") {
    // Sanity: per-OP solver must populate the Lo1/Lo2 stash so that
    // get_extra_components_inputs can emit two outputInductor MAS Inputs.
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, /*Q*/2.0);
    j["rectifierType"] = "currentDoubler";
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);

    auto extras = src.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    // Expected layout: { Cr (CAS::Inputs), Lr (MAS Inputs), Lo1 (MAS), Lo2 (MAS) }.
    REQUIRE(extras.size() == 4);
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[0]));   // Cr
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[1]));  // Lr
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[2]));  // Lo1
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[3]));  // Lo2

    auto& lo1 = std::get<OpenMagnetics::Inputs>(extras[2]);
    auto& lo2 = std::get<OpenMagnetics::Inputs>(extras[3]);
    REQUIRE(lo1.get_design_requirements().get_name().value() == "outputInductor1");
    REQUIRE(lo2.get_design_requirements().get_name().value() == "outputInductor2");
    // Lo sizing convention (Vout / (4·fs·0.30·Iout) = 48 / (4·100k·0.3·10)).
    const double Lo_design = 48.0 / (4.0 * 100e3 * 0.30 * 10.0);
    const double Lo1_nom = lo1.get_design_requirements()
                              .get_magnetizing_inductance()
                              .get_nominal().value();
    REQUIRE_THAT(Lo1_nom, WithinRel(Lo_design, 1e-6));
    // KCL: at any sample, ILo1 + ILo2 should be ≈ Iout_dc.
    auto& lo1Ops = lo1.get_operating_points();
    auto& lo2Ops = lo2.get_operating_points();
    REQUIRE(!lo1Ops.empty()); REQUIRE(!lo2Ops.empty());
    auto i1 = lo1Ops[0].get_excitations_per_winding()[0].get_current();
    auto i2 = lo2Ops[0].get_excitations_per_winding()[0].get_current();
    REQUIRE(i1.has_value()); REQUIRE(i2.has_value());
    auto i1Wfm = i1->get_waveform();
    auto i2Wfm = i2->get_waveform();
    REQUIRE(i1Wfm.has_value()); REQUIRE(i2Wfm.has_value());
    // Average the sum over the period; should land at Iout_dc = 10 A.
    const auto& i1Data = i1Wfm->get_data();
    const auto& i2Data = i2Wfm->get_data();
    REQUIRE(i1Data.size() == i2Data.size());
    REQUIRE(!i1Data.empty());
    double sumMean = 0.0;
    for (size_t k = 0; k < i1Data.size(); ++k) sumMean += (i1Data[k] + i2Data[k]);
    sumMean /= i1Data.size();
    REQUIRE_THAT(sumMean, WithinAbs(10.0, 0.5));   // ≈ Iout_dc, ripple-free DC
}


// =====================================================================
// CENTER_TAPPED_DIODE rectifier (analytical, Phase-2 extension)
// =====================================================================

TEST_CASE("SRC: centerTappedDiode design emits 2 windings per output", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["rectifierType"] = "centerTappedDiode";
    Src src(j);

    REQUIRE(src.get_effective_rectifier_type() == SrcRectifierType::CENTER_TAPPED_DIODE);
    REQUIRE(src.is_center_tapped());
    REQUIRE(src.windings_per_output() == 2);

    auto dr = src.process_design_requirements();
    // turns_ratios: wpo=2 entries per output (both equal to the same n_design).
    REQUIRE(dr.get_turns_ratios().size() == 2);
    double n0 = dr.get_turns_ratios()[0].get_nominal().value();
    double n1 = dr.get_turns_ratios()[1].get_nominal().value();
    REQUIRE_THAT(n0, WithinRel(0.5 * 400.0 / 48.0, 0.02));   // same as FB
    REQUIRE(n0 == n1);
    // isolation_sides: 1 primary + 2 (wpo) secondary entries per output.
    REQUIRE(dr.get_isolation_sides().has_value());
    REQUIRE(dr.get_isolation_sides()->size() == 3);
}


TEST_CASE("SRC: centerTappedDiode emits two half-windings with alternating conduction",
          "[src-topology]") {
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, /*Q*/2.0);
    j["rectifierType"] = "centerTappedDiode";
    Src src(j);
    auto ops = src.process_operating_points(std::vector<double>{}, 0);
    REQUIRE(!ops.empty());

    auto& windings = ops[0].get_excitations_per_winding();
    // primary + 2 secondary halves
    REQUIRE(windings.size() == 3);

    auto halfA_i = windings[1].get_current()->get_waveform()->get_data();
    auto halfA_v = windings[1].get_voltage()->get_waveform()->get_data();
    auto halfB_i = windings[2].get_current()->get_waveform()->get_data();
    auto halfB_v = windings[2].get_voltage()->get_waveform()->get_data();
    REQUIRE(!halfA_i.empty());
    REQUIRE(halfA_i.size() == halfB_i.size());

    // Each half-winding current must be non-negative everywhere (CT rectifier
    // only allows one direction per half via the diode).
    double maxA = 0, maxB = 0;
    for (size_t k = 0; k < halfA_i.size(); ++k) {
        REQUIRE(halfA_i[k] >= -1e-9);
        REQUIRE(halfB_i[k] >= -1e-9);
        maxA = std::max(maxA, halfA_i[k]);
        maxB = std::max(maxB, halfB_i[k]);
    }
    REQUIRE(maxA > 0);
    REQUIRE(maxB > 0);
    // Symmetric peaks (within numerical tolerance from polarity-step interpolation).
    REQUIRE_THAT(maxA, WithinRel(maxB, 0.05));

    // Voltage waveforms: bipolar ±Vout square wave, with opposite polarity
    // between halves. Verify peak magnitudes match Vout (samples near the
    // polarity step interpolate to lower magnitudes — accept the peak only).
    double peakA = 0, peakB = 0;
    for (size_t k = 0; k < halfA_v.size(); ++k) {
        peakA = std::max(peakA, std::abs(halfA_v[k]));
        peakB = std::max(peakB, std::abs(halfB_v[k]));
    }
    REQUIRE_THAT(peakA, WithinRel(48.0, 1e-3));
    REQUIRE_THAT(peakB, WithinRel(48.0, 1e-3));

    // Half-windings 180° out of phase: when one is +Vout, the other is -Vout
    // (away from the polarity-step interpolation region). Cross-correlation
    // at lag 0 should be negative ⇒ Σ vA·vB < 0.
    double dot = 0.0;
    for (size_t k = 0; k < halfA_v.size(); ++k)
        dot += halfA_v[k] * halfB_v[k];
    REQUIRE(dot < 0.0);
}


TEST_CASE("SRC: centerTappedDiode tank waveforms match fullBridgeDiode", "[src-topology]") {
    // CT and FB have the same FHA model — only the winding decomposition differs.
    // Primary current peak and Vcr peak must match within 1 %.
    auto jFB = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src srcFB(jFB);
    srcFB.process_operating_points(std::vector<double>{}, 0);

    auto jCT = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    jCT["rectifierType"] = "centerTappedDiode";
    Src srcCT(jCT);
    srcCT.process_operating_points(std::vector<double>{}, 0);

    REQUIRE_THAT(srcCT.get_last_ir_peak(),
                 WithinRel(srcFB.get_last_ir_peak(), 1e-3));
    REQUIRE_THAT(srcCT.get_last_vcr_peak(),
                 WithinRel(srcFB.get_last_vcr_peak(), 1e-3));
    REQUIRE_THAT(srcCT.get_last_gain(),
                 WithinRel(srcFB.get_last_gain(), 1e-9));
}


TEST_CASE("SRC: currentDoubler SPICE codegen succeeds (CD branch lands)",
          "[src-topology]") {
    // CD is now fully supported in both analytical and SPICE paths.
    // Sanity-check that the netlist contains the canonical CD signatures
    // (Lo1/Lo2 returning to vout_neg, per-Lo current sense, no freewheel
    // diodes). A full PtP regression lives in
    // tests/TestSrcReferenceDesignsPtp.cpp.
    auto j = make_src_json(400, 48, 10, 100e3, 100e3);
    j["rectifierType"] = "currentDoubler";
    Src src(j);
    src.process_design_requirements();   // populate Lr/Cr before codegen
    std::string netlist = src.generate_ngspice_circuit(std::vector<double>{4.17}, 0);
    REQUIRE(netlist.find("CD diode rectifier") != std::string::npos);
    REQUIRE(netlist.find("Lo1_o1 ") != std::string::npos);
    REQUIRE(netlist.find("Lo2_o1 ") != std::string::npos);
    REQUIRE(netlist.find("VLo1_sense_o1") != std::string::npos);
    REQUIRE(netlist.find("VLo2_sense_o1") != std::string::npos);
    // Canonical CD topology: no freewheel diodes (each Lo freewheels through
    // the opposite diode + load), and Lo terminates at vout_neg, not vout_pos.
    REQUIRE(netlist.find("Dfw") == std::string::npos);
    REQUIRE(netlist.find(" vout_neg_o1 ") != std::string::npos);
}


// =====================================================================
// get_extra_components_inputs — Cr (CAS) + Lr (MAS Inputs)
// =====================================================================

TEST_CASE("SRC: get_extra_components_inputs throws before processing", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    REQUIRE_THROWS_AS(
        src.get_extra_components_inputs(ExtraComponentsMode::IDEAL),
        std::runtime_error);
}


TEST_CASE("SRC: get_extra_components_inputs REAL requires magnetic", "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);
    REQUIRE_THROWS_AS(
        src.get_extra_components_inputs(ExtraComponentsMode::REAL, std::nullopt),
        std::invalid_argument);
}


TEST_CASE("SRC: get_extra_components_inputs IDEAL returns Cr + Lr", "[src-topology]") {
    // HB, fr=100kHz, Q=2, Vin=400, Vout=48, Iout=10A
    auto j = make_src_json(400, 48, 10, /*fsw*/100e3, /*fr*/100e3, /*Q*/2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);

    auto extras = src.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 2);

    // [0] = Cr as CAS::Inputs (RESONANT application)
    REQUIRE(std::holds_alternative<CAS::Inputs>(extras[0]));
    auto& crInputs = std::get<CAS::Inputs>(extras[0]);
    auto& crDr = crInputs.get_design_requirements();
    REQUIRE(crDr.get_capacitance().get_nominal().has_value());
    double Cr_emit = crDr.get_capacitance().get_nominal().value();
    REQUIRE_THAT(Cr_emit, WithinRel(src.get_computed_resonant_capacitance(), 1e-9));
    REQUIRE(crDr.get_rated_voltage() > 0);
    REQUIRE(crDr.get_role().has_value());
    REQUIRE(crDr.get_role().value() == CAS::Application::RESONANT);
    REQUIRE(!crInputs.get_operating_points().empty());
    {
        auto& exc = crInputs.get_operating_points()[0].get_excitation();
        REQUIRE(exc.get_voltage().has_value());
        REQUIRE(exc.get_current().has_value());
        REQUIRE(!exc.get_voltage()->get_waveform()->get_data().empty());
        REQUIRE(!exc.get_current()->get_waveform()->get_data().empty());
    }

    // [1] = Lr as MAS Inputs (always external in IDEAL)
    REQUIRE(std::holds_alternative<OpenMagnetics::Inputs>(extras[1]));
    auto& lrInputs = std::get<OpenMagnetics::Inputs>(extras[1]);
    auto& lrDr = lrInputs.get_design_requirements();
    REQUIRE(lrDr.get_magnetizing_inductance().get_nominal().has_value());
    double Lr_emit = lrDr.get_magnetizing_inductance().get_nominal().value();
    REQUIRE_THAT(Lr_emit, WithinRel(src.get_computed_resonant_inductance(), 1e-9));
    REQUIRE(lrDr.get_topology().has_value());
    REQUIRE(lrDr.get_topology().value() == Topologies::SERIES_RESONANT_CONVERTER);
    REQUIRE(!lrInputs.get_operating_points().empty());
    {
        auto& exc = lrInputs.get_operating_points()[0].get_excitations_per_winding()[0];
        REQUIRE(exc.get_current().has_value());
        REQUIRE(exc.get_voltage().has_value());
    }
}


TEST_CASE("SRC: get_extra_components_inputs round-trips fr = 1/(2π√(LrCr))",
          "[src-topology]") {
    auto j = make_src_json(400, 48, 10, 100e3, 100e3, 2.0);
    Src src(j);
    src.process_operating_points(std::vector<double>{}, 0);
    auto extras = src.get_extra_components_inputs(ExtraComponentsMode::IDEAL);
    REQUIRE(extras.size() == 2);
    double Cr = std::get<CAS::Inputs>(extras[0])
                    .get_design_requirements().get_capacitance()
                    .get_nominal().value();
    double Lr = std::get<OpenMagnetics::Inputs>(extras[1])
                    .get_design_requirements().get_magnetizing_inductance()
                    .get_nominal().value();
    double fr = 1.0 / (2.0 * M_PI * std::sqrt(Lr * Cr));
    REQUIRE_THAT(fr, WithinRel(100e3, 1e-3));
}

// =========================================================================
// SPICE §8a.5 probe contract — see TestTopologyLlc.cpp counterpart.
// =========================================================================
TEST_CASE("Test_Src_Spice_Probe_Contract",
          "[converter-model][src-topology][spice][probe-contract]") {
    auto runOne = [](const std::string& bridge,
                     OpenMagnetics::BridgeSimulationMode mode,
                     const std::string& label) {
        auto j = make_src_json(/*Vin*/400, /*Vout*/48, /*Iout*/10,
                               /*fsw*/100e3, /*fr*/100e3, /*Q*/2.0, bridge);
        OpenMagnetics::Src src(j);
        src.set_bridge_simulation_mode(mode);
        auto req = src.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios()) turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        std::string netlist = src.generate_ngspice_circuit(turnsRatios, Lm, 0, 0);

        INFO("scenario=" << label);
        REQUIRE(netlist.find("Vq1_sense leg_a_src") == std::string::npos);
        REQUIRE(netlist.find("Vq1_sense bridge_a")  == std::string::npos);
        REQUIRE(netlist.find("Vq3_sense leg_a_src") == std::string::npos);
        REQUIRE(netlist.find("Vq3_sense bridge_a")  == std::string::npos);

        auto savePos = netlist.find(".save");
        REQUIRE(savePos != std::string::npos);
        auto saveLineEnd = netlist.find('\n', savePos);
        std::string saveLine = netlist.substr(savePos, saveLineEnd - savePos);
        if (saveLine.find("i(Vq") != std::string::npos) {
            REQUIRE(netlist.find("Vq1_sense vin_dc") != std::string::npos);
        }
    };

    runOne("fullBridge", OpenMagnetics::BridgeSimulationMode::BEHAVIORAL_PULSE,         "FB_BEHAVIORAL");
    runOne("fullBridge", OpenMagnetics::BridgeSimulationMode::VOLTAGE_CONTROLLED_SWITCH, "FB_SW1");
    runOne("halfBridge", OpenMagnetics::BridgeSimulationMode::BEHAVIORAL_PULSE,         "HB_BEHAVIORAL");
    runOne("halfBridge", OpenMagnetics::BridgeSimulationMode::VOLTAGE_CONTROLLED_SWITCH, "HB_SW1");
}
