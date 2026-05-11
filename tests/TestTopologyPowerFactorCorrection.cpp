// ─────────────────────────────────────────────────────────────────────────────
// TestTopologyPowerFactorCorrection.cpp
//
// Per CONVERTER_MODELS_GOLDEN_GUIDE.md §5.0 / §5.1, every converter topology
// must expose `simulate_and_extract_topology_waveforms()` returning
// `ConverterWaveforms` whose external-port signals satisfy the appropriate
// converter-port gate.
//
// PFC is the one topology in the family whose input port is fundamentally AC
// (rectified line voltage), so the standard `check_dc_ports` helper does NOT
// apply to its input. A PFC-specific helper, `check_pfc_ports`, validates:
//   - input_voltage line-cycle MEAN ≈ 0.9·Vrms and RMS ≈ Vrms (rectified |sin|)
//   - output_voltage MEAN ≈ Vbus_nominal (DC bus)
//   - output_current MEAN ≈ Iload_nominal (DC load)
//
// This file mirrors the structure of TestTopologyBuck.cpp / TestTopologyBoost.cpp
// for the standard golden-tier set: design-requirements, operating-point
// generation, and the §5.1 converter-port gate.
// ─────────────────────────────────────────────────────────────────────────────

#include "converter_models/PowerFactorCorrection.h"
#include "ConverterPortChecks.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <vector>

using namespace OpenMagnetics;

namespace {

// Universal PFC nameplate (230 Vrms / 50 Hz / 400 V bus / 500 W).
constexpr double kVinRmsNom = 230.0;
constexpr double kVinRmsMin = 185.0;
constexpr double kVinRmsMax = 265.0;
constexpr double kVbus      = 400.0;
constexpr double kPout      = 500.0;
constexpr double kFsw       = 100e3;
constexpr double kFline     = 50.0;
constexpr double kEta       = 0.95;
constexpr double kRipple    = 0.3;
constexpr double kVdiode    = 0.6;
// 500 µF bulk cap → ΔVbus_pp = Pout / (π·fline·Cbus·Vbus)
//                  = 500 / (π·50·500e-6·400) ≈ 15.9 V_pp (~4 % of Vbus).
constexpr double kCbus      = 500e-6;

OpenMagnetics::PowerFactorCorrection make_default_pfc(bool sweepInputVoltage = true) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(kVinRmsNom);
    if (sweepInputVoltage) {
        iv.set_minimum(kVinRmsMin);
        iv.set_maximum(kVinRmsMax);
    }
    pfc.set_input_voltage(iv);
    pfc.set_output_voltage(kVbus);
    pfc.set_output_power(kPout);
    pfc.set_switching_frequency(kFsw);
    pfc.set_line_frequency(kFline);
    pfc.set_efficiency(kEta);
    pfc.set_current_ripple_ratio(kRipple);
    pfc.set_diode_voltage_drop(kVdiode);
    pfc.set_bulk_capacitance(kCbus);
    pfc.set_mode(PfcModes::CONTINUOUS_CONDUCTION_MODE);
    pfc.set_ambient_temperature(25.0);
    return pfc;
}

// ─────────────────────────────────────────────────────────────────────────────
// Design requirements: process_design_requirements() must produce a finite,
// positive boost-inductor inductance window (and a single primary-side winding).
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_Design_Requirements",
          "[converter-model][pfc-topology][design-requirements]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    auto dr  = pfc.process_design_requirements();

    REQUIRE((dr.get_magnetizing_inductance().get_nominal().has_value() ||
             dr.get_magnetizing_inductance().get_minimum().has_value()));

    const double L = dr.get_magnetizing_inductance().get_nominal().value_or(
                     dr.get_magnetizing_inductance().get_minimum().value_or(0.0));
    INFO("PFC computed boost inductance: " << L * 1e6 << " µH");
    CHECK(L > 0.0);
    CHECK(std::isfinite(L));
}

// ─────────────────────────────────────────────────────────────────────────────
// Operating point generation: process_operating_points() must return at least
// one OperatingPoint with non-empty current/voltage waveforms on a valid time
// grid spanning at least one mains period.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_OperatingPoints_Generation",
          "[converter-model][pfc-topology][operating-points]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_num_periods_to_extract(1);
    const double L = pfc.calculate_inductance_ccm();
    REQUIRE(L > 0.0);

    auto ops = pfc.process_operating_points({}, L);
    REQUIRE(!ops.empty());

    const auto& exc = ops[0].get_excitations_per_winding();
    REQUIRE(exc.size() == 1);

    auto iWfOpt = exc[0].get_current().value().get_waveform();
    auto vWfOpt = exc[0].get_voltage().value().get_waveform();
    REQUIRE(iWfOpt.has_value());
    REQUIRE(vWfOpt.has_value());
    const auto& iWf = iWfOpt.value();
    const auto& vWf = vWfOpt.value();
    REQUIRE(!iWf.get_data().empty());
    REQUIRE(!vWf.get_data().empty());
    REQUIRE(iWf.get_time().has_value());
    REQUIRE(vWf.get_time().has_value());

    const auto iTime = iWf.get_time().value();
    const double tSpan = iTime.back() - iTime.front();
    INFO("Operating-point time span: " << tSpan * 1e3 << " ms (≥ "
         << (1.0 / kFline) * 1e3 << " ms expected for 1 line cycle)");
    CHECK(tSpan >= 0.99 / kFline);
}

// ─────────────────────────────────────────────────────────────────────────────
// §5.1 converter-port gate (PFC-flavoured).
//
// Sweeps nominal/min/max input RMS; validates that on each ConverterWaveforms:
//   - rectified input voltage line-cycle MEAN ≈ 0.9·Vrms, RMS ≈ Vrms
//   - DC bus output voltage MEAN ≈ Vbus
//   - DC load current MEAN ≈ Pout/Vbus
// per ConverterPortChecks::check_pfc_ports.
//
// This is the analytical PFC waveform synthesis (no ngspice call needed); the
// underlying SPICE netlist is exercised separately by TestNgspiceRunner.cpp.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_ConverterPortWaveforms",
          "[converter-port-waveforms][pfc-topology]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/true);
    const double L = pfc.calculate_inductance_ccm();
    REQUIRE(L > 0.0);

    auto wfs = pfc.simulate_and_extract_topology_waveforms(L, /*numberOfCycles=*/1);
    REQUIRE(!wfs.empty());
    REQUIRE(wfs.size() == 3); // nominal, minimum, maximum (collect_input_voltages order)

    const std::vector<double> vinSweep = {kVinRmsNom, kVinRmsMin, kVinRmsMax};
    const double iLoadNom = kPout / kVbus; // 1.25 A

    for (size_t i = 0; i < wfs.size(); ++i) {
        ConverterPortChecks::check_pfc_ports(wfs[i], "PFC", i,
                                             vinSweep[i], kVbus, iLoadNom);
    }

    // ───────────────────────────────────────────────────────────────────
    // Cbus state-variable model: bus voltage should exhibit a twice-line-
    // frequency ripple of amplitude
    //   ΔVbus_pk = Pin_avg / (2·ω_line·Cbus·Vbus)
    //            = (Pout/η) / (2·2π·fline·Cbus·Vbus)
    // For the parameters above (500 W / 0.95 / 50 Hz / 500 µF / 400 V):
    //   ΔVbus_pk ≈ 8.4 V    (peak-to-peak ≈ 16.8 V, ~4 % of Vbus)
    // ───────────────────────────────────────────────────────────────────
    const double expectedRipplePk = (kPout / kEta) /
                                    (2.0 * 2.0 * M_PI * kFline * kCbus * kVbus);
    for (size_t i = 0; i < wfs.size(); ++i) {
        const auto& d = wfs[i].get_output_voltages()[0].get_data();
        const double vMax = *std::max_element(d.begin(), d.end());
        const double vMin = *std::min_element(d.begin(), d.end());
        const double actualPk = 0.5 * (vMax - vMin);
        INFO("PFC OP " << i << " Vbus ripple pk: actual=" << actualPk
             << " V, expected=" << expectedRipplePk << " V");
        CHECK(actualPk > 0.7 * expectedRipplePk);
        CHECK(actualPk < 1.3 * expectedRipplePk);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Topology variant: TOTEM_POLE PFC (bridgeless, GaN/SiC).
//
// Magnetic-design difference from BOOST:
//   - inductor sees AC bipolar voltage (not rectified |sin|), so the
//     synthesized inductor-current waveform should swing positive and negative
//     across a line cycle (true sine envelope).
//   - duty-cycle calc uses Vd=0 (sync rectification, no boost diode), so the
//     CCM inductance is slightly smaller than for an equivalent boost.
// CCM totemPole REQUIRES wideBandgapSwitch=true; otherwise validate throws.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_TotemPole_RequiresWideBandgapInCcm",
          "[converter-model][pfc-topology][totem-pole]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_topology_variant(PfcTopologyVariants::TOTEM_POLE);
    // Don't set wide_bandgap_switch — defaults to false (Si MOSFET).
    REQUIRE_THROWS_WITH(pfc.calculate_inductance_ccm(),
                        Catch::Matchers::ContainsSubstring("wideBandgapSwitch"));
}

TEST_CASE("Test_Pfc_TotemPole_BipolarInductorCurrent",
          "[converter-model][pfc-topology][totem-pole]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_topology_variant(PfcTopologyVariants::TOTEM_POLE);
    pfc.set_wide_bandgap_switch(true);
    pfc.set_num_periods_to_extract(1);

    // CCM inductance should be slightly smaller than equivalent boost
    // because totem-pole duty calc uses Vd=0 (sync rectification).
    const double L_totem = pfc.calculate_inductance_ccm();
    REQUIRE(L_totem > 0.0);

    auto ops = pfc.process_operating_points({}, L_totem);
    REQUIRE(!ops.empty());
    auto iWf = ops[0].get_excitations_per_winding()[0]
                     .get_current().value().get_waveform().value();
    const auto& iData = iWf.get_data();
    REQUIRE(!iData.empty());

    // Bipolar swing: max>0 AND min<0 with comparable magnitudes (true sine
    // envelope). For boost (rectified) min would be ≥ −ripple/2, ≪ |max|.
    const double iMax = *std::max_element(iData.begin(), iData.end());
    const double iMin = *std::min_element(iData.begin(), iData.end());
    INFO("Totem-pole inductor current: min=" << iMin << " A, max=" << iMax << " A");
    CHECK(iMax > 0.0);
    CHECK(iMin < 0.0);
    // Symmetry: |min| ≈ |max| (sine envelope, ripple is symmetric).
    CHECK(std::abs(iMin + iMax) < 0.2 * std::abs(iMax));
}

// ─────────────────────────────────────────────────────────────────────────────
// Topology variant: INTERLEAVED_BOOST PFC (N parallel boost cells).
//
// Magnetic-design difference from BOOST:
//   - per-phase inductor carries Pout/N → smaller per-phase Ipeak
//   - per-phase ΔI = Ipeak_phase·ripple_ratio is also 1/N of single-phase
//   - per-phase L = Vpk·D/(ΔI·fsw) is therefore N× LARGER than single-phase
//     (total stored magnetic energy across N cores remains equal).
// numberOfPhases < 2 throws; > 3 throws.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_InterleavedBoost_RequiresNumberOfPhases",
          "[converter-model][pfc-topology][interleaved-boost]") {
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_topology_variant(PfcTopologyVariants::INTERLEAVED_BOOST);
    // numberOfPhases unset → defaults to 1, which is invalid for interleaved.
    REQUIRE_THROWS_WITH(pfc.calculate_inductance_ccm(),
                        Catch::Matchers::ContainsSubstring("numberOfPhases"));
    pfc.set_number_of_phases(4);
    REQUIRE_THROWS_WITH(pfc.calculate_inductance_ccm(),
                        Catch::Matchers::ContainsSubstring("must be 2 or 3"));
}

TEST_CASE("Test_Pfc_InterleavedBoost_PerPhaseInductanceScalesWithN",
          "[converter-model][pfc-topology][interleaved-boost]") {
    // Single-phase boost reference.
    auto pfcBoost = make_default_pfc(/*sweepInputVoltage=*/false);
    const double L_single = pfcBoost.calculate_inductance_ccm();
    REQUIRE(L_single > 0.0);

    // Interleaved boost N=2: per-phase L should be ~2× single-phase L.
    auto pfcN2 = make_default_pfc(/*sweepInputVoltage=*/false);
    pfcN2.set_topology_variant(PfcTopologyVariants::INTERLEAVED_BOOST);
    pfcN2.set_number_of_phases(2);
    const double L_n2 = pfcN2.calculate_inductance_ccm();
    INFO("Single-phase boost L: " << L_single * 1e6 << " µH; N=2 per-phase L: "
         << L_n2 * 1e6 << " µH (expected ~2×)");
    CHECK(L_n2 > 1.8 * L_single);
    CHECK(L_n2 < 2.2 * L_single);

    // N=3.
    auto pfcN3 = make_default_pfc(/*sweepInputVoltage=*/false);
    pfcN3.set_topology_variant(PfcTopologyVariants::INTERLEAVED_BOOST);
    pfcN3.set_number_of_phases(3);
    const double L_n3 = pfcN3.calculate_inductance_ccm();
    CHECK(L_n3 > 2.7 * L_single);
    CHECK(L_n3 < 3.3 * L_single);
}

// ─────────────────────────────────────────────────────────────────────────────
// Unsupported variants throw at every engineering entry point.
// ─────────────────────────────────────────────────────────────────────────────
TEST_CASE("Test_Pfc_UnsupportedVariants_Throw",
          "[converter-model][pfc-topology][variants]") {
    for (auto variant : {PfcTopologyVariants::BRIDGELESS,
                         PfcTopologyVariants::SEMI_BRIDGELESS,
                         PfcTopologyVariants::BUCK,
                         PfcTopologyVariants::BUCK_BOOST,
                         PfcTopologyVariants::SEPIC,
                         PfcTopologyVariants::CUK}) {
        auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
        pfc.set_topology_variant(variant);
        REQUIRE_THROWS_WITH(pfc.calculate_inductance_ccm(),
                            Catch::Matchers::ContainsSubstring("not yet implemented"));
    }
    // Vienna gets a dedicated error message redirecting to VIENNA_PLAN.md.
    auto pfc = make_default_pfc(/*sweepInputVoltage=*/false);
    pfc.set_topology_variant(PfcTopologyVariants::VIENNA);
    REQUIRE_THROWS_WITH(pfc.calculate_inductance_ccm(),
                        Catch::Matchers::ContainsSubstring("VIENNA_PLAN.md"));
}

// ─────────────────────────────────────────────────────────────────────────────
// §3.D Phase 6 — industry reference-design Values trios.
//
// Three commercial PFC controller EVMs across the power range (BOOST CCM,
// 230 Vrms / 50 Hz / single-phase):
//
//   1. NCP1654  100 W (low-power consumer / LED driver)
//   2. UCC28180 360 W (mid-range adapter / appliance)
//   3. L4981   1000 W (industrial / server PSU)
//
// Each spec is checked against the closed-form analytical model with
// η = 1 and Vd = 0, so the implementation formulas in
// PowerFactorCorrection.cpp reduce cleanly to:
//
//   D_atVpk   = 1 − Vin_pk / Vbus
//   I_in_avg  = Pout / Vrms
//   I_Lpeak_env = √2 · I_in_avg            (envelope peak, mid-cycle)
//   ΔI        = ripple · I_Lpeak_env       (peak-to-peak switching ripple
//                                            at envelope peak)
//   L_CCM     = Vin_pk · D / (ΔI · Fsw)
//   I_Lpeak   = I_Lpeak_env + ΔI / 2       (envelope + half ripple)
//
// We anchor calculate_inductance_ccm, calculate_duty_cycle,
// calculate_peak_current, and determine_actual_mode to those numbers.
//
// ── Why no analytical-vs-ngspice PtP gate for PFC ─────────────────────────
// The peers in §3.D Phase 6 (Buck / Boost / FlyBack / IsoBuck / IBB) pair an
// "anchor on closed-form" Values test with a "compare analytical waveform to
// switching-circuit ngspice output" PtP test. PFC is structurally different:
//
//   * Its `simulate_and_extract_*` family is purely analytical (rectified-
//     sine envelope + triangular ripple synthesised in C++; no ngspice call).
//   * Its `generate_ngspice_circuit` emits a behavioural-source netlist
//     mathematically equivalent to that synthesis — running ngspice on it
//     does not provide independent confirmation.
//   * A real switching boost-PFC SPICE netlist would require a production-
//     grade average-current-mode controller (multiplier + voltage PI loop +
//     current PI loop, per TI SNVA408B, NXP AN5257, Plexim PFC tutorial).
//     Open-loop feed-forward duty (D = 1 − Vin/Vbus) is unstable when Vbus
//     is the live node, and feed-forward against Vbus_nom needs many line
//     cycles for the bulk-cap dynamic to settle (Rload·Cbus ≫ Tline). A
//     hysteretic peak-current scheme oscillates infinitely at line zero
//     crossings. Building, validating, and converging such a controller is
//     several days of SPICE engineering and is out of scope for a magnetic-
//     component test harness.
//
// The PFC analytical model is exercised end-to-end by the §5.1 converter-port
// gate (Test_Pfc_ConverterPortWaveforms above), which validates the synthesis
// on three input-voltage operating points against the PFC port contract
// (rectified-sine input mean/RMS, DC bus voltage mean, DC load current mean,
// twice-line-frequency Cbus ripple amplitude). That, plus the closed-form
// Values gates below, fully covers the magnetic-design surface of the PFC.
// ─────────────────────────────────────────────────────────────────────────────

struct PfcRefDesignSpec {
    const char* name;
    double      vRms;       // line-cycle RMS input voltage
    double      vBus;       // regulated DC bus voltage
    double      pOut;       // output power
    double      fSw;        // switching frequency
    double      cBus;       // bulk-capacitance (per-design sizing, 1–2 µF/W)
    double      ripple;     // peak-to-peak switching ripple ratio at envelope peak
};

constexpr PfcRefDesignSpec kPfcRefDesign1{
    "NCP1654-100W",  230.0, 400.0,  100.0, 100e3,  100e-6, 0.3 };
constexpr PfcRefDesignSpec kPfcRefDesign2{
    "UCC28180-360W", 230.0, 390.0,  360.0,  65e3,  470e-6, 0.3 };
constexpr PfcRefDesignSpec kPfcRefDesign3{
    "L4981-1000W",   230.0, 400.0, 1000.0,  50e3, 1500e-6, 0.3 };

OpenMagnetics::PowerFactorCorrection build_pfc_from_spec(const PfcRefDesignSpec& s) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(s.vRms);
    // Use min = nom so the worst-case CCM formula (which evaluates at
    // Vrms_min) gives the same answer as the closed-form table at nominal.
    iv.set_minimum(s.vRms);
    iv.set_maximum(s.vRms);
    pfc.set_input_voltage(iv);
    pfc.set_output_voltage(s.vBus);
    pfc.set_output_power(s.pOut);
    pfc.set_switching_frequency(s.fSw);
    pfc.set_line_frequency(50.0);
    pfc.set_efficiency(1.0);              // ideal, so analytical formulas reduce
    pfc.set_diode_voltage_drop(0.0);      // ideal boost diode
    pfc.set_current_ripple_ratio(s.ripple);
    pfc.set_bulk_capacitance(s.cBus);
    pfc.set_mode(PfcModes::CONTINUOUS_CONDUCTION_MODE);
    pfc.set_ambient_temperature(25.0);
    return pfc;
}

void assert_pfc_refdesign_values(const PfcRefDesignSpec& s) {
    auto pfc = build_pfc_from_spec(s);

    const double vinPk        = s.vRms * std::sqrt(2.0);
    const double dExpected    = 1.0 - vinPk / s.vBus;                  // η=1, Vd=0
    const double iInAvg       = s.pOut / s.vRms;                       // η=1
    const double iLpeakEnv    = std::sqrt(2.0) * iInAvg;               // envelope peak
    const double deltaI       = s.ripple * iLpeakEnv;                  // p-p ripple
    const double lCcmExpected = vinPk * dExpected / (deltaI * s.fSw);
    const double iPkExpected  = iLpeakEnv + 0.5 * deltaI;              // envelope + ripple/2

    INFO("Ref design: " << s.name);

    // Inductance — closed-form match within 1 % (formulas algebraically equal).
    const double lCcm = pfc.calculate_inductance_ccm();
    INFO("L_CCM: actual=" << lCcm * 1e6 << " µH, expected=" << lCcmExpected * 1e6 << " µH");
    CHECK_THAT(lCcm, Catch::Matchers::WithinRel(lCcmExpected, 0.01));

    // Duty cycle at line peak — exact reduction with Vd=0.
    const double D = pfc.calculate_duty_cycle(vinPk, s.vBus);
    INFO("D@Vpk: actual=" << D << ", expected=" << dExpected);
    CHECK_THAT(D, Catch::Matchers::WithinRel(dExpected, 0.01));

    // Peak inductor current at sized inductance — should match
    // I_Lpeak_env + ΔI/2 within 1 %.
    const double iPk = pfc.calculate_peak_current(vinPk, lCcm);
    INFO("I_Lpeak: actual=" << iPk << " A, expected=" << iPkExpected << " A");
    CHECK_THAT(iPk, Catch::Matchers::WithinRel(iPkExpected, 0.01));

    // Mode classification: with ripple=0.3, L_CCM = (2/ripple)·L_CrCM = 6.67·L_CrCM,
    // safely above the 5 % tolerance band → must classify as CCM.
    const std::string mode = pfc.determine_actual_mode(lCcm);
    INFO("Mode: " << mode);
    CHECK(mode == "Continuous Conduction Mode");
}

TEST_CASE("Test_Pfc_RefDesign_NCP1654_100W_Values",
          "[converter-model][pfc-topology][refdesign][values]") {
    assert_pfc_refdesign_values(kPfcRefDesign1);
}

TEST_CASE("Test_Pfc_RefDesign_UCC28180_360W_Values",
          "[converter-model][pfc-topology][refdesign][values]") {
    assert_pfc_refdesign_values(kPfcRefDesign2);
}

TEST_CASE("Test_Pfc_RefDesign_L4981_1000W_Values",
          "[converter-model][pfc-topology][refdesign][values]") {
    assert_pfc_refdesign_values(kPfcRefDesign3);
}

} // namespace
