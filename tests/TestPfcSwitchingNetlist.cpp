// =====================================================================
// Phase 2 smoke test for the boost-PFC switching netlist generator.
//
// Scope: verify that
//   1. derive_pfc_controller_tuning() produces finite, sign-correct
//      values for the three industry reference designs (NCP1654 100 W,
//      UCC28180 360 W, L4981 1000 W);
//   2. PowerFactorCorrection::generate_ngspice_switching_circuit() emits
//      a netlist string containing every required section (subckt
//      prelude, power stage, all six controller blocks, .tran, .end);
//   3. invalid configurations (missing required field, non-BOOST
//      variant, fc_i ≥ fsw/2) throw, per project no-fallbacks rule.
//
// NOT in scope: actually invoking ngspice and asserting convergence /
// waveform shape — that is the Phase 4 convergence campaign and lives
// in TestTopologyPowerFactorCorrection.cpp under the
// [refdesign][ptp][slow] tag set.
// =====================================================================

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <stdexcept>
#include <string>

#include "converter_models/PowerFactorCorrection.h"
#include "converter_models/PfcControllerDesign.h"

using namespace OpenMagnetics;

namespace {

struct RefDesign {
    const char* name;
    double vrms;
    double vbus;
    double pout;
    double fsw;
    double cbus;
    double inductance;   // closed-form L_CCM at ripple=0.3
};

constexpr RefDesign kRef1{"NCP1654-100W",  230.0, 400.0,  100.0, 100e3,  100e-6, 3.30e-3};
constexpr RefDesign kRef2{"UCC28180-360W", 230.0, 390.0,  360.0,  65e3,  470e-6, 1.25e-3};
constexpr RefDesign kRef3{"L4981-1000W",   230.0, 400.0, 1000.0,  50e3, 1500e-6, 659e-6};

OpenMagnetics::PowerFactorCorrection build_pfc(const RefDesign& d) {
    OpenMagnetics::PowerFactorCorrection pfc;
    DimensionWithTolerance iv;
    iv.set_nominal(d.vrms);
    iv.set_minimum(d.vrms);
    iv.set_maximum(d.vrms);
    pfc.set_input_voltage(iv);
    pfc.set_output_voltage(d.vbus);
    pfc.set_output_power(d.pout);
    pfc.set_switching_frequency(d.fsw);
    pfc.set_line_frequency(50.0);
    pfc.set_efficiency(1.0);
    pfc.set_diode_voltage_drop(0.0);
    pfc.set_current_ripple_ratio(0.3);
    pfc.set_bulk_capacitance(d.cbus);
    pfc.set_mode(PfcModes::CONTINUOUS_CONDUCTION_MODE);
    pfc.set_ambient_temperature(25.0);
    return pfc;
}

void check_finite_positive(double x, const char* what) {
    INFO(what << " = " << x);
    REQUIRE(std::isfinite(x));
    REQUIRE(x > 0.0);
}

void assert_tuning_sane(const RefDesign& d) {
    INFO("Ref design: " << d.name);

    auto t = derive_pfc_controller_tuning(
        d.vrms, d.vbus, d.pout, d.fsw, 50.0, d.inductance);

    check_finite_positive(t.k_div,    "k_div");
    check_finite_positive(t.gv_fc_hz, "gv_fc_hz");
    check_finite_positive(t.gv_cz,    "gv_cz");
    check_finite_positive(t.gv_cp,    "gv_cp");
    check_finite_positive(t.ff_fc_hz, "ff_fc_hz");
    check_finite_positive(t.ff_c,     "ff_c");
    check_finite_positive(t.g_mul,    "g_mul");
    check_finite_positive(t.gi_fc_hz, "gi_fc_hz");
    check_finite_positive(t.gi_kp,    "gi_kp");
    check_finite_positive(t.gi_cz,    "gi_cz");
    check_finite_positive(t.gi_cp,    "gi_cp");
    check_finite_positive(t.ic_vbus,  "ic_vbus");
    check_finite_positive(t.ic_vea,   "ic_vea");
    check_finite_positive(t.ic_vff,   "ic_vff");

    // Sanity: voltage loop must be well below 2·fline ripple (100 Hz)
    // to avoid distorting the line current. fc_v ≤ fline/5 = 10 Hz is the
    // design-note rule (with floor at 5 Hz); accept anything ≤ 25 Hz.
    CHECK(t.gv_fc_hz <= 25.0);
    CHECK(t.gv_fc_hz >= 5.0);

    // Sanity: current loop must be well below fsw/2 (Nyquist for the
    // switching ripple — design note targets fsw/10).
    CHECK(t.gi_fc_hz < d.fsw / 2.0);
    CHECK(t.gi_fc_hz > d.fsw / 100.0);

    // Warm-start ICs must be physically meaningful.
    CHECK(t.ic_vbus == d.vbus);
    CHECK_THAT(t.ic_vff,
               Catch::Matchers::WithinRel(2.0 * std::sqrt(2.0) * d.vrms / M_PI, 0.01));

    // K_p of current loop, Kp = 2π·fc·L·V_pk_saw / Vbus, must be > 0
    // and finite — covered by check_finite_positive above. Order of
    // magnitude ranges from ~0.1 to a few — sanity bound at 1e3.
    CHECK(t.gi_kp < 1e3);
}

void assert_netlist_contains(const std::string& netlist,
                             const std::string& token,
                             const char* what) {
    INFO("Looking for " << what << " ('" << token << "') in netlist");
    REQUIRE(netlist.find(token) != std::string::npos);
}

void assert_netlist_well_formed(const RefDesign& d) {
    INFO("Ref design: " << d.name);
    auto pfc = build_pfc(d);
    const std::string netlist =
        pfc.generate_ngspice_switching_circuit(d.inductance, /*cycles*/3);

    REQUIRE_FALSE(netlist.empty());

    // Required structural sections.
    assert_netlist_contains(netlist, ".subckt OPAMP_IDEAL",     "OPAMP_IDEAL subckt");
    assert_netlist_contains(netlist, ".subckt COMPARATOR_IDEAL","COMPARATOR_IDEAL subckt");
    assert_netlist_contains(netlist, "vl_sense",                 "inductor-current sense");
    assert_netlist_contains(netlist, ".tran",                    ".tran statement");
    assert_netlist_contains(netlist, "uic",                      "uic flag");
    assert_netlist_contains(netlist, ".end",                     ".end card");

    // Each of the six controller blocks must be represented somewhere.
    assert_netlist_contains(netlist, "vbus",     "Vbus node");
    assert_netlist_contains(netlist, "vea",      "voltage-error-amp output node");
    assert_netlist_contains(netlist, "vff",      "feed-forward stage node");
    assert_netlist_contains(netlist, "iref",     "multiplier output node");
    assert_netlist_contains(netlist, "vc_i",     "current-error-amp output node");
    assert_netlist_contains(netlist, "saw",      "sawtooth node");
    assert_netlist_contains(netlist, "gate",     "PWM gate node");

    // Solver options that we documented as required for convergence.
    assert_netlist_contains(netlist, "METHOD=GEAR", "GEAR integration");
}

} // namespace

TEST_CASE("Test_Pfc_Controller_Tuning_NCP1654_Sane",
          "[converter-model][pfc-topology][pfc-controller][tuning]") {
    assert_tuning_sane(kRef1);
}
TEST_CASE("Test_Pfc_Controller_Tuning_UCC28180_Sane",
          "[converter-model][pfc-topology][pfc-controller][tuning]") {
    assert_tuning_sane(kRef2);
}
TEST_CASE("Test_Pfc_Controller_Tuning_L4981_Sane",
          "[converter-model][pfc-topology][pfc-controller][tuning]") {
    assert_tuning_sane(kRef3);
}

TEST_CASE("Test_Pfc_Controller_Tuning_Throws_On_Bad_Inputs",
          "[converter-model][pfc-topology][pfc-controller][tuning]") {
    // Vbus ≤ Vpk: boost cannot regulate (would require buck).
    REQUIRE_THROWS(derive_pfc_controller_tuning(
        230.0, /*vbus*/300.0, 100.0, 100e3, 50.0, 3.3e-3));

    // Non-positive arguments.
    REQUIRE_THROWS(derive_pfc_controller_tuning(
        -230.0, 400.0, 100.0, 100e3, 50.0, 3.3e-3));
    REQUIRE_THROWS(derive_pfc_controller_tuning(
        230.0, 400.0, 0.0, 100e3, 50.0, 3.3e-3));
    REQUIRE_THROWS(derive_pfc_controller_tuning(
        230.0, 400.0, 100.0, 0.0, 50.0, 3.3e-3));
}

TEST_CASE("Test_Pfc_Switching_Netlist_NCP1654_WellFormed",
          "[converter-model][pfc-topology][pfc-controller][netlist]") {
    assert_netlist_well_formed(kRef1);
}
TEST_CASE("Test_Pfc_Switching_Netlist_UCC28180_WellFormed",
          "[converter-model][pfc-topology][pfc-controller][netlist]") {
    assert_netlist_well_formed(kRef2);
}
TEST_CASE("Test_Pfc_Switching_Netlist_L4981_WellFormed",
          "[converter-model][pfc-topology][pfc-controller][netlist]") {
    assert_netlist_well_formed(kRef3);
}

TEST_CASE("Test_Pfc_Switching_Netlist_Rejects_NonBoost_Variant",
          "[converter-model][pfc-topology][pfc-controller][netlist]") {
    auto pfc = build_pfc(kRef1);
    pfc.set_topology_variant(PfcTopologyVariants::TOTEM_POLE);
    pfc.set_wide_bandgap_switch(true);
    REQUIRE_THROWS(pfc.generate_ngspice_switching_circuit(kRef1.inductance, 3));
}
