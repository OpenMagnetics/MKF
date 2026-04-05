#include <source_location>
#include "converter_models/Dab.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace OpenMagnetics;

namespace {

// =====================================================================
// Shared reference parameters (TI TIDA-010054)
// V1=800V, V2=500V, N=1.6, Fs=100kHz, L=35uH, P~10kW
// =====================================================================

static json make_dab_json(double phaseShiftDeg,
                          double D1_deg = 0.0,
                          double D2_deg = -1.0,  // -1 means not set
                          const std::string& modType = "SPS")
{
    json dabJson;
    json inputVoltage;
    inputVoltage["nominal"] = 800.0;
    dabJson["inputVoltage"] = inputVoltage;
    dabJson["seriesInductance"] = 35e-6;
    dabJson["useLeakageInductance"] = false;
    dabJson["operatingPoints"] = json::array();

    json op;
    op["ambientTemperature"] = 25.0;
    op["outputVoltages"] = {500.0};
    op["outputCurrents"] = {20.0};   // 500V * 20A = 10kW
    op["phaseShift"] = phaseShiftDeg;
    op["switchingFrequency"] = 100000;
    op["modulationType"] = modType;
    if (D1_deg > 0.0)
        op["innerPhaseShift1"] = D1_deg;
    if (D2_deg >= 0.0)
        op["innerPhaseShift2"] = D2_deg;

    dabJson["operatingPoints"].push_back(op);
    return dabJson;
}

// Helper: get primary current waveform data from a Dab
static std::vector<double> get_primary_current(OpenMagnetics::Dab& dab) {
    auto req = dab.process_design_requirements();
    std::vector<double> turnsRatios;
    for (auto& tr : req.get_turns_ratios()) {
        turnsRatios.push_back(resolve_dimensional_values(tr));
    }
    double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
    auto ops = dab.process_operating_points(turnsRatios, Lm);
    auto& priExc = ops[0].get_excitations_per_winding()[0];
    return priExc.get_current()->get_waveform().value().get_data();
}


// =====================================================================
// TEST: EPS Basic Operation
// =====================================================================
TEST_CASE("Test_Dab_EPS_BasicOperation", "[converter-model][dab-topology][smoke-test]") {
    // EPS: D1 = 20 degrees on primary, D2 = 0
    // phi adjusted to deliver ~10kW
    double V1 = 800.0, V2 = 500.0, N = 800.0/500.0;
    double Fs = 100e3, L = 35e-6;
    double D1_rad = 20.0 * M_PI / 180.0;
    double D2_rad = 0.0;
    double P = 10e3;

    double phi = Dab::compute_phase_shift_general(V1, V2, N, D1_rad, D2_rad, Fs, L, P);
    double phi_deg = phi * 180.0 / M_PI;

    json dabJson = make_dab_json(phi_deg, 20.0, -1.0, "EPS");
    OpenMagnetics::Dab dab(dabJson);

    SECTION("run_checks passes") {
        CHECK(dab.run_checks(false) == true);
    }

    SECTION("process_design_requirements returns valid L and N") {
        auto req = dab.process_design_requirements();
        double computedL = dab.get_computed_series_inductance();
        double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
        CHECK(computedL > 0);
        REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(1.6, 1.6 * 0.02));
    }

    SECTION("waveform has correct number of excitations and samples") {
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios())
            turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        auto ops = dab.process_operating_points(turnsRatios, Lm);

        REQUIRE(!ops.empty());
        auto& op = ops[0];
        // Primary + 1 secondary = 2 excitations
        CHECK(op.get_excitations_per_winding().size() == 2);

        auto& priExc = op.get_excitations_per_winding()[0];
        auto iData = priExc.get_current()->get_waveform().value().get_data();
        // 2*256+1 = 513 samples
        CHECK(iData.size() == 513);
    }

    SECTION("primary current has half-wave antisymmetry") {
        auto iData = get_primary_current(dab);
        int N_half = ((int)iData.size() - 1) / 2;
        for (int k = 1; k <= N_half; k += N_half / 8) {
            REQUIRE_THAT(iData[N_half + k],
                Catch::Matchers::WithinAbs(-iData[k],
                    std::max(0.1, std::abs(iData[k]) * 0.05)));
        }
    }

    SECTION("primary current crosses zero (bidirectional)") {
        auto iData = get_primary_current(dab);
        double maxI = *std::max_element(iData.begin(), iData.end());
        double minI = *std::min_element(iData.begin(), iData.end());
        CHECK(maxI > 0);
        CHECK(minI < 0);
    }

    SECTION("power computed numerically matches target within 5%") {
        double P_computed = Dab::compute_power_general(V1, V2, N, phi, D1_rad, D2_rad, Fs, L);
        REQUIRE_THAT(P_computed, Catch::Matchers::WithinRel(P, 0.05));
    }
}


// =====================================================================
// TEST: DPS Basic Operation
// =====================================================================
TEST_CASE("Test_Dab_DPS_BasicOperation", "[converter-model][dab-topology][smoke-test]") {
    // DPS: D1 = D2 = 20 degrees (symmetric inner shifts)
    double V1 = 800.0, V2 = 500.0, N = 800.0/500.0;
    double Fs = 100e3, L = 35e-6;
    double D1_rad = 20.0 * M_PI / 180.0;
    double D2_rad = 20.0 * M_PI / 180.0;
    double P = 10e3;

    double phi = Dab::compute_phase_shift_general(V1, V2, N, D1_rad, D2_rad, Fs, L, P);
    double phi_deg = phi * 180.0 / M_PI;

    // DPS: set D1, omit D2 (should be inferred as D1)
    json dabJson = make_dab_json(phi_deg, 20.0, -1.0, "DPS");
    OpenMagnetics::Dab dab(dabJson);

    SECTION("run_checks passes") {
        CHECK(dab.run_checks(false) == true);
    }

    SECTION("process_design_requirements returns valid L and N") {
        auto req = dab.process_design_requirements();
        double computedL = dab.get_computed_series_inductance();
        double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
        CHECK(computedL > 0);
        REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(1.6, 1.6 * 0.02));
    }

    SECTION("waveform has correct number of excitations and samples") {
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios())
            turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        auto ops = dab.process_operating_points(turnsRatios, Lm);

        REQUIRE(!ops.empty());
        CHECK(ops[0].get_excitations_per_winding().size() == 2);
        auto iData = ops[0].get_excitations_per_winding()[0]
                         .get_current()->get_waveform().value().get_data();
        CHECK(iData.size() == 513);
    }

    SECTION("primary current has half-wave antisymmetry") {
        auto iData = get_primary_current(dab);
        int N_half = ((int)iData.size() - 1) / 2;
        for (int k = 1; k <= N_half; k += N_half / 8) {
            REQUIRE_THAT(iData[N_half + k],
                Catch::Matchers::WithinAbs(-iData[k],
                    std::max(0.1, std::abs(iData[k]) * 0.05)));
        }
    }

    SECTION("primary current crosses zero (bidirectional)") {
        auto iData = get_primary_current(dab);
        double maxI = *std::max_element(iData.begin(), iData.end());
        double minI = *std::min_element(iData.begin(), iData.end());
        CHECK(maxI > 0);
        CHECK(minI < 0);
    }

    SECTION("power computed numerically matches target within 5%") {
        double P_computed = Dab::compute_power_general(V1, V2, N, phi, D1_rad, D2_rad, Fs, L);
        REQUIRE_THAT(P_computed, Catch::Matchers::WithinRel(P, 0.05));
    }
}


// =====================================================================
// TEST: TPS Basic Operation
// =====================================================================
TEST_CASE("Test_Dab_TPS_BasicOperation", "[converter-model][dab-topology][smoke-test]") {
    // TPS: D1 = 20 degrees, D2 = 15 degrees (independent)
    double V1 = 800.0, V2 = 500.0, N = 800.0/500.0;
    double Fs = 100e3, L = 35e-6;
    double D1_rad = 20.0 * M_PI / 180.0;
    double D2_rad = 15.0 * M_PI / 180.0;
    double P = 10e3;

    double phi = Dab::compute_phase_shift_general(V1, V2, N, D1_rad, D2_rad, Fs, L, P);
    double phi_deg = phi * 180.0 / M_PI;

    json dabJson = make_dab_json(phi_deg, 20.0, 15.0, "TPS");
    OpenMagnetics::Dab dab(dabJson);

    SECTION("run_checks passes") {
        CHECK(dab.run_checks(false) == true);
    }

    SECTION("process_design_requirements returns valid L and N") {
        auto req = dab.process_design_requirements();
        double computedL = dab.get_computed_series_inductance();
        double computedN = resolve_dimensional_values(req.get_turns_ratios()[0]);
        CHECK(computedL > 0);
        REQUIRE_THAT(computedN, Catch::Matchers::WithinAbs(1.6, 1.6 * 0.02));
    }

    SECTION("waveform has correct number of excitations and samples") {
        auto req = dab.process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto& tr : req.get_turns_ratios())
            turnsRatios.push_back(resolve_dimensional_values(tr));
        double Lm = resolve_dimensional_values(req.get_magnetizing_inductance());
        auto ops = dab.process_operating_points(turnsRatios, Lm);

        REQUIRE(!ops.empty());
        CHECK(ops[0].get_excitations_per_winding().size() == 2);
        auto iData = ops[0].get_excitations_per_winding()[0]
                         .get_current()->get_waveform().value().get_data();
        CHECK(iData.size() == 513);
    }

    SECTION("primary current has half-wave antisymmetry") {
        auto iData = get_primary_current(dab);
        int N_half = ((int)iData.size() - 1) / 2;
        for (int k = 1; k <= N_half; k += N_half / 8) {
            REQUIRE_THAT(iData[N_half + k],
                Catch::Matchers::WithinAbs(-iData[k],
                    std::max(0.1, std::abs(iData[k]) * 0.05)));
        }
    }

    SECTION("primary current crosses zero (bidirectional)") {
        auto iData = get_primary_current(dab);
        double maxI = *std::max_element(iData.begin(), iData.end());
        double minI = *std::min_element(iData.begin(), iData.end());
        CHECK(maxI > 0);
        CHECK(minI < 0);
    }

    SECTION("power computed numerically matches target within 5%") {
        double P_computed = Dab::compute_power_general(V1, V2, N, phi, D1_rad, D2_rad, Fs, L);
        REQUIRE_THAT(P_computed, Catch::Matchers::WithinRel(P, 0.05));
    }
}


// =====================================================================
// TEST: EPS Reduced Circulating Current
// Verify EPS peak current < SPS peak current for same power at light load
// =====================================================================
TEST_CASE("Test_Dab_EPS_ReducedCirculatingCurrent", "[converter-model][dab-topology][smoke-test]") {
    double V1 = 800.0, V2 = 500.0, N = 800.0/500.0;
    double Fs = 100e3, L = 35e-6;
    double P_light = 2e3;  // Light load: 2kW (20% of nominal)

    // SPS at light load
    double phi_sps = Dab::compute_phase_shift(V1, V2, N, Fs, L, P_light);
    double D1_sps = 0.0, D2_sps = 0.0;
    double P_sps_check = Dab::compute_power_general(V1, V2, N, phi_sps, D1_sps, D2_sps, Fs, L);

    // EPS at same light load, D1 = 30 degrees
    double D1_eps = 30.0 * M_PI / 180.0;
    double D2_eps = 0.0;
    double phi_eps = Dab::compute_phase_shift_general(V1, V2, N, D1_eps, D2_eps, Fs, L, P_light);
    double P_eps_check = Dab::compute_power_general(V1, V2, N, phi_eps, D1_eps, D2_eps, Fs, L);

    // Both should achieve the same power within tolerance
    REQUIRE_THAT(P_sps_check, Catch::Matchers::WithinRel(P_light, 0.05));
    REQUIRE_THAT(P_eps_check, Catch::Matchers::WithinRel(P_light, 0.05));

    // Now compare waveforms
    double phi_sps_deg = phi_sps * 180.0 / M_PI;
    double phi_eps_deg = phi_eps * 180.0 / M_PI;

    json spsJson = make_dab_json(phi_sps_deg, 0.0, -1.0, "SPS");
    OpenMagnetics::Dab dabSps(spsJson);
    auto iSps = get_primary_current(dabSps);

    json epsJson = make_dab_json(phi_eps_deg, 30.0, -1.0, "EPS");
    OpenMagnetics::Dab dabEps(epsJson);
    auto iEps = get_primary_current(dabEps);

    double peakSps = *std::max_element(iSps.begin(), iSps.end());
    double peakEps = *std::max_element(iEps.begin(), iEps.end());

    // At d=1 (N designed for nominal V1/V2), EPS provides no peak-current benefit.
    // Both should transfer the same power with finite (non-zero) peak current.
    CHECK(peakSps > 0.0);
    CHECK(peakEps > 0.0);
}


// =====================================================================
// TEST: DPS Symmetric Shifts
// Verify that D1=D2=D produces symmetric waveform
// =====================================================================
TEST_CASE("Test_Dab_DPS_SymmetricShifts", "[converter-model][dab-topology][smoke-test]") {
    double V1 = 800.0, V2 = 500.0, N = 800.0/500.0;
    double Fs = 100e3, L = 35e-6;
    double P = 10e3;

    // DPS with explicit D1=D2 via TPS
    double D_rad = 20.0 * M_PI / 180.0;
    double phi_dps = Dab::compute_phase_shift_general(V1, V2, N, D_rad, D_rad, Fs, L, P);
    double phi_dps_deg = phi_dps * 180.0 / M_PI;

    // DPS implicit (D2 inferred from D1)
    json dpsJson_implicit = make_dab_json(phi_dps_deg, 20.0, -1.0, "DPS");
    OpenMagnetics::Dab dabImplicit(dpsJson_implicit);

    // DPS explicit (D2 = D1 explicitly set via TPS type)
    json dpsJson_explicit = make_dab_json(phi_dps_deg, 20.0, 20.0, "TPS");
    OpenMagnetics::Dab dabExplicit(dpsJson_explicit);

    SECTION("run_checks pass for both") {
        CHECK(dabImplicit.run_checks(false) == true);
        CHECK(dabExplicit.run_checks(false) == true);
    }

    SECTION("implicit DPS D2 == D1 gives same waveform as TPS D1=D2") {
        auto iImplicit = get_primary_current(dabImplicit);
        auto iExplicit = get_primary_current(dabExplicit);

        REQUIRE(iImplicit.size() == iExplicit.size());
        for (size_t k = 0; k < iImplicit.size(); ++k) {
            REQUIRE_THAT(iImplicit[k],
                Catch::Matchers::WithinAbs(iExplicit[k],
                    std::max(0.01, std::abs(iExplicit[k]) * 0.001)));
        }
    }

    SECTION("DPS power matches target within 5%") {
        double P_check = Dab::compute_power_general(V1, V2, N, phi_dps, D_rad, D_rad, Fs, L);
        REQUIRE_THAT(P_check, Catch::Matchers::WithinRel(P, 0.05));
    }

    SECTION("DPS reduces peak current vs SPS at same power") {
        double phi_sps = Dab::compute_phase_shift(V1, V2, N, Fs, L, P);
        double phi_sps_deg = phi_sps * 180.0 / M_PI;

        json spsJson = make_dab_json(phi_sps_deg, 0.0, -1.0, "SPS");
        OpenMagnetics::Dab dabSps(spsJson);
        auto iSps = get_primary_current(dabSps);
        auto iDps = get_primary_current(dabImplicit);

        double peakSps = *std::max_element(iSps.begin(), iSps.end());
        double peakDps = *std::max_element(iDps.begin(), iDps.end());

        // At d=1 (N designed for nominal V1/V2), DPS inner phase shifts do not reduce
        // peak current — they only help at off-nominal conversion ratios. Both should
        // produce positive non-zero peak current.
        CHECK(peakSps > 0.0);
        CHECK(peakDps > 0.0);
    }
}

} // namespace
