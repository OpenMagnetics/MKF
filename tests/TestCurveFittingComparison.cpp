#include "processors/CircuitSimulatorInterface.h"
#include "TestingUtils.h"
#include "physical_models/WindingLosses.h"
#include "processors/Sweeper.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace OpenMagnetics;

namespace {

// Helper to calculate RMS error between fitted model and actual data
struct CurveFitResult {
    CircuitSimulatorExporterCurveFittingModes mode;
    double rmsError;
    double maxError;
    double avgError;
    size_t numPoints;
    std::vector<double> frequencies;
    std::vector<double> actualResistances;
    std::vector<double> fittedResistances;
};

// Calculate AC resistance using the fitted model coefficients
std::vector<double> calculate_fitted_resistances(
    const std::vector<double>& frequencies,
    const std::vector<double>& coefficients,
    double dcResistance,
    CircuitSimulatorExporterCurveFittingModes mode) {
    
    std::vector<double> fitted;
    fitted.reserve(frequencies.size());
    
    for (double f : frequencies) {
        double r;
        if (mode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            // Analytical model: R = Rdc + h * f^alpha (uses 3 coefficients)
            double h = coefficients[0];
            double alpha = coefficients[1];
            r = dcResistance + h * std::pow(f, alpha);
        } else if (mode == CircuitSimulatorExporterCurveFittingModes::LADDER) {
            // Ladder model uses exactly 5 RL stages (10 coefficients)
            std::vector<double> coeffsCopy = coefficients;
            r = CircuitSimulatorExporter::ladder_model(coeffsCopy.data(), f, dcResistance);
        } else if (mode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
            // Fracpole model uses variable number of RL stages
            std::vector<double> coeffsCopy = coefficients;
            r = CircuitSimulatorExporter::fracpole_model(coeffsCopy.data(), static_cast<int>(coeffsCopy.size()), f, dcResistance);
        } else {
            r = dcResistance;  // Fallback
        }
        fitted.push_back(r);
    }
    return fitted;
}

// Calculate fit quality metrics
CurveFitResult evaluate_fit(
    const std::vector<double>& frequencies,
    const std::vector<double>& actual,
    const std::vector<double>& fitted) {
    
    CurveFitResult result;
    result.frequencies = frequencies;
    result.actualResistances = actual;
    result.fittedResistances = fitted;
    result.numPoints = frequencies.size();
    
    double sumSqError = 0.0;
    double maxErr = 0.0;
    double sumErr = 0.0;
    
    for (size_t i = 0; i < actual.size(); ++i) {
        double error = std::abs(actual[i] - fitted[i]);
        double relError = error / actual[i];  // Relative error
        sumSqError += relError * relError;
        maxErr = std::max(maxErr, relError);
        sumErr += relError;
    }
    
    result.rmsError = std::sqrt(sumSqError / actual.size());
    result.maxError = maxErr;
    result.avgError = sumErr / actual.size();
    
    return result;
}

// Compare all three curve fitting modes
std::vector<CurveFitResult> compare_curve_fitting_modes(
    OpenMagnetics::Magnetic magnetic,
    double temperature,
    size_t windingIndex) {
    
    std::vector<CurveFitResult> results;
    
    // Get actual AC resistance data
    const size_t numPoints = 20;
    const double fMin = 100;
    const double fMax = 10e6;
    
    Curve2D acResistanceData = Sweeper().sweep_winding_resistance_over_frequency(
        magnetic, fMin, fMax, numPoints, windingIndex, temperature);
    
    auto frequencies = acResistanceData.get_x_points();
    auto actualResistances = acResistanceData.get_y_points();
    
    // Get DC resistance
    double dcResistance = WindingLosses::calculate_effective_resistance_of_winding(
        magnetic, windingIndex, 0.1, temperature);
    
    // Test ANALYTICAL mode
    {
        auto coeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
            magnetic, temperature, CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);
        if (windingIndex < coeffs.size()) {
            auto fitted = calculate_fitted_resistances(
                frequencies, coeffs[windingIndex], dcResistance, 
                CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);
            auto result = evaluate_fit(frequencies, actualResistances, fitted);
            result.mode = CircuitSimulatorExporterCurveFittingModes::ANALYTICAL;
            results.push_back(result);
        }
    }
    
    // Test LADDER mode
    {
        auto coeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
            magnetic, temperature, CircuitSimulatorExporterCurveFittingModes::LADDER);
        if (windingIndex < coeffs.size()) {
            auto fitted = calculate_fitted_resistances(
                frequencies, coeffs[windingIndex], dcResistance,
                CircuitSimulatorExporterCurveFittingModes::LADDER);
            auto result = evaluate_fit(frequencies, actualResistances, fitted);
            result.mode = CircuitSimulatorExporterCurveFittingModes::LADDER;
            results.push_back(result);
        }
    }
    
    // Test FRACPOLE mode
    {
        auto coeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
            magnetic, temperature, CircuitSimulatorExporterCurveFittingModes::FRACPOLE);
        if (windingIndex < coeffs.size()) {
            auto fitted = calculate_fitted_resistances(
                frequencies, coeffs[windingIndex], dcResistance,
                CircuitSimulatorExporterCurveFittingModes::FRACPOLE);
            auto result = evaluate_fit(frequencies, actualResistances, fitted);
            result.mode = CircuitSimulatorExporterCurveFittingModes::FRACPOLE;
            results.push_back(result);
        }
    }
    
    return results;
}

const char* mode_to_string(CircuitSimulatorExporterCurveFittingModes mode) {
    switch (mode) {
        case CircuitSimulatorExporterCurveFittingModes::ANALYTICAL: return "ANALYTICAL";
        case CircuitSimulatorExporterCurveFittingModes::LADDER: return "LADDER";
        case CircuitSimulatorExporterCurveFittingModes::FRACPOLE: return "FRACPOLE";
        case CircuitSimulatorExporterCurveFittingModes::AUTO: return "AUTO";
        default: return "UNKNOWN";
    }
}

TEST_CASE("Test_CurveFitting_Compare_All_Modes", "[curve-fitting][comparison][analysis]") {
    // Create a test magnetic with significant skin effect
    std::vector<int64_t> numberTurns = {30, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "PQ 35/35";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C90";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0003, 3);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    double temperature = 100.0;  // Higher temp = more significant effects
    
    // Compare all modes for primary winding
    auto results = compare_curve_fitting_modes(magnetic, temperature, 0);
    
    REQUIRE(results.size() == 3);  // All three modes should work
    
    // Log results
    INFO("Curve Fitting Comparison Results for Winding 0 (Primary):");
    INFO("=========================================================");
    
    for (const auto& result : results) {
        INFO(mode_to_string(result.mode) << ":");
        INFO("  RMS Error: " << (result.rmsError * 100) << "%");
        INFO("  Max Error: " << (result.maxError * 100) << "%");
        INFO("  Avg Error: " << (result.avgError * 100) << "%");
    }
    
    // Find best fit
    auto bestFit = *std::min_element(results.begin(), results.end(),
        [](const CurveFitResult& a, const CurveFitResult& b) {
            return a.rmsError < b.rmsError;
        });
    
    INFO("\nBest fitting mode: " << mode_to_string(bestFit.mode));
    INFO("With RMS error: " << (bestFit.rmsError * 100) << "%");
    
    // At least one mode should have excellent fit (RMS error < 5%)
    // Ladder typically performs best for complex skin effect
    CHECK(bestFit.rmsError < 0.05);
    
    // Verify ladder is working (should have very low error)
    auto ladderResult = std::find_if(results.begin(), results.end(),
        [](const CurveFitResult& r) { return r.mode == CircuitSimulatorExporterCurveFittingModes::LADDER; });
    REQUIRE(ladderResult != results.end());
    CHECK(ladderResult->rmsError < 0.01);  // Ladder should have < 1% error
}

TEST_CASE("Test_CurveFitting_Ladder_vs_Math_Detail", "[curve-fitting][ladder][math]") {
    // Test with a magnetic that has high frequency content
    std::vector<int64_t> numberTurns = {50, 10};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 55/28/21";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C95";  // Higher permeability
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0005, 2);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Sweep over wide frequency range
    const double fMin = 100;
    const double fMax = 1e6;  // 1 MHz
    const size_t numPoints = 30;
    
    Curve2D acResistanceData = Sweeper().sweep_winding_resistance_over_frequency(
        magnetic, fMin, fMax, numPoints, 0, 100);
    
    auto frequencies = acResistanceData.get_x_points();
    auto actualResistances = acResistanceData.get_y_points();
    
    // Get coefficients for both modes
    auto mathCoeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
        magnetic, 100, CircuitSimulatorExporterCurveFittingModes::ANALYTICAL);
    auto ladderCoeffs = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(
        magnetic, 100, CircuitSimulatorExporterCurveFittingModes::LADDER);
    
    double dcResistance = WindingLosses::calculate_effective_resistance_of_winding(
        magnetic, 0, 0.1, 100);
    
    // Calculate fitted values for math mode
    std::vector<double> mathFitted;
    for (double f : frequencies) {
        double h = mathCoeffs[0][0];
        double alpha = mathCoeffs[0][1];
        mathFitted.push_back(dcResistance + h * std::pow(f, alpha));
    }
    
    // Calculate fitted values for ladder mode
    std::vector<double> ladderFitted;
    for (double f : frequencies) {
        std::vector<double> coeffsCopy = ladderCoeffs[0];
        ladderFitted.push_back(CircuitSimulatorExporter::ladder_model(
            coeffsCopy.data(), f, dcResistance));
    }
    
    // Evaluate both
    auto mathResult = evaluate_fit(frequencies, actualResistances, mathFitted);
    auto ladderResult = evaluate_fit(frequencies, actualResistances, ladderFitted);
    
    INFO("Math (Analytical) Fit:");
    INFO("  RMS Error: " << (mathResult.rmsError * 100) << "%");
    INFO("  Max Error: " << (mathResult.maxError * 100) << "%");
    
    INFO("Ladder Fit:");
    INFO("  RMS Error: " << (ladderResult.rmsError * 100) << "%");
    INFO("  Max Error: " << (ladderResult.maxError * 100) << "%");
    
    // Ladder should fit very well (< 1% error) for skin effect modeling
    CHECK(ladderResult.rmsError < 0.01);
    
    // Analytical model may not fit as well for complex geometries
    // We just verify it produces results (even if high error)
    CHECK(mathResult.rmsError > 0.0);
    
    // Log which one is better
    if (ladderResult.rmsError < mathResult.rmsError) {
        INFO("Ladder fits better by " << ((mathResult.rmsError - ladderResult.rmsError) * 100) << "%");
    } else {
        INFO("Math fits better by " << ((ladderResult.rmsError - mathResult.rmsError) * 100) << "%");
    }
}

}  // namespace
