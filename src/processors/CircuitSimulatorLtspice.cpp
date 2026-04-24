#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "Defaults.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

// Forward declarations for free functions defined in CircuitSimulatorInterface.cpp.
std::string to_string(double d, size_t precision);
std::vector<std::string> to_string(std::vector<double> v, size_t precision);

// Emit saturating magnetizing inductance for LTspice using Flux= syntax.
// LTspice supports nonlinear inductors with Flux={expression} directly.
static std::string emit_saturating_inductor_ltspice(
    const SaturationParameters& sat,
    const std::string& windingIndex,
    const std::string& nodeIn,
    const std::string& nodeOut) {

    if (!sat.valid) {
        return "";  // Fall back to linear inductor
    }

    std::string s;
    const double mu0 = 4e-7 * M_PI;

    // Calculate saturation current and flux linkage
    double Isat = sat.Isat();
    double lambdaSat = sat.fluxLinkageSat();

    // Calculate gap contribution
    double effectiveMuR = sat.Lmag / (mu0 * sat.primaryTurns * sat.primaryTurns * sat.Ae / sat.le);
    double gapFactor = std::min(0.99, 1.0 / effectiveMuR);
    double Lgap = sat.Lmag * gapFactor;

    // Ensure reasonable values
    if (Isat < 1e-6) Isat = 1e-6;
    if (lambdaSat < 1e-9) lambdaSat = 1e-9;

    s += "* Saturating magnetizing inductance (winding " + windingIndex + ")\n";
    s += "* Bsat=" + to_string(sat.Bsat, 4) + "T, Isat=" + to_string(Isat, 4) + "A\n";

    // LTspice flux linkage model: Lambda = Lambda_sat * tanh(I/I_sat) + L_gap * I
    // The tanh term models core saturation, L_gap*I models the linear air gap
    s += ".param Lambda_sat_" + windingIndex + "=" + to_string(lambdaSat, 12) + "\n";
    s += ".param I_sat_" + windingIndex + "=" + to_string(Isat, 6) + "\n";
    s += ".param L_gap_" + windingIndex + "=" + to_string(Lgap, 12) + "\n";

    // Use LTspice Flux= syntax for nonlinear inductor
    // x is the current through the inductor
    s += "Lmag_" + windingIndex + " " + nodeIn + " " + nodeOut;
    s += " Flux={Lambda_sat_" + windingIndex + "*tanh(x/I_sat_" + windingIndex + ")+L_gap_" + windingIndex + "*x}\n";

    return s;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());

    // Resolve AUTO mode from settings before computing coefficients so we use the right model
    auto resolvedMode_lt = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature, resolvedMode_lt);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance_all_windings(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    std::vector<FractionalPoleNetwork> fracpoleNets_lt;
    std::optional<FractionalPoleNetwork> coreFracNet_lt;
    if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets_lt = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet_lt = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
        catch (...) {}
    }

    parametersString += ".param MagnetizingInductance_Value=" + std::to_string(magnetizingInductance) + "\n";
    parametersString += ".param Permeance=MagnetizingInductance_Value/NumberTurns_1**2\n";

    // Check if saturation modeling is enabled
    auto& settings = Settings::GetInstance();
    bool includeSaturation = settings.get_circuit_simulator_include_saturation();
    SaturationParameters satParams;
    if (includeSaturation) {
        satParams = get_saturation_parameters(magnetic, temperature);
        // Saturation only supported for single winding (inductor) - K coupling doesn't work with nonlinear inductors
        if (satParams.valid && coil.get_functional_description().size() > 1) {
            // For transformers, fall back to linear model with warning comment
            circuitString += "* WARNING: Saturation modeling not supported for multi-winding transformers in LTspice\n";
            circuitString += "* K coupling requires linear inductors. Using linear model.\n";
            includeSaturation = false;
        }
    }

    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, temperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index]);
            if (leakageInductance >= magnetizingInductance) {
                leakageInductance = magnetizingInductance * 0.1;
            }
            double couplingCoefficient = sqrt((magnetizingInductance - leakageInductance) / magnetizingInductance);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient) + "\n";
        }

        std::vector<std::string> c = to_string(acResistanceCoefficientsPerWinding[index], 12);

        // FRACPOLE mode: emit fracpole-based winding subcircuit
        if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
            if (index < fracpoleNets_lt.size()) {
                circuitString += emit_fracpole_winding_spice(fracpoleNets_lt[index], index, is, effectiveResistanceThisWinding);
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            circuitString += "E" + is + " P" + is + "+ Node_R_Lmag_" + is + " P" + is + "+ Node_R_Lmag_" + is + " Laplace = 1 /(" + c[0] + " + " + c[1] + " * sqrt(abs(s)/(2*pi)) + " + c[2] + " * abs(s)/(2*pi))\n";
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "P" + is + "-", "Node_R_Lmag_" + is);
            } else {
                circuitString += "Lmag_" + is + " P" + is + "- Node_R_Lmag_" + is + " {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else if (resolvedMode_lt == CircuitSimulatorExporterCurveFittingModes::ROSANO) {
            const auto& rosanoCoeffs = acResistanceCoefficientsPerWinding[index];
            if (!rosanoCoeffs.empty()) {
                circuitString += emit_winding_rosano_spice(rosanoCoeffs, is, effectiveResistanceThisWinding);
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else {
            // LADDER mode (default)
            // Check if ladder coefficients are valid (inductance values should be small, < 0.1H)
            // If fitting failed, coefficients can be very large (1.0H or more) which breaks the circuit
            bool validLadderCoeffs = true;
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                double inductanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex + 1];
                double resistanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex];
                // Sanity check: inductance should be < 0.1H (100mH), resistance should be < 100 Ohm per ladder element
                if (inductanceVal > 0.1 || inductanceVal < 1e-7 || resistanceVal > 100.0 || resistanceVal < 1e-3) {
                    validLadderCoeffs = false;
                    break;
                }
            }

            if (validLadderCoeffs && acResistanceCoefficientsPerWinding[index].size() >= 2) {
                // Use full ladder network for AC resistance modeling
                for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                    std::string ladderIndexs = std::to_string(ladderIndex);
                    circuitString += "Lladder" + is + "_" + ladderIndexs + " P" + is + "+ Node_Lladder_" + is + "_" + ladderIndexs + " " + c[ladderIndex + 1] + "\n";
                    if (ladderIndex == 0) {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_R_Lmag_" + is + " " + c[ladderIndex] + "\n";
                    }
                    else {
                        circuitString += "Rladder" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + ladderIndexs + " Node_Lladder_" + is + "_" + std::to_string(ladderIndex - 2) + " " + c[ladderIndex] + "\n";
                    }
                }
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            } else {
                // Ladder fitting failed - use simplified model with just DC resistance
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ltspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        if (index > 0) {
            // Each K statement needs a unique name in LTspice (K1, K2, etc.)
            circuitString += "K" + is + " Lmag_1 Lmag_" + is + " {CouplingCoefficient_1" + is + "_Value}\n";
        }

        headerString += " P" + is + "+ P" + is + "-";

    }

    // Core losses: frequency-dependent resistance network in parallel with Lmag
    auto coreLossTopology_lt = static_cast<CoreLossTopology>(settings.get_circuit_simulator_core_loss_topology());
    auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature, coreLossTopology_lt);
    if (!coreResistanceCoefficients.empty()) {
        if (coreLossTopology_lt == CoreLossTopology::ROSANO) {
            circuitString += emit_core_rosano_spice(coreResistanceCoefficients, coil.get_functional_description().size());
        } else {
            circuitString += emit_core_ladder_spice(coreResistanceCoefficients, coil.get_functional_description().size());
        }
    }

    // Mutual resistance: auxiliary winding network for cross-coupling losses (Hesterman 2020)
    // Only for multi-winding transformers when setting is enabled
    size_t numWindings = coil.get_functional_description().size();
    bool includeMutualResistance = settings.get_circuit_simulator_include_mutual_resistance();
    if (includeMutualResistance && numWindings >= 2) {
        auto mutualResistanceCoeffs = CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(magnetic, temperature);
        if (!mutualResistanceCoeffs.empty()) {
            circuitString += emit_mutual_resistance_network_spice(mutualResistanceCoeffs, magnetizingInductance, numWindings);
        }
    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

std::string CircuitSimulatorExporterLtspiceModel::export_magnetic_as_symbol(Magnetic magnetic, std::optional<std::string> filePathOrFile) {
    std::string symbolString = "Version 4\n";
    symbolString += "SymbolType BLOCK\n";

    auto coil = magnetic.get_coil();

    int rectangleSemiWidth = 72;

    int leftSideSize = 16;
    int rightSideSize = 16;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {
        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            leftSideSize += 64;
        }
        else {
            rightSideSize += 64;
        }
    }

    int rectangleHeight = std::max(leftSideSize, rightSideSize);


    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(-rectangleHeight / 2 + 8) + " Left 0 " + magnetic.get_reference() + "\n";
    symbolString += "TEXT " + std::to_string(-rectangleSemiWidth + 8) + " " + std::to_string(rectangleHeight / 2 - 8) + " Left 0 Made with OpenMagnetics\n";


    symbolString += "RECTANGLE Normal " + std::to_string(-rectangleSemiWidth) + " -" + std::to_string(rectangleHeight / 2) + " " + std::to_string(rectangleSemiWidth) + " " + std::to_string(rectangleHeight / 2) + "\n";
    symbolString += "SYMATTR Prefix X\n";
    symbolString += "SYMATTR Value " + fix_filename(magnetic.get_reference()) + "\n";
    symbolString += "SYMATTR ModelFile " + fix_filename(magnetic.get_reference()) + ".cir\n";

    int currentSpiceOrder = 1;
    int currentRectangleLeftSideHeight = -leftSideSize / 2 + 24;
    int currentRectangleRightSideHeight = -rightSideSize / 2 + 24;
    for (size_t index = 0; index < coil.get_functional_description().size(); index++) {

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else{
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "+\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;

        if (coil.get_functional_description()[index].get_isolation_side() == IsolationSide::PRIMARY) {
            symbolString += "PIN " + std::to_string(-rectangleSemiWidth) + " " + std::to_string(currentRectangleLeftSideHeight) + " LEFT 8\n";
            currentRectangleLeftSideHeight += 32;
        }
        else {
            symbolString += "PIN " + std::to_string(rectangleSemiWidth) + " " + std::to_string(currentRectangleRightSideHeight) + " RIGHT 8\n";
            currentRectangleRightSideHeight += 32;
        }
        symbolString += "PINATTR PinName P" + std::to_string(index + 1) + "-\n";
        symbolString += "PINATTR SpiceOrder " + std::to_string(currentSpiceOrder) + "\n";
        currentSpiceOrder++;
    }


    return symbolString;
}

} // namespace OpenMagnetics
