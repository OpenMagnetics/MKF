#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "Defaults.h"

#include <cmath>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

// Forward declarations for free functions defined in CircuitSimulatorInterface.cpp.
std::string to_string(double d, size_t precision);
std::vector<std::string> to_string(std::vector<double> v, size_t precision);

// Emit saturating magnetizing inductance subcircuit fragment for ngspice.
// This behavioral-source form is ngspice-specific; LTspice uses Flux={expr} instead.
static std::string emit_saturating_inductor_ngspice(
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

    // The model uses: V = d(Lambda)/dt where Lambda = Lambda_sat * tanh(I/I_sat) + L_gap * I
    // L_gap is the linear component (from air gap), approximately L_mag * (1 - 1/effective_permeability)
    // For a typical gapped core, most inductance comes from the gap, so L_gap ≈ L_mag
    double effectiveMuR = sat.Lmag / (mu0 * sat.primaryTurns * sat.primaryTurns * sat.Ae / sat.le);
    double gapFactor = 1.0 / effectiveMuR;  // Fraction of reluctance from gap
    double Lgap = sat.Lmag * gapFactor;  // Linear (gap) component

    // Ensure reasonable values
    if (Isat < 1e-6) Isat = 1e-6;
    if (lambdaSat < 1e-9) lambdaSat = 1e-9;
    if (Lgap < 1e-12) Lgap = sat.Lmag * 0.5;  // Fallback

    s += "* Saturating magnetizing inductance (winding " + windingIndex + ")\n";
    s += "* Bsat=" + to_string(sat.Bsat, 4) + "T, Isat=" + to_string(Isat, 4) + "A\n";

    // ngspice behavioral inductor using flux linkage
    // V = dLambda/dt = dLambda/dI * dI/dt
    // For Lambda = Lsat*tanh(I/Isat) + Lgap*I:
    // dLambda/dI = Lsat/Isat * sech²(I/Isat) + Lgap
    //            = Lsat/Isat * (1 - tanh²(I/Isat)) + Lgap
    // We use the GVALUE element with Laplace to model this

    // Simpler approach: Use the inductor with polynomial current dependence
    // L(I) = L0 / (1 + (I/Isat)^2) gives smooth saturation
    double L0 = sat.Lmag;

    s += ".param Lmag_" + windingIndex + "_L0=" + to_string(L0, 12) + "\n";
    s += ".param Lmag_" + windingIndex + "_Isat=" + to_string(Isat, 6) + "\n";

    // Use behavioral source: V = L(I) * dI/dt
    // L(I) = L0 / (1 + (abs(I)/Isat)^2)
    // Note: ngspice behavioral sources use I() for current through voltage source
    std::string senseName = "Vmag_sense_" + windingIndex;
    std::string fluxName = "Lmag_flux_" + windingIndex;

    // Current sense (zero volt source)
    s += senseName + " " + nodeIn + " " + fluxName + " 0\n";

    // Behavioral inductor: V = L(I) * dI/dt
    // Using EL (behavioral voltage source) with ddt() function
    s += "BLmag_" + windingIndex + " " + fluxName + " " + nodeOut;
    s += " V='Lmag_" + windingIndex + "_L0/(1+pow(abs(I(" + senseName + "))/Lmag_" + windingIndex + "_Isat,2))*ddt(I(" + senseName + "))'\n";

    return s;
}

std::string CircuitSimulatorExporterNgspiceModel::export_magnetic_as_subcircuit(Magnetic magnetic, double frequency, double temperature, std::optional<std::string> filePathOrFile, CircuitSimulatorExporterCurveFittingModes mode) {
    std::string headerString = "* Magnetic model made with OpenMagnetics\n";
    headerString += "* " + magnetic.get_reference() + "\n\n";
    headerString += ".subckt " + fix_filename(magnetic.get_reference());
    std::string circuitString = "";
    std::string parametersString = "";
    std::string footerString = ".ends " + fix_filename(magnetic.get_reference());

    auto coil = magnetic.get_coil();

    double magnetizingInductance = resolve_dimensional_values(MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());

    // Resolve AUTO mode from settings before computing coefficients so we use the right model
    auto resolvedMode = CircuitSimulatorExporter::resolve_curve_fitting_mode(mode);
    auto acResistanceCoefficientsPerWinding = CircuitSimulatorExporter::calculate_ac_resistance_coefficients_per_winding(magnetic, temperature, resolvedMode);
    auto leakageInductances = LeakageInductance().calculate_leakage_inductance_all_windings(magnetic, Defaults().measurementFrequency).get_leakage_inductance_per_winding();

    std::vector<FractionalPoleNetwork> fracpoleNets;
    std::optional<FractionalPoleNetwork> coreFracNet;
    if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
        fracpoleNets = CircuitSimulatorExporter::calculate_fracpole_networks_per_winding(magnetic, temperature);
        try { coreFracNet = CircuitSimulatorExporter::calculate_core_fracpole_network(magnetic, temperature); }
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
        // Saturation only supported for single winding (inductor) - K coupling doesn't work with behavioral sources
        if (satParams.valid && coil.get_functional_description().size() > 1) {
            // For transformers, fall back to linear model with warning comment
            circuitString += "* WARNING: Saturation modeling not supported for multi-winding transformers in ngspice\n";
            circuitString += "* K coupling requires linear inductors. Using linear model.\n";
            includeSaturation = false;
        }
    }

    // Store coupling coefficients for each pair (indexed from winding 1)
    std::vector<double> couplingCoeffs;
    size_t numWindings = coil.get_functional_description().size();

    for (size_t index = 0; index < numWindings; index++) {
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
            couplingCoeffs.push_back(couplingCoefficient);
            parametersString += ".param Llk_" + is + "_Value=" + std::to_string(leakageInductance) + "\n";
            parametersString += ".param CouplingCoefficient_1" + is + "_Value=" + std::to_string(couplingCoefficient) + "\n";
        }

        std::vector<std::string> c = to_string(acResistanceCoefficientsPerWinding[index], 12);

        // FRACPOLE mode: emit fracpole-based winding subcircuit
        if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::FRACPOLE) {
            if (index < fracpoleNets.size()) {
                circuitString += emit_fracpole_winding_spice(fracpoleNets[index], index, is, effectiveResistanceThisWinding);
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            throw std::invalid_argument("Analytical mode not supported in NgSpice");
        }
        else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::ROSANO) {
            // Rosano winding model: Rdc + (R1||L1) + ... + (R6||L6) flat series
            const auto& rosanoCoeffs = acResistanceCoefficientsPerWinding[index];
            if (!rosanoCoeffs.empty()) {
                circuitString += emit_winding_rosano_spice(rosanoCoeffs, is, effectiveResistanceThisWinding);
            } else {
                // Fitting failed — fall back to DC resistance only
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }
        else {
            // LADDER mode (default)
            bool validLadderCoeffs = true;
            for (size_t ladderIndex = 0; ladderIndex < acResistanceCoefficientsPerWinding[index].size(); ladderIndex+=2) {
                double inductanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex + 1];
                double resistanceVal = acResistanceCoefficientsPerWinding[index][ladderIndex];
                if (inductanceVal > 0.1 || inductanceVal < 1e-7 || resistanceVal > 100.0 || resistanceVal < 1e-3) {
                    validLadderCoeffs = false;
                    break;
                }
            }

            if (validLadderCoeffs && acResistanceCoefficientsPerWinding[index].size() >= 2) {
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
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {NumberTurns_" + is + "**2*Permeance}\n";
            }
        }

        headerString += " P" + is + "+ P" + is + "-";
    }

    // Generate pairwise K statements for magnetic coupling
    // ngspice has a bug with 3+ inductor K statements inside subcircuits - the third inductor
    // doesn't get the proper hierarchical path prefix. Use pairwise K statements as workaround.
    // Each K statement gets a unique name (K12, K13, K23, etc.)
    // Use per-pair leakage inductance calculation for accurate coupling coefficients
    if (numWindings == 2) {
        // Simple 2-winding case - use already calculated coupling, capped at 0.98 for stability
        double k12 = couplingCoeffs.size() > 0 ? std::min(0.98, couplingCoeffs[0]) : 0.98;
        circuitString += "K Lmag_1 Lmag_2 " + std::to_string(k12) + "\n";
    } else if (numWindings >= 3) {
        // For 3+ windings, calculate coupling for each pair
        for (size_t i = 0; i < numWindings; i++) {
            for (size_t j = i + 1; j < numWindings; j++) {
                // Calculate leakage inductance between winding i and j
                double leakageIJ = 0.0;
                try {
                    auto leakageResult = LeakageInductance().calculate_leakage_inductance(
                        magnetic, Defaults().measurementFrequency, i, j);
                    auto leakagePerWinding = leakageResult.get_leakage_inductance_per_winding();
                    if (leakagePerWinding.size() > 0) {
                        leakageIJ = resolve_dimensional_values(leakagePerWinding[0]);
                    }
                } catch (...) {
                    // Fallback: use high coupling if calculation fails
                    leakageIJ = magnetizingInductance * 0.02;  // ~99% coupling
                }

                // Clamp leakage inductance to avoid negative or very low coupling
                if (leakageIJ >= magnetizingInductance) {
                    leakageIJ = magnetizingInductance * 0.1;  // Limit to ~95% coupling minimum
                }
                if (leakageIJ < 0) {
                    leakageIJ = magnetizingInductance * 0.02;  // ~99% coupling
                }

                double kij = sqrt((magnetizingInductance - leakageIJ) / magnetizingInductance);
                // Ensure coupling is in valid range
                // Cap at 0.98 for numerical stability (matches IDEAL mode default)
                kij = std::max(0.5, std::min(0.98, kij));

                std::string kName = "K" + std::to_string(i + 1) + std::to_string(j + 1);
                circuitString += kName + " Lmag_" + std::to_string(i + 1) + " Lmag_" + std::to_string(j + 1) + " " + std::to_string(kij) + "\n";
            }
        }
    }

    // Core losses: frequency-dependent resistance network in parallel with Lmag
    // The network impedance must be >> Lmag impedance, otherwise it shorts the magnetizing inductance
    auto coreLossTopology = static_cast<CoreLossTopology>(settings.get_circuit_simulator_core_loss_topology());
    auto coreResistanceCoefficients = CircuitSimulatorExporter::calculate_core_resistance_coefficients(magnetic, temperature, coreLossTopology);
    if (!coreResistanceCoefficients.empty()) {
        // Sanity check: core loss impedance at mid-frequency must be > 10x Lmag impedance
        // Otherwise the network shorts the magnetizing inductance instead of modeling losses
        double fMid = std::sqrt(1000.0 * 300000.0);
        double wMid = 2 * std::numbers::pi * fMid;
        double zLmag = wMid * magnetizingInductance;
        double zCoreLoss;
        if (coreLossTopology == CoreLossTopology::ROSANO) {
            zCoreLoss = CircuitSimulatorExporter::core_rosano_model(coreResistanceCoefficients.data(), fMid);
        } else {
            zCoreLoss = CircuitSimulatorExporter::core_ladder_model(coreResistanceCoefficients.data(), fMid, coreResistanceCoefficients[0]);
        }
        if (zCoreLoss > zLmag * 10) {
            if (coreLossTopology == CoreLossTopology::ROSANO) {
                circuitString += emit_core_rosano_spice(coreResistanceCoefficients, numWindings);
            } else {
                circuitString += emit_core_ladder_spice(coreResistanceCoefficients, numWindings);
            }
        }
        // If zCoreLoss <= 10*zLmag, skip core loss network to avoid shorting Lmag
    }

    // Mutual resistance: auxiliary winding network for cross-coupling losses (Hesterman 2020)
    // Only for multi-winding transformers when setting is enabled
    bool includeMutualResistance = settings.get_circuit_simulator_include_mutual_resistance();
    if (includeMutualResistance && numWindings >= 2) {
        auto mutualResistanceCoeffs = CircuitSimulatorExporter::calculate_mutual_resistance_coefficients(magnetic, temperature);
        if (!mutualResistanceCoeffs.empty()) {
            circuitString += emit_mutual_resistance_network_spice(mutualResistanceCoeffs, magnetizingInductance, numWindings);
        }
    }

    return headerString + "\n" + circuitString + "\n" + parametersString + "\n" + footerString;
}

} // namespace OpenMagnetics
