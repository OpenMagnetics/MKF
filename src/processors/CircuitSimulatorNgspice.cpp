#include "processors/CircuitSimulatorInterface.h"
#include "processors/CircuitSimulatorExporterHelpers.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/ExtendedCantilever.h"
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

    double Isat = sat.Isat;
    if (Isat < 1e-6) Isat = 1e-6;            // numerical floor

    s += "* Saturating magnetizing inductance (winding " + windingIndex + ")\n";
    s += "* Bsat=" + to_string(sat.Bsat, 4) + "T, Isat=" + to_string(Isat, 4) + "A\n";

    // Smooth-saturation behavioural inductor: L(I) = L0 / (1 + (I/Isat)^2), realized below in a numerically
    // stable flux-integral form (NOT V = L*dI/dt, which is ddt-unstable in ngspice). Needs only L0 and Isat. (The flux-linkage / air-gap quantities
    // the LTspice and NL5 exporters use — Lambda_sat = N*Ae*Bsat, L_gap = Lmag/effective_mu_r — belong to
    // THEIR tanh model; ngspice does not use them, so they are not computed here.)
    double L0 = sat.Lmag;

    s += ".param Lmag_" + windingIndex + "_L0=" + to_string(L0, 12) + "\n";
    s += ".param Lmag_" + windingIndex + "_Isat=" + to_string(Isat, 6) + "\n";

    // Flux-integral realization (numerically stable). The flux linkage of
    // L(I)=L0/(1+(I/Isat)^2) is lambda = integral_0^I L(i) di = L0*Isat*atan(I/Isat),
    // so the constitutive winding current is I(lambda) = Isat*tan(lambda/(L0*Isat)).
    // We drive the current from the INTEGRATED winding voltage (a flux state node)
    // instead of a behavioural V = L(I)*dI/dt source: a numerical current derivative
    // ddt(I) inside a B-source is numerically unstable in ngspice (Gear cannot start
    // it -> "timestep too small" at t~0; trapezoidal rings and runs away to nonsense
    // currents). Only the SPICE realization changes here; the physics (L0, Isat) are
    // the authoritative MKF values. Single-winding only (saturation is disabled for
    // coupled transformers above, whose K coupling needs linear inductors).
    std::string fluxName = "Lmag_flux_" + windingIndex;
    std::string L0p = "Lmag_" + windingIndex + "_L0";
    std::string Isp = "Lmag_" + windingIndex + "_Isat";

    // lambda = integral(V(nodeIn,nodeOut)) dt : a VCCS of the winding voltage into a
    // 1 F capacitor makes V(fluxName) the time-integral of the winding voltage.
    s += "Gflux_" + windingIndex + " 0 " + fluxName + " " + nodeIn + " " + nodeOut + " 1\n";
    s += "Cflux_" + windingIndex + " " + fluxName + " 0 1\n";
    s += "Rflux_" + windingIndex + " " + fluxName + " 0 1e12\n";  // DC leak: keep the flux node non-floating

    // winding current I(lambda) = Isat * tan(lambda / (L0*Isat))
    s += "Bind_" + windingIndex + " " + nodeIn + " " + nodeOut +
         " I='" + Isp + "*tan(V(" + fluxName + ")/(" + L0p + "*" + Isp + "))'\n";

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

    // For 3+ windings the historical pairwise-K approach computed each coupling from an
    // independently-measured leakage and then CAPPED it at 0.98 — a fallback that is neither
    // self-consistent nor physical. Instead build the full inductance matrix L = M + Λ once
    // (Erickson/Maksimovic extended-cantilever foundation); it is positive-definite, so each
    // self-inductance L_ii and coupling k_ij = L_ij/sqrt(L_ii·L_jj) is consistent and |k_ij| < 1
    // with no cap. The 2-winding path is left on its validated formula.
    std::vector<std::vector<double>> inductanceMatrix;
    if (numWindings >= 3) {
        inductanceMatrix = ExtendedCantilever::calculate_inductance_matrix(magnetic, Defaults().measurementFrequency);
    }
    auto magnetizingInductorValue = [&](size_t windingIdx, const std::string& is) -> std::string {
        if (numWindings >= 3) {
            return to_string(inductanceMatrix[windingIdx][windingIdx], 12);  // full self-inductance L_ii
        }
        return "NumberTurns_" + is + "**2*Permeance";
    };

    for (size_t index = 0; index < numWindings; index++) {
        auto effectiveResistanceThisWinding = WindingLosses::calculate_effective_resistance_of_winding(magnetic, index, 0.1, temperature);
        std::string is = std::to_string(index + 1);
        parametersString += ".param Rdc_" + is + "_Value=" + std::to_string(effectiveResistanceThisWinding) + "\n";
        parametersString += ".param NumberTurns_" + is + "=" + std::to_string(coil.get_functional_description()[index].get_number_turns()) + "\n";
        if (index > 0) {
            double leakageInductance = resolve_dimensional_values(leakageInductances[index]);
            if (leakageInductance < 0 || leakageInductance >= magnetizingInductance) {
                throw std::runtime_error("Unphysical leakage inductance (" + std::to_string(leakageInductance) + " H vs Lmag " + std::to_string(magnetizingInductance) + " H) for winding " + is);
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
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {" + magnetizingInductorValue(index, is) + "}\n";
            }
        }
        else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::ANALYTICAL) {
            throw std::invalid_argument("Analytical mode not supported in NgSpice");
        }
        else if (resolvedMode == CircuitSimulatorExporterCurveFittingModes::ROSANO ||
                 resolvedMode == CircuitSimulatorExporterCurveFittingModes::ROSANO_RLC) {
            // Rosano winding model: Rdc + (R1||L1) + ... + (Rn||Ln) flat series.
            // ROSANO_RLC additionally allows the last stage to be R||(L+C
            // series). The fitter encodes which topology was selected via a
            // trailing marker flag in the coefficient vector;
            // emit_winding_rosano_spice decodes and emits the right network.
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
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {" + magnetizingInductorValue(index, is) + "}\n";
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
                // Fitted model (ladder_model): Z = Rdc + L1||(R1 + L2||(R2 + ...))
                // Rdc in series first, then nested L||(R + rest) stages.
                size_t numStages = acResistanceCoefficientsPerWinding[index].size() / 2;
                circuitString += "Rdc" + is + " P" + is + "+ Node_Lladder_" + is + "_0 {Rdc_" + is + "_Value}\n";
                for (size_t stage = 0; stage < numStages; ++stage) {
                    std::string stageNode = "Node_Lladder_" + is + "_" + std::to_string(stage);
                    std::string nextNode = (stage + 1 == numStages) ? ("Node_R_Lmag_" + is)
                                                                    : ("Node_Lladder_" + is + "_" + std::to_string(stage + 1));
                    // L_stage in parallel with (R_stage + rest of the ladder)
                    circuitString += "Lladder" + is + "_" + std::to_string(stage) + " " + stageNode + " Node_R_Lmag_" + is + " " + c[stage * 2 + 1] + "\n";
                    circuitString += "Rladder" + is + "_" + std::to_string(stage) + " " + stageNode + " " + nextNode + " " + c[stage * 2] + "\n";
                }
            } else {
                circuitString += "Rdc" + is + " P" + is + "+ Node_R_Lmag_" + is + " {Rdc_" + is + "_Value}\n";
            }
            // Emit magnetizing inductance (saturating or linear)
            if (includeSaturation && satParams.valid) {
                circuitString += emit_saturating_inductor_ngspice(satParams, is, "Node_R_Lmag_" + is, "P" + is + "-");
            } else {
                circuitString += "Lmag_" + is + " Node_R_Lmag_" + is + " P" + is + "- {" + magnetizingInductorValue(index, is) + "}\n";
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
        // Consistent coupling from the full inductance matrix L = M + Λ (positive-definite):
        //   k_ij = L_ij / sqrt(L_ii · L_jj),  with the inductors emitted as L_ii above.
        // This is self-consistent across all pairs and yields |k_ij| < 1 inherently, so no
        // stability cap is applied. (ngspice still needs pairwise K statements rather than a
        // single coupling matrix because of its 3+-inductor-in-subckt path-prefix bug.)
        for (size_t i = 0; i < numWindings; i++) {
            for (size_t j = i + 1; j < numWindings; j++) {
                double Lii = inductanceMatrix[i][i];
                double Ljj = inductanceMatrix[j][j];
                double Mij = inductanceMatrix[i][j];
                if (Lii <= 0 || Ljj <= 0) {
                    throw std::runtime_error("Non-positive self-inductance building coupling for windings " + std::to_string(i) + "-" + std::to_string(j));
                }
                double kij = Mij / std::sqrt(Lii * Ljj);
                if (std::abs(kij) >= 1.0) {
                    throw std::runtime_error("Inconsistent coupling coefficient (|k|=" + std::to_string(std::abs(kij)) + " >= 1) for windings " + std::to_string(i) + "-" + std::to_string(j) + "; inductance matrix is not positive-definite");
                }
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
            // Odd size => trailing element is the measured DC core resistance R0
            double dcCoreResistance = (coreResistanceCoefficients.size() % 2 == 1) ? coreResistanceCoefficients.back() : 0.0;
            zCoreLoss = CircuitSimulatorExporter::core_ladder_model(coreResistanceCoefficients.data(), fMid, dcCoreResistance);
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
