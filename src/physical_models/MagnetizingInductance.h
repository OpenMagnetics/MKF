#pragma once
#include "Constants.h"
#include "Defaults.h"
#include "support/Utils.h"

#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include <MAS.hpp>
#include "constructive_models/Coil.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <streambuf>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

class MagnetizingInductance {
  private:
    std::map<std::string, std::string> _models;

  protected:
  public:
    MagnetizingInductance() {
        _models["gapReluctance"] = to_string(Defaults().reluctanceModelDefault);
    }

    MagnetizingInductance(ReluctanceModels model) {
        _models["gapReluctance"] = to_string(model);
    }

    MagnetizingInductance(std::string model) {
        _models["gapReluctance"] = model;
    }
    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(Core core,
                                                                                   Coil coil,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    MagnetizingInductanceOutput calculate_inductance_from_number_turns_and_gapping(Magnetic magnetic,
                                                                                   OperatingPoint* operatingPoint = nullptr);

    std::vector<CoreGap> calculate_gapping_from_number_turns_and_inductance(Core core,
                                                                      Coil coil,
                                                                      Inputs* inputs,
                                                                      GappingType gappingType,
                                                                      size_t decimals = 4);

    int calculate_number_turns_from_gapping_and_inductance(Core core, Inputs* inputs, DimensionalValues preferredValue = DimensionalValues::NOMINAL);

    /**
     * @brief Calculate gap length to prevent saturation given current and target B-field
     * 
     * Iteratively determines the gap length required to keep B-field below target
     * for given operating current. Similar structure to calculate_number_turns_from_gapping_and_inductance.
     * 
     * @param core The core to gap
     * @param inputs Design inputs with operating points
     * @param targetMagneticFluxDensity Maximum allowed B-field (e.g., 0.9 * Bsat)
     * @param magnetizingCurrentPeak Peak magnetizing current from operating point
     * @return double The calculated gap length in meters
     */
    double calculate_gap_from_saturation_constraint(Core core,
                                                   Inputs* inputs,
                                                   double targetMagneticFluxDensity,
                                                   double magnetizingCurrentPeak);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        Core core,
        Coil coil,
        OperatingPoint* operatingPoint = nullptr);

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> calculate_inductance_and_magnetic_flux_density(
        Magnetic magnetic,
        OperatingPoint* operatingPoint = nullptr);

    double calculate_inductance_air_solenoid(
        Magnetic magnetic);

    double calculate_inductance_air_solenoid(
        Core core,
        Coil coil);

    /**
     * @brief Calculate optimal gap and turns simultaneously (Option 1: Direct Analytical)
     * 
     * Solves the coupled gap/turns equations analytically to find the optimal
     * combination that meets inductance and saturation constraints.
     * 
     * The key insight is:
     *   - From saturation: N_min = (L * I_peak) / (B_max * A_core)
     *   - From inductance: N² = L * R_total(gap)
     * 
     * These are solved together to find (gap, turns) in one pass.
     * 
     * @param core The core (not modified)
     * @param inputs Design inputs with inductance requirement and operating points
     * @param maxFluxDensity Maximum allowed B-field (typically 0.7 * Bsat)
     * @param peakCurrent Peak magnetizing current
     * @return std::pair<double, double> {gap_length, number_turns}, or {-1, -1} if no valid solution
     */
    std::pair<double, double> calculate_optimal_gap_and_turns(
        const Core& core,
        Inputs* inputs,
        double maxFluxDensity,
        double peakCurrent);

    /**
     * @brief Calculate turns required for a given gap to achieve target inductance
     * 
     * Simple analytical calculation: N = sqrt(L * R_total)
     * Accounts for permeability variation with DC bias if operating point is provided.
     * 
     * @param core Core with gap already set
     * @param targetInductance Target magnetizing inductance in Henries
     * @param temperature Operating temperature for permeability calculation
     * @param frequency Operating frequency for permeability calculation
     * @return double Number of turns (rounded to nearest integer)
     */
    double calculate_turns_for_gap(
        Core& core,
        double targetInductance,
        double temperature = 25.0,
        double frequency = 100000.0);

    /**
     * @brief Calculate peak flux density for given turns, current, and core
     * 
     * B_peak = (N * I_peak) / (R_total * A_effective)
     *        = L * I_peak / (N * A_effective)
     * 
     * @param core Core with gap set
     * @param numberTurns Number of primary turns
     * @param peakCurrent Peak magnetizing current
     * @param temperature Operating temperature
     * @param frequency Operating frequency
     * @return double Peak magnetic flux density in Tesla
     */
    double calculate_flux_density_peak(
        Core& core,
        double numberTurns,
        double peakCurrent,
        double temperature = 25.0,
        double frequency = 100000.0);

    /**
     * @brief Calculate peak flux density from voltage using Faraday's Law
     * 
     * For transformers where B is determined by applied voltage, not current.
     * B_peak = V_peak / (N * Ae * ω)
     * 
     * This method does NOT iterate on permeability since transformers
     * ideally have minimal magnetizing current and no DC bias.
     * 
     * @param core Core with effective area
     * @param numberTurns Number of primary turns
     * @param voltagePeak Peak voltage applied to primary
     * @param frequency Operating frequency
     * @return double Peak magnetic flux density in Tesla
     */
    double calculate_flux_density_peak_from_voltage(
        Core& core,
        double numberTurns,
        double voltagePeak,
        double frequency = 100000.0);

    /**
     * @brief Calculate minimum turns to avoid saturation from voltage
     * 
     * For transformers, calculates minimum turns needed to keep B below max.
     * N_min = V_peak / (Ae * ω * B_max)
     * 
     * @param core Core with effective area
     * @param voltagePeak Peak voltage applied to primary
     * @param frequency Operating frequency
     * @param maxFluxDensity Maximum allowed flux density (typically 0.7 * Bsat)
     * @return double Minimum number of turns (rounded up)
     */
    double calculate_turns_from_voltage_and_max_flux_density(
        Core& core,
        double voltagePeak,
        double frequency,
        double maxFluxDensity);
};

} // namespace OpenMagnetics