#pragma once
#include <MAS.hpp>
#include "support/Utils.h"
#include "Defaults.h"
#include "physical_models/Impedance.h"
#include "physical_models/WindingLosses.h"
#include "constructive_models/Magnetic.h"
#include <complex>

using namespace MAS;

namespace OpenMagnetics {

class Sweeper {
    private:
    protected:
    public:
    // fast=true (default): the series leakage term uses DC + skin-effect winding
    // resistance (analytic, cheap). fast=false additionally includes proximity
    // effect (field-based, slower) for a more accurate leakage-resonance damping.
    // fastCapacitance selects the stray-capacitance model, and the caller must choose it --
    // there is no silent fallback between the two. false (default) keeps the full per-turn
    // energy sum, which is the accurate answer for a design being ANALYSED. true selects the
    // Massarini-based one-layer model, which is what the advisers already use to RANK
    // candidates, and it is the only option available when the coil cannot be wound at all:
    // the full model winds the coil to sum turn energies and throws when the winder returns
    // no turns (ABT #850), so a part whose wire does not fit its bore has no full-model
    // capacitance to compute. Its impedance is still meaningful -- for a bead |Z| is set by
    // the core, and the capacitance only moves the self-resonance -- but the caller has to
    // say it is accepting the one-layer estimate.
    static Curve2D sweep_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Impedance over frequency", bool fast = true, bool fastCapacitance = false);
    // Common-mode impedance (all windings driven in parallel — the CMC datasheet /
    // REDEXPERT CM measurement): the magnetizing tank alone, no leakage resonance.
    // Use this, not sweep_impedance_over_frequency, for a choke's CM curve (ABT #167).
    // fastCapacitance as above: the common-mode model shunts the magnetizing tank with the
    // winding self-capacitance, which the full model can only get from wound turns.
    static Curve2D sweep_common_mode_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Common-mode impedance over frequency", bool fastCapacitance = false);
    static Curve2D sweep_differential_mode_impedance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Differential-mode impedance over frequency");
    static Curve2D sweep_q_factor_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, std::string mode="log", std::string title = "Impedance over frequency");
    static Curve2D sweep_magnetizing_inductance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Magnetizing Inductance over frequency");
    static Curve2D sweep_magnetizing_inductance_over_temperature(Magnetic magnetic, double start, double stop, size_t numberElements, double frequency = defaults.measurementFrequency, std::string mode="linear", std::string title = "Magnetizing Inductance over temperature");
    static Curve2D sweep_magnetizing_inductance_over_dc_bias(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="linear", std::string title = "Magnetizing Inductance over DC bias");
    static Curve2D sweep_winding_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, size_t windingIndex, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Winding Resistance over frequency");
    static Curve2D sweep_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Resistance over frequency");
    static Curve2D sweep_core_resistance_over_frequency(Magnetic magnetic, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Core Resistance over frequency");
    static Curve2D sweep_core_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Core Losses over frequency");
    static Curve2D sweep_winding_losses_over_frequency(Magnetic magnetic, OperatingPoint operatingPoint, double start, double stop, size_t numberElements, double temperature = defaults.ambientTemperature, std::string mode="log", std::string title = "Winding Losses over frequency");
    
    /**
     * @brief Sweep the full resistance matrix (including mutual resistance) over frequency.
     * 
     * Returns a vector of ScalarMatrixAtFrequency objects, one per frequency point.
     * Each matrix contains:
     * - Diagonal elements R_ii: self-resistance of winding i
     * - Off-diagonal elements R_ij: mutual resistance between windings i and j
     * 
     * Based on Hesterman (2020) "Mutual Resistance" and Spreen (1990).
     * 
     * @param magnetic The magnetic component
     * @param start Starting frequency in Hz
     * @param stop Ending frequency in Hz
     * @param numberElements Number of frequency points
     * @param temperature Temperature in Celsius
     * @param mode "log" or "linear" spacing
     * @return Vector of resistance matrices at each frequency
     */
    static std::vector<ScalarMatrixAtFrequency> sweep_resistance_matrix_over_frequency(
        Magnetic magnetic, 
        double start, 
        double stop, 
        size_t numberElements, 
        double temperature = defaults.ambientTemperature, 
        std::string mode = "log");

};

} // namespace OpenMagnetics