#include "advisers/CoreAdviser.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // is_energy_storing_topology
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/Impedance.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace OpenMagnetics {

void add_gapping(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    MagneticEnergy magneticEnergy;
    
    // Gapping Decision Logic (same as add_gapping_standard_cores):
    // - Skip gapping only if we have minimum-only inductance (want maximum L, no specific target)
    // - Energy-storing topologies or specific inductance targets need gap calculation
    auto inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
    bool hasNominalInductance = inductanceReq.get_nominal().has_value();
    bool hasMaxInductance = inductanceReq.get_maximum().has_value();
    bool needsSpecificInductance = hasNominalInductance || hasMaxInductance;
    
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    bool isTransformer = topology.has_value() ? !isEnergyStoring : 
        (inductanceReq.get_minimum() && !hasNominalInductance && !hasMaxInductance);
    
    bool skipGapping = isTransformer && !needsSpecificInductance;
    
    if (skipGapping)
    {
        for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
            Core core = (*magneticsWithScoring)[i].first.get_core();
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
        }

        return;
    }
    auto requiredMagneticEnergy = resolve_dimensional_values(magneticEnergy.calculate_required_magnetic_energy(inputs), DimensionalValues::MAXIMUM);
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        Core core = (*magneticsWithScoring)[i].first.get_core();

        if (core.get_material_name() == DUMMY_SENTINEL_NAME) {
            core.set_material_initial_permeability(defaults.ferriteInitialPermeability);
        }
        if (!core.get_processed_description()) {
            core.process_data();
        }
        if (core.get_shape_family() != CoreShapeFamily::T) {
            // Use realistic ferrite Bsat (~350 mT) instead of dummy material's 500 mT
// CA-LOGIC-3 FIX: Use actual material Bsat instead of hardcoded ferrite default
            // Phase 1 fix: split named-sentinel (Dummy material from
            // shape-only pre-filter stage) from real-material lookup. Previously
            // both shared a silent `catch (...)` that swallowed real-material
            // failures and contaminated the gap calc with a generic ferrite
            // default Bsat — letting unsuitable candidates pass.
            double realisticBsat;
            if (core.get_material_name() == DUMMY_SENTINEL_NAME) {
                // Shape-only pre-filter stage (see CoreAdviser.cpp:2344 fan-out):
                // material isn't bound yet, use the ferrite reference Bsat as
                // the documented contract of this stage.
                realisticBsat = defaults.ferriteSaturationFluxDensity;
            }
            else {
                try {
                    auto coreMat = core.resolve_material();
                    double opTemp = 25.0;
                    for (auto& op : inputs.get_operating_points()) {
                        opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
                    }
                    realisticBsat = Core::get_magnetic_flux_density_saturation(coreMat, opTemp);
                }
                catch (const std::exception& e) {
                    // Real material's Bsat lookup failed → cannot compute a
                    // valid gap for this candidate. Skip the gap mutation; the
                    // core stays in its incoming (likely ungapped) state and
                    // downstream Bsat / energy filters reject it on its own
                    // merits. Logged so the underlying material-data gap is
                    // visible (loud).
                    logEntry(std::string("Bsat lookup failed for material '")
                             + core.get_material_name() + "' on core "
                             + core.get_name().value_or("unnamed") + ": "
                             + e.what() + " — leaving core ungapped (downstream filters will reject)",
                             "CoreAdviser");
                    continue;
                }
            }
            double bSatTarget = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation; // Use configured proportion of Bsat for safety margin
            
            // Calculate gap based on energy storage requirement
            double gapLength = roundFloat(magneticEnergy.calculate_gap_length_by_magnetic_energy(core.get_gapping()[0], realisticBsat, requiredMagneticEnergy), 5);
            double gapEnergy = gapLength;
            
            // Iterative refinement to find gap that prevents saturation
            double magnetizingCurrentPeak = 0;
            for (auto& op : inputs.get_operating_points()) {
                auto excitation = Inputs::get_primary_excitation(op);
                if (excitation.get_magnetizing_current() && excitation.get_magnetizing_current()->get_processed() && excitation.get_magnetizing_current()->get_processed()->get_peak()) {
                    magnetizingCurrentPeak = std::max(magnetizingCurrentPeak, excitation.get_magnetizing_current()->get_processed()->get_peak().value());
                }
            }
            
            // Get the actual current from primary excitation (this is what saturation filter will use)
            double actualCurrentPeak = 0;
            for (auto& op : inputs.get_operating_points()) {
                auto excitation = Inputs::get_primary_excitation(op);
                if (excitation.get_current() && excitation.get_current()->get_processed() && excitation.get_current()->get_processed()->get_peak()) {
                    actualCurrentPeak = std::max(actualCurrentPeak, excitation.get_current()->get_processed()->get_peak().value());
                }
            }
            
            // Use the larger of magnetizing current or actual current
            double currentForCalculation = std::max(magnetizingCurrentPeak, actualCurrentPeak);
            
            if (currentForCalculation > 0) {
                // Use the proper method to calculate gap from saturation constraint
                MagnetizingInductance magnetizingInductance;
                double gapFromSaturation = magnetizingInductance.calculate_gap_from_saturation_constraint(
                    core, &inputs, bSatTarget, currentForCalculation);
                
                // Use the maximum of energy gap and saturation gap
                gapLength = std::max(gapEnergy, gapFromSaturation);
            }

            core.set_ground_gapping(gapLength);
            core.process_gap();
            std::stringstream ss;
            ss << roundFloat(gapLength * 1000, 2);
            if (gapLength > 0) {
                core.set_name(core.get_name().value_or("unnamed") + " gapped " + ss.str()  + " mm");
            }
            else {
                core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            }
        }

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

// ============================================================================
// New Gapping Algorithm for STANDARD_CORES Mode
// ============================================================================

CoreAdviser::GappingConstraints CoreAdviser::calculate_gapping_constraints(Inputs inputs, Core core) {
    GappingConstraints constraints;

    // Get required magnetic energy
    MagneticEnergy magneticEnergy;
    auto requiredMagneticEnergy = resolve_dimensional_values(
        magneticEnergy.calculate_required_magnetic_energy(inputs),
        DimensionalValues::MAXIMUM);

    // Account for stacked cores - each stack shares the energy
    auto numStacksOpt = core.get_functional_description().get_number_stacks();
    int64_t numStacks = numStacksOpt ? numStacksOpt.value() : 1;
    if (numStacks > 1) {
        requiredMagneticEnergy /= numStacks;
    }

    // Get realistic Bsat for the material at operating temperature
    double opTemp = 25.0;
    for (auto& op : inputs.get_operating_points()) {
        opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
    }
    double realisticBsat = core.get_magnetic_flux_density_saturation(opTemp, false);

    // 1. Calculate minimum gap: energy storage requirement
    // Use the maximum allowed flux density (proportional to Bsat) for energy calculation
    double maxAllowedB = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double minGap = magneticEnergy.calculate_gap_length_by_magnetic_energy(
        core.get_gapping()[0], maxAllowedB, requiredMagneticEnergy);
    constraints.minGap = minGap;

    // 2. Calculate maximum gap constraints
    // 2a. Maximum gap based on column width (30%)
    double columnWidth = core.get_columns()[0].get_width();
    double maxGapByColumn = 0.30 * columnWidth;

    // 2b. Maximum gap based on fringing factor (25%)
    double maxGapByFringing = calculate_gap_for_fringing_factor(defaults.coreAdviserMaxFringingFactor, core);

    // Use the more restrictive limit
    constraints.maxGap = std::min(maxGapByColumn, maxGapByFringing);

    // Ensure maxGap >= minGap (sanity check)
    if (constraints.maxGap < constraints.minGap) {
        constraints.maxGap = constraints.minGap;
    }

    // 3. Calculate saturation constraint gap (ensure we don't exceed 70% Bsat)
    MagnetizingInductance magnetizingInductance;
    double targetB = realisticBsat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double peakCurrent = get_peak_current(inputs);
    
    // Account for stacked cores - current is divided among stacks
    if (numStacks > 1) {
        peakCurrent /= numStacks;
    }
    
    double saturationGap = magnetizingInductance.calculate_gap_from_saturation_constraint(
        core, &inputs, targetB, peakCurrent);

    // The minimum gap must satisfy BOTH energy storage AND saturation constraints
    double effectiveMinGap = std::max(constraints.minGap, saturationGap);
    
    // 4. Check if we can satisfy both saturation and fringing constraints
    if (effectiveMinGap > constraints.maxGap) {
        // Cannot satisfy both - prioritize saturation constraint
        // This core will have higher fringing but won't saturate
        constraints.maxGap = effectiveMinGap * 1.5;  // Allow some margin for optimization
    }
    
    // 5. Find optimal gap
    if (settings.get_gapping_strategy() == GappingOptimizationStrategy::GOLDEN_SECTION) {
        try {
            constraints.optimalGap = optimize_gap_golden_section(
                effectiveMinGap, constraints.maxGap, inputs, core);
        }
        catch (const std::exception& e) {
            // Phase 1 fix: golden-section needs computable losses across the
            // bracket. If the loss model fails for this core (e.g. missing
            // Steinmetz coefficients at the operating frequency), fall back
            // explicitly to the SIMPLE strategy rather than letting the
            // optimizer work on DBL_MAX poison values (which it used to
            // silently receive from calculate_core_losses_for_gap). This
            // fallback is named, logged, and bounded — not a hidden default.
            logEntry(std::string("Golden-section gap optimization failed (")
                     + e.what() + "), falling back to SIMPLE strategy for core "
                     + core.get_name().value_or("unnamed"), "CoreAdviser");
            constraints.optimalGap = effectiveMinGap;
        }
    } else {
        // SIMPLE strategy: use effective minimum gap
        constraints.optimalGap = effectiveMinGap;
    }

    return constraints;
}

double CoreAdviser::calculate_gap_for_fringing_factor(double targetFringingFactor, Core core) {
    // Use the reluctance model to find gap length for target fringing factor
    auto reluctanceModel = ReluctanceModel::factory(_models);
    return reluctanceModel->get_gapping_by_fringing_factor(core, targetFringingFactor);
}

double CoreAdviser::get_peak_current(Inputs inputs) {
    double peakCurrent = 0.0;

    // For transformer topologies (forward converters), core saturation is driven by
    // magnetizing current only. The reflected secondary current is balanced and does
    // not contribute to net flux. Using actual current would oversize the core.
    // When topology is unset, fall back to inductance-based heuristic:
    // minimum-only inductance = transformer; nominal/max inductance = inductor.
    auto topology = inputs.get_design_requirements().get_topology();
    bool isTransformerTopology;
    if (topology.has_value()) {
        isTransformerTopology = !is_energy_storing_topology(topology);
    } else {
        auto& inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
        isTransformerTopology = inductanceReq.get_minimum() &&
                                !inductanceReq.get_nominal() &&
                                !inductanceReq.get_maximum();
    }

    for (auto& op : inputs.get_operating_points()) {
        auto excitation = Inputs::get_primary_excitation(op);

        // Check magnetizing current
        if (excitation.get_magnetizing_current() &&
            excitation.get_magnetizing_current()->get_processed() &&
            excitation.get_magnetizing_current()->get_processed()->get_peak()) {
            peakCurrent = std::max(peakCurrent,
                excitation.get_magnetizing_current()->get_processed()->get_peak().value());
        }

        // Check actual current only for energy-storing topologies (e.g. flyback, boost)
        // and inductors (identified by having a nominal inductance).
        // For pure transformer topologies the actual current includes reflected secondary
        // load current which does not bias the core.
        if (!isTransformerTopology &&
            excitation.get_current() &&
            excitation.get_current()->get_processed() &&
            excitation.get_current()->get_processed()->get_peak()) {
            peakCurrent = std::max(peakCurrent,
                excitation.get_current()->get_processed()->get_peak().value());
        }
    }

    return peakCurrent;
}

double CoreAdviser::calculate_core_losses_for_gap(double gap, Inputs inputs, Core core) {
    // Phase 1 fix: previously wrapped in `try { ... } catch (const std::exception&)
    // { return DBL_MAX; }`, which silently turned uncomputable losses into a
    // "very high loss" signal that contaminated the golden-section optimizer
    // with no way for it to know the optimum was meaningless. Per the
    // no-silent-fallbacks policy, let the exception propagate; the single
    // caller in add_gapping_standard_cores wraps the optimizer call and
    // falls back to the SIMPLE gap strategy with an explicit logEntry on
    // failure (named fallback, not a hidden default).
    core.set_ground_gapping(gap);
    core.process_gap();

    auto coreLossesModel = CoreLossesModel::factory(
        std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));

    double totalLosses = 0.0;
    for (auto& op : inputs.get_operating_points()) {
        auto excitation = Inputs::get_primary_excitation(op);
        double temperature = op.get_conditions().get_ambient_temperature();
        auto losses = coreLossesModel->get_core_losses(core, excitation, temperature);
        totalLosses += losses.get_core_losses();
    }

    return totalLosses;
}

double CoreAdviser::optimize_gap_golden_section(double minGap, double maxGap, Inputs inputs, Core core) {
    const double PHI = 1.618033988749895;  // Golden ratio
    const double RESPHI = 2.0 - PHI;       // 1/phi^2
    const int MAX_ITERATIONS = 10;

    double a = minGap;
    double b = maxGap;

    // Initial points
    double c = a + RESPHI * (b - a);
    double d = b - RESPHI * (b - a);

    double fc = calculate_core_losses_for_gap(c, inputs, core);
    double fd = calculate_core_losses_for_gap(d, inputs, core);

    for (int i = 0; i < MAX_ITERATIONS && (b - a) > 1e-6; ++i) {
        if (fc < fd) {
            // Minimum is in [a, d]
            b = d;
            d = c;
            fd = fc;
            c = a + RESPHI * (b - a);
            fc = calculate_core_losses_for_gap(c, inputs, core);
        } else {
            // Minimum is in [c, b]
            a = c;
            c = d;
            fc = fd;
            d = b - RESPHI * (b - a);
            fd = calculate_core_losses_for_gap(d, inputs, core);
        }
    }

    // Return the midpoint of the final interval
    return (a + b) / 2.0;
}

void CoreAdviser::add_gapping_standard_cores(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                              Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }

    // Gapping Decision Logic:
    //
    // The decision to gap depends on TWO factors:
    // 1. Whether the topology stores energy (needs gap for DC bias)
    // 2. Whether a specific inductance value is required (needs gap to control AL)
    //
    // Cases:
    // - Energy-storing (Flyback, Buck, etc.): ALWAYS needs gap for energy storage
    // - Transformer with nominal/max inductance: MAY need gap to hit specific L value
    //   (e.g., LLC resonant converter needs controlled magnetizing inductance)
    // - Transformer with minimum-only inductance: NO gap needed (want maximum L)
    //
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    
    auto inductanceReq = inputs.get_design_requirements().get_magnetizing_inductance();
    bool hasNominalInductance = inductanceReq.get_nominal().has_value();
    bool hasMaxInductance = inductanceReq.get_maximum().has_value();

    // Transformer topologies skip gapping UNLESS a specific Lm value is required.
    // When Lm is specified (nominal or maximum), the transformer core may need a gap:
    // - If AL_ungapped > Lm (e.g. CLLLC Lm=1.32mH << ferrite AL), a gap reduces AL
    //   to achieve the target Lm at the volt-seconds-derived turn count.
    // - The gap is computed directly: gap = (N² / Lm - R_core) × μ₀ × Ae,
    //   where N comes from the same volt-seconds estimate as add_initial_turns_by_inductance.
    // When only a minimum inductance is specified (want high Lm), no gap is needed.
    bool isTransformer = topology.has_value() ? !isEnergyStoring :
        (inductanceReq.get_minimum() && !hasNominalInductance && !hasMaxInductance);

    bool skipGapping = isTransformer && !hasNominalInductance && !hasMaxInductance;

    // Two transformer paths share the same outcome HERE — leave the core
    // ungapped — but for opposite reasons:
    //
    //   1. Pure transformer (no specific Lm target): we want Lm
    //      maximised, so no gap, ever.
    //   2. Transformer WITH specific Lm: the proper gap depends on the
    //      core's real material (μᵢ) and the final turn count N. Both
    //      are unknown at this stage — `add_ferrite_materials_by_losses`
    //      assigns the material later, and `add_initial_turns_by_inductance`
    //      sets N after that. So we defer the gap to that stage, where
    //      it can call `MagnetizingInductance::calculate_gapping_from_
    //      number_turns_and_inductance` against the real reluctance
    //      model rather than the hand-rolled
    //      `gap = (N²/Lm − R_core) · μ₀ · Aₑ` approximation an earlier
    //      draft used here.
    if (isTransformer) {
        for (size_t i = 0; i < magneticsWithScoring->size(); ++i) {
            Core core = (*magneticsWithScoring)[i].first.get_core();
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
        }
        return;
    }

    for (size_t i = 0; i < magneticsWithScoring->size(); ++i) {
        Core core = (*magneticsWithScoring)[i].first.get_core();

        if (core.get_shape_family() == CoreShapeFamily::T) {
            // Toroidal cores don't have gaps
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
            continue;
        }

        // Process core data if needed
        if (!core.get_processed_description()) {
            core.process_data();
        }

        // Calculate gapping constraints
        auto constraints = calculate_gapping_constraints(inputs, core);

        // Apply the optimal gap
        core.set_ground_gapping(constraints.optimalGap);
        core.process_gap();

        // Update name with gap info (avoid duplicates)
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << (constraints.optimalGap * 1000);
        auto currentName = core.get_name().value_or("unnamed");
        // Remove existing gap info if present to avoid duplicates
        size_t gappedPos = currentName.find(" gapped ");
        if (gappedPos != std::string::npos) {
            currentName = currentName.substr(0, gappedPos);
        }
        size_t ungappedPos = currentName.find(" ungapped");
        if (ungappedPos != std::string::npos) {
            currentName = currentName.substr(0, ungappedPos);
        }
        if (constraints.optimalGap > 0) {
            core.set_name(currentName + " gapped " + ss.str() + " mm");
        } else {
            core.set_name(currentName + " ungapped");
        }

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}

void CoreAdviser::refine_gaps_for_saturation(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                              Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }

    MagnetizingInductance magnetizingInductance;
    MagneticSimulator magneticSimulator;
    const int MAX_ITERATIONS = 5;  // Reduced from 10 to improve performance

    for (size_t i = 0; i < magneticsWithScoring->size();) {
        Core core = (*magneticsWithScoring)[i].first.get_core();
        
        if (core.get_shape_family() == CoreShapeFamily::T) {
            // Toroidal cores don't have gaps to refine
            ++i;
            continue;
        }

        // Skip powder cores - they have distributed gaps and are designed for high DC current
        // Powder cores are identified by their material type
        auto materialName = core.get_material_name();
        std::string mat = materialName;
        if (mat.find("Kool M") != std::string::npos || 
            mat.find("High Flux") != std::string::npos ||
            mat.find("XFlux") != std::string::npos ||
            mat.find("Edge") != std::string::npos ||
            mat.find("Fluxsan") != std::string::npos ||
            mat.find("FS ") == 0 ||  // Fluxsan prefix
            mat.find("HF ") == 0 ||  // High Flux prefix
            mat.find("GX ") == 0) {  // Another powder series
            // Powder core - skip gap refinement, keep as-is
            ++i;
            continue;
        }

        // Target safe operating flux for this core (material Bsat curve at opTemp,
        // capped by defaults.maximumProportionMagneticFluxDensitySaturation).
        double opTemp = 25.0;
        for (auto& op : inputs.get_operating_points()) {
            opTemp = std::max(opTemp, op.get_conditions().get_ambient_temperature());
        }
        double maxB = core.get_magnetic_flux_density_saturation(opTemp, /*proportion=*/true);

        // Get current gap
        double currentGap = 0;
        if (!core.get_functional_description().get_gapping().empty()) {
            currentGap = core.get_functional_description().get_gapping()[0].get_length();
        }

        // Iterative refinement to achieve target saturation
        int iteration = 0;
        double bPeak = 0;
        bool converged = false;
        bool shouldRemove = false;

        while (iteration < MAX_ITERATIONS && !converged && !shouldRemove) {
            // Apply current gap
            core.set_ground_gapping(currentGap);
            core.process_gap();


            // Simulate to get actual Bpeak
            bool gotBpeak = false;
            try {
                // Wind the coil first
                auto magnetic = (*magneticsWithScoring)[i].first;
                magnetic.get_mutable_coil().wind();
                
                auto tempMas = Mas();
                tempMas.set_magnetic(magnetic);
                tempMas.set_inputs(inputs);
                auto simulatedMas = magneticSimulator.simulate(tempMas, true);  // fastMode = true
                
                // Get Bpeak from the excitation after simulation
                if (!simulatedMas.get_inputs().get_operating_points().empty()) {
                    auto& op = simulatedMas.get_inputs().get_operating_points()[0];
                    auto excitation = Inputs::get_primary_excitation(op);
                    if (excitation.get_magnetic_flux_density().has_value()) {
                        auto bField = excitation.get_magnetic_flux_density().value();
                        if (bField.get_processed().has_value() && 
                            bField.get_processed()->get_peak().has_value()) {
                            bPeak = bField.get_processed()->get_peak().value();
                            gotBpeak = true;
                        }
                    }
                }
            } catch (const std::exception& e) {
                // Phase 1 fix: previously silently set `converged = true`
                // (mislabelled an unhandled exception as convergence with no
                // logging). We cannot certify this candidate's Bpeak stays
                // below Bsat for the iterated gap — the conservative response
                // is to drop the candidate rather than ship an unverified
                // design. Logged so the underlying simulator failure is
                // visible.
                logEntry(std::string("Gap-iteration simulator failed for core ")
                         + core.get_name().value_or("unnamed") + ": " + e.what()
                         + " — removing candidate (Bsat could not be verified)",
                         "CoreAdviser");
                shouldRemove = true;
                break;
            }
            
            if (!gotBpeak) {
                converged = true;
                break;
            }

            // Check if we're at or below the safe operating flux.
            if (bPeak <= maxB * 1.01) {  // Allow 1% tolerance
                converged = true;
            } else {
                // Increase gap proportionally to the overshoot, with a 10% margin.
                double newGap = currentGap * (bPeak / maxB) * 1.1;
                
                // Check practical limits
                double columnWidth = core.get_columns()[0].get_width();
                double maxPracticalGap = columnWidth * defaults.coreAdviserMaxPracticalGapColumnWidthFraction;  // fraction of column width
                
                if (newGap > maxPracticalGap) {
                    shouldRemove = true;
                } else {
                    currentGap = newGap;
                    // Recalculate turns for the new gap to maintain inductance
                    double newTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(
                        core, &inputs, DimensionalValues::MINIMUM);
                    (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(newTurns);
                }
            }

            iteration++;
        }

        if (shouldRemove) {
            // Remove this core from the list
            magneticsWithScoring->erase(magneticsWithScoring->begin() + i);
        } else {
            // Apply final gap and update name (avoid duplicates)
            core.set_ground_gapping(currentGap);
            core.process_gap();
            
            std::stringstream ss;
            ss << std::fixed << std::setprecision(2) << (currentGap * 1000);
            auto currentName = core.get_name().value_or("unnamed");
            // Remove existing gap info if present to avoid duplicates
            size_t gappedPos = currentName.find(" gapped ");
            if (gappedPos != std::string::npos) {
                currentName = currentName.substr(0, gappedPos);
            }
            size_t ungappedPos = currentName.find(" ungapped");
            if (ungappedPos != std::string::npos) {
                currentName = currentName.substr(0, ungappedPos);
            }
            if (currentGap > 0) {
                core.set_name(currentName + " gapped " + ss.str() + " mm");
            } else {
                core.set_name(currentName + " ungapped");
            }
            
            (*magneticsWithScoring)[i].first.set_core(core);
            ++i;
        }
    }
}

// ============================================================================
// Option 2: Binary Search Gap Optimization with Analytical Cost Function
// ============================================================================

double CoreAdviser::calculate_gap_cost_analytical(double gap, Inputs& inputs, Core& core, 
                                                   MagnetizingInductance& magnetizingInductance) {
    auto constants = OpenMagnetics::Constants();
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    // Set gap on core
    core.set_ground_gapping(gap);
    core.process_gap();
    
    // Get target inductance
    double targetInductance = resolve_dimensional_values(
        inputs.get_design_requirements().get_magnetizing_inductance(),
        DimensionalValues::NOMINAL);
    
    // Calculate turns for this gap using MagnetizingInductance method
    double turns = magnetizingInductance.calculate_turns_for_gap(core, targetInductance, temperature, frequency);
    
    // Get peak current
    double peakCurrent = get_peak_current(inputs);
    
    // Calculate flux density using MagnetizingInductance method
    double bPeak = magnetizingInductance.calculate_flux_density_peak(core, turns, peakCurrent, temperature, frequency);
    
    // Safe operating flux cap for this material at the operating temperature.
    double maxAllowedB = core.get_magnetic_flux_density_saturation(temperature, /*proportion=*/true);

    // Check saturation constraint
    if (bPeak > maxAllowedB) {
        return std::numeric_limits<double>::max();  // Constraint violation
    }
    
    // Calculate fringing factor penalty using ReluctanceModel
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double fringingFactor = 0.0;
    if (!core.get_gapping().empty()) {
        auto gapReluctance = reluctanceModel->get_gap_reluctance(core.get_gapping()[0]);
        fringingFactor = gapReluctance.get_fringing_factor();
    }
    
    // Estimate core losses (simplified - just use volume and frequency as proxy)
    // The actual losses depend on B, but we're already constraining B
    double coreVolume = core.get_processed_description()->get_effective_parameters().get_effective_volume();
    double lossProxy = coreVolume * frequency * bPeak * bPeak;  // Simplified Steinmetz-like estimate
    
    // Cost function: minimize losses while penalizing high fringing
    // Fringing factor above defaults.coreAdviserMaxFringingFactor is considered excessive
    double fringingPenalty = 0.0;
    if (fringingFactor > defaults.coreAdviserMaxFringingFactor) {
        fringingPenalty = pow((fringingFactor - defaults.coreAdviserMaxFringingFactor) * 10, 2);  // Quadratic penalty
    }
    
    return lossProxy + fringingPenalty;
}

std::pair<double, double> CoreAdviser::optimize_gap_and_turns_binary_search(Inputs& inputs, Core& core) {
    MagnetizingInductance magnetizingInductance(_models["gapReluctance"]);
    auto constants = OpenMagnetics::Constants();
    
    // Ensure core is processed
    if (!core.get_processed_description()) {
        core.process_data();
    }
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs.get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    double targetInductance = resolve_dimensional_values(
        inputs.get_design_requirements().get_magnetizing_inductance(),
        DimensionalValues::NOMINAL);
    double peakCurrent = get_peak_current(inputs);
    double bSat = core.get_magnetic_flux_density_saturation(temperature, false);
    double maxAllowedB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Define gap search bounds
    double columnWidth = core.get_columns()[0].get_width();
    double gMin = constants.residualGap;
    double gMax = 0.4 * columnWidth;  // 40% of column width
    
    // Quick check: can this core possibly work?
    // Minimum turns from saturation: N_min = L * I / (B_max * A)
    double minTurnsFromSaturation = (targetInductance * peakCurrent) / (maxAllowedB * effectiveArea);
    if (minTurnsFromSaturation > 1000) {
        // Unreasonable number of turns - core too small
        return {-1.0, -1.0};
    }
    
    // Golden section search for optimal gap
    const double PHI = 1.618033988749895;
    const double RESPHI = 2.0 - PHI;
    const int MAX_ITERATIONS = 15;
    
    double a = gMin;
    double b = gMax;
    
    // Initial points
    double c = a + RESPHI * (b - a);
    double d = b - RESPHI * (b - a);
    
    Core coreC = core;
    Core coreD = core;
    double fc = calculate_gap_cost_analytical(c, inputs, coreC, magnetizingInductance);
    double fd = calculate_gap_cost_analytical(d, inputs, coreD, magnetizingInductance);
    
    for (int i = 0; i < MAX_ITERATIONS && (b - a) > 1e-6; ++i) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = a + RESPHI * (b - a);
            coreC = core;
            fc = calculate_gap_cost_analytical(c, inputs, coreC, magnetizingInductance);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = b - RESPHI * (b - a);
            coreD = core;
            fd = calculate_gap_cost_analytical(d, inputs, coreD, magnetizingInductance);
        }
    }
    
    // Get optimal gap (midpoint of final interval)
    double optimalGap = (a + b) / 2.0;
    
    // Verify this gap is valid
    Core finalCore = core;
    double finalCost = calculate_gap_cost_analytical(optimalGap, inputs, finalCore, magnetizingInductance);
    
    if (finalCost >= std::numeric_limits<double>::max() / 2) {
        // No valid solution found in search range
        // Try using the analytical method from Option 1 as fallback
        double maxB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
        auto [fallbackGap, fallbackTurns] = magnetizingInductance.calculate_optimal_gap_and_turns(
            core, &inputs, maxB, peakCurrent);
        return {fallbackGap, fallbackTurns};
    }
    
    // Calculate final turns for optimal gap
    core.set_ground_gapping(optimalGap);
    core.process_gap();
    double optimalTurns = magnetizingInductance.calculate_turns_for_gap(core, targetInductance, temperature, frequency);
    
    return {optimalGap, optimalTurns};
}

void CoreAdviser::add_gapping_and_turns_analytical(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring,
                                                    Inputs inputs) {
    if (magneticsWithScoring->empty()) {
        return;
    }
    
    MagnetizingInductance magnetizingInductance(_models["gapReluctance"]);
    
    // Get operating conditions
    double temperature = Defaults().ambientTemperature;
    if (inputs.get_operating_points().size() > 0) {
        temperature = inputs.get_operating_point(0).get_conditions().get_ambient_temperature();
    }
    double peakCurrent = get_peak_current(inputs);
    
    for (size_t i = 0; i < magneticsWithScoring->size();) {
        Core core = (*magneticsWithScoring)[i].first.get_core();
        
        // Skip toroidal cores - they don't have gaps
        if (core.get_shape_family() == CoreShapeFamily::T) {
            core.set_name(core.get_name().value_or("unnamed") + " ungapped");
            (*magneticsWithScoring)[i].first.set_core(core);
            ++i;
            continue;
        }
        
        // Ensure core is processed
        if (!core.get_processed_description()) {
            core.process_data();
        }
        
        // Get Bsat for this material
        double bSat = core.get_magnetic_flux_density_saturation(temperature, false);
        double maxAllowedB = bSat * defaults.maximumProportionMagneticFluxDensitySaturation;
        
        // Use Option 1: Direct analytical calculation (O(1) per core)
        auto [optimalGap, optimalTurns] = magnetizingInductance.calculate_optimal_gap_and_turns(
            core, &inputs, maxAllowedB, peakCurrent);
        
        if (optimalGap < 0 || optimalTurns < 0) {
            // No valid solution - remove this core
            magneticsWithScoring->erase(magneticsWithScoring->begin() + i);
            continue;
        }
        
        // Apply optimal gap
        core.set_ground_gapping(optimalGap);
        core.process_gap();
        
        // Update name with gap info
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << (optimalGap * 1000);
        auto currentName = core.get_name().value_or("unnamed");
        // Remove existing gap info if present
        size_t gappedPos = currentName.find(" gapped ");
        if (gappedPos != std::string::npos) {
            currentName = currentName.substr(0, gappedPos);
        }
        size_t ungappedPos = currentName.find(" ungapped");
        if (ungappedPos != std::string::npos) {
            currentName = currentName.substr(0, ungappedPos);
        }
        if (optimalGap > 0.0001) {
            core.set_name(currentName + " gapped " + ss.str() + " mm");
        } else {
            core.set_name(currentName + " ungapped");
        }
        
        // Apply optimal turns
        (*magneticsWithScoring)[i].first.set_core(core);
        (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(optimalTurns);
        
        ++i;
    }
}

} // namespace OpenMagnetics
