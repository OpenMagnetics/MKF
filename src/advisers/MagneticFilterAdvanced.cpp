#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/CoreLosses.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Temperature.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace OpenMagnetics {

std::pair<bool, double> MagnetomotiveForce::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto coil = magnetic->get_coil();
    double maximumMagnetomotiveForce = 0;
    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        std::vector<double> currentRmsPerParallelPerWinding;
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = inputs->get_operating_points()[operatingPointIndex].get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException("Current is missing in excitation");
            }
            if (!excitation.get_current()->get_processed()) {
                throw InvalidInputException("Current is not processed");
            }
            if (!excitation.get_current()->get_processed()->get_rms()) {
                throw InvalidInputException("Current RMS is not processed");
            }
            auto currentRms = excitation.get_current()->get_processed()->get_rms().value();
            currentRmsPerParallelPerWinding.push_back(currentRms / coil.get_functional_description()[windingIndex].get_number_parallels());
            if (!coil.get_layers_description()) {
                throw CoilNotProcessedException("Coil not wound");
            }
        }
        std::vector<double> magnetomotiveForcePerLayer;
        magnetomotiveForcePerLayer.push_back(0);
        auto layers = coil.get_layers_description().value();
        for (auto layer : layers) {
            double magnetomotiveForceThisLayer = magnetomotiveForcePerLayer.back();
            if (layer.get_type() == ElectricalType::CONDUCTION) {

                auto windingIndex = coil.get_winding_index_by_name(layer.get_partial_windings()[0].get_winding());
                auto numberTurns = coil.get_functional_description()[windingIndex].get_number_turns();  
                auto numberPhysicalTurnsInLayer = 0;
                for (auto parallelProportion : layer.get_partial_windings()[0].get_parallels_proportion()) {
                    numberPhysicalTurnsInLayer += round(numberTurns * parallelProportion);
                }
                numberPhysicalTurnsInLayer *= layer.get_partial_windings().size();
                if (coil.get_functional_description()[windingIndex].get_isolation_side() == IsolationSide::PRIMARY) {
                    magnetomotiveForceThisLayer += numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
                else {
                    magnetomotiveForceThisLayer -= numberPhysicalTurnsInLayer * currentRmsPerParallelPerWinding[windingIndex];
                }
            }
            // Phase 1 fix: previously push_back happened only in the
            // INSULATION (else) branch, so the CONDUCTION branch's update
            // to magnetomotiveForceThisLayer was dropped (dead local) and
            // the per-layer MMF trace was always zero for conductors.
            magnetomotiveForcePerLayer.push_back(magnetomotiveForceThisLayer);
        }

        double maximumMagnetomotiveForceThisOperatingPoint = *max_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        double minimumMagnetomotiveForceThisOperatingPoint = *min_element(magnetomotiveForcePerLayer.begin(), magnetomotiveForcePerLayer.end());
        maximumMagnetomotiveForceThisOperatingPoint = std::max(fabs(maximumMagnetomotiveForceThisOperatingPoint), fabs(minimumMagnetomotiveForceThisOperatingPoint));
        maximumMagnetomotiveForce = std::max(maximumMagnetomotiveForce, maximumMagnetomotiveForceThisOperatingPoint);

    }
    return {true, maximumMagnetomotiveForce};
}

std::pair<bool, double> MagneticFilterLeakageInductance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    // Leakage inductance filter for CMC optimization
    // For Common Mode Chokes, we want to minimize leakage inductance to maximize coupling coefficient
    // Coupling coefficient k = 1 - (Lk / Lm), where Lk is leakage and Lm is magnetizing inductance
    // Lower leakage means higher coupling, which means better common mode rejection
    
    if (magnetic->get_coil().get_functional_description().size() < 2) {
        // Single winding - no leakage calculation possible
        return {true, 0.0};
    }

    // Use the first operating point frequency for leakage calculation
    double frequency = defaults.measurementFrequency;
    if (inputs->get_operating_points().size() > 0) {
        frequency = inputs->get_operating_points()[0].get_excitations_per_winding()[0].get_frequency();
    }

    // Phase 1 fix: previously this method swallowed every exception and
    // returned {false, DBL_MAX} as an in-band sentinel, and used the same
    // sentinel when magnetizingInductance was zero. Both hid real bugs
    // (the DBL_MAX value also poisoned any downstream normalization).
    // Let exceptions propagate and throw explicitly on zero Lm.
    LeakageInductance leakageModel;
    auto leakageOutput = leakageModel.calculate_leakage_inductance(*magnetic, frequency, 0, 1);
    double leakageInductance = resolve_dimensional_values(leakageOutput.get_leakage_inductance_per_winding()[0]);

    // Get magnetizing inductance for normalization
    OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");
    OperatingPoint* operatingPoint = nullptr;
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = &inputs->get_mutable_operating_points()[0];
    }
    auto magnetizingOutput = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
        magnetic->get_mutable_core(),
        magnetic->get_mutable_coil(),
        operatingPoint
    );
    double magnetizingInductance = resolve_dimensional_values(magnetizingOutput.get_magnetizing_inductance());

    if (magnetizingInductance <= 0) {
        throw InvalidInputException(
            "MagneticFilterLeakageInductance: magnetizing inductance is non-positive ("
            + std::to_string(magnetizingInductance) + " H), cannot compute leakage ratio");
    }

    // Score is leakage ratio Lk/Lm - lower is better for CMC.
    // A perfect CMC has k ≈ 1, meaning Lk/Lm ≈ 0.
    double scoring = leakageInductance / magnetizingInductance;

    if (outputs != nullptr && inputs->get_operating_points().size() > 0) {
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            while (outputs->size() < operatingPointIndex + 1) {
                outputs->push_back(Outputs());
            }
            InductanceOutput inductanceOutput;
            if ((*outputs)[operatingPointIndex].get_inductance()) {
                inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
            }
            inductanceOutput.set_leakage_inductance(leakageOutput);
            (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
        }
    }

    return {true, scoring};
}

MagneticFilterTemperature::MagneticFilterTemperature(Inputs inputs, double maximumTemperature)
    : _maximumTemperature(maximumTemperature)
{
    _coreLossesModel = CoreLossesModel::factory(
        std::map<std::string, std::string>{{"coreLosses", "Steinmetz"}});
}

std::pair<bool, double> MagneticFilterTemperature::evaluate_magnetic(
    Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs)
{
    // Phase 1 fix: previously wrapped in `try {} catch (...) { return {true,
    // _maximumTemperature}; }` which silently treated EVERY failure as a
    // pass at the threshold. That hid model bugs (e.g. null _coreLossesModel)
    // and produced "always-acceptable" temperature scores. Let exceptions
    // propagate so callers can either handle them or fail loudly.
    const auto& core = magnetic->get_core();
    double coreLosses = 0.0;
    double ambientTemperature = 25.0;

    const auto& coil = magnetic->get_coil();
    const std::string magneticRef = magnetic->get_reference();
    size_t opIndex = 0;
    for (auto& op : inputs->get_operating_points()) {
        ambientTemperature = op.get_conditions().get_ambient_temperature();
        auto excitation = op.get_excitations_per_winding()[0];

        // Phase 8 (perf): if a prior filter (typically SATURATION on the
        // inductor branch) already ran calculate_inductance_and_magnetic_
        // flux_density for this (magnetic, OP) with the default reluctance
        // model, reuse its peak B. Both filters default-construct
        // MagnetizingInductance and use it on the same fixed coil, so the
        // peak B matches by construction. On cache miss, fall through to
        // the iterative path.
        SignalDescriptor magneticFluxDensity;
        if (auto cached = get_cached_inductance_flux(magneticRef, opIndex)) {
            magneticFluxDensity = cached->magneticFluxDensity;
        } else {
            auto opCopy = op;
            auto [L, B] = _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &opCopy);
            (void)L;
            magneticFluxDensity = B;
        }
        excitation.set_magnetic_flux_density(magneticFluxDensity);
        auto coreLossesMethods = core.get_available_core_losses_methods();
        CoreLossesOutput cl;
        if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(),
                      VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
            cl = _coreLossesModel->get_core_losses(core, excitation, ambientTemperature);
        } else {
            auto proprietaryModel = CoreLossesModel::factory(
                std::map<std::string, std::string>{{"coreLosses", "Proprietary"}});
            cl = proprietaryModel->get_core_losses(core, excitation, ambientTemperature);
        }
        coreLosses += cl.get_core_losses();
        ++opIndex;
    }
    if (!inputs->get_operating_points().empty())
        coreLosses /= inputs->get_operating_points().size();

    TemperatureConfig config;
    config.coreOnly = true;
    config.coreLosses = coreLosses;
    config.ambientTemperature = ambientTemperature;
    config.plotSchematic = false;
    if (!inputs->get_operating_points().empty()) {
        auto& cond = inputs->get_operating_points()[0].get_conditions();
        if (cond.get_cooling()) config.masCooling = cond.get_cooling();
    }

    Temperature temp(*magnetic, config);
    auto result = temp.calculateTemperatures();

    return {result.maximumTemperature <= _maximumTemperature, result.maximumTemperature};
}

} // namespace OpenMagnetics
