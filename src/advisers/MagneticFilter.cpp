#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "physical_models/Temperature.h"
#include "constructive_models/NumberTurns.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/Impedance.h"
#include "physical_models/LeakageInductance.h"
#include <magic_enum.hpp>

#include <cfloat>
#include <cmath>
#include <limits>
#include <algorithm>
#include "support/Exceptions.h"
#include "support/Logger.h"

namespace OpenMagnetics {

// Phase 5: helpers previously defined in this file (is_pqi_or_ui_shape,
// default_loss_filter_models, compute_maximum_power_mean_and_maybe_force_steinmetz,
// prepare_bobbin_for_non_pqi, CoreLossesPick/compute_core_losses_with_negative_guard,
// finalize_losses_scoring, cached_or_compute_scoring) now live in
// advisers/MagneticFilterInternal.h so future split-out filter .cpp files
// can share them without duplicating code or losing the documentation.

// Phase 5: is_energy_storing_topology() now lives as an `inline` function
// in advisers/MagneticFilterInternal.h so the extracted physical/saturation
// filter TU can share it without duplicating the topology switch.

std::shared_ptr<MagneticFilter> MagneticFilter::factory(MagneticFilters filterName, std::optional<Inputs> inputs) {
    switch(filterName) {
        case MagneticFilters::AREA_PRODUCT:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter AREA_PRODUCT");
            }
            return std::make_shared<MagneticFilterAreaProduct>(inputs.value());
        case MagneticFilters::ENERGY_STORED:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter ENERGY_STORED");
            }
            return std::make_shared<MagneticFilterEnergyStored>(inputs.value());
        case MagneticFilters::ESTIMATED_COST:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter ESTIMATED_COST");
            }
            return std::make_shared<MagneticFilterEstimatedCost>(inputs.value());
        case MagneticFilters::COST:
            return std::make_shared<MagneticFilterCost>();
        case MagneticFilters::CORE_AND_DC_LOSSES:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter CORE_AND_DC_LOSSES");
            }
            return std::make_shared<MagneticFilterCoreAndDcLosses>(inputs.value());
        case MagneticFilters::CORE_DC_AND_SKIN_LOSSES:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter CORE_DC_AND_SKIN_LOSSES");
            }
            return std::make_shared<MagneticFilterCoreDcAndSkinLosses>(inputs.value());
        case MagneticFilters::LOSSES:
            return std::make_shared<MagneticFilterLosses>();
        case MagneticFilters::LOSSES_NO_PROXIMITY:
            return std::make_shared<MagneticFilterLossesNoProximity>();
        case MagneticFilters::DIMENSIONS:
            return std::make_shared<MagneticFilterDimensions>();
        case MagneticFilters::CORE_MINIMUM_IMPEDANCE:
            return std::make_shared<MagneticFilterCoreMinimumImpedance>();
        case MagneticFilters::AREA_NO_PARALLELS:
            return std::make_shared<MagneticFilterAreaNoParallels>();
        case MagneticFilters::AREA_WITH_PARALLELS:
            return std::make_shared<MagneticFilterAreaWithParallels>();
        case MagneticFilters::EFFECTIVE_RESISTANCE:
            return std::make_shared<MagneticFilterEffectiveResistance>();
        case MagneticFilters::PROXIMITY_FACTOR:
            return std::make_shared<MagneticFilterProximityFactor>();
        case MagneticFilters::SOLID_INSULATION_REQUIREMENTS:
            return std::make_shared<MagneticFilterSolidInsulationRequirements>();
        case MagneticFilters::TURNS_RATIOS:
            return std::make_shared<MagneticFilterTurnsRatios>();
        case MagneticFilters::MAXIMUM_DIMENSIONS:
            return std::make_shared<MagneticFilterMaximumDimensions>();
        case MagneticFilters::SATURATION:
            return std::make_shared<MagneticFilterSaturation>();
        case MagneticFilters::DC_CURRENT_DENSITY:
            return std::make_shared<MagneticFilterDcCurrentDensity>();
        case MagneticFilters::EFFECTIVE_CURRENT_DENSITY:
            return std::make_shared<MagneticFilterEffectiveCurrentDensity>();
        case MagneticFilters::IMPEDANCE:
            return std::make_shared<MagneticFilterImpedance>();
        case MagneticFilters::MAGNETIZING_INDUCTANCE:
            return std::make_shared<MagneticFilterMagnetizingInductance>();
        case MagneticFilters::SKIN_LOSSES_DENSITY:
            return std::make_shared<MagneticFilterSkinLossesDensity>();
        case MagneticFilters::FRINGING_FACTOR:
            // Phase 1 fix: previously used default ctor, leaving
            // _reluctanceModel null → SIGSEGV in evaluate_magnetic.
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter FRINGING_FACTOR");
            }
            return std::make_shared<MagneticFilterFringingFactor>(inputs.value());
        case MagneticFilters::VOLUME:
            return std::make_shared<MagneticFilterVolume>();
        case MagneticFilters::AREA:
            return std::make_shared<MagneticFilterArea>();
        case MagneticFilters::HEIGHT:
            return std::make_shared<MagneticFilterHeight>();
        case MagneticFilters::TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterTemperatureRise>();
        case MagneticFilters::LOSSES_TIMES_VOLUME:
            return std::make_shared<MagneticFilterLossesTimesVolume>();
        case MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterVolumeTimesTemperatureRise>();
        case MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE:
            return std::make_shared<MagneticFilterLossesTimesVolumeTimesTemperatureRise>();
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter LOSSES_NO_PROXIMITY_TIMES_VOLUME");
            }
            return std::make_shared<MagneticFilterLossesNoProximityTimesVolume>();
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE:
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE");
            }
            return std::make_shared<MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise>();
        case MagneticFilters::MAGNETOMOTIVE_FORCE:
            return std::make_shared<MagnetomotiveForce>();
        case MagneticFilters::LEAKAGE_INDUCTANCE:
            return std::make_shared<MagneticFilterLeakageInductance>();
        case MagneticFilters::TEMPERATURE:
            // Phase 1 fix: previously used default ctor, leaving
            // _coreLossesModel null → NPE in evaluate_magnetic.
            // _maximumTemperature defaults to 130 °C as before.
            if (!inputs) {
                throw InvalidInputException("Inputs needed for filter TEMPERATURE");
            }
            return std::make_shared<MagneticFilterTemperature>(inputs.value(), 130.0);
        case MagneticFilters::TURN_COUNT:
            return std::make_shared<MagneticFilterTurnCount>();
        default:
            throw ModelNotAvailableException("Unknown filter, available options are: {AREA_PRODUCT, ENERGY_STORED, ESTIMATED_COST, COST, CORE_AND_DC_LOSSES, CORE_DC_AND_SKIN_LOSSES, LOSSES, LOSSES_NO_PROXIMITY, DIMENSIONS, CORE_MINIMUM_IMPEDANCE, AREA_NO_PARALLELS, AREA_WITH_PARALLELS, EFFECTIVE_RESISTANCE, PROXIMITY_FACTOR, SOLID_INSULATION_REQUIREMENTS, TURNS_RATIOS, MAXIMUM_DIMENSIONS, SATURATION, DC_CURRENT_DENSITY, EFFECTIVE_CURRENT_DENSITY, IMPEDANCE, MAGNETIZING_INDUCTANCE, FRINGING_FACTOR, SKIN_LOSSES_DENSITY, VOLUME, AREA, HEIGHT, TEMPERATURE_RISE, LOSSES_TIMES_VOLUME, VOLUME_TIMES_TEMPERATURE_RISE, LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE, LOSSES_NO_PROXIMITY_TIMES_VOLUME, LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE, LEAKAGE_INDUCTANCE, TEMPERATURE}");
    }
}









std::pair<bool, double> MagneticFilterSolidInsulationRequirements::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = false;
    double scoring = 0;
    auto coreType = magnetic->get_core().get_functional_description().get_type();
    auto patterns = Coil::get_patterns(*inputs, coreType);
    auto repetitions = Coil::get_repetitions(*inputs, coreType);

    for (auto repetition : repetitions) {
        for (auto pattern : patterns) {
            auto aux = magnetic->get_mutable_coil().check_pattern_and_repetitions_integrity(pattern, repetition);
            pattern = aux.first;
            repetition = aux.second;
            auto combinationsSolidInsulationRequirementsForWires = InsulationCoordinator::get_solid_insulation_requirements_for_wires(*inputs, pattern, repetition);
            for(size_t insulationIndex = 0; insulationIndex < combinationsSolidInsulationRequirementsForWires.size(); ++insulationIndex) {
                auto solidInsulationRequirementsForWires = combinationsSolidInsulationRequirementsForWires[insulationIndex];
                bool combinationValid = true;
                double combinationScoring = 0;
                for (size_t windingIndex = 0; windingIndex < magnetic->get_coil().get_functional_description().size(); ++windingIndex) {
                    auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
                    auto [auxValid, auxScoring] = evaluate_magnetic(winding, solidInsulationRequirementsForWires[windingIndex]);
                    combinationValid &= auxValid;
                    combinationScoring += auxScoring;
                }

                valid |= combinationValid;
                if (valid) {
                    scoring = std::max(scoring, combinationScoring);
                    return {valid, scoring};
                }
            }
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterSolidInsulationRequirements::evaluate_magnetic(Winding winding, WireSolidInsulationRequirements wireSolidInsulationRequirements) {
    auto wire = Coil::resolve_wire(winding);

    if (wire.get_type() == WireType::FOIL || wire.get_type() == WireType::PLANAR) {
        return {true, 0.0};
    }

    if (!wire.resolve_coating()) {
        return {false, 0.0};
    }

    auto coating = wire.resolve_coating().value();

    if (wire.get_type() == WireType::LITZ) {
        auto strand = wire.resolve_strand();
        coating = Wire::resolve_coating(strand).value();
    }

    if (!coating.get_breakdown_voltage()) {
        return {false, 0.0};
        // throw std::runtime_error("Wire " + wire.get_name().value() + " is missing breakdown voltage");
    }

    if (coating.get_breakdown_voltage().value() < wireSolidInsulationRequirements.get_minimum_breakdown_voltage()) {
        return {false, 0.0};
    }

    if (wireSolidInsulationRequirements.get_minimum_grade() && coating.get_grade()) {
        if (coating.get_grade().value() < wireSolidInsulationRequirements.get_minimum_grade().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_minimum_number_layers() && coating.get_number_layers()) {
        if (coating.get_number_layers().value() < wireSolidInsulationRequirements.get_minimum_number_layers().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_minimum_number_layers() || wireSolidInsulationRequirements.get_minimum_grade()) {
        return {false, 0.0};
    }
    
    if (wireSolidInsulationRequirements.get_maximum_grade() && coating.get_grade()) {
        if (coating.get_grade().value() > wireSolidInsulationRequirements.get_maximum_grade().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_maximum_number_layers() && coating.get_number_layers()) {
        if (coating.get_number_layers().value() > wireSolidInsulationRequirements.get_maximum_number_layers().value()) {
            return {false, 0.0};
        }
    }
    else if (wireSolidInsulationRequirements.get_maximum_number_layers() || wireSolidInsulationRequirements.get_maximum_grade()) {
        return {false, 0.0};
    }

    double scoring = 0;
    if (wireSolidInsulationRequirements.get_minimum_breakdown_voltage() > 0) {
        scoring = wireSolidInsulationRequirements.get_minimum_breakdown_voltage() - coating.get_breakdown_voltage().value();
    }

    return {true, scoring};
}



std::pair<bool, double> MagneticFilterMagnetizingInductance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    // If no inductance requirement is set (e.g. CMC/DMC with only impedance requirements),
    // this filter has nothing to check — pass with neutral score.
    {
        auto& lmReq = inputs->get_design_requirements().get_magnetizing_inductance();
        if (!lmReq.get_minimum() && !lmReq.get_nominal() && !lmReq.get_maximum()) {
            return {true, 0};
        }
        // Lm = 0 means "not specified" (same convention as pre_process_inputs
        // in CoreAdviserPipeline.cpp): LLC inputs export nominal 0 when Lm is
        // left free. Zero inductance is not a designable target, and a literal
        // check would build an empty [0, 0] acceptance band and reject every
        // core (Test_FastAdviser_LLC_Lm_Zero).
        if (lmReq.get_nominal() && lmReq.get_nominal().value() <= 0 &&
            !lmReq.get_minimum() && !lmReq.get_maximum()) {
            return {true, 0};
        }
    }

    // For interference-suppression topologies (CMC, DMC), the inductance requirement was
    // derived from Z/(2πf). The correct gate is the impedance filter, which already checked
    // the complex permeability at the noise frequency. Checking inductance here would use
    // DC-biased real permeability (which rolls off significantly), producing a false rejection
    // of physically valid cores. The impedance check is authoritative for these topologies.
    {
        auto topology = inputs->get_design_requirements().get_topology();
        if (topology.has_value() &&
            (topology.value() == Topologies::COMMON_MODE_CHOKE ||
             topology.value() == Topologies::DIFFERENTIAL_MODE_CHOKE)) {
            return {true, 0};
        }
    }

    // Transformer vs Inductor Detection
    // ==================================
    // Transformers: Use initial permeability without DC bias iteration
    //   - Skip operating point to avoid feedback loop where magnetizing current causes permeability collapse
    //   - Ungapped high-permeability ferrite will naturally exceed minimum inductance requirement
    //
    // Inductors/Energy-storing: Use full iterative calculation with DC bias
    //   - Account for permeability rolloff with magnetizing current DC offset
    //
    // Detection priority:
    // 1. Use topology if specified (most reliable)
    // 2. Fall back to inductance field heuristic (minimum-only = transformer)
    //
    auto topology = inputs->get_design_requirements().get_topology();
    bool isTransformer;
    if (windings_on_single_isolation_side(inputs->get_design_requirements().get_isolation_sides())) {
        // All windings on one isolation side -> (coupled) inductor, never a transformer,
        // regardless of the converter topology (e.g. Weinberg L1 input coupled inductor).
        isTransformer = false;
    } else if (topology.has_value()) {
        isTransformer = !is_energy_storing_topology(topology);
    } else {
        // Legacy heuristic: minimum-only inductance = transformer
        isTransformer = inputs->get_design_requirements().get_magnetizing_inductance().get_minimum() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_nominal() &&
                         !inputs->get_design_requirements().get_magnetizing_inductance().get_maximum();
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
        auto operatingPoint = inputs->get_operating_points()[operatingPointIndex];
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");

        // For transformers, calculate inductance without operating point to use initial permeability
        MagnetizingInductanceOutput aux;
        try {
            if (isTransformer) {
                aux = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
                    magnetic->get_mutable_core(), magnetic->get_mutable_coil(), nullptr);
            } else {
                aux = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(
                    magnetic->get_mutable_core(), magnetic->get_mutable_coil(), &operatingPoint);
            }
        } catch (const std::exception& e) {
            return {false, 0};
        }
        double magnetizingInductance = resolve_dimensional_values(aux.get_magnetizing_inductance());
        double requiredInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance());
        scoring += fabs(requiredInductance - magnetizingInductance);

        bool reqOk = check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance);

        if (!reqOk) {
            valid = false;
        }
        else {
            if (outputs != nullptr) {
                while (outputs->size() < operatingPointIndex + 1) {
                    outputs->push_back(Outputs());
                }
                InductanceOutput inductanceOutput;
                if ((*outputs)[operatingPointIndex].get_inductance()) {
                    inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
                }
                inductanceOutput.set_magnetizing_inductance(aux);
                (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
            }
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

MagneticFilterFringingFactor::MagneticFilterFringingFactor(Inputs inputs, std::map<std::string, std::string> models) {
    _models = models;
    _magneticEnergy = MagneticEnergy(models);
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
    ReluctanceModels reluctanceModel;
    from_json(_models["gapReluctance"], reluctanceModel);
    _reluctanceModel = ReluctanceModel::factory(reluctanceModel);
}

MagneticFilterFringingFactor::MagneticFilterFringingFactor(Inputs inputs) {
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
    _reluctanceModel = ReluctanceModel::factory();
}

std::pair<bool, double> MagneticFilterFringingFactor::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    if (core.get_shape_family() == CoreShapeFamily::T) {
        return {true, 1};
    }
    else if (core.get_gapping().size() == 0) {
        return {true, 1};
    }
    else {
        double maximumGapLength = _reluctanceModel->get_gapping_by_fringing_factor(core, _fringingFactorLitmit);
        double gapLength = core.get_gapping()[0].get_length();
        bool valid = gapLength <= maximumGapLength;
        if (!valid) {
            return {false, 1};
        }
        else {
            return {true, 1};
        }
    }
}

std::pair<bool, double> MagneticFilterSkinLossesDensity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    auto temperature = inputs->get_maximum_temperature();

    if (inputs->get_operating_points().size() > 0 && magnetic->get_mutable_coil().get_functional_description().size() != inputs->get_operating_points()[0].get_excitations_per_winding().size()) {
        return {false, 0.0};        
    }

    for (const auto& operatingPoint : inputs->get_operating_points()) {
        for (size_t windingIndex = 0; windingIndex < magnetic->get_mutable_coil().get_functional_description().size(); ++windingIndex) {
            auto excitation = operatingPoint.get_excitations_per_winding()[windingIndex];
            if (!excitation.get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
            }
            auto current = excitation.get_current().value();
            auto winding = magnetic->get_coil().get_functional_description()[windingIndex];
            auto [auxValid, auxScoring] = evaluate_magnetic(winding, current, temperature);
            valid &= auxValid;
            scoring += auxScoring;
        }
    }

    scoring /= inputs->get_operating_points().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterSkinLossesDensity::evaluate_magnetic(Winding winding, SignalDescriptor current, double temperature) {
    auto wire = Coil::resolve_wire(winding);

    double skinEffectLossesPerMeter = WindingSkinEffectLosses::calculate_skin_effect_losses_per_meter(wire, current, temperature).first;
    double valid = true;
    return {valid, skinEffectLossesPerMeter};
}












} // namespace OpenMagnetics
