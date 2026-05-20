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







std::pair<bool, double> MagneticFilterDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto core = magnetic->get_core();

    auto depth = core.get_depth();

    if (magnetic->get_coil().get_layers_description()) {
        auto coilDepth = magnetic->get_mutable_core().get_columns()[0].get_depth();
        auto layers = magnetic->get_coil().get_layers_description().value();
        for (auto layer : layers) {
            coilDepth += layer.get_dimensions()[0] * 2;
        }
        depth = std::max(depth, coilDepth);
    }

    double volume = core.get_width() * core.get_height() * depth;

    // Hard-reject planar-family cores when wiring is wound wire. PLANAR_E / PLANAR_EL /
    // PLANAR_ER / LP cores assume PCB/planar conductors; a wound coil cannot make
    // efficient use of the wide-shallow window and the resulting design is always
    // dominated by a properly-shaped core. (Note: winding-window width/height are not
    // populated for planar shapes, so we key off the shape family enum.)
    if (inputs && inputs->get_wiring_technology() == WiringTechnology::WOUND) {
        auto family = magnetic->get_mutable_core().get_shape_family();
        if (family == CoreShapeFamily::PLANAR_E ||
            family == CoreShapeFamily::PLANAR_EL ||
            family == CoreShapeFamily::PLANAR_ER ||
            family == CoreShapeFamily::LP) {
            return {false, volume};
        }
    }

    return {true, volume};
}

std::pair<bool, double> MagneticFilterTurnCount::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto coil = magnetic->get_coil();

    // Sum N across all windings. For a CMC the windings are equal so this
    // is 2 × N; for a single inductor it's just N.
    double totalTurns = 0;
    for (const auto& winding : coil.get_functional_description()) {
        totalTurns += static_cast<double>(winding.get_number_turns());
    }

    if (totalTurns <= 0) {
        // Turns not yet assigned — return a neutral, non-rejecting score.
        // Caller is responsible for ordering (run after add_initial_turns_*).
        return {true, 0.0};
    }

    // Manufacturability proxy = wire length used ≈ N × (core size). We must
    // multiply by a size proxy: scoring N alone is monotonically minimised
    // by the largest core in the catalogue (because for a fixed |Z| target
    // Z = µ_complex · ω · N² · Aₑ / lₑ, N ∝ 1/√(µ·Aₑ/lₑ), so bigger core
    // → fewer turns). That made the CMC adviser pick T 140-class toroids
    // for plain 230 V / 5 A / 1 kΩ@150 kHz wizard defaults — see
    // TestTopologyCmc::Test_Cmc_AdviserMustNotPickOversizedToroid_WizardDefaults.
    //
    // get_width() returns the OD for toroids and the A-dimension for
    // two-piece sets; both scale linearly with overall core size and are
    // a reasonable stand-in for mean turn length without forcing the
    // filter to instantiate a bobbin (which the cheap pre-filter pipeline
    // explicitly avoids). Falls back to raw N if width is unavailable.
    // Avoid get_width() here: it throws CoreNotProcessedException on
    // unprocessed cores, and exception throwing in the hot pre-filter loop
    // (thousands of cores × multiple passes) costs seconds. Read the
    // processed description directly instead.
    auto& core = magnetic->get_mutable_core();
    auto processed = core.get_processed_description();
    if (!processed) {
        return {true, totalTurns};
    }
    double widthProxy = processed->get_width();
    if (!std::isfinite(widthProxy) || widthProxy <= 0.0) {
        return {true, totalTurns};
    }

    return {true, totalTurns * widthProxy};
}



std::pair<bool, double> MagneticFilterSolidInsulationRequirements::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = false;
    double scoring = 0;
    auto core = magnetic->get_core();
    auto coreType = core.get_functional_description().get_type();
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

std::pair<bool, double> MagneticFilterTurnsRatios::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_turns_ratios().size() > 0) {
        auto magneticTurnsRatios = magnetic->get_turns_ratios();
        if (magneticTurnsRatios.size() != inputs->get_design_requirements().get_turns_ratios().size()) {
            return {false, 0.0};
        }
        for (size_t i = 0; i < inputs->get_design_requirements().get_turns_ratios().size(); ++i) {
            auto turnsRatioRequirement = inputs->get_design_requirements().get_turns_ratios()[i];
            // TODO: Try all combinations of turns ratios, not just the default order
            if (!check_requirement(turnsRatioRequirement, magneticTurnsRatios[i])) {
                return {false, 0.0};
            }
            scoring += abs(resolve_dimensional_values(turnsRatioRequirement) - resolve_dimensional_values(magneticTurnsRatios[i]));
        }
    }
    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterMaximumDimensions::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    if (inputs->get_design_requirements().get_maximum_dimensions()) {
        auto maximumDimensions = inputs->get_design_requirements().get_maximum_dimensions().value();
        auto magneticDimensions = magnetic->get_maximum_dimensions();
        scoring = sqrt(pow(maximumDimensions.get_width().value() - magneticDimensions[0], 2) + pow(maximumDimensions.get_height().value() - magneticDimensions[1], 2)+ pow(maximumDimensions.get_depth().value() - magneticDimensions[2], 2));
        if (!magnetic->fits(maximumDimensions, true)) {
            valid = false;
        }
    }
    return {valid, scoring};
}



std::pair<bool, double> MagneticFilterMagnetizingInductance::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

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
    if (topology.has_value()) {
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

        if (!check_requirement(inputs->get_design_requirements().get_magnetizing_inductance(), magnetizingInductance)) {
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
    auto core = magnetic->get_core();

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

    for (auto operatingPoint : inputs->get_operating_points()) {
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

std::pair<bool, double> MagneticFilterVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double volume = maximumDimensions[0] * maximumDimensions[1] * maximumDimensions[2];
    return {true, volume};
}

std::pair<bool, double> MagneticFilterArea::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double area = maximumDimensions[0] * maximumDimensions[2];
    return {true, area};
}

std::pair<bool, double> MagneticFilterHeight::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto maximumDimensions = magnetic->get_maximum_dimensions();
    double height = maximumDimensions[1];
    return {true, height};
}


std::pair<bool, double> MagneticFilterTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);

    auto coreTemperatureModel = CoreTemperatureModel::factory(defaults.coreTemperatureModelDefault);

    auto coreTemperature = coreTemperatureModel->get_core_temperature(magnetic->get_core(), losses, defaults.ambientTemperature);
    double calculatedTemperature = coreTemperature.get_maximum_temperature();

    return {true, calculatedTemperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    MagneticFilterLosses fallback;
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES, fallback);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}

std::pair<bool, double> MagneticFilterVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    MagneticFilterLosses lossesFallback;
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES, lossesFallback);
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolume::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring};
}

std::pair<bool, double> MagneticFilterLossesNoProximityTimesVolumeTimesTemperatureRise::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    double losses = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::LOSSES_NO_PROXIMITY, _magneticFilterLossesNoProximity);
    double temperature = cached_or_compute_scoring(magnetic, inputs, outputs, MagneticFilters::TEMPERATURE_RISE, _magneticFilterTemperatureRise);
    auto [volumeValid, volumeScoring] = MagneticFilterVolume().evaluate_magnetic(magnetic, inputs, outputs);
    return {true, losses * volumeScoring * temperature};
}

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
    auto core = magnetic->get_core();
    double coreLosses = 0.0;
    double ambientTemperature = 25.0;

    auto coil = magnetic->get_coil();
    for (auto& op : inputs->get_operating_points()) {
        ambientTemperature = op.get_conditions().get_ambient_temperature();
        auto excitation = op.get_excitations_per_winding()[0];
        auto opCopy = op;
        auto [magnetizingInductance, magneticFluxDensity] =
            _magnetizingInductance.calculate_inductance_and_magnetic_flux_density(core, coil, &opCopy);
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
