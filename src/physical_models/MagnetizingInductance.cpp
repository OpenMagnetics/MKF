#include "physical_models/MagnetizingInductance.h"

#include "processors/Inputs.h"
#include "physical_models/MagneticField.h"
#include "physical_models/Reluctance.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

namespace OpenMagnetics {

double calculate_air_inductance(int64_t numberTurnsPrimary, Core core) {

    auto bobbin = Bobbin::create_quick_bobbin(core);
    auto bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    auto bobbinWindingWindowWidth = bobbin.get_winding_window_dimensions()[0];
    auto bobbinWindingWindowHeight = bobbin.get_winding_window_dimensions()[1];
    auto meanLengthRadius = (bobbinColumnDepth + bobbinColumnWidth) / 2 + bobbinWindingWindowWidth / 4;

    double coilInternalArea = std::numbers::pi * pow(meanLengthRadius, 2);

    double coreColumnArea = core.get_processed_description()->get_columns()[0].get_area();

    double airAreaProportion = (coilInternalArea - coreColumnArea) / coilInternalArea;

    return Constants().vacuumPermeability * pow(numberTurnsPrimary, 2) * (coilInternalArea * airAreaProportion * 2) / bobbinWindingWindowHeight;
}

std::pair<MagnetizingInductanceOutput, SignalDescriptor> MagnetizingInductance::calculate_inductance_and_magnetic_flux_density(Magnetic magnetic, OperatingPoint* operatingPoint) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_inductance_and_magnetic_flux_density(core, coil, operatingPoint);
}


std::pair<MagnetizingInductanceOutput, SignalDescriptor> MagnetizingInductance::calculate_inductance_and_magnetic_flux_density(Core core, Coil coil, OperatingPoint* operatingPoint) {

    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;

    if (operatingPoint) {

        temperature = operatingPoint->get_conditions().get_ambient_temperature();
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            frequency = operatingPoint->get_mutable_excitations_per_winding()[0].get_frequency();
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);

            if (excitation.get_current()) {
                excitation.set_current(OpenMagnetics::Inputs::standardize_waveform(excitation.get_current().value(), frequency));
            }
            if (excitation.get_voltage()) {
                excitation.set_voltage(OpenMagnetics::Inputs::standardize_waveform(excitation.get_voltage().value(), frequency));
            }
            operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;

            Inputs::make_waveform_size_power_of_two(operatingPoint);
        }
    }

    std::pair<MagnetizingInductanceOutput, SignalDescriptor> result;
    double numberWindings = coil.get_functional_description().size();
    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    double currentInitialPermeability;

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double currentTotalReluctance;
    double modifiedTotalReluctance = 0;
    double modifiedMagnetizingInductance = 5e-3;
    double currentMagnetizingInductance;


    if (operatingPoint) {
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
            if (!excitation.get_voltage()) {
                auto current = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value();
                auto currentWaveform = current.get_waveform().value();
                if (!is_size_power_of_2(currentWaveform.get_data())) {
                    auto currentSampledWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency);
                    current.set_waveform(currentSampledWaveform);
                    operatingPoint->get_mutable_excitations_per_winding()[0].set_current(current);
                }
            }

            if (excitation.get_voltage()) {
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("voltage_data vector size is not a power of 2");
                }
            }
            if (excitation.get_current()) {
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_current().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("current_data vector size is not a power of 2");
                }
            }
            if (!excitation.get_voltage()) {
                if (!excitation.get_magnetizing_current()) {
                    Inputs::set_current_as_magnetizing_current(operatingPoint);
                }
                auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
                if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                    throw std::invalid_argument("magnetizing_current_data vector size from current is not a power of 2");
                }

            }
        }
    }

    MagnetizingInductanceOutput magnetizingInductanceOutput;
    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    modifiedTotalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / modifiedTotalReluctance;

    size_t externalTimeout = 1;
    do {
        currentMagnetizingInductance = modifiedMagnetizingInductance;
        size_t internalTimeout = 1;

        do {
            currentTotalReluctance = modifiedTotalReluctance;
            modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / currentTotalReluctance;


            if (operatingPoint) {
                if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
                    OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);

                    // If the converter model already computed the magnetizing current
                    // (e.g. DMC summing all winding currents), respect it — don't overwrite.
                    if (excitation.get_magnetizing_current()) {
                        // Already set — skip derivation
                    }
                    else if (numberWindings == 1 && excitation.get_current()) {
                        Inputs::set_current_as_magnetizing_current(operatingPoint);
                    }
                    // CMC check must come BEFORE is_multiport_inductor. CMCs
                    // have every winding on the same isolation side (L, N, PE
                    // all primary-side mains), which makes is_multiport_inductor
                    // return true — wrongly routing to the multiport path that
                    // uses the primary winding's raw current (DM + CM) as
                    // magnetizing current. That pumps the DM line current into
                    // the core flux calculation and overstates B by the ratio
                    // I_line / I_cm. The correct CMC path averages the winding
                    // currents and drops the DM DC offset.
                    else if (Inputs::can_be_common_mode_choke(*operatingPoint) && core.get_type() == CoreType::TOROIDAL) {
                        auto magnetizingCurrent = Inputs::get_common_mode_choke_magnetizing_current(*operatingPoint);
                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }
                    else if (Inputs::is_multiport_inductor(*operatingPoint, coil.get_isolation_sides())) {
                        auto magnetizingCurrent = Inputs::get_multiport_inductor_magnetizing_current(*operatingPoint);
                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }
                    else if (excitation.get_voltage()) {
                        auto voltage = operatingPoint->get_mutable_excitations_per_winding()[0].get_voltage().value();
                        auto sampledVoltageWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);

                        auto turnsRatios = coil.get_turns_ratios();
                        bool addOffset = Inputs::include_dc_offset_into_magnetizing_current(*operatingPoint, turnsRatios);

                        auto magnetizingCurrent = Inputs::calculate_magnetizing_current(excitation,
                                                                                                sampledVoltageWaveform,
                                                                                                modifiedMagnetizingInductance,
                                                                                                false,
                                                                                                addOffset);

                        auto sampledMagnetizingCurrentWaveform = Inputs::calculate_sampled_waveform(magnetizingCurrent.get_waveform().value(), excitation.get_frequency());
                        magnetizingCurrent.set_harmonics(Inputs::calculate_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency()));
                        magnetizingCurrent.set_processed(Inputs::calculate_processed_data(magnetizingCurrent, sampledMagnetizingCurrentWaveform, false));

                        excitation.set_magnetizing_current(magnetizingCurrent);
                        operatingPoint->get_mutable_excitations_per_winding()[0] = excitation;
                    }

                    auto aux = operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value().get_waveform().value();
                    if (aux.get_data().size() > 0 && ((aux.get_data().size() & (aux.get_data().size() - 1)) != 0)) {
                        throw std::invalid_argument("magnetizing_current_data vector size from voltage is not a power of 2");
                    }

                    if (!operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current()->get_waveform()->get_time()) {
                        auto magnetizingCurrent = Inputs::standardize_waveform(operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), excitation.get_frequency());
                        operatingPoint->get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                    }

                    auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint->get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), currentTotalReluctance, numberTurnsPrimary);
                    auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
                    result.second = magneticFluxDensity;
                    auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);
                    double switchingFrequency = Inputs::get_switching_frequency(operatingPoint->get_mutable_excitations_per_winding()[0]);

                    double hFieldDcBias = magneticFieldStrength.get_processed().value().get_offset();
                    if (magneticFieldStrength.get_harmonics().value().get_frequencies()[1] < switchingFrequency) {
                        for (size_t i = 0; i < magneticFieldStrength.get_harmonics().value().get_frequencies().size() - 1; ++i) {
                            if (magneticFieldStrength.get_harmonics().value().get_frequencies()[i] >= switchingFrequency) {
                                break;
                            }
                            hFieldDcBias = std::max(hFieldDcBias, magneticFieldStrength.get_harmonics().value().get_amplitudes()[i]);
                        }
                    }

                    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, hFieldDcBias, frequency);

                    magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
                    modifiedTotalReluctance = magnetizingInductanceOutput.get_core_reluctance();
                    modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / modifiedTotalReluctance;
                }
            }

            internalTimeout--;
            if (internalTimeout == 0) {
                break;
            }
        } while (fabs(currentTotalReluctance - modifiedTotalReluctance) / modifiedTotalReluctance >= 0.1);

        externalTimeout--;
        if (externalTimeout == 0) {
            break;
        }
    } while (fabs(currentMagnetizingInductance - modifiedMagnetizingInductance) / modifiedMagnetizingInductance >= 0.1);

    if (operatingPoint) {
        if (operatingPoint->get_mutable_excitations_per_winding().size() > 0) {
            OperatingPointExcitation excitation = Inputs::get_primary_excitation(*operatingPoint);
            if (!excitation.get_voltage()) {
                operatingPoint->get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, modifiedMagnetizingInductance));
            }
        }
    }

    auto& settings = Settings::GetInstance();
    if (settings.get_magnetizing_inductance_include_air_inductance()) {
        modifiedMagnetizingInductance += calculate_air_inductance(numberTurnsPrimary, core);
    }

    DimensionWithTolerance magnetizingInductanceWithTolerance;
    magnetizingInductanceWithTolerance.set_nominal(modifiedMagnetizingInductance);
    magnetizingInductanceOutput.set_magnetizing_inductance(magnetizingInductanceWithTolerance);

    result.first = magnetizingInductanceOutput;
    return result;
}

double MagnetizingInductance::calculate_inductance_air_solenoid(Magnetic magnetic) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_inductance_air_solenoid(core, coil);
}


double MagnetizingInductance::calculate_inductance_air_solenoid(Core core, Coil coil) {
    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double airCoreReluctance = reluctanceModel->get_air_cored_reluctance(coil.resolve_bobbin());
    auto modifiedMagnetizingInductance = pow(numberTurnsPrimary, 2) / airCoreReluctance;
    return modifiedMagnetizingInductance;
}

MagnetizingInductanceOutput MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping(Core core, Coil coil, OperatingPoint* operatingPoint) {
    auto inductance_and_magnetic_flux_density = calculate_inductance_and_magnetic_flux_density(core, coil, operatingPoint);

    return inductance_and_magnetic_flux_density.first;
}

MagnetizingInductanceOutput MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping(Magnetic magnetic, OperatingPoint* operatingPoint) {
    auto inductance_and_magnetic_flux_density = calculate_inductance_and_magnetic_flux_density(magnetic, operatingPoint);

    return inductance_and_magnetic_flux_density.first;
}

int MagnetizingInductance::calculate_number_turns_from_gapping_and_inductance(Core core, Inputs* inputs, DimensionalValues preferredValue) {
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance(), preferredValue);
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPoint operatingPoint;

    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    }
    OpenMagnetics::InitialPermeability initialPermeability;
    int numberTurnsPrimary;
    size_t timeout = 10;
    double currentInitialPermeability;
    double modifiedInitialPermeability;

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double totalReluctance;

    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    if (inputs->get_operating_points().size() > 0) {
        OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);
        if (!excitation.get_magnetizing_current() && !excitation.get_voltage()) {
            Inputs::set_current_as_magnetizing_current(&operatingPoint);
            inputs->set_operating_point_by_index(operatingPoint, 0);
        }
    }

    while (true) {

        auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
        totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
        numberTurnsPrimary = std::round(sqrt(desiredMagnetizingInductance * totalReluctance));

        if (inputs->get_operating_points().size() > 0) {
            auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), totalReluctance, numberTurnsPrimary);
            auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
            auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);

            auto Hoffset = magneticFieldStrength.get_processed().value().get_mutable_offset();
            modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, Hoffset, frequency);

            if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
                break;
            }
            else {
                currentInitialPermeability = modifiedInitialPermeability;
                timeout--;
            }
        }
    }

    if (inputs->get_operating_points().size() > 0) {
        OperatingPointExcitation excitation = Inputs::get_primary_excitation(operatingPoint);
        if (!excitation.get_voltage()) {
            operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, desiredMagnetizingInductance));
            inputs->set_operating_point_by_index(operatingPoint, 0);
        }
    }

    return std::max(1, numberTurnsPrimary);
}

double MagnetizingInductance::calculate_gap_from_saturation_constraint(Core core,
                                                                       Inputs* inputs,
                                                                       double targetMagneticFluxDensity,
                                                                       double magnetizingCurrentPeak) {
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance(), DimensionalValues::NOMINAL);
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    // Start with energy-based gap as initial guess
    auto constants = OpenMagnetics::Constants();
    double initialGap = (2 * desiredMagnetizingInductance * pow(magnetizingCurrentPeak, 2) * constants.vacuumPermeability) / 
                        (effectiveArea * pow(targetMagneticFluxDensity, 2));
    
    double gapLength = initialGap;
    const double gapStepFactor = 1.5;  // Increase gap by 50% per iteration
    const int maxIterations = 15;
    
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        // Create core with current gap
        Core testCore = core;
        auto columns = testCore.get_processed_description().value().get_columns();
        
        auto basicCentralGap = CoreGap();
        basicCentralGap.set_type(GapType::SUBTRACTIVE);
        basicCentralGap.set_length(gapLength);
        basicCentralGap.set_area(columns[0].get_area());
        
        auto basicLateralGap = CoreGap();
        basicLateralGap.set_type(GapType::RESIDUAL);
        basicLateralGap.set_length(constants.residualGap);
        
        std::vector<CoreGap> gapping;
        gapping.push_back(basicCentralGap);
        for (size_t i = 1; i < columns.size(); ++i) {
            auto lateralGap = basicLateralGap;
            lateralGap.set_area(columns[i].get_area());
            gapping.push_back(lateralGap);
        }
        testCore.get_mutable_functional_description().set_gapping(gapping);
        testCore.process_gap();
        
        // Calculate total reluctance manually (classic formula)
        // Core reluctance: R_core = l_e / (μ₀ * μ_r * A_e)
        double effectiveLength = testCore.get_processed_description()->get_effective_parameters().get_effective_length();
        double coreReluctance = effectiveLength / (constants.vacuumPermeability * 2000 * effectiveArea);
        
        // Gap reluctance: R_gap = l_gap / (μ₀ * A_gap)
        double gapArea = columns[0].get_area();
        double gapReluctance = gapLength / (constants.vacuumPermeability * gapArea);
        
        double totalReluctance = coreReluctance + gapReluctance;
        
        // Calculate turns for this gap to achieve target inductance
        int numberTurns = std::round(sqrt(desiredMagnetizingInductance * totalReluctance));
        numberTurns = std::max(1, numberTurns);
        
        // Calculate B-field: B = (N * I) / (R_total * A_e)
        double magneticFlux = numberTurns * magnetizingCurrentPeak / totalReluctance;
        double bPeak = magneticFlux / effectiveArea;
        
        // Check convergence
        if (bPeak <= targetMagneticFluxDensity) {
            return gapLength;
        }
        
        // Increase gap for next iteration
        gapLength *= gapStepFactor;
    }
    
    // Return last calculated gap even if not converged
    return gapLength;
}

Core get_core_with_ground_gapping(Core core, double gapLength) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

Core get_core_with_distributed_gapping(Core core, double gapLength, size_t numberDistributedGaps) {
    auto constants = OpenMagnetics::Constants();
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::SUBTRACTIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::RESIDUAL);
    basicLateralGap.set_length(constants.residualGap);
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberDistributedGaps; ++i) {
        gapping.push_back(basicCentralGap);
    }
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

Core get_core_with_spacer_gapping(Core core, double gapLength) {
    auto basicCentralGap = CoreGap();
    basicCentralGap.set_type(GapType::ADDITIVE);
    basicCentralGap.set_length(gapLength);
    auto basicLateralGap = CoreGap();
    basicLateralGap.set_type(GapType::ADDITIVE);
    basicLateralGap.set_length(gapLength);
    std::vector<CoreGap> gapping;
    gapping.push_back(basicCentralGap);
    for (size_t i = 0; i < core.get_processed_description().value().get_columns().size() - 1; ++i) {
        gapping.push_back(basicLateralGap);
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    core.process_gap();
    return core;
}

std::vector<CoreGap> MagnetizingInductance::calculate_gapping_from_number_turns_and_inductance(Core core,
                                                                                               Coil coil,
                                                                                               Inputs* inputs,
                                                                                               GappingType gappingType,
                                                                                               size_t decimals) {
    double frequency = Defaults().coreAdviserFrequencyReference;
    double temperature = Defaults().ambientTemperature;
    OperatingPointExcitation excitation;
    OperatingPoint operatingPoint;

    auto constants = OpenMagnetics::Constants();
    if (inputs->get_operating_points().size() > 0) {
        operatingPoint = inputs->get_operating_point(0);
        excitation = Inputs::get_primary_excitation(operatingPoint);
        temperature = operatingPoint.get_conditions().get_ambient_temperature();
        frequency = operatingPoint.get_mutable_excitations_per_winding()[0].get_frequency();
    }

    double numberTurnsPrimary = coil.get_functional_description()[0].get_number_turns();
    double desiredMagnetizingInductance = resolve_dimensional_values(inputs->get_design_requirements().get_magnetizing_inductance());
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    OpenMagnetics::InitialPermeability initialPermeability;
    size_t timeout = 10;
    double modifiedInitialPermeability;
    double currentInitialPermeability;

    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    double neededTotalReluctance = pow(numberTurnsPrimary, 2) / desiredMagnetizingInductance;

    currentInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);

    if (!excitation.get_voltage() && excitation.get_current()) {
        Inputs::set_current_as_magnetizing_current(&operatingPoint);
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    while (true) {
        if (!excitation.get_magnetizing_current()) {
            break;
        }
        auto magneticFlux = OpenMagnetics::MagneticField::calculate_magnetic_flux(operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value(), neededTotalReluctance, numberTurnsPrimary);
        auto magneticFluxDensity = OpenMagnetics::MagneticField::calculate_magnetic_flux_density(magneticFlux, effectiveArea);
        auto magneticFieldStrength = OpenMagnetics::MagneticField::calculate_magnetic_field_strength(magneticFluxDensity, currentInitialPermeability);

        modifiedInitialPermeability = initialPermeability.get_initial_permeability(core.resolve_material(), temperature, magneticFieldStrength.get_processed().value().get_mutable_offset(), frequency);

        if (fabs(currentInitialPermeability - modifiedInitialPermeability) < 1 || timeout == 0) {
            break;
        }
        else {
            currentInitialPermeability = modifiedInitialPermeability;
            timeout--;
        }
    }

    if (!excitation.get_voltage() && excitation.get_current()) {
        operatingPoint.get_mutable_excitations_per_winding()[0].set_voltage(Inputs::calculate_induced_voltage(excitation, desiredMagnetizingInductance));
        inputs->set_operating_point_by_index(operatingPoint, 0);
    }

    double gapLength = constants.residualGap;
    double gapLengthModification = constants.initialGapLengthForSearching;
    bool increasingGap = true;
    double fringingFactorOneGap = 0;

    double reluctance = 0;
    size_t numberDistributedGaps = 3;
    timeout = 100;
    Core gappedCore;

    while (true) {
        reluctance = 0;
        switch (gappingType) {
            case GappingType::GROUND:
                gappedCore = get_core_with_ground_gapping(core, gapLength);
                break;
            case GappingType::SPACER:
                gappedCore = get_core_with_spacer_gapping(core, gapLength);
                break;
            case GappingType::RESIDUAL:
                throw GapException("Residual type cannot be chosen to calculate the needed gapping");
                break;
            case GappingType::DISTRIBUTED:
                while (numberDistributedGaps > 3) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap < constants.minimumDistributedFringingFactor &&
                        numberDistributedGaps > 1) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps -= 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                while (true) {
                    gappedCore = get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps);
                    fringingFactorOneGap = reluctanceModel->get_gap_reluctance(gappedCore.get_gapping()[0]).get_fringing_factor();
                    if (fringingFactorOneGap > constants.maximumDistributedFringingFactor) {
                        gapLength *= numberDistributedGaps;
                        numberDistributedGaps += 2;
                        gapLength /= numberDistributedGaps;
                    }
                    else {
                        break;
                    }
                }
                break;
            default:
                throw GapException("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
        }

        auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(gappedCore, currentInitialPermeability);
        reluctance = magnetizingInductanceOutput.get_core_reluctance();

        if (fabs(neededTotalReluctance - reluctance) / neededTotalReluctance < 0.001 || timeout == 0) {
            break;
        }
        else {
            if (neededTotalReluctance < reluctance && increasingGap) {
                increasingGap = false;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (neededTotalReluctance > reluctance && !increasingGap) {
                increasingGap = true;
                gapLengthModification = std::max(gapLengthModification / 2., constants.residualGap);
            }
            if (increasingGap) {
                gapLength += gapLengthModification;
            }
            else {
                gapLength = std::max(constants.residualGap, gapLength - gapLengthModification);
            }

            timeout--;
        }
    }

    gapLength = std::max(constants.residualGap, roundFloat(gapLength, decimals));

    switch (gappingType) {
        case GappingType::GROUND:
            return get_core_with_ground_gapping(core, gapLength).get_gapping();
        case GappingType::SPACER:
            return get_core_with_spacer_gapping(core, gapLength).get_gapping();
        case GappingType::RESIDUAL:
            throw GapException("Residual type cannot be chosen to calculate the needed gapping");
            break;
        case GappingType::DISTRIBUTED:
            return get_core_with_distributed_gapping(core, gapLength, numberDistributedGaps).get_gapping();
            break;
        default:
            throw GapException("Unknown type of gap, options are {GROUND, SPACER, RESIDUAL, DISTRIBUTED}");
    }
}

double MagnetizingInductance::calculate_turns_for_gap(
    Core& core,
    double targetInductance,
    double temperature,
    double frequency) 
{
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        core.resolve_material(), temperature, std::nullopt, frequency);
    
    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    double totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    
    double turns = sqrt(targetInductance * totalReluctance);
    return std::max(1.0, std::round(turns));
}

double MagnetizingInductance::calculate_flux_density_peak(
    Core& core,
    double numberTurns,
    double peakCurrent,
    double temperature,
    double frequency) 
{
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        core.resolve_material(), temperature, std::nullopt, frequency);
    
    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(core, currentInitialPermeability);
    double totalReluctance = magnetizingInductanceOutput.get_core_reluctance();
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // B = Φ / A = (N * I / R) / A
    double magneticFlux = numberTurns * peakCurrent / totalReluctance;
    double bPeak = magneticFlux / effectiveArea;
    
    return bPeak;
}

double MagnetizingInductance::calculate_flux_density_peak_from_voltage(
    Core& core,
    double numberTurns,
    double voltagePeak,
    double frequency) 
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Faraday's Law: V = N * Ae * dB/dt
    // For sinusoidal: V_peak = N * Ae * ω * B_peak
    // Therefore: B_peak = V_peak / (N * Ae * ω)
    double omega = 2 * std::numbers::pi * frequency;
    
    if (numberTurns <= 0 || effectiveArea <= 0 || omega <= 0) {
        return 0.0;
    }
    
    double bPeak = voltagePeak / (numberTurns * effectiveArea * omega);
    return bPeak;
}

double MagnetizingInductance::calculate_turns_from_voltage_and_max_flux_density(
    Core& core,
    double voltagePeak,
    double frequency,
    double maxFluxDensity)
{
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // From Faraday's Law: B_peak = V_peak / (N * Ae * ω)
    // Therefore: N = V_peak / (Ae * ω * B_max)
    double omega = 2 * std::numbers::pi * frequency;
    
    if (effectiveArea <= 0 || omega <= 0 || maxFluxDensity <= 0) {
        return 1.0;
    }
    
    double turns = voltagePeak / (effectiveArea * omega * maxFluxDensity);
    return std::max(1.0, std::ceil(turns));
}

std::pair<double, double> MagnetizingInductance::calculate_optimal_gap_and_turns(
    const Core& core,
    Inputs* inputs,
    double maxFluxDensity,
    double peakCurrent) 
{
    auto constants = OpenMagnetics::Constants();
    double targetInductance = resolve_dimensional_values(
        inputs->get_design_requirements().get_magnetizing_inductance(), 
        DimensionalValues::NOMINAL);
    
    double temperature = Defaults().ambientTemperature;
    double frequency = Defaults().coreAdviserFrequencyReference;
    if (inputs->get_operating_points().size() > 0) {
        temperature = inputs->get_operating_point(0).get_conditions().get_ambient_temperature();
        frequency = inputs->get_operating_point(0).get_excitations_per_winding()[0].get_frequency();
    }
    
    // Create a mutable copy for processing
    Core tempCore = core;
    if (!tempCore.get_processed_description()) {
        tempCore.process_data();
    }
    
    double effectiveArea = tempCore.get_processed_description()->get_effective_parameters().get_effective_area();
    
    // Step 1: Calculate minimum turns from saturation constraint
    // B_max = L * I_peak / (N * A)  =>  N_min = L * I_peak / (B_max * A)
    double turnsFromSaturation = (targetInductance * peakCurrent) / (maxFluxDensity * effectiveArea);
    double minTurns = std::max(1.0, std::ceil(turnsFromSaturation));
    
    // Step 2: Calculate required reluctance for these turns
    // L = N² / R  =>  R = N² / L
    double requiredReluctance = (minTurns * minTurns) / targetInductance;
    
    // Step 3: Get core reluctance (ungapped) to determine needed gap reluctance
    OpenMagnetics::InitialPermeability initialPermeability;
    ReluctanceModels reluctanceModelEnum;
    from_json(_models["gapReluctance"], reluctanceModelEnum);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory(reluctanceModelEnum);
    
    double currentInitialPermeability = initialPermeability.get_initial_permeability(
        tempCore.resolve_material(), temperature, std::nullopt, frequency);
    
    // Estimate core reluctance from effective parameters
    double effectiveLength = tempCore.get_processed_description()->get_effective_parameters().get_effective_length();
    double coreReluctance = effectiveLength / (constants.vacuumPermeability * currentInitialPermeability * effectiveArea);
    
    // Step 4: Calculate needed gap reluctance
    double neededGapReluctance = requiredReluctance - coreReluctance;
    
    if (neededGapReluctance <= 0) {
        // No gap needed - core reluctance alone is enough
        // But verify saturation constraint is met
        double actualTurns = sqrt(targetInductance * coreReluctance);
        return {0.0, std::max(1.0, std::round(actualTurns))};
    }
    
    // Step 5: Calculate gap length from gap reluctance
    // R_gap = g / (μ₀ * A_gap)  =>  g = R_gap * μ₀ * A_gap
    // Note: This is simplified - doesn't account for fringing
    double gapArea = tempCore.get_processed_description()->get_columns()[0].get_area();
    double gapLength = neededGapReluctance * constants.vacuumPermeability * gapArea;
    
    // Step 6: Validate gap is practical
    double columnWidth = tempCore.get_columns()[0].get_width();
    double maxGap = 0.4 * columnWidth;  // 40% of column width max
    double minGap = constants.residualGap;
    
    if (gapLength > maxGap) {
        // Gap too large - core is too small for this application
        return {-1.0, -1.0};
    }
    
    gapLength = std::max(minGap, gapLength);
    
    // Step 7: Refine turns with actual gap (accounts for fringing via reluctance model)
    Core gappedCore = tempCore;
    gappedCore.set_ground_gapping(gapLength);
    gappedCore.process_gap();
    
    double finalTurns = calculate_turns_for_gap(gappedCore, targetInductance, temperature, frequency);
    
    // Step 8: Verify saturation constraint with final values
    double finalBpeak = calculate_flux_density_peak(gappedCore, finalTurns, peakCurrent, temperature, frequency);
    
    if (finalBpeak > maxFluxDensity * 1.05) {  // Allow 5% tolerance
        // Need to increase turns to reduce B
        finalTurns = std::ceil((targetInductance * peakCurrent) / (maxFluxDensity * effectiveArea));
        // Recalculate gap for new turns
        double newRequiredReluctance = (finalTurns * finalTurns) / targetInductance;
        double newNeededGapReluctance = newRequiredReluctance - coreReluctance;
        if (newNeededGapReluctance > 0) {
            gapLength = newNeededGapReluctance * constants.vacuumPermeability * gapArea;
            gapLength = std::min(maxGap, std::max(minGap, gapLength));
        }
    }
    
    return {gapLength, finalTurns};
}

} // namespace OpenMagnetics
