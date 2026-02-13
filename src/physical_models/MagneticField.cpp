#include "physical_models/MagneticField.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Magnetic.h"
#include "support/CoilMesher.h"
#include "Models.h"
#include "physical_models/Reluctance.h"
#include "processors/MagneticSimulator.h"
#include "physical_models/InitialPermeability.h"
#include "support/Settings.h"
#include "support/Utils.h"

#include <cmath>

// Some platforms (Emscripten, Apple libc++) don't have std::comp_ellint_1/2
// Use custom implementations from Utils.h via namespace injection
#if defined(__EMSCRIPTEN__) || defined(__APPLE__)
namespace std {
    using OpenMagnetics::comp_ellint_1;
    using OpenMagnetics::comp_ellint_2;
}
#endif
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"

namespace OpenMagnetics {

SignalDescriptor MagneticField::calculate_magnetic_flux(SignalDescriptor magnetizingCurrent,
                                                          double reluctance,
                                                          double numberTurns) {
    SignalDescriptor magneticFlux;
    Waveform magneticFluxWaveform;
    std::vector<double> magneticFluxData;
    auto compressedMagnetizingCurrentWaveform = magnetizingCurrent.get_waveform().value();

    if (Inputs::is_waveform_sampled(compressedMagnetizingCurrentWaveform)) {
        compressedMagnetizingCurrentWaveform = Inputs::compress_waveform(compressedMagnetizingCurrentWaveform);
    }

    for (auto& datum : compressedMagnetizingCurrentWaveform.get_data()) {
        magneticFluxData.push_back(datum * numberTurns / reluctance);
    }

    if (compressedMagnetizingCurrentWaveform.get_time()) {
        magneticFluxWaveform.set_time(compressedMagnetizingCurrentWaveform.get_time());
    }

    magneticFluxWaveform.set_data(magneticFluxData);
    magneticFlux.set_waveform(magneticFluxWaveform);
    if (magnetizingCurrent.get_harmonics()) {
        auto harmonics = magnetizingCurrent.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] *= numberTurns / reluctance;
        }
        magneticFlux.set_harmonics(harmonics);
    }
    return magneticFlux;
}

SignalDescriptor MagneticField::calculate_magnetic_flux_density(SignalDescriptor magneticFlux,
                                                                  double area) {
    SignalDescriptor magneticFluxDensity;
    Waveform magneticFluxDensityWaveform;
    std::vector<double> magneticFluxDensityData;
    auto magneticFluxWaveform = magneticFlux.get_waveform().value();

    if (magneticFluxWaveform.get_time()) {
        magneticFluxDensityWaveform.set_time(magneticFluxWaveform.get_time());
    }

    for (auto& datum : magneticFluxWaveform.get_data()) {
        magneticFluxDensityData.push_back(datum / area);
    }

    magneticFluxDensityWaveform.set_data(magneticFluxDensityData);
    magneticFluxDensity.set_waveform(magneticFluxDensityWaveform);
    if (magneticFlux.get_harmonics()) {
        auto harmonics = magneticFlux.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= area;
        }
        magneticFluxDensity.set_harmonics(harmonics);
    }
    magneticFluxDensity.set_processed(
        Inputs::calculate_basic_processed_data(magneticFluxDensityWaveform));

    return magneticFluxDensity;
}

SignalDescriptor MagneticField::calculate_magnetic_field_strength(SignalDescriptor magneticFluxDensity,
                                                                    double initialPermeability) {
    SignalDescriptor magneticFieldStrength;
    Waveform magneticFieldStrengthWaveform;
    std::vector<double> magneticFieldStrengthData;
    auto constants = Constants();
    auto magneticFluxDensityWaveform = magneticFluxDensity.get_waveform().value();

    if (magneticFluxDensityWaveform.get_time()) {
        magneticFieldStrengthWaveform.set_time(magneticFluxDensityWaveform.get_time());
    }

    for (auto& datum : magneticFluxDensityWaveform.get_data()) {
        magneticFieldStrengthData.push_back(datum / (initialPermeability * constants.vacuumPermeability));
    }

    magneticFieldStrengthWaveform.set_data(magneticFieldStrengthData);
    magneticFieldStrength.set_waveform(magneticFieldStrengthWaveform);
    if (magneticFluxDensity.get_harmonics()) {
        auto harmonics = magneticFluxDensity.get_harmonics().value();
        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            harmonics.get_mutable_amplitudes()[harmonicIndex] /= (initialPermeability * constants.vacuumPermeability);
        }
        magneticFieldStrength.set_harmonics(harmonics);
    }
    magneticFieldStrength.set_processed(
        Inputs::calculate_basic_processed_data(magneticFieldStrengthWaveform));

    return magneticFieldStrength;
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory(MagneticFieldStrengthModels modelName) {
    if (modelName == MagneticFieldStrengthModels::BINNS_LAWRENSON) {
        return std::make_shared<MagneticFieldStrengthBinnsLawrensonModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::LAMMERANER) {
        return std::make_shared<MagneticFieldStrengthLammeranerModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::WANG) {
        return std::make_shared<MagneticFieldStrengthWangModel>();
    }
    else if (modelName == MagneticFieldStrengthModels::ALBACH) {
        return std::make_shared<MagneticFieldStrengthAlbach2DModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Magnetic Field Strength model, available options are: {BINNS_LAWRENSON, LAMMERANER, WANG, ALBACH}");
}

std::shared_ptr<MagneticFieldStrengthFringingEffectModel> MagneticField::factory(MagneticFieldStrengthFringingEffectModels modelName) {
    if (modelName == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        return std::make_shared<MagneticFieldStrengthAlbachModel>();
    }
    else if (modelName == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
        return std::make_shared<MagneticFieldStrengthRoshenModel>();
    }
    else if (modelName == MagneticFieldStrengthFringingEffectModels::SULLIVAN) {
        return std::make_shared<MagneticFieldStrengthSullivanModel>();
    }
    else
        throw ModelNotAvailableException("Unknown Magnetic Field Strength Fringing Effect model, available options are: {ALBACH, ROSHEN, SULLIVAN}");
}

std::shared_ptr<MagneticFieldStrengthModel> MagneticField::factory() {
    auto defaults = Defaults();
    return factory(defaults.magneticFieldStrengthModelDefault);
}

bool is_inside_inducing_turns(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, Wire* inducingWire) {
    double distanceX = fabs(inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0]);
    double distanceY = fabs(inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1]);
    if (inducingWire->get_type() == WireType::ROUND || inducingWire->get_type() == WireType::LITZ) {
        if (hypot(distanceX, distanceY) < inducingWire->get_maximum_outer_width() / 2) {
            return true;
        }
    }
    else {
        if (distanceX < inducingWire->get_maximum_outer_width() / 2 && distanceY < inducingWire->get_maximum_outer_height() / 2) {
            return true;
        }
    }

    return false;
}


bool is_inside_turns(std::vector<Turn> turns, FieldPoint inducedFieldPoint, std::vector<Wire> wires, Magnetic magnetic) {
    for (auto turn : turns) {
        auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turn.get_winding());
        double distanceX = fabs(turn.get_coordinates()[0] - inducedFieldPoint.get_point()[0]);
        double distanceY = fabs(turn.get_coordinates()[1] - inducedFieldPoint.get_point()[1]);
        if (wires[windingIndex].get_type() == WireType::ROUND || wires[windingIndex].get_type() == WireType::LITZ) {
            if (hypot(distanceX, distanceY) < wires[windingIndex].get_maximum_outer_width() / 2) {
                return true;
            }
        }
        else {
            if (distanceX < wires[windingIndex].get_maximum_outer_width() / 2 && distanceY < wires[windingIndex].get_maximum_outer_height() / 2) {
                return true;
            }
        }
    }

    return false;
}

bool is_inside_core(FieldPoint inducedFieldPoint, double coreColumnWidth, double coreWidth, CoreShapeFamily coreShapeFamily) {
    if (coreShapeFamily != CoreShapeFamily::T) {
        return false;
    }
    double radius = sqrt(pow(inducedFieldPoint.get_point()[0], 2) + pow(inducedFieldPoint.get_point()[1], 2));

    if (radius * 1.05 > coreWidth / 2) {
        return false;
    }
    if (radius * 0.95 < (coreWidth / 2 - coreColumnWidth)) {
        return false;
    }
    return true;
}

double get_magnetic_field_strength_gap(OperatingPoint& operatingPoint, Magnetic magnetic, double frequency) {
    auto numberTurns = magnetic.get_mutable_coil().get_number_turns(0);
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    OpenMagnetics::InitialPermeability initial_permeability;
    double initialPermeability = initial_permeability.get_initial_permeability(magnetic.get_mutable_core().resolve_material(), std::nullopt, std::nullopt, frequency);
    double reluctance = reluctanceModel->get_core_reluctance(magnetic.get_core(), initialPermeability).get_core_reluctance();
    
    // Calculate magnetizing current if missing
    if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
        auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
        auto includeDcCurrent = Inputs::include_dc_offset_into_magnetizing_current(operatingPoint, magnetic.get_turns_ratios());
        auto magnetizingCurrent = Inputs::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                               resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()),
                                                                               true, includeDcCurrent);
        operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
    }

    auto magnetizingCurrent = operatingPoint.get_mutable_excitations_per_winding()[0].get_magnetizing_current().value();
    if (!magnetizingCurrent.get_waveform()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Magnetizing current is missing waveform");
    }
    if (!magnetizingCurrent.get_waveform()->get_time()) {
        magnetizingCurrent = Inputs::standardize_waveform(magnetizingCurrent, frequency);
    }
    auto magneticFlux = MagneticField::calculate_magnetic_flux(magnetizingCurrent, reluctance, numberTurns);
    auto magneticFluxDensity = MagneticField::calculate_magnetic_flux_density(magneticFlux, magnetic.get_core().get_processed_description()->get_effective_parameters().get_effective_area());

    double magneticFieldStrengthGap = magneticFluxDensity.get_processed()->get_peak().value() / Constants().vacuumPermeability;
    return magneticFieldStrengthGap;
}

WindingWindowMagneticStrengthFieldOutput MagneticField::calculate_magnetic_field_strength_field(OperatingPoint operatingPoint, Magnetic magnetic, std::optional<Field> externalInducedField, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding, std::optional<CoilMesherModels> coilMesherModel) {
    auto& settings = OpenMagnetics::Settings::GetInstance();
    auto includeFringing = settings.get_magnetic_field_include_fringing();
    CoilMesher coilMesher; 
    std::vector<Field> inducingFields;
    auto core = magnetic.get_core();

    if (!core.is_gap_processed()) {
        core.process_gap();
    }
    auto gapping = core.get_functional_description().get_gapping();
    double coreColumnWidth = core.get_columns()[0].get_width();
    auto processedDescription = core.get_processed_description().value();
    double coreWidth = processedDescription.get_width();
    auto coreShapeFamily = core.get_shape_family();

    std::vector<int8_t> currentDirectionPerWinding;
    if (!customCurrentDirectionPerWinding) {
        currentDirectionPerWinding.push_back(1);
        for (size_t windingIndex = 1; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
            currentDirectionPerWinding.push_back(-1);
        }
    }
    else {
        currentDirectionPerWinding = customCurrentDirectionPerWinding.value();
    }

    if (externalInducedField){
        auto aux = coilMesher.generate_mesh_inducing_coil(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold(), currentDirectionPerWinding, coilMesherModel);
        // We only process the harmonic that comes from the external field
        for (auto field : aux) {
            if (field.get_frequency() == externalInducedField.value().get_frequency()) {
                inducingFields.push_back(field);
                break;
            }
        }
    }
    else {
        inducingFields = coilMesher.generate_mesh_inducing_coil(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold(), currentDirectionPerWinding);
    }

    if (!magnetic.get_coil().get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description in coil");
    }
    _wirePerWinding = magnetic.get_mutable_coil().get_wires();
    _model->_wirePerWinding = _wirePerWinding;
    auto turns = magnetic.get_coil().get_turns_description().value();

    // Set up ALBACH model if being used
    if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
        auto albach2DModel = std::dynamic_pointer_cast<MagneticFieldStrengthAlbach2DModel>(_model);
        if (albach2DModel) {
            albach2DModel->setupFromMagnetic(magnetic, _wirePerWinding);
        }
    }

    std::vector<ComplexField> complexFieldPerHarmonic;
    for (auto& fieldPerHarmonic : inducingFields) {
        ComplexField field;
        field.set_frequency(fieldPerHarmonic.get_frequency());
        complexFieldPerHarmonic.push_back(field);
    }

    std::vector<Field> inducedFields;
    if (externalInducedField) {
        inducedFields.push_back(externalInducedField.value());
    }
    else {
        inducedFields = coilMesher.generate_mesh_induced_coil(magnetic, operatingPoint, settings.get_harmonic_amplitude_threshold());
    }

    // For ALBACH model, we only compute the turns contribution (air coil field).
    // The gap fringing field is handled by whichever fringing effect model is configured
    // (ALBACH or ROSHEN fringing). This allows users to choose their preferred fringing model.
    // ALBACH supports any number of gaps since the fringing is computed separately.
    bool isAlbach = (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH);
    
    // Use the configured fringing model (default is ROSHEN, can be overridden to ALBACH)
    std::shared_ptr<MagneticFieldStrengthFringingEffectModel> fringingModel = _fringingEffectModel;

    if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH) {
        if (includeFringing) {
            if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {

                auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
                auto includeDcCurrent = Inputs::include_dc_offset_into_magnetizing_current(operatingPoint, magnetic.get_turns_ratios());
                auto magnetizingCurrent = Inputs::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                                       resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()),
                                                                                       true, includeDcCurrent);

                operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                // throw std::runtime_error("Operating point is missing magnetizing current");
            }
            if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()) {
                auto excitations = operatingPoint.get_excitations_per_winding();
                auto magnetizingCurrent = excitations[0].get_magnetizing_current().value();
                auto processed = Inputs::calculate_basic_processed_data(magnetizingCurrent.get_waveform().value());
                magnetizingCurrent.set_processed(processed);
                excitations[0].set_magnetizing_current(magnetizingCurrent);
                operatingPoint.set_excitations_per_winding(excitations);
                // throw std::runtime_error("Operating point is missing magnetizing current processed data");

            }

            for (size_t harmonicIndex = 0; harmonicIndex < inducingFields.size(); ++harmonicIndex){

                if (inducingFields[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {

                    double frequency = inducingFields[harmonicIndex].get_frequency();
                    double magneticFieldStrengthGap = get_magnetic_field_strength_gap(operatingPoint, magnetic, frequency);
                    for (auto& gap : gapping) {
                        if (gap.get_coordinates().value()[0] < 0) {
                            continue;
                        }
                        auto fieldPoint = fringingModel->get_equivalent_inducing_point_for_gap(gap, magneticFieldStrengthGap);
                        
                        inducingFields[harmonicIndex].get_mutable_data().push_back(fieldPoint);
                    }
                }
            }
        }
    }

    for (size_t harmonicIndex = 0; harmonicIndex < inducingFields.size(); ++harmonicIndex){
        std::vector<ComplexFieldPoint> fieldPoints;

        if (inducedFields[harmonicIndex].get_data().size() == 0) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Empty complexField");
        }

        // For ALBACH model, use a more efficient approach that calculates
        // the total field from all turns at once for each induced point.
        if (_magneticFieldStrengthModel == MagneticFieldStrengthModels::ALBACH) {
            auto albach2DModel = std::dynamic_pointer_cast<MagneticFieldStrengthAlbach2DModel>(_model);
            if (albach2DModel) {
                // Update turn currents based on harmonic data
                auto& harmonicData = inducingFields[harmonicIndex].get_data();
                std::vector<double> turnCurrents(turns.size(), 0.0);
                
                // Collect fringing field inducing points (those without turn_index)
                std::vector<FieldPoint> fringingPoints;
                
                for (auto& inducingPoint : harmonicData) {
                    if (inducingPoint.get_turn_index()) {
                        size_t turnIdx = inducingPoint.get_turn_index().value();
                        if (turnIdx < turnCurrents.size()) {
                            turnCurrents[turnIdx] = inducingPoint.get_value();
                        }
                    } else {
                        // This is a fringing field equivalent point (from ALBACH fringing model)
                        fringingPoints.push_back(inducingPoint);
                    }
                }
                albach2DModel->updateTurnCurrents(turnCurrents);

                // Update skin depths for frequency-dependent current distribution (Wang 2018)
                // At high frequency, current concentrates at conductor edges
                double frequency = inducingFields[harmonicIndex].get_frequency();
                if (frequency > 0 && !_wirePerWinding.empty()) {
                    // Use the first wire to calculate a representative skin depth
                    // In practice, all wires in a winding should have similar material
                    double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(
                        _wirePerWinding[0], frequency, operatingPoint.get_conditions().get_ambient_temperature());
                    albach2DModel->updateSkinDepths(skinDepth);
                }
                
                // Create a model for computing fringing field contribution from equivalent current loops (ALBACH model)
                // Using BINNS_LAWRENSON to compute the field from the equivalent current loops
                auto fringingFieldModel = factory(MagneticFieldStrengthModels::BINNS_LAWRENSON);
                
                // For ROSHEN fringing, we need the gap field strength
                double magneticFieldStrengthGap = 0;
                if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN && includeFringing) {
                    double frequency = inducingFields[harmonicIndex].get_frequency();
                    if (frequency == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {
                        magneticFieldStrengthGap = get_magnetic_field_strength_gap(operatingPoint, magnetic, frequency);
                    }
                }
                
                // Calculate field at each induced point directly from all turns
                for (auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
                    // Skip points inside the core
                    if (is_inside_core(inducedFieldPoint, coreColumnWidth, coreWidth, coreShapeFamily)) {
                        continue;
                    }
                    
                    auto complexFieldPoint = albach2DModel->calculateTotalFieldAtPoint(inducedFieldPoint);
                    
                    // Add fringing field contribution based on configured fringing model
                    if (includeFringing && inducingFields[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {
                        if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ALBACH) {
                            // ALBACH fringing: use equivalent current loops
                            for (auto& fringingPoint : fringingPoints) {
                                auto fringingContrib = fringingFieldModel->get_magnetic_field_strength_between_two_points(fringingPoint, inducedFieldPoint);
                                complexFieldPoint.set_real(complexFieldPoint.get_real() + fringingContrib.get_real());
                                complexFieldPoint.set_imaginary(complexFieldPoint.get_imaginary() + fringingContrib.get_imaginary());
                            }
                        } else if (_magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
                            // ROSHEN fringing: compute field directly from each gap
                            for (auto& gap : gapping) {
                                if (gap.get_coordinates().value()[0] < 0) {
                                    continue;
                                }
                                auto fringingContrib = _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, inducedFieldPoint);
                                complexFieldPoint.set_real(complexFieldPoint.get_real() + fringingContrib.get_real());
                                complexFieldPoint.set_imaginary(complexFieldPoint.get_imaginary() + fringingContrib.get_imaginary());
                            }
                        }
                    }
                    
                    if (std::isnan(complexFieldPoint.get_real()) || std::isnan(complexFieldPoint.get_imaginary())) {
                        throw NaNResultException("NaN found in ALBACH magnetic field calculation");
                    }
                    
                    fieldPoints.push_back(complexFieldPoint);
                }
                complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
                continue; // Skip the standard per-turn-pair loop for this harmonic
            }
        }

        for (auto& inducedFieldPoint : inducedFields[harmonicIndex].get_data()) {
            double totalInducedFieldX = 0;
            double totalInducedFieldY = 0;

            // ROSHEN fringing is computed per-point in this loop (not via equivalent current loops)
            // Skip if using ALBACH since ROSHEN fringing is already added in the ALBACH branch above
            if (!isAlbach && _magneticFieldStrengthFringingEffectModel == MagneticFieldStrengthFringingEffectModels::ROSHEN) {
                // For the main harmonic we calculate the fringing effect for each gap
                if (includeFringing && inducedFields[harmonicIndex].get_frequency() == operatingPoint.get_excitations_per_winding()[0].get_frequency()) {
                    if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()) {
                        auto magnetizingInductance = MagneticSimulator().calculate_magnetizing_inductance(operatingPoint, magnetic);
                        auto includeDcCurrent = Inputs::include_dc_offset_into_magnetizing_current(operatingPoint, magnetic.get_turns_ratios());
                        auto magnetizingCurrent = Inputs::calculate_magnetizing_current(operatingPoint.get_mutable_excitations_per_winding()[0],
                                                                                               resolve_dimensional_values(magnetizingInductance.get_magnetizing_inductance()),
                                                                                               true, includeDcCurrent);

                        operatingPoint.get_mutable_excitations_per_winding()[0].set_magnetizing_current(magnetizingCurrent);
                        // throw std::runtime_error("Operating point is missing magnetizing current");
                    }
                    if (!operatingPoint.get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()) {
                        auto excitations = operatingPoint.get_excitations_per_winding();
                        auto magnetizingCurrent = excitations[0].get_magnetizing_current().value();
                        auto processed = Inputs::calculate_basic_processed_data(magnetizingCurrent.get_waveform().value());
                        magnetizingCurrent.set_processed(processed);
                        excitations[0].set_magnetizing_current(magnetizingCurrent);
                        operatingPoint.set_excitations_per_winding(excitations);
                        // throw std::runtime_error("Operating point is missing magnetizing current processed data");
                    }
                    double frequency = inducingFields[harmonicIndex].get_frequency();
                    double magneticFieldStrengthGap = get_magnetic_field_strength_gap(operatingPoint, magnetic, frequency);

                    for (auto& gap : gapping) {
                        if (gap.get_coordinates().value()[0] < 0) {
                            continue;
                        }
                        auto complexFieldPoint = _fringingEffectModel->get_magnetic_field_strength_between_gap_and_point(gap, magneticFieldStrengthGap, inducedFieldPoint);

                        totalInducedFieldX += complexFieldPoint.get_real();
                        totalInducedFieldY += complexFieldPoint.get_imaginary();
                        if (std::isnan(complexFieldPoint.get_real())) {
                            throw NaNResultException("NaN found in Roshen's fringing field");
                        }
                        if (std::isnan(complexFieldPoint.get_imaginary())) {
                            throw NaNResultException("NaN found in Roshen's fringing field");
                        }
                    }
                }
            }

            for (auto& inducingFieldPoint : inducingFields[harmonicIndex].get_data()) {
                std::optional<size_t> windingIndex = std::nullopt;
                if (inducingFieldPoint.get_turn_index()) {
                    windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turns[inducingFieldPoint.get_turn_index().value()].get_winding());
                }
                if (inducingFieldPoint.get_turn_index()) {
                    if (inducedFieldPoint.get_turn_index()) {
                        if (inducedFieldPoint.get_turn_index().value() == inducingFieldPoint.get_turn_index().value()) {
                            continue;
                        }
                    }
                    // else if (is_inside_inducing_turns(inducingFieldPoint, inducedFieldPoint, &_wirePerWinding[windingIndex])) {
                    //     continue;
                    // }
                    // else if (is_inside_turns(turns, inducedFieldPoint, _wirePerWinding, magnetic)) {
                    //     continue;
                    // }
                    else if (is_inside_core(inducedFieldPoint, coreColumnWidth, coreWidth, coreShapeFamily)) {
                        continue;
                    }
                }

                auto complexFieldPoint = _model->get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint, windingIndex);

                totalInducedFieldX += complexFieldPoint.get_real();
                totalInducedFieldY += complexFieldPoint.get_imaginary();
                if (std::isnan(complexFieldPoint.get_real())) {
                    throw NaNResultException("NaN found in magnetic field calculation");
                }
                if (std::isnan(complexFieldPoint.get_imaginary())) {
                    throw NaNResultException("NaN found in magnetic field calculation");
                }
            }
            ComplexFieldPoint complexFieldPoint;
            complexFieldPoint.set_point(inducedFieldPoint.get_point());
            complexFieldPoint.set_real(totalInducedFieldX);
            complexFieldPoint.set_imaginary(totalInducedFieldY);
            if (inducedFieldPoint.get_turn_index()) {
                complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
            }
            if (inducedFieldPoint.get_label()) {
                complexFieldPoint.set_label(inducedFieldPoint.get_label().value());
            }
            fieldPoints.push_back(complexFieldPoint);
        }
        complexFieldPerHarmonic[harmonicIndex].set_data(fieldPoints);
    }

    WindingWindowMagneticStrengthFieldOutput windingWindowMagneticStrengthFieldOutput;
    windingWindowMagneticStrengthFieldOutput.set_field_per_frequency(complexFieldPerHarmonic);
    windingWindowMagneticStrengthFieldOutput.set_method_used(to_string(_magneticFieldStrengthModel));
    windingWindowMagneticStrengthFieldOutput.set_origin(ResultOrigin::SIMULATION);
    return windingWindowMagneticStrengthFieldOutput;
}

ComplexFieldPoint MagneticFieldStrengthWangModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    double Hx = 0;
    double Hy = 0;
    if (!inducingWireIndex) {
        return MagneticFieldStrengthLammeranerModel().get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint);
    }
    else {
        double c = 0;
        double h = 0;
        if (_wirePerWinding[inducingWireIndex.value()].get_type() == WireType::FOIL) {
            c = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value());
            h = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value());
        }
        else {
            h = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value());
            c = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value());
        }
        double k = c / h;
        double lambda = 0.01 * k + 0.66;

        if (inducedFieldPoint.get_label() && inducingFieldPoint.get_label()) {
            auto inducingLabel = inducingFieldPoint.get_label().value();
            auto inducedLabel = inducedFieldPoint.get_label().value();
            auto current = inducingFieldPoint.get_value();
            double distanceX = inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0];
            double distanceY = inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1];
            double distance = hypot(distanceX, distanceY);
            if (inducingLabel == "left") {
                if (inducedLabel == "left") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "right") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "top") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "bottom") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
                else {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducedLabel: " + inducedLabel);
                }
            }
            else if (inducingLabel == "right") {
                if (inducedLabel == "right") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "left") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "bottom") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "top") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
            }
            else if (inducingLabel == "bottom") {
                if (inducedLabel == "bottom") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "top") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "right") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "left") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
                else {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducedLabel: " + inducedLabel);
                }
            }
            else if (inducingLabel == "top") {
                if (inducedLabel == "top") {
                    double tetha1 = asin(fabs(distanceY) / distance);
                    Hy = 0.5 * current / (2 * std::numbers::pi * lambda * h) + 0.5 * current * cos(tetha1) / (2 * std::numbers::pi * sqrt(pow(lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "bottom") {
                    double tetha2 = asin(fabs(distanceY) / distance);
                    Hy = -0.5 * current / (2 * std::numbers::pi * (c - lambda * h)) - 0.5 * current * cos(tetha2) / (2 * std::numbers::pi * sqrt(pow(c - lambda * h, 2) + pow(distanceY, 2)));
                }
                else if (inducedLabel == "left") {
                    if (distanceY > 0) {
                        Hx = 0;
                    }
                    else {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                }
                else if (inducedLabel == "right") {
                    if (distanceY > 0) {
                        Hx = (current - Hy * 2 * h) / (2 * c);
                    }
                    else {
                        Hx = 0;
                    }
                }
            }
            else {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wrong inducingLabel: " + inducingLabel);
            }
        }
        else {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wang Magnetic Field model must be used woth his CoilMesher model");
        }
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;   
}

ComplexFieldPoint MagneticFieldStrengthBinnsLawrensonModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    double Hx;
    double Hy;

    double distanceX = inducingFieldPoint.get_point()[0] - inducedFieldPoint.get_point()[0];
    double distanceY = inducingFieldPoint.get_point()[1] - inducedFieldPoint.get_point()[1];

    bool computeRectangularFiled = false;

    if (inducingWireIndex) {
        if (_wirePerWinding[inducingWireIndex.value()].get_type() != WireType::ROUND && _wirePerWinding[inducingWireIndex.value()].get_type() != WireType::LITZ) {
            computeRectangularFiled = true;
        }
        else {
            if (hypot(distanceX, distanceY) < _wirePerWinding[inducingWireIndex.value()].get_maximum_outer_width() / 2) {
                Hx = 0;
                Hy = 0;
            }
            else { 
                double divisor = 2 * std::numbers::pi * (pow(distanceY, 2) + pow(distanceX, 2));
                Hx = -inducingFieldPoint.get_value() * (distanceY) / divisor;
                Hy = inducingFieldPoint.get_value() * (distanceX) / divisor;
            }
            if (std::isnan(Hx) || std::isnan(Hy)) {
                throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
            }
        }
    }
    else {
        double divisor = 2 * std::numbers::pi * (pow(distanceY, 2) + pow(distanceX, 2));
        Hx = -inducingFieldPoint.get_value() * (distanceY) / divisor;
        Hy = inducingFieldPoint.get_value() * (distanceX) / divisor;
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }


    if (computeRectangularFiled) {
        double a = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_width().value()) / 2;
        double b = resolve_dimensional_values(_wirePerWinding[inducingWireIndex.value()].get_conducting_height().value()) / 2;
        double x = inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0];
        double y = inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1];

        if (inducingFieldPoint.get_rotation()) {
            double modulo = hypot(x, y);
            double currentAngle = atan2(y, x);
            double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
            if (currentAngle < 0) {
                currentAngle += 2 * std::numbers::pi;
            }
            double totalAngle = currentAngle - turnAngle;

            x = modulo * cos(totalAngle);
            y = modulo * sin(totalAngle);
        }

        double r1 = hypot(y + b, x - a);
        double r2 = hypot(y + b, x + a);
        double r3 = hypot(y - b, x + a);
        double r4 = hypot(y - b, x - a);

        double tetha1 = atan((y + b) / (x - a));
        double tetha2 = atan((y + b) / (x + a));
        double tetha3 = atan((y - b) / (x + a));
        double tetha4 = atan((y - b) / (x - a));
        if (fabs(x) < a && fabs(y) < b) {
            Hx = 0;
            Hy = 0;
        }
        else if (std::isnan(tetha1) || std::isnan(tetha2) || std::isnan(tetha3) || std::isnan(tetha4)) {
            Hx = 0;
            Hy = 0;
        } 
        else {
            if (x == a) {
                if ((y + b) > 0) {
                    tetha1 = std::numbers::pi / 2;
                }
                else {
                    tetha1 = -std::numbers::pi / 2;
                }
                if ((y - b) > 0) {
                    tetha4 = std::numbers::pi / 2;
                }
                else {
                    tetha4 = -std::numbers::pi / 2;
                }
            }

            if (x > a && -b < y && y < b) {

            }
            else {
                if (x > a && y < -b) {
                    tetha1 += 2 * std::numbers::pi;
                }
                else if (x < a || y < -b) {
                    tetha1 += std::numbers::pi;
                }

                if (x > -a && y < -b) {
                    tetha2 += 2 * std::numbers::pi;
                }
                else if (x < -a || y < -b) {
                    tetha2 += std::numbers::pi;
                }

                if (x > -a && y < b) {
                    tetha3 += 2 * std::numbers::pi;
                }
                else if (x < -a || y < b) {
                    tetha3 += std::numbers::pi;
                }

                if (x > a && y < b) {
                    tetha4 += 2 * std::numbers::pi;
                }
                else if (x < a || y < b) {
                    tetha4 += std::numbers::pi;
                }
            }
        }

        double common_part = inducingFieldPoint.get_value() / (8.0 * std::numbers::pi * a * b);
        Hx = common_part * ((y + b) * (tetha1 - tetha2) - (y - b) * (tetha4 - tetha3) + (x + a) * log(r2 / r3) - (x - a) * log(r1 / r4));

        Hy = -common_part * ((x + a) * (tetha2 - tetha3) - (x - a) * (tetha1 - tetha4) + (y + b) * log(r2 / r1) - (y - b) * log(r3 / r4));
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }


    if (inducingFieldPoint.get_rotation()) {
        double modulo = hypot(Hx, Hy);
        double currentAngle = atan2(Hy, Hx);
        if (currentAngle < 0) {
            currentAngle += 2 * std::numbers::pi;
        }
        double turnAngle = inducingFieldPoint.get_rotation().value() / 180 * std::numbers::pi;
        double totalAngle = currentAngle + turnAngle;
        Hx = modulo * cos(totalAngle);
        Hy = modulo * sin(totalAngle);
        if (std::isnan(Hx) || std::isnan(Hy)) {
            throw NaNResultException("NaN found in Binns Lawrenson's model for magnetic field");
        }
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;   
}

ComplexFieldPoint MagneticFieldStrengthLammeranerModel::get_magnetic_field_strength_between_two_points(FieldPoint inducingFieldPoint, FieldPoint inducedFieldPoint, std::optional<size_t> inducingWireIndex) {
    double Hx;
    double Hy;

    double turnLength = 1;
    if (inducingFieldPoint.get_turn_length()) {
        turnLength = inducingFieldPoint.get_turn_length().value();
    }
    double distance = hypot(inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1], inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0]);
    double angle = atan2(inducedFieldPoint.get_point()[0] - inducingFieldPoint.get_point()[0], inducedFieldPoint.get_point()[1] - inducingFieldPoint.get_point()[1]);

    if (inducingWireIndex) {
        if (_wirePerWinding[inducingWireIndex.value()].get_type() != WireType::ROUND && _wirePerWinding[inducingWireIndex.value()].get_type() != WireType::LITZ) {
            return MagneticFieldStrengthBinnsLawrensonModel().get_magnetic_field_strength_between_two_points(inducingFieldPoint, inducedFieldPoint);
        }
        if (distance < _wirePerWinding[inducingWireIndex.value()].get_maximum_outer_width() / 2) {
            Hx = 0;
            Hy = 0;
        }
        else {
            double ex = cos(angle - std::numbers::pi / 2);
            double ey = sin(angle - std::numbers::pi / 2);
            double magneticFiledStrengthModule = -inducingFieldPoint.get_value() / 2 / std::numbers::pi / distance * turnLength / hypot(turnLength, distance);
            Hx = magneticFiledStrengthModule * ex;
            Hy = magneticFiledStrengthModule * ey;
        }
    }
    else {
        double ex = cos(angle - std::numbers::pi / 2);
        double ey = sin(angle - std::numbers::pi / 2);
        double magneticFiledStrengthModule = -inducingFieldPoint.get_value() / 2 / std::numbers::pi / distance * turnLength / hypot(turnLength, distance);
        Hx = magneticFiledStrengthModule * ex;
        Hy = magneticFiledStrengthModule * ey;
    }

    if (std::isnan(Hx) || std::isnan(Hy)) {
        throw NaNResultException("NaN found in Lammeraner's model for magnetic field");
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_turn_length()) {
        complexFieldPoint.set_turn_length(inducedFieldPoint.get_turn_length().value());
    }
    return complexFieldPoint;
}

FieldPoint MagneticFieldStrengthAlbachModel::get_equivalent_inducing_point_for_gap(CoreGap gap, double magneticFieldStrengthGap) {
    if (!gap.get_section_dimensions()) {
        throw GapException("Gap is missing section dimensions");
    }
    if (!gap.get_coordinates()) {
        throw GapException("Gap is missing coordinates");
    }
    double rc = gap.get_section_dimensions().value()[0] / 2;
    double xi = gap.get_length() / (2 * rc);
    double x = 1 - 1.05 * xi - 2.88 * pow(xi, 2) - 8.8 * pow(xi, 3);
    if (x < 0) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something went wrong with Albach method with x");
    }
    double current = (magneticFieldStrengthGap * gap.get_length()) / (0.25 - 1.569 * xi + 4.34 * pow(xi, 2) - 7.042 * pow(xi, 3));
    double eta = x * rc;

    if (eta > gap.get_section_dimensions().value()[0] / 2) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something went wrong with Albach method with eta");
    }
    FieldPoint fieldPoint;
    // Position the equivalent wire at distance eta from the center leg axis
    // For center leg gaps (x=0), the wire is placed at radius eta
    // For lateral gaps (x>0 or x<0), adjust based on gap position
    if (gap.get_coordinates().value()[0] > 0) {
        // Gap on positive x side - wire is at gap_x - eta (closer to center)
        fieldPoint.set_point({gap.get_coordinates().value()[0] - eta, gap.get_coordinates().value()[1]});
    }
    else if (gap.get_coordinates().value()[0] < 0) {
        // Gap on negative x side - wire is at gap_x + eta (closer to center)
        fieldPoint.set_point({gap.get_coordinates().value()[0] + eta, gap.get_coordinates().value()[1]});
    }
    else {
        // Center leg gap (x=0) - wire is placed at radius eta from axis
        fieldPoint.set_point({eta, gap.get_coordinates().value()[1]});
    }
    fieldPoint.set_value(current);
    return fieldPoint;
}

ComplexFieldPoint MagneticFieldStrengthRoshenModel::get_magnetic_field_strength_between_gap_and_point(CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) {
    double distanceFromCenterEdgeGapX;
    if (gap.get_coordinates().value()[0] == 0) {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] + gap.get_section_dimensions().value()[0] / 2);
    }
    else {
        distanceFromCenterEdgeGapX = inducedFieldPoint.get_point()[0] - (gap.get_coordinates().value()[0] - gap.get_section_dimensions().value()[0] / 2);
    }
    double distanceFromCenterEdgeGapY = inducedFieldPoint.get_point()[1] - gap.get_coordinates().value()[1];
    double halfGapLength = gap.get_length() / 2;

    double magneticIntensityXDividend = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY - halfGapLength, 2);
    double magneticIntensityXDivisor = pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY + halfGapLength, 2);
    double Hx = -0.9 * magneticFieldStrengthGap / 2 / std::numbers::pi * log(magneticIntensityXDividend / magneticIntensityXDivisor);

    double m;
    if (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) > pow(halfGapLength, 2)) {
        m = 0;
    }
    else {
        m = 1;
    }

    double x = distanceFromCenterEdgeGapX * halfGapLength / (pow(distanceFromCenterEdgeGapX, 2) + pow(distanceFromCenterEdgeGapY, 2) - pow(halfGapLength, 2));
    double Hy = -0.9 * magneticFieldStrengthGap / std::numbers::pi * (atan(x) + m * std::numbers::pi);

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    complexFieldPoint.set_real(Hx);
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}


// ============================================================================
// MagneticFieldStrengthAlbach2DModel Implementation (Air Coil / Biot-Savart)
// ============================================================================

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::get_magnetic_field_strength_between_two_points(
    FieldPoint inducingFieldPoint, 
    FieldPoint inducedFieldPoint, 
    std::optional<size_t> inducingWireIndex
) {
    // ALBACH model calculates field from all turns at once via calculateTotalFieldAtPoint()
    // This per-turn-pair method should not be called
    (void)inducingFieldPoint;
    (void)inducedFieldPoint;
    (void)inducingWireIndex;
    throw std::runtime_error("ALBACH model does not support per-turn-pair field calculation. Use calculateTotalFieldAtPoint() instead.");
}

std::pair<double, double> MagneticFieldStrengthAlbach2DModel::calculateMagneticField(double r, double z) const {
    // Calculate H directly using analytical Biot-Savart formulas for circular filaments
    // Uses complete elliptic integrals of the first and second kind
    
    // Handle r near zero to avoid division by zero
    if (r < 1e-10) {
        return {0.0, 0.0};
    }
    
    double H_r_total = 0.0;
    double H_z_total = 0.0;
    
    for (const auto& turn : _turns) {
        double I = turn.current;
        if (std::abs(I) < 1e-15) continue;
        
        if (turn.isRectangular()) {
            // Rectangular conductor: use filamentary subdivision
            double width = turn.width;
            double height = turn.height;
            int numR = 3, numZ = 3;
            double dI = I / (numR * numZ);
            
            for (int ir = 0; ir < numR; ++ir) {
                for (int iz = 0; iz < numZ; ++iz) {
                    double fr = (ir + 0.5) / numR;
                    double fz = (iz + 0.5) / numZ;
                    
                    double rf = turn.r - width/2 + width * fr;
                    double zf = turn.z - height/2 + height * fz;
                    
                    if (rf < 1e-10) continue;
                    
                    // Direct Biot-Savart calculation for circular filament
                    double deltaZ = z - zf;
                    double sumR = r + rf;
                    double diffR = r - rf;
                    
                    double denom = sumR * sumR + deltaZ * deltaZ;
                    if (denom < 1e-20) continue;
                    
                    double k2 = 4 * r * rf / denom;
                    double k = std::sqrt(k2);
                    if (k > 0.999999) k = 0.999999;
                    
                    if (k > 1e-10) {
                        double K_k = std::comp_ellint_1(k);
                        double E_k = std::comp_ellint_2(k);
                        
                        double sqrtDenom = std::sqrt(denom);
                        double denomDiffR = diffR * diffR + deltaZ * deltaZ;
                        
                        if (denomDiffR > 1e-20) {
                            double prefactor = dI / (2 * std::numbers::pi);
                            
                            H_r_total += prefactor * deltaZ / (r * sqrtDenom) * 
                                  (-K_k + E_k * (rf*rf + r*r + deltaZ*deltaZ) / denomDiffR);
                            
                            H_z_total += prefactor / sqrtDenom * 
                                  (K_k + E_k * (rf*rf - r*r - deltaZ*deltaZ) / denomDiffR);
                        }
                    }
                }
            }
        } else {
            // Round wire: single filament
            double r0 = turn.r;
            double z0 = turn.z;
            
            double deltaZ = z - z0;
            double sumR = r + r0;
            double diffR = r - r0;
            
            double denom = sumR * sumR + deltaZ * deltaZ;
            if (denom < 1e-20) continue;
            
            double k2 = 4 * r * r0 / denom;
            double k = std::sqrt(k2);
            if (k > 0.999999) k = 0.999999;
            
            if (k > 1e-10) {
                double K_k = std::comp_ellint_1(k);
                double E_k = std::comp_ellint_2(k);
                
                double sqrtDenom = std::sqrt(denom);
                double denomDiffR = diffR * diffR + deltaZ * deltaZ;
                
                if (denomDiffR > 1e-20) {
                    double prefactor = I / (2 * std::numbers::pi);
                    
                    H_r_total += prefactor * deltaZ / (r * sqrtDenom) * 
                          (-K_k + E_k * (r0*r0 + r*r + deltaZ*deltaZ) / denomDiffR);
                    
                    H_z_total += prefactor / sqrtDenom * 
                          (K_k + E_k * (r0*r0 - r*r - deltaZ*deltaZ) / denomDiffR);
                }
            }
        }
    }
    
    return {H_r_total, H_z_total};
}

ComplexFieldPoint MagneticFieldStrengthAlbach2DModel::calculateTotalFieldAtPoint(FieldPoint inducedFieldPoint) {
    // Extract induced point coordinates - in 2D cross section: [0] = x (radial), [1] = y (axial)
    double r = std::abs(inducedFieldPoint.get_point()[0]);
    double z = inducedFieldPoint.get_point()[1];
    
    // Calculate H field from all turns
    auto [H_r, H_z] = calculateMagneticField(r, z);
    
    // Convert to 2D Cartesian: real = radial (Hx), imaginary = axial (Hy)
    ComplexFieldPoint result;
    result.set_real(H_r);
    result.set_imaginary(H_z);
    result.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_turn_index()) {
        result.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    if (inducedFieldPoint.get_label()) {
        result.set_label(inducedFieldPoint.get_label().value());
    }
    
    return result;
}

void MagneticFieldStrengthAlbach2DModel::setupFromMagnetic(
    Magnetic magnetic, 
    const std::vector<Wire>& wirePerWinding
) {
    // Get turns description
    if (!magnetic.get_coil().get_turns_description()) {
        throw std::runtime_error("Missing turns description in coil");
    }
    auto turns = magnetic.get_coil().get_turns_description().value();
    
    // Set up all turns from the coil
    _turns.clear();
    for (size_t turnIdx = 0; turnIdx < turns.size(); ++turnIdx) {
        auto& turn = turns[turnIdx];
        AlbachTurnPosition albachTurn;
        
        // In 2D cross-section, x = radial, y = axial
        albachTurn.r = std::abs(turn.get_coordinates()[0]);
        albachTurn.z = turn.get_coordinates()[1];
        albachTurn.current = 1.0; // Will be scaled per harmonic
        albachTurn.turnIndex = turnIdx;

        // Get wire info for this turn to set dimensions for rectangular wires
        auto windingIndex = magnetic.get_mutable_coil().get_winding_index_by_name(turn.get_winding());
        if (windingIndex < wirePerWinding.size()) {
            auto& wire = wirePerWinding[windingIndex];
            if (wire.get_type() != WireType::ROUND && wire.get_type() != WireType::LITZ) {
                // Rectangular, foil, or planar wire - set dimensions for subdivision
                if (wire.get_conducting_width()) {
                    albachTurn.width = resolve_dimensional_values(wire.get_conducting_width().value());
                }
                if (wire.get_conducting_height()) {
                    albachTurn.height = resolve_dimensional_values(wire.get_conducting_height().value());
                }
            }
            // For round/litz wires, width and height stay at 0 (point filament)
        }
        _turns.push_back(albachTurn);
    }
}



// ============================================================================
// MagneticFieldStrengthSullivanModel Implementation
// (2D Image Method / Biot-Savart for gap fringing field)
// ============================================================================
//
// THEORY (from Sullivan's shapeopt MATLAB code):
// -------
// The air gap is modeled as a set of current filaments distributed along the
// gap length. For each filament at position R_gap:
//   - A "cross" current (+I_per_div, into page) is placed at the gap face
//   - A "dot" current (-I_per_div, out of page) is placed at the mirror
//     position (reflected about x=0 for center gaps)
//
// The winding window (width bw, height hw) is the fundamental unit cell.
// Image copies are tiled in both x and y:
//   x: at x_center + n * 2*hw,  n in [-imageUnitsX, +imageUnitsX]
//   y: at y_center + m * bw,    m in [-imageUnitsY, +imageUnitsY]
//
// Total field at point P is superposition of all image contributions:
//   B(P) = sum (mu_0*I)/(2*pi) * (P - R_fil) / |P - R_fil|^2
// Then H = B / mu_0
//
// A 0.9 attenuation factor is applied (same as Roshen) to account for
// the fraction of MMF that produces external fringing field.
//
// MAPPING FROM MATLAB CODE:
// -------------------------
// In the original shapeopt code (function Bfinite):
//   - pvec(1:2) = unit_of_X, unit_of_Y  -> _imageUnitsX, _imageUnitsY
//   - pvec(3:4) = bw, hw                -> estimated from gap geometry
//   - pvec(7:8) = gw, gap_div           -> gap.get_length(), _gapDivisions
//   - pvec(9)   = I_per_gap_div         -> I_total / _gapDivisions
//   - pvec(10)  = Rgbase(k)             -> filament position (complex)
//   - const1 = u0*(j)*I/(2*pi)          -> Biot-Savart coefficient
//   - R1_1/R1_2: cross/dot current pair
//   - center_matrix: image unit centers
//

ComplexFieldPoint MagneticFieldStrengthSullivanModel::get_magnetic_field_strength_between_gap_and_point(
    CoreGap gap, double magneticFieldStrengthGap, FieldPoint inducedFieldPoint) {

    if (!gap.get_section_dimensions()) {
        throw GapException("Gap is missing section dimensions");
    }
    if (!gap.get_coordinates()) {
        throw GapException("Gap is missing coordinates");
    }

    // ---- Extract gap geometry ----
    double gapLength = gap.get_length();
    double gapX = gap.get_coordinates().value()[0];
    double gapY = gap.get_coordinates().value()[1];
    double columnWidth = gap.get_section_dimensions().value()[0];

    // ---- Estimate winding window dimensions ----
    // bw: window breadth in y-direction (along the gap), approx = column width
    // hw: window height in x-direction (perpendicular to gap)
    // For center gap (gapX=0): hw ~ columnWidth (the window extends from
    //   the centerpost edge outward)
    // For lateral gap: hw ~ 2*|gapX|
    double bw = columnWidth;
    double hw;
    if (std::abs(gapX) > 1e-10) {
        hw = 2.0 * std::abs(gapX);
    } else {
        hw = columnWidth;
    }

    // ---- Compute total gap current (Ampere's law: NI = H_gap * gapLength) ----
    double I_total = magneticFieldStrengthGap * gapLength;
    double I_per_div = I_total / _gapDivisions;

    // ---- Empirical attenuation (consistent with Roshen) ----
    double attenuationFactor = 0.9;

    // ---- Gap filament spacing ----
    double gapGrid = gapLength / _gapDivisions;

    // ---- Point of interest ----
    double xP = inducedFieldPoint.get_point()[0];
    double yP = inducedFieldPoint.get_point()[1];

    // Accumulate B field components
    double Bx_total = 0.0;
    double By_total = 0.0;

    double u0 = Constants().vacuumPermeability;

    // ---- For each gap filament ----
    for (int gapIdx = 0; gapIdx < _gapDivisions; ++gapIdx) {
        // Y-position of this filament relative to gap center
        double filY;
        if (_gapDivisions == 1) {
            filY = 0.0;
        } else {
            filY = -(gapLength - gapGrid) / 2.0 + gapIdx * gapGrid;
        }

        // Absolute position of the "cross" current (into page)
        // For center-leg gap (gapX~0): place at the centerpost edge
        double crossX, crossY;
        if (std::abs(gapX) < 1e-10) {
            crossX = 0.0;
            crossY = gapY + filY;
        } else {
            crossX = gapX;
            crossY = gapY + filY;
        }

        // Mirror image "dot" current (out of page): reflected about x=0
        double dotX = -crossX;
        double dotY = crossY;

        // ---- Sum over all image units ----
        // This is the core of the method of images from Sullivan's code:
        // center_matrix = ones(size(y_temp))' * x_temp + j*y_temp' * ones(size(x_temp))
        // where x_temp = 2*hw * linspace(-unit_of_X, unit_of_X, ...)
        //       y_temp = bw  * linspace(-unit_of_Y, unit_of_Y, ...)
        for (int nx = -_imageUnitsX; nx <= _imageUnitsX; ++nx) {
            for (int ny = -_imageUnitsY; ny <= _imageUnitsY; ++ny) {
                double unitCenterX = nx * 2.0 * hw;
                double unitCenterY = ny * bw;

                // "Cross" current position in this image unit
                double srcCrossX = crossX + unitCenterX;
                double srcCrossY = crossY + unitCenterY;

                // "Dot" current position in this image unit
                double srcDotX = dotX + unitCenterX;
                double srcDotY = dotY + unitCenterY;

                // Biot-Savart for "cross" current (INTO page, +z direction)
                // B_x = +(mu_0*I)/(2*pi) * dy/r^2
                // B_y = -(mu_0*I)/(2*pi) * dx/r^2
                {
                    double dx = xP - srcCrossX;
                    double dy = yP - srcCrossY;
                    double r2 = dx * dx + dy * dy;
                    // Avoid division by zero (same approach as MATLAB code:
                    // Rp_abs = ((Rp_abs == 0) + Rp_abs); )
                    if (r2 < 1e-30) r2 = 1.0;

                    double coeff = u0 * I_per_div / (2.0 * std::numbers::pi * r2);
                    Bx_total += coeff * dy;
                    By_total -= coeff * dx;
                }

                // Biot-Savart for "dot" current (OUT OF page, -z direction)
                // Opposite sign current
                {
                    double dx = xP - srcDotX;
                    double dy = yP - srcDotY;
                    double r2 = dx * dx + dy * dy;
                    if (r2 < 1e-30) r2 = 1.0;

                    double coeff = u0 * (-I_per_div) / (2.0 * std::numbers::pi * r2);
                    Bx_total += coeff * dy;
                    By_total -= coeff * dx;
                }
            }
        }
    }

    // Convert B to H: H = B / mu_0
    double Hx = attenuationFactor * Bx_total / u0;
    double Hy = attenuationFactor * By_total / u0;

    if (std::isnan(Hx) || std::isnan(Hy)) {
        throw NaNResultException("NaN found in Sullivan's fringing field model");
    }

    ComplexFieldPoint complexFieldPoint;
    complexFieldPoint.set_real(Hx);
    complexFieldPoint.set_imaginary(Hy);
    complexFieldPoint.set_point(inducedFieldPoint.get_point());
    if (inducedFieldPoint.get_turn_index()) {
        complexFieldPoint.set_turn_index(inducedFieldPoint.get_turn_index().value());
    }
    return complexFieldPoint;
}


} // namespace OpenMagnetics
