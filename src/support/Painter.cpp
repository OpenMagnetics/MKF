#include "physical_models/MagneticField.h"
#include <numbers>
#include "physical_models/StrayCapacitance.h"
#include "physical_models/LeakageInductance.h"
#include "support/Painter.h"
#include "support/CoilMesher.h"
#include "constructive_models/Coil.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <cfloat>
#include <chrono>
#include <thread>
#include <filesystem>
#include <set>
#include <cmath>
#include <iostream>
#include "support/Exceptions.h"


namespace OpenMagnetics {

static size_t find_dominant_non_dc_harmonic_index(const Harmonics& harmonics) {
    size_t dominantIndex = 1;
    double maxMetric = 0.0;
    for (size_t i = 1; i < harmonics.get_amplitudes().size(); ++i) {
        double metric = harmonics.get_amplitudes()[i] * std::sqrt(harmonics.get_frequencies()[i]);
        if (metric > maxMetric) {
            maxMetric = metric;
            dominantIndex = i;
        }
    }
    return dominantIndex;
}

static void resolve_harmonic_index_for_painting(const Harmonics& harmonics, size_t& harmonicIndex) {
    if (harmonicIndex != 1 || harmonics.get_amplitudes().size() <= 2) {
        return;
    }
    double metric1 = harmonics.get_amplitudes()[1] * std::sqrt(harmonics.get_frequencies()[1]);
    size_t dominant = find_dominant_non_dc_harmonic_index(harmonics);
    double metricDominant = harmonics.get_amplitudes()[dominant] * std::sqrt(harmonics.get_frequencies()[dominant]);
    if (metricDominant > 0 && metric1 / metricDominant < 0.01) {
        harmonicIndex = dominant;
    }
}

ComplexField PainterInterface::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    if (magnetic.get_core().get_shape_family() == MAS::CoreShapeFamily::T) {
        return calculate_toroidal_magnetic_field(operatingPoint, magnetic, harmonicIndex);
    }
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics()) {
            auto current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();
            if (!current.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from current");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(current.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(current.get_waveform().value(), sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            current.set_harmonics(harmonics);
            if (!current.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                current.set_processed(processed);
            }
            else {
                if (!current.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    current.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_current(current);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    resolve_harmonic_index_for_painting(harmonics, harmonicIndex);
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    int mirroringDimension = settings.get_painter_mirroring_dimension();  // int (0/1/2/3 mirroring planes), was truncated through bool

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, true, true, true).first;

    auto modelOverride = settings.get_painter_magnetic_field_strength_model();
    // Use painter override if set, otherwise use the simulation magnetic field strength model
    auto magneticFieldModel = modelOverride.value_or(settings.get_magnetic_field_strength_model());
    // Use the configured fringing effect model from settings
    auto fringingEffectModel = settings.get_magnetic_field_strength_fringing_effect_model();
    MagneticField magneticField(magneticFieldModel, fringingEffectModel);
    // RAII: restore the global magnetic-field settings on scope exit. They used to be
    // overwritten permanently, leaking the painter's fringing/mirroring choice into
    // every later physics computation in the process.
    SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing, &Settings::set_magnetic_field_include_fringing, includeFringing);
    SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension, &Settings::set_magnetic_field_mirroring_dimension, mirroringDimension);
    ComplexField field;
    {
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

    }
    auto turns = magnetic.get_coil().get_turns_description().value();

    if (turns[0].get_additional_coordinates()) {
        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
            if (turns[turnIndex].get_additional_coordinates()) {
                turns[turnIndex].set_coordinates(turns[turnIndex].get_additional_coordinates().value()[0]);
            }
        }
        magnetic.get_mutable_coil().set_turns_description(turns);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        auto additionalField = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        for (size_t pointIndex = 0; pointIndex < field.get_data().size(); ++pointIndex) {
            field.get_mutable_data()[pointIndex].set_real(field.get_mutable_data()[pointIndex].get_real() + additionalField.get_mutable_data()[pointIndex].get_real());
            field.get_mutable_data()[pointIndex].set_imaginary(field.get_mutable_data()[pointIndex].get_imaginary() + additionalField.get_mutable_data()[pointIndex].get_imaginary());
        }
    }


    return field;
}

ComplexField PainterInterface::calculate_toroidal_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    // The ring-plane field of the toroidal leakage model (LeakageInductance.h, "Toroidal cores"; white
    // paper section 9.3): per-region Kelvin images with the net-current term, built and evaluated by the
    // same code that gives the leakage energy, so that the picture is the field whose energy is the
    // model's L_DM. Each winding is driven by its current phasor at the painted harmonic, signed by the
    // MAS direction convention (CoilMesher::calculate_current_direction_per_winding). Like the other cores'
    // painter field (MagneticField's in-phase output), each point holds the in-phase field, i.e. the field of
    // the peak currents' components in phase with the gauge winding (CoilMesher::calculate_current_phase_per_winding),
    // as real = Hx and imaginary = Hy.
    if (magnetic.get_core().get_shape_family() != MAS::CoreShapeFamily::T) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "calculate_toroidal_magnetic_field needs a toroidal core");
    }
    auto coil = magnetic.get_coil();
    size_t numberWindings = coil.get_functional_description().size();
    if (operatingPoint.get_excitations_per_winding().size() != numberWindings) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Toroidal field: " + std::to_string(operatingPoint.get_excitations_per_winding().size()) +
                                    " excitations for " + std::to_string(numberWindings) + " windings");
    }
    auto excitationFrequency = operatingPoint.get_excitations_per_winding()[0].get_frequency();
    for (const auto& excitation : operatingPoint.get_excitations_per_winding()) {
        if (!excitation.get_current() || !excitation.get_current()->get_waveform()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Toroidal field: every winding needs a current waveform");
        }
        if (excitation.get_frequency() != excitationFrequency) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal field: the windings are excited at different frequencies");
        }
    }
    // Every winding's harmonics, from its waveform when MAS does not list them, so that the phases come from
    // the same gauge as every other core's painted field (CoilMesher::calculate_current_phase_per_winding).
    auto excitations = operatingPoint.get_excitations_per_winding();
    for (auto& excitation : excitations) {
        auto current = excitation.get_current().value();
        if (!current.get_harmonics()) {
            auto sampledWaveform = Inputs::calculate_sampled_waveform(current.get_waveform().value(), excitationFrequency);
            current.set_harmonics(Inputs::calculate_harmonics_data(current.get_waveform().value(), sampledWaveform, excitationFrequency));
            excitation.set_current(current);
        }
    }
    operatingPoint.set_excitations_per_winding(excitations);
    Harmonics harmonics = excitations[0].get_current()->get_harmonics().value();
    resolve_harmonic_index_for_painting(harmonics, harmonicIndex);
    if (harmonicIndex >= harmonics.get_frequencies().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal field: harmonic " + std::to_string(harmonicIndex) + " is not in the current");
    }
    double frequency = harmonics.get_frequencies()[harmonicIndex];

    // In-phase current of winding w: A_w cos(phi_w - phi_gauge), signed by its MAS direction.
    auto directions = CoilMesher::calculate_current_direction_per_winding(coil);
    auto phases = CoilMesher::calculate_current_phase_per_winding(coil, operatingPoint, {harmonicIndex})[0];
    std::vector<double> realCurrents(numberWindings);
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        const auto amplitudes = excitations[windingIndex].get_current()->get_harmonics()->get_amplitudes();
        double amplitude = harmonicIndex < amplitudes.size() ? amplitudes[harmonicIndex] : 0.0;  // a winding without this harmonic carries none of it
        realCurrents[windingIndex] = amplitude * std::cos(phases[windingIndex]) * static_cast<double>(directions[windingIndex]);
    }
    auto conductors = LeakageInductance::calculate_toroidal_ring_plane_conductors(magnetic, realCurrents);

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field grid = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, true, true, true).first;
    ComplexField field;
    field.set_frequency(frequency);
    std::vector<ComplexFieldPoint> data;
    data.reserve(grid.get_data().size());
    for (const auto& gridPoint : grid.get_data()) {
        const auto& point = gridPoint.get_point();
        auto h = LeakageInductance::calculate_toroidal_ring_plane_field(conductors, point[0], point[1]);
        ComplexFieldPoint datum;
        datum.set_point(point);
        datum.set_real(h[0]);
        datum.set_imaginary(h[1]);
        data.push_back(datum);
    }
    field.set_data(data);
    return field;
}

// ==================== Electric energy painting ====================

PainterInterface::ElectricEnergyPainting PainterInterface::calculate_electric_energy_painting(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    if (!magnetic.get_core().get_processed_description()) {
        magnetic.get_mutable_core().process_data();
    }
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Electric energy painting: the coil has no turns description");
    }
    auto turns = coil.get_turns_description().value();
    auto windings = coil.get_functional_description();
    auto excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.size() < windings.size()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Electric energy painting: the operating point has " + std::to_string(excitations.size()) +
                                    " excitations for a coil with " + std::to_string(windings.size()) + " windings");
    }

    // The harmonic, from the first winding's voltage.
    if (!excitations[0].get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Electric energy painting: the first winding's excitation has no voltage");
    }
    auto firstVoltage = excitations[0].get_voltage().value();
    Harmonics firstHarmonics;
    if (firstVoltage.get_harmonics()) {
        firstHarmonics = firstVoltage.get_harmonics().value();
    }
    else if (firstVoltage.get_waveform()) {
        auto sampledWaveform = Inputs::calculate_sampled_waveform(firstVoltage.get_waveform().value(), excitations[0].get_frequency());
        firstHarmonics = Inputs::calculate_harmonics_data(firstVoltage.get_waveform().value(), sampledWaveform, excitations[0].get_frequency());
    }
    else {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Electric energy painting: the first winding's voltage has neither waveform nor harmonics");
    }
    if (harmonicIndex >= firstHarmonics.get_frequencies().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Electric energy painting: harmonic " + std::to_string(harmonicIndex) +
                                    " requested, the voltage has " + std::to_string(firstHarmonics.get_frequencies().size()));
    }
    const double harmonicFrequency = firstHarmonics.get_frequencies()[harmonicIndex];

    // Rms phasor of every winding's voltage at that harmonic.
    std::vector<std::complex<double>> phasorPerWinding;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        if (!excitations[windingIndex].get_voltage()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Electric energy painting: winding " + windings[windingIndex].get_name() + " has no voltage");
        }
        auto voltage = excitations[windingIndex].get_voltage().value();
        if (voltage.get_waveform()) {
            phasorPerWinding.push_back(CoilMesher::calculate_harmonic_phasor(voltage.get_waveform().value(), excitations[windingIndex].get_frequency(), harmonicFrequency) / std::sqrt(2.0));
        }
        else if (windings.size() == 1 && voltage.get_harmonics()) {
            // One winding: its phase is the reference, the amplitude is all that matters.
            phasorPerWinding.push_back(std::complex<double>(firstHarmonics.get_amplitudes()[harmonicIndex] / std::sqrt(2.0), 0.0));
        }
        else {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Electric energy painting: winding " + windings[windingIndex].get_name() +
                                        " has no voltage waveform, and with more than one winding the relative phases are needed");
        }
    }

    // Each winding's linear divider: 1 at its start terminal, 0 at its end.
    std::map<std::string, double> unitVoltagePerWinding;
    for (auto& winding : windings) {
        unitVoltagePerWinding[winding.get_name()] = 1.0;
    }
    auto dividerOutput = StrayCapacitance::calculate_voltages_per_turn(coil, unitVoltagePerWinding);
    if (!dividerOutput.get_voltage_per_turn()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Electric energy painting: no voltage per turn");
    }
    auto dividerPerTurn = dividerOutput.get_voltage_per_turn().value();
    if (dividerPerTurn.size() != turns.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Electric energy painting: voltage divider and turns differ in size");
    }
    std::vector<std::complex<double>> turnPotentials(turns.size());
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        turnPotentials[turnIndex] = phasorPerWinding[coil.get_winding_index_by_name(turns[turnIndex].get_winding())] * dividerPerTurn[turnIndex];
    }

    // The core node: held by its electrical reference, or floating.
    std::optional<std::complex<double>> fixedCorePotential;
    auto coreElectricalReference = magnetic.get_core_electrical_reference();
    auto unitCorePotential = StrayCapacitance::resolve_core_reference_potential(coil, coreElectricalReference, dividerOutput, unitVoltagePerWinding);
    if (unitCorePotential) {
        if (coreElectricalReference->get_type() == CoreElectricalReferenceType::TIED_TO_WINDING) {
            fixedCorePotential = phasorPerWinding[coil.get_winding_index_by_name(coreElectricalReference->get_winding().value())] * unitCorePotential.value();
        }
        else {
            fixedCorePotential = std::complex<double>(unitCorePotential.value(), 0.0);
        }
    }

    StrayCapacitance strayCapacitance(settings.get_stray_capacitance_model());
    ElectricEnergyPainting painting;
    auto core = magnetic.get_core();
    painting.elements = strayCapacitance.calculate_electric_energy_elements(coil, core, harmonicFrequency);
    auto distribution = StrayCapacitance::calculate_electric_energy_per_element(painting.elements, coil, core, turnPotentials, fixedCorePotential);
    painting.gapEnergy = distribution.gapEnergy;
    painting.totalEnergy = distribution.gapEnergy;

    // The painter's pixel tiling, exactly: the cells paint_electric_field draws.
    auto [cellWidth, cellHeight] = Painter::get_pixel_dimensions(magnetic);
    if (!(cellWidth > 0) || !(cellHeight > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Electric energy painting: the core gives no pixel dimensions");
    }
    painting.cellWidth = cellWidth;
    painting.cellHeight = cellHeight;
    const size_t numberPointsX = settings.get_painter_number_points_x();
    const size_t numberPointsY = settings.get_painter_number_points_y();
    double originX;
    const double originY = -0.5 * numberPointsY * cellHeight;
    if (magnetic.get_core().get_shape_family() == CoreShapeFamily::T) {
        originX = -0.5 * numberPointsX * cellWidth;
    }
    else {
        originX = magnetic.get_mutable_core().get_columns()[0].get_width() / 2;
    }
    std::vector<double> energyPerLengthPerCell(numberPointsX * numberPointsY, 0.0);
    const double samplingSpacing = std::min(cellWidth, cellHeight);

    for (size_t elementIndex = 0; elementIndex < painting.elements.size(); ++elementIndex) {
        const auto& element = painting.elements[elementIndex];
        if (!(element.length > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Electric energy painting: element " + std::to_string(elementIndex) + " has no length");
        }
        double energy = distribution.energyPerElement[elementIndex];
        painting.totalEnergy += energy;
        double energyPerLength = energy / element.length;
        painting.energyPerLength.push_back(energyPerLength);
        if (!element.inPlane) {
            painting.outOfPlaneEnergyPerLength += energyPerLength;
            continue;
        }
        painting.totalEnergyPerLengthInPlane += energyPerLength;
        if (energyPerLength == 0) {
            continue;
        }
        painting.unplacedEnergyPerLength += StrayCapacitance::sample_electric_energy_element(element, energyPerLength, samplingSpacing,
            [&](double x, double y, double sampleEnergyPerLength) {
                double column = std::floor((x - originX) / cellWidth);
                double row = std::floor((y - originY) / cellHeight);
                if (column < 0 || row < 0 || column >= static_cast<double>(numberPointsX) || row >= static_cast<double>(numberPointsY)) {
                    painting.outsideGridEnergyPerLength += sampleEnergyPerLength;
                    return;
                }
                energyPerLengthPerCell[static_cast<size_t>(row) * numberPointsX + static_cast<size_t>(column)] += sampleEnergyPerLength;
                painting.placedEnergyPerLength += sampleEnergyPerLength;
            });
    }

    std::vector<FieldPoint> data;
    data.reserve(energyPerLengthPerCell.size());
    const double cellArea = cellWidth * cellHeight;
    for (size_t row = 0; row < numberPointsY; ++row) {
        for (size_t column = 0; column < numberPointsX; ++column) {
            FieldPoint fieldPoint;
            fieldPoint.set_point({originX + (column + 0.5) * cellWidth, originY + (row + 0.5) * cellHeight});
            fieldPoint.set_value(energyPerLengthPerCell[row * numberPointsX + column] / cellArea);
            data.push_back(fieldPoint);
        }
    }
    painting.field.set_data(data);
    painting.field.set_frequency(harmonicFrequency);
    return painting;
}

Field PainterInterface::calculate_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, ElectricFieldVisualizationModel model) {
    if (model == ElectricFieldVisualizationModel::SDF_PHYSICS) {
        return calculate_electric_energy_painting(operatingPoint, magnetic, harmonicIndex).field;
    }
    // === LEGACY implementation below ===
    if (!operatingPoint.get_excitations_per_winding()[0].get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "voltage is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage()->get_harmonics()) {
            auto voltage = operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage().value();
            if (!voltage.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from voltage");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(voltage.get_waveform().value(), sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            voltage.set_harmonics(harmonics);
            if (!voltage.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                voltage.set_processed(processed);
            }
            else {
                if (!voltage.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    voltage.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_voltage(voltage);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    int mirroringDimension = settings.get_painter_mirroring_dimension();  // int (0/1/2/3 mirroring planes), was truncated through bool

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field inducedField;
    {
        // RAII (ABT #113 sweep): exception-safe replacement for the manual
        // save/set/restore around the mesh generation.
        SettingsGuard<double> mesherFactorGuard(settings, &Settings::get_coil_mesher_inside_turns_factor, &Settings::set_coil_mesher_inside_turns_factor, 1.2);
        inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, false, false, true).first;
    }

    auto strayCapacitanceModel = settings.get_stray_capacitance_model();
    StrayCapacitance strayCapacitance(strayCapacitanceModel);
    // RAII: restore the global magnetic-field settings on scope exit. They used to be
    // overwritten permanently, leaking the painter's fringing/mirroring choice into
    // every later physics computation in the process.
    SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing, &Settings::set_magnetic_field_include_fringing, includeFringing);
    SettingsGuard<int> mirroringGuard(settings, &Settings::get_magnetic_field_mirroring_dimension, &Settings::set_magnetic_field_mirroring_dimension, mirroringDimension);

    auto [pixelXDimension, pixelYDimension] = Painter::get_pixel_dimensions(magnetic);

    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();

    auto capacitanceOutput = strayCapacitance.calculate_capacitance(coil);
    auto electricEnergyAmongTurns = capacitanceOutput.get_electric_energy_among_turns().value();
    auto voltageDropAmongTurns = capacitanceOutput.get_voltage_drop_among_turns().value();

    std::set<std::pair<size_t, size_t>> turnsCombinations;
    for (size_t pointIndex = 0; pointIndex < inducedField.get_data().size(); ++pointIndex) {
        auto inducedFieldPoint = inducedField.get_data()[pointIndex];
        double fieldValue = 0;
        std::set<std::pair<std::string, std::string>> turnsCombinations;

        for (auto [firstTurnName, aux] : electricEnergyAmongTurns) {
            auto firstTurn = coil.get_turn_by_name(firstTurnName);
            for (auto [secondTurnName, energy] : aux) {
                auto key = std::make_pair(firstTurnName, secondTurnName);

                if (turnsCombinations.contains(key) || turnsCombinations.contains(std::make_pair(secondTurnName, firstTurnName))) {
                    continue;
                }
                turnsCombinations.insert(key);

                auto secondTurn = coil.get_turn_by_name(secondTurnName);
                double pixelArea = Painter::get_pixel_area_between_turns(firstTurn.get_coordinates(), firstTurn.get_dimensions().value(), firstTurn.get_cross_sectional_shape().value(), secondTurn.get_coordinates(), secondTurn.get_dimensions().value(), secondTurn.get_cross_sectional_shape().value(), inducedFieldPoint.get_point(), std::max(pixelXDimension, pixelYDimension));
                if (pixelArea > 0) {
                    // Calculate gap between turns
                    double dx = firstTurn.get_coordinates()[0] - secondTurn.get_coordinates()[0];
                    double dy = firstTurn.get_coordinates()[1] - secondTurn.get_coordinates()[1];
                    double centerDist = std::sqrt(dx*dx + dy*dy);
                    double r1 = firstTurn.get_dimensions().value()[0] / 2.0;
                    double r2 = secondTurn.get_dimensions().value()[0] / 2.0;
                    double gap = centerDist - r1 - r2;
                    if (gap < 1e-9) gap = 1e-9;  // Minimum gap to avoid division by zero
                    
                    // Get voltage drop between these turns
                    double voltageDrop = 0.0;
                    if (voltageDropAmongTurns.contains(firstTurnName) && 
                        voltageDropAmongTurns[firstTurnName].contains(secondTurnName)) {
                        voltageDrop = fabs(voltageDropAmongTurns[firstTurnName][secondTurnName]);
                    }
                    
                    if (voltageDrop > 0) {
                        // Compute energy density directly from E-field: u = 0.5 * ε * (V/gap)²
                        // This matches the SDF physics approach
                        constexpr double epsilon0 = 8.854187817e-12;
                        double E = voltageDrop / gap;
                        double energyDensityVolume = 0.5 * epsilon0 * E * E;
                        fieldValue += energyDensityVolume * pixelArea;
                    }
                }
            }
        }
        inducedField.get_mutable_data()[pointIndex].set_value(fieldValue);
    }

    return inducedField;
}


double Painter::get_pixel_area_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    return dimension * dimension * get_pixel_proportion_between_turns(firstTurnCoordinates, firstTurnDimensions, firstTurncrossSectionalShape, secondTurnCoordinates, secondTurnDimensions, secondTurncrossSectionalShape, pixelCoordinates, dimension);
}

double Painter::get_pixel_proportion_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    // auto factor = Defaults().overlappingFactorSurroundingTurns;
    auto x1 = firstTurnCoordinates[0];
    auto y1 = firstTurnCoordinates[1];
    auto x2 = secondTurnCoordinates[0];
    auto y2 = secondTurnCoordinates[1];

    if (y2 == y1 && x2 == x1) {
        return 0;
    }


    double firstTurnMaximumDimension = 0;
    double secondTurnMaximumDimension = 0;

    if (firstTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        firstTurnMaximumDimension = hypot(firstTurnDimensions[0], firstTurnDimensions[1]);
    }
    else {
        firstTurnMaximumDimension = firstTurnDimensions[0];
    }
    if (secondTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        secondTurnMaximumDimension = hypot(secondTurnDimensions[0], secondTurnDimensions[1]);
    }
    else {
        secondTurnMaximumDimension = secondTurnDimensions[0];
    }

    double semiAverageDimensionOf12 = (firstTurnMaximumDimension + secondTurnMaximumDimension) / 4;

    auto x0 = pixelCoordinates[0];
    auto y0 = pixelCoordinates[1];

    auto distanceFrom0toLine12 = fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(pow(y2 - y1, 2) + pow(x2 - x1, 2));
    auto distanceFrom0toCenter1 = hypot(x1 - x0, y1 - y0);
    auto distanceFrom0toCenter2 = hypot(x2 - x0, y2 - y0);
    auto distanceFromCenter1toCenter2 = hypot(x2 - x1, y2 - y1);

    if (distanceFrom0toCenter1 > distanceFromCenter1toCenter2 || distanceFrom0toCenter2 > distanceFromCenter1toCenter2) {
        return 0;   
    }

    double proportion;
    if (distanceFrom0toLine12 - dimension / 2 > semiAverageDimensionOf12) {
        proportion = 0;
    }
    else if (distanceFrom0toLine12 + dimension / 2 < semiAverageDimensionOf12) {
        proportion = 1;
    }
    else {
        proportion = (semiAverageDimensionOf12 - (distanceFrom0toLine12 - dimension / 2)) / dimension;
    }

    proportion *= (semiAverageDimensionOf12 - distanceFrom0toLine12) / semiAverageDimensionOf12;

    return proportion;
}

std::pair<double, double> Painter::get_pixel_dimensions(Magnetic magnetic) {

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();

    auto columns = magnetic.get_mutable_core().get_columns();
    if (columns.empty()) {
        return {0.0, 0.0};
    }
    double coreColumnHeight = columns[0].get_height();

    auto windingWindow = magnetic.get_mutable_core().get_winding_window();
    double coreWindingWindowWidth;
    double extraDimension = 1.0;
    
    if (windingWindow.get_width()) {
        // Rectangular winding window - use width directly
        coreWindingWindowWidth = windingWindow.get_width().value();
    } else {
        // Round winding window (toroidal) - use core width (diameter) as effective width
        auto processedDesc = magnetic.get_mutable_core().get_processed_description();
        if (processedDesc && processedDesc->get_width() > 0) {
            coreWindingWindowWidth = processedDesc->get_width();
            // For toroidal cores, get extraDimension to match grid spacing in CoilMesher
            extraDimension = Coil::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());
        } else {
            // Fallback: estimate from column dimensions
            coreWindingWindowWidth = columns[0].get_width() * 2.5;  // Approximate toroid diameter
        }
    }

    // Apply extraDimension for toroidal cores to match grid spacing
    double pixelXDimension = (coreWindingWindowWidth * extraDimension) / numberPointsX;
    double pixelYDimension = (coreColumnHeight * extraDimension) / numberPointsY;

    return {pixelXDimension, pixelYDimension};
}

} // namespace OpenMagnetics
