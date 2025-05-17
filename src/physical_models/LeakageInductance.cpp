#include "support/CoilMesher.h"
#include "physical_models/MagneticField.h"
#include "physical_models/LeakageInductance.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <cfloat>

namespace OpenMagnetics {

std::pair<ComplexField, double> LeakageInductance::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];
    auto numberPointsX = settings->get_magnetic_field_number_points_x();
    auto numberPointsY = settings->get_magnetic_field_number_points_y();

    auto aux = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY);
    Field inducedField = aux.first;
    double dA = aux.second;

    MagneticField magneticField;

    std::vector<int8_t> currentDirectionPerWinding;
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (windingIndex == sourceIndex) {
            currentDirectionPerWinding.push_back(1);
        }
        else if (windingIndex == destinationIndex) {
            currentDirectionPerWinding.push_back(-1);
        }
        else {
            currentDirectionPerWinding.push_back(0);
        }
    }

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
        auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
        auto windingWindows = bobbin.get_processed_description()->get_winding_windows();

        double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

        magnetic.get_mutable_coil().set_turns_description(turns);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        auto additionalField = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        for (size_t pointIndex = 0; pointIndex < field.get_data().size(); ++pointIndex) {
            if (hypot(field.get_data()[pointIndex].get_point()[0], field.get_data()[pointIndex].get_point()[1]) > windingWindowRadialHeight) {
                field.get_mutable_data()[pointIndex].set_real(additionalField.get_data()[pointIndex].get_real());
                field.get_mutable_data()[pointIndex].set_imaginary(additionalField.get_data()[pointIndex].get_imaginary());
            }
        }
    }

    return {field, dA};
}

LeakageInductanceOutput LeakageInductance::calculate_leakage_inductance(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    auto originallyIncludeFringing = settings->get_magnetic_field_include_fringing();
    settings->set_magnetic_field_include_fringing(false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }


    Processed sourceProcessed;
    sourceProcessed.set_peak_to_peak(2);
    sourceProcessed.set_duty_cycle(0.5);
    sourceProcessed.set_offset(0);
    sourceProcessed.set_rms(1 / sqrt(2));
    sourceProcessed.set_label(WaveformLabel::SINUSOIDAL);
    auto sourceWaveform = Inputs::create_waveform(sourceProcessed, frequency);
    SignalDescriptor sourceCurrent;
    sourceCurrent.set_waveform(sourceWaveform);
    sourceCurrent.set_processed(sourceProcessed);
    sourceCurrent.set_harmonics(Inputs::calculate_harmonics_data(sourceWaveform, frequency));
    OperatingPointExcitation sourceExcitation;
    sourceExcitation.set_current(sourceCurrent);


    Processed destinationProcessed = sourceProcessed;
    double sourceDestinationTurnsRatio = double(magnetic.get_mutable_coil().get_number_turns(sourceIndex)) / magnetic.get_mutable_coil().get_number_turns(destinationIndex);
    destinationProcessed.set_peak_to_peak(2.0 * sourceDestinationTurnsRatio);
    auto destinationWaveform = Inputs::create_waveform(destinationProcessed, frequency);
    SignalDescriptor destinationCurrent;
    destinationProcessed.set_rms(sourceDestinationTurnsRatio / sqrt(2.0));
    destinationCurrent.set_waveform(destinationWaveform);
    destinationCurrent.set_processed(destinationProcessed);
    destinationCurrent.set_harmonics(Inputs::calculate_harmonics_data(destinationWaveform, frequency));
    OperatingPointExcitation destinationExcitation;
    destinationExcitation.set_current(destinationCurrent);


    Processed restProcessed = sourceProcessed;
    restProcessed.set_peak_to_peak(1e-9);
    restProcessed.set_rms(1e-9);
    auto restWaveform = Inputs::create_waveform(restProcessed, frequency);
    SignalDescriptor restCurrent;
    restCurrent.set_waveform(restWaveform);
    restCurrent.set_processed(restProcessed);
    restCurrent.set_harmonics(Inputs::calculate_harmonics_data(restWaveform, frequency));
    OperatingPointExcitation restExcitation;
    restExcitation.set_current(restCurrent);

    std::vector<OperatingPointExcitation> excitationPerWinding;
    for (size_t windingIndex = 0; windingIndex <  magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (windingIndex == sourceIndex) {
            excitationPerWinding.push_back(sourceExcitation);
        }
        else if (windingIndex == destinationIndex) {
            excitationPerWinding.push_back(destinationExcitation);
        }
        else {
            excitationPerWinding.push_back(restExcitation);
        }
    }

    OperatingPoint operatingPoint;
    operatingPoint.set_excitations_per_winding(excitationPerWinding);

    auto aux = calculate_magnetic_field(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex);
    ComplexField field = aux.first;
    double dA = aux.second;

    double bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    double bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    auto coreColumnShape = magnetic.get_mutable_core().get_columns()[0].get_shape();

    double windingWindowRadialHeight = 0;
    if (bobbinWindingWindowShape == WindingWindowShape::ROUND) {
        auto windingWindows = bobbin.get_processed_description()->get_winding_windows();
        windingWindowRadialHeight = windingWindows[0].get_radial_height().value();
    }


    double vacuumPermeability = Constants().vacuumPermeability;

    double energy = 0;
    for (auto datum : field.get_data()){
        double length = 0;

        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            if (coreColumnShape == ColumnShape::ROUND) {
                length = 2 * std::numbers::pi * datum.get_point()[0];
            }
            else {
                length = 2 * std::numbers::pi * (datum.get_point()[0] - bobbinColumnWidth) + 2 * 2 * bobbinColumnWidth + 2 * 2 * bobbinColumnDepth;
            }
        }
        else {
            auto cartesianCoordinates = Coil::cartesian_to_polar(datum.get_point(), windingWindowRadialHeight);
            // Only half turn, as we are integrating also field from outside
            double radialHeightFromCenter = fabs(cartesianCoordinates[0] + bobbinColumnWidth);
            if (coreColumnShape == ColumnShape::ROUND) {
                length = std::numbers::pi * radialHeightFromCenter;
            }
            else {
                length = std::numbers::pi * (radialHeightFromCenter - bobbinColumnWidth) + 2 * bobbinColumnWidth + 2 * bobbinColumnDepth;
            }
        }

        double HeSquared = pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
        energy += 0.5 * vacuumPermeability * HeSquared * dA * length;
    }

    double currentRms = operatingPoint.get_excitations_per_winding()[sourceIndex].get_current()->get_processed()->get_rms().value();
    double leakageInductance = 2.0 / pow(currentRms, 2) * energy;
    LeakageInductanceOutput leakageInductanceOutput;

    leakageInductanceOutput.set_method_used("Energy");
    leakageInductanceOutput.set_origin(ResultOrigin::SIMULATION);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(leakageInductance);
    leakageInductanceOutput.set_leakage_inductance_per_winding({dimensionWithTolerance});

    settings->set_magnetic_field_include_fringing(originallyIncludeFringing);

    return leakageInductanceOutput;
}

ComplexField LeakageInductance::calculate_leakage_magnetic_field(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    settings->set_magnetic_field_include_fringing(false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }


    Processed sourceProcessed;
    sourceProcessed.set_peak_to_peak(2);
    sourceProcessed.set_duty_cycle(0.5);
    sourceProcessed.set_offset(0);
    sourceProcessed.set_rms(1 / sqrt(2));
    sourceProcessed.set_label(WaveformLabel::SINUSOIDAL);
    auto sourceWaveform = Inputs::create_waveform(sourceProcessed, frequency);
    SignalDescriptor sourceCurrent;
    sourceCurrent.set_waveform(sourceWaveform);
    sourceCurrent.set_processed(sourceProcessed);
    sourceCurrent.set_harmonics(Inputs::calculate_harmonics_data(sourceWaveform, frequency));
    OperatingPointExcitation sourceExcitation;
    sourceExcitation.set_current(sourceCurrent);


    Processed destinationProcessed = sourceProcessed;
    double sourceDestinationTurnsRatio = double(magnetic.get_mutable_coil().get_number_turns(sourceIndex)) / magnetic.get_mutable_coil().get_number_turns(destinationIndex);
    destinationProcessed.set_peak_to_peak(2.0 * sourceDestinationTurnsRatio);
    auto destinationWaveform = Inputs::create_waveform(destinationProcessed, frequency);
    SignalDescriptor destinationCurrent;
    destinationProcessed.set_rms(sourceDestinationTurnsRatio / sqrt(2.0));
    destinationCurrent.set_waveform(destinationWaveform);
    destinationCurrent.set_processed(destinationProcessed);
    destinationCurrent.set_harmonics(Inputs::calculate_harmonics_data(destinationWaveform, frequency));
    OperatingPointExcitation destinationExcitation;
    destinationExcitation.set_current(destinationCurrent);


    Processed restProcessed = sourceProcessed;
    restProcessed.set_peak_to_peak(1e-9);
    restProcessed.set_rms(1e-9);
    auto restWaveform = Inputs::create_waveform(restProcessed, frequency);
    SignalDescriptor restCurrent;
    restCurrent.set_waveform(restWaveform);
    restCurrent.set_processed(restProcessed);
    restCurrent.set_harmonics(Inputs::calculate_harmonics_data(restWaveform, frequency));
    OperatingPointExcitation restExcitation;
    restExcitation.set_current(restCurrent);

    std::vector<OperatingPointExcitation> excitationPerWinding;
    for (size_t windingIndex = 0; windingIndex <  magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (windingIndex == sourceIndex) {
            excitationPerWinding.push_back(sourceExcitation);
        }
        else if (windingIndex == destinationIndex) {
            excitationPerWinding.push_back(destinationExcitation);
        }
        else {
            excitationPerWinding.push_back(restExcitation);
        }
    }

    OperatingPoint operatingPoint;
    operatingPoint.set_excitations_per_winding(excitationPerWinding);

    auto aux = calculate_magnetic_field(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex);
    return aux.first;
}


LeakageInductanceOutput LeakageInductance::calculate_leakage_inductance_all_windings(Magnetic magnetic, double frequency, size_t sourceIndex, size_t harmonicIndex) {
    LeakageInductanceOutput leakageInductanceOutput;

    leakageInductanceOutput.set_method_used("Energy");
    leakageInductanceOutput.set_origin(ResultOrigin::SIMULATION);
    std::vector<DimensionWithTolerance> leakageInductancePerWinding;

    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        double leakageInductance = 0;
        if (windingIndex != sourceIndex) {
            leakageInductance = calculate_leakage_inductance(magnetic, frequency, sourceIndex, windingIndex, harmonicIndex).get_leakage_inductance_per_winding()[0].get_nominal().value();
        }

        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_nominal(leakageInductance);
        leakageInductancePerWinding.push_back(dimensionWithTolerance);
    }

    leakageInductanceOutput.set_leakage_inductance_per_winding(leakageInductancePerWinding);

    return leakageInductanceOutput;
}

} // namespace OpenMagnetics
