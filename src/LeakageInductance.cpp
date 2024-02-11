#include "MagneticField.h"
#include "LeakageInductance.h"
#include "Settings.h"
#include "MAS.hpp"
#include "Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>

namespace OpenMagnetics {

ComplexField LeakageInductance::calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto numberPointsX = settings->get_magnetic_field_number_points_x();
    auto numberPointsY = settings->get_magnetic_field_number_points_y();

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }
    double bobbinWidthStart = bobbin.get_processed_description()->get_winding_windows()[0].get_coordinates().value()[0] - bobbin.get_processed_description()->get_winding_windows()[0].get_width().value() / 2;
    double bobbinWidth = bobbin.get_processed_description()->get_winding_windows()[0].get_width().value();

    double coreColumnWidth = magnetic.get_mutable_core().get_columns()[0].get_width();
    double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    std::vector<double> bobbinPointsX = matplot::linspace(coreColumnWidth / 2, bobbinWidthStart + bobbinWidth, numberPointsX);
    std::vector<double> bobbinPointsY = matplot::linspace(-coreColumnHeight / 2, coreColumnHeight / 2, numberPointsY);
    std::vector<OpenMagnetics::FieldPoint> points;
    for (size_t j = 0; j < bobbinPointsY.size(); ++j) {
        for (size_t i = 0; i < bobbinPointsX.size(); ++i) {
            OpenMagnetics::FieldPoint fieldPoint;
            fieldPoint.set_point(std::vector<double>{bobbinPointsX[i], bobbinPointsY[j]});
            points.push_back(fieldPoint);
        }
    }

    Field inducedField;
    inducedField.set_data(points);
    inducedField.set_frequency(frequency);

    OpenMagnetics::MagneticField magneticField;

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

    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField, currentDirectionPerWinding);

    auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
    return field;
}

LeakageInductanceOutput LeakageInductance::calculate_leakage_inductance(MagneticWrapper magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto numberPointsX = settings->get_magnetic_field_number_points_x();
    auto numberPointsY = settings->get_magnetic_field_number_points_y();
    auto originallyIncludeFringing = settings->get_magnetic_field_include_fringing();
    settings->set_magnetic_field_include_fringing(false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }
    double bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    double bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto coreColumnShape = magnetic.get_mutable_core().get_columns()[0].get_shape();
    double coreWindingWindowWidth = magnetic.get_mutable_core().get_winding_windows()[0].get_width().value();
    double coreWindingWindowHeight = magnetic.get_mutable_core().get_winding_windows()[0].get_height().value();


    OpenMagnetics::Processed sourceProcessed;
    sourceProcessed.set_peak_to_peak(2);
    sourceProcessed.set_duty_cycle(0.5);
    sourceProcessed.set_offset(0);
    sourceProcessed.set_rms(1 / sqrt(2));
    sourceProcessed.set_label(WaveformLabel::SINUSOIDAL);
    auto sourceWaveform = OpenMagnetics::InputsWrapper::create_waveform(sourceProcessed, frequency);
    SignalDescriptor sourceCurrent;
    sourceCurrent.set_waveform(sourceWaveform);
    sourceCurrent.set_processed(sourceProcessed);
    sourceCurrent.set_harmonics(InputsWrapper::calculate_harmonics_data(sourceWaveform, frequency));
    OperatingPointExcitation sourceExcitation;
    sourceExcitation.set_current(sourceCurrent);


    OpenMagnetics::Processed destinationProcessed = sourceProcessed;
    double sourceDestinationTurnsRatio = double(magnetic.get_mutable_coil().get_number_turns(sourceIndex)) / magnetic.get_mutable_coil().get_number_turns(destinationIndex);
    destinationProcessed.set_peak_to_peak(2.0 * sourceDestinationTurnsRatio);
    auto destinationWaveform = OpenMagnetics::InputsWrapper::create_waveform(destinationProcessed, frequency);
    SignalDescriptor destinationCurrent;
    destinationProcessed.set_rms(sourceDestinationTurnsRatio / sqrt(2.0));
    destinationCurrent.set_waveform(destinationWaveform);
    destinationCurrent.set_processed(destinationProcessed);
    destinationCurrent.set_harmonics(InputsWrapper::calculate_harmonics_data(destinationWaveform, frequency));
    OperatingPointExcitation destinationExcitation;
    destinationExcitation.set_current(destinationCurrent);


    OpenMagnetics::Processed restProcessed = sourceProcessed;
    restProcessed.set_peak_to_peak(1e-9);
    restProcessed.set_rms(1e-9);
    auto restWaveform = OpenMagnetics::InputsWrapper::create_waveform(restProcessed, frequency);
    SignalDescriptor restCurrent;
    restCurrent.set_waveform(restWaveform);
    restCurrent.set_processed(restProcessed);
    restCurrent.set_harmonics(InputsWrapper::calculate_harmonics_data(restWaveform, frequency));
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

    ComplexField field = calculate_magnetic_field(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex);

    double dx = coreWindingWindowWidth / numberPointsX;
    double dy = coreWindingWindowHeight / numberPointsY;
    double dA = dx * dy;
    double vacuumPermeability = Constants().vacuumPermeability;

    double energy = 0;
    for (auto datum : field.get_data()){
        double length;

        if (coreColumnShape == ColumnShape::ROUND) {
            length = 2 * std::numbers::pi * datum.get_point()[0];
        }
        else {
            length = 2 * std::numbers::pi * (datum.get_point()[0] - bobbinColumnWidth) + 2 * 2 * bobbinColumnWidth + 2 * 2 * bobbinColumnDepth;
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


LeakageInductanceOutput LeakageInductance::calculate_leakage_inductance_all_windings(MagneticWrapper magnetic, double frequency, size_t sourceIndex, size_t harmonicIndex) {
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
