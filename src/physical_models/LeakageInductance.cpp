#include "support/CoilMesher.h"
#include "physical_models/MagneticField.h"
#include "physical_models/LeakageInductance.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <cfloat>
#include "support/Exceptions.h"

namespace OpenMagnetics {

std::pair<size_t, size_t> LeakageInductance::calculate_number_points_needed_for_leakage(Coil coil) {
    if (!coil.get_layers_description()) {
        throw CoilNotProcessedException("Layers description is missing");
    }
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Turns description is missing");
    }

    double minimumDistanceHorizontallyOrRadially = DBL_MAX;
    double minimumDistanceVerticallyOrAngular = DBL_MAX;

    auto layers = coil.get_layers_description().value();
    auto turns = coil.get_turns_description().value();

    for (auto layer : layers) {
        if (layer.get_dimensions()[0] > 0) {
            minimumDistanceHorizontallyOrRadially = std::min(minimumDistanceHorizontallyOrRadially, layer.get_dimensions()[0]);
        }
        if (layer.get_dimensions()[1] > 0) {
            minimumDistanceVerticallyOrAngular = std::min(minimumDistanceVerticallyOrAngular, layer.get_dimensions()[1]);
        }
    }
    for (auto turn : turns) {
        if (turn.get_dimensions().value()[0] > 0) {
            minimumDistanceHorizontallyOrRadially = std::min(minimumDistanceHorizontallyOrRadially, turn.get_dimensions().value()[0]);
        }
        if (turn.get_dimensions().value()[1] > 0) {
            minimumDistanceVerticallyOrAngular = std::min(minimumDistanceVerticallyOrAngular, turn.get_dimensions().value()[1]);
        }
    }

    auto bobbin = coil.resolve_bobbin();
    auto windingWindowsDimensions = bobbin.get_winding_window_dimensions();
    auto windingWindowShape = bobbin.get_winding_window_shape();
    size_t numberPointsX;
    size_t numberPointsY;
    if (windingWindowShape == WindingWindowShape::ROUND) {
        double angleInDistance = 2 * std::numbers::pi * windingWindowsDimensions[0] * (windingWindowsDimensions[1] / DEGREES_IN_CIRCLE);
        numberPointsX = size_t(ceil(windingWindowsDimensions[0] / minimumDistanceHorizontallyOrRadially));
        numberPointsY = size_t(ceil(angleInDistance / minimumDistanceVerticallyOrAngular));
    }
    else {
        numberPointsX = size_t(ceil(windingWindowsDimensions[0] / minimumDistanceHorizontallyOrRadially));
        numberPointsY = size_t(ceil(windingWindowsDimensions[1] / minimumDistanceVerticallyOrAngular));
    }

    return {numberPointsX, numberPointsY};
}

std::pair<ComplexField, double> LeakageInductance::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];
    auto [numberPointsX, numberPointsY] = calculate_grid_points(magnetic, frequency);

    auto meshResult = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY);
    Field inducedField = meshResult.first;

        if (inducedField.get_data().size() == 0) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Mesh generation failed: induced field data is empty");
        }
    double dA = meshResult.second;

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
        CoilMesherModels modelToUse = select_mesh_model(magnetic);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField, std::nullopt, modelToUse);
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
    auto originallyIncludeFringing = settings.get_magnetic_field_include_fringing();
    settings.set_magnetic_field_include_fringing(false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw CoilNotProcessedException("Cannot calculate leakage inductance: bobbin description has not been processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Cannot calculate leakage inductance: bobbin column width is not defined");
    }


    OperatingPoint operatingPoint = create_leakage_operating_point(magnetic, sourceIndex, destinationIndex, frequency);

    auto magneticFieldResult = calculate_magnetic_field(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex);
    ComplexField field = magneticFieldResult.first;
    double dA = magneticFieldResult.second;

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

        double magneticFieldStrengthSquared = pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
        energy += 0.5 * vacuumPermeability * magneticFieldStrengthSquared * dA * length;
    }

    double currentRms = operatingPoint.get_excitations_per_winding()[sourceIndex].get_current()->get_processed()->get_rms().value();
    double leakageInductance = 2.0 / pow(currentRms, 2) * energy;
    LeakageInductanceOutput leakageInductanceOutput;

    leakageInductanceOutput.set_method_used("Energy");
    leakageInductanceOutput.set_origin(ResultOrigin::SIMULATION);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(leakageInductance);
    leakageInductanceOutput.set_leakage_inductance_per_winding({dimensionWithTolerance});

    settings.set_magnetic_field_include_fringing(originallyIncludeFringing);

    return leakageInductanceOutput;
}

ComplexField LeakageInductance::calculate_leakage_magnetic_field(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    settings.set_magnetic_field_include_fringing(false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw CoilNotProcessedException("Cannot calculate leakage inductance: bobbin description has not been processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Cannot calculate leakage inductance: bobbin column width is not defined");
    }


    OperatingPoint operatingPoint = create_leakage_operating_point(magnetic, sourceIndex, destinationIndex, frequency);

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

OperatingPoint LeakageInductance::create_leakage_operating_point(Magnetic& magnetic, size_t sourceIndex, size_t destinationIndex, double frequency) {
    // Create source excitation (1 A RMS sinusoidal)
    Processed sourceProcessed;
    sourceProcessed.set_peak_to_peak(SINUSOIDAL_PEAK_TO_PEAK);
    sourceProcessed.set_duty_cycle(SINUSOIDAL_DUTY_CYCLE);
    sourceProcessed.set_offset(SINUSOIDAL_OFFSET);
    sourceProcessed.set_rms(1.0 / std::sqrt(SINUSOIDAL_PEAK_TO_PEAK));
    sourceProcessed.set_label(WaveformLabel::SINUSOIDAL);
    auto sourceWaveform = Inputs::create_waveform(sourceProcessed, frequency);
    SignalDescriptor sourceCurrent;
    sourceCurrent.set_waveform(sourceWaveform);
    sourceCurrent.set_processed(sourceProcessed);
    sourceCurrent.set_harmonics(Inputs::calculate_harmonics_data(sourceWaveform, frequency));
    OperatingPointExcitation sourceExcitation;
    sourceExcitation.set_current(sourceCurrent);

    // Create destination excitation with turns ratio
    Processed destinationProcessed = sourceProcessed;
    double sourceDestinationTurnsRatio = double(magnetic.get_mutable_coil().get_number_turns(sourceIndex)) /
                                         magnetic.get_mutable_coil().get_number_turns(destinationIndex);
    destinationProcessed.set_peak_to_peak(SINUSOIDAL_PEAK_TO_PEAK * sourceDestinationTurnsRatio);
    destinationProcessed.set_rms(sourceDestinationTurnsRatio / std::sqrt(SINUSOIDAL_PEAK_TO_PEAK));
    auto destinationWaveform = Inputs::create_waveform(destinationProcessed, frequency);
    SignalDescriptor destinationCurrent;
    destinationCurrent.set_waveform(destinationWaveform);
    destinationCurrent.set_processed(destinationProcessed);
    destinationCurrent.set_harmonics(Inputs::calculate_harmonics_data(destinationWaveform, frequency));
    OperatingPointExcitation destinationExcitation;
    destinationExcitation.set_current(destinationCurrent);

    // Create rest excitation with negligible current
    Processed restProcessed = sourceProcessed;
    restProcessed.set_peak_to_peak(NEGLIGIBLE_CURRENT);
    restProcessed.set_rms(NEGLIGIBLE_CURRENT);
    auto restWaveform = Inputs::create_waveform(restProcessed, frequency);
    SignalDescriptor restCurrent;
    restCurrent.set_waveform(restWaveform);
    restCurrent.set_processed(restProcessed);
    restCurrent.set_harmonics(Inputs::calculate_harmonics_data(restWaveform, frequency));
    OperatingPointExcitation restExcitation;
    restExcitation.set_current(restCurrent);

    // Build excitation vector for all windings
    std::vector<OperatingPointExcitation> excitationPerWinding;
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
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
    return operatingPoint;
}

CoilMesherModels LeakageInductance::select_mesh_model(Magnetic& magnetic) {
    auto isPlanar = magnetic.get_wires()[0].get_type() == WireType::PLANAR;

    if (isPlanar) {
        double minimumRatio = DBL_MAX;
        for (auto wire : magnetic.get_wires()) {
            minimumRatio = std::min(minimumRatio, wire.get_maximum_conducting_width() / wire.get_maximum_conducting_height());
        }
        auto isThickPlanar = minimumRatio < PLANAR_THICKNESS_RATIO_THRESHOLD;
        if (!isThickPlanar) {
            return CoilMesherModels::WANG;
        }
    }

    return CoilMesherModels::CENTER;
}

std::pair<size_t, size_t> LeakageInductance::calculate_grid_points(Magnetic& magnetic, double frequency) {
    size_t numberPointsX;
    size_t numberPointsY;
    auto isPlanar = magnetic.get_wires()[0].get_type() == WireType::PLANAR;

    
    if (settings.get_leakage_inductance_grid_auto_scaling()) {
        auto gridPoints = calculate_number_points_needed_for_leakage(magnetic.get_coil());
        numberPointsX = gridPoints.first;
        numberPointsY = gridPoints.second;
        double precisionLevel = 1;
        if (isPlanar)  {
            precisionLevel = settings.get_leakage_inductance_grid_precision_level_planar();
        }
        else {
            precisionLevel = settings.get_leakage_inductance_grid_precision_level_wound();
        }
        precisionLevel = std::max(MINIMUM_PRECISION_LEVEL, precisionLevel);
        numberPointsX *= precisionLevel;
        numberPointsY *= precisionLevel;
    }
    else {
        numberPointsX = settings.get_magnetic_field_number_points_x();
        numberPointsY = settings.get_magnetic_field_number_points_y();
        if (isPlanar) {
            // If planar, we swap the number of points, as Y is larger by default
            numberPointsX = settings.get_magnetic_field_number_points_y();
            numberPointsY = settings.get_magnetic_field_number_points_x();
        }
    }

    return {numberPointsX, numberPointsY};
}

} // namespace OpenMagnetics
