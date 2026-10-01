#include "support/CoilMesher.h"
#include "physical_models/MagneticField.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/Reluctance.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <cfloat>
#include <algorithm>
#include <array>
#include <complex>
#include <map>
#include <mutex>
#include <tuple>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
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

std::pair<ComplexField, double> LeakageInductance::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding) {
    auto phasorField = calculate_magnetic_field_phasor(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex, customCurrentDirectionPerWinding);
    return {phasorField.inPhase, phasorField.dA};
}

LeakageInductance::PhasorField LeakageInductance::calculate_magnetic_field_phasor(OperatingPoint operatingPoint, Magnetic magnetic, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding) {

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];
    auto [numberPointsX, numberPointsY] = calculate_grid_points(magnetic, frequency);

    // ABT #1240. The window field is the exact 2-D field of the conductors (true cross-section, true
    // current, uniform density) inside the high-permeability window, built from a complete image lattice,
    // and its energy is integrated over the WHOLE core winding window. Measured against OMFEM 2D, each of
    // the removed shortcuts under-read the leakage: copper excluded from the grid (planar), points far from
    // every turn dropped (tall windows, separated sections), the grid stopping at the bobbin window edge,
    // the Wang two-filament planar mesh (no images), the truncated [-M, M] lattice (side-by-side sections
    // 0.68-0.85) and the zeroed in-conductor field.
    SettingsGuard<bool> completeCellsGuard(settings, &Settings::get_magnetic_field_mirroring_complete_cells, &Settings::set_magnetic_field_mirroring_complete_cells, true);
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    bool singleRectangularWindow = bobbin.get_winding_window_shape() == WindingWindowShape::RECTANGULAR &&
                                   magnetic.get_mutable_core().get_shape_family() != CoreShapeFamily::T &&
                                   magnetic.get_mutable_core().get_winding_windows().size() == 1;
    std::pair<Field, double> meshResult;
    if (singleRectangularWindow) {
        // The point counts resolve the bobbin window (window 0, a single chamber on a sectioned bobbin); keep
        // that pitch over the larger core window.
        auto bobbinWindowDimensions = bobbin.get_winding_window_dimensions();
        auto coreWindingWindow = magnetic.get_mutable_core().get_winding_windows()[0];
        if (!coreWindingWindow.get_width() || !coreWindingWindow.get_height()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Cannot calculate leakage inductance: the core winding window has no width or height");
        }
        double coreWindowWidth = coreWindingWindow.get_width().value();
        double coreWindowHeight = coreWindingWindow.get_height().value();
        size_t coreNumberPointsX = static_cast<size_t>(std::ceil(static_cast<double>(numberPointsX) * coreWindowWidth / bobbinWindowDimensions[0]));
        size_t coreNumberPointsY = static_cast<size_t>(std::ceil(static_cast<double>(numberPointsY) * coreWindowHeight / bobbinWindowDimensions[1]));
        meshResult = CoilMesher::generate_mesh_core_winding_window_grid(magnetic, frequency, coreNumberPointsX, coreNumberPointsY);
    }
    else {
        meshResult = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY);
    }
    Field inducedField = meshResult.first;

        if (inducedField.get_data().size() == 0) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Mesh generation failed: induced field data is empty");
        }
    double dA = meshResult.second;

    // Use dedicated leakage inductance H-field model (default: BINNS_LAWRENSON, which works best for air-only calculations)
    MagneticField magneticField(settings.get_leakage_inductance_magnetic_field_strength_model(), settings.get_magnetic_field_strength_fringing_effect_model());

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
    ComplexField quadratureField;
    {
        // Every conductor as its true cross-section with its images (CENTER); the Wang mesh splits a planar
        // track into two full-current filaments without images, a proximity-loss construction.
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField, customCurrentDirectionPerWinding, CoilMesherModels::CENTER);
        field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        quadratureField = windingWindowMagneticStrengthFieldOutput.get_quadrature_field_per_frequency()[0];
    }
    auto turns = magnetic.get_coil().get_turns_description().value();

    // Toroid-only outer-half stitching: concentric multi-window turns now also carry
    // additional coordinates (their second crossing), but this radius-gated overwrite
    // is derived for the round bore geometry — rectangular windows must skip it.
    if (turns[0].get_additional_coordinates() &&
        magnetic.get_mutable_coil().resolve_bobbin().get_winding_window_shape() == WindingWindowShape::ROUND) {
        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
            if (turns[turnIndex].get_additional_coordinates()) {
                turns[turnIndex].set_coordinates(turns[turnIndex].get_additional_coordinates().value()[0]);
            }
        }
        auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
        auto windingWindows = bobbin.get_processed_description()->get_winding_windows();

        double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

        magnetic.get_mutable_coil().set_turns_description(turns);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField, customCurrentDirectionPerWinding);
        auto additionalField = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        auto additionalQuadratureField = windingWindowMagneticStrengthFieldOutput.get_quadrature_field_per_frequency()[0];
        for (size_t pointIndex = 0; pointIndex < field.get_data().size(); ++pointIndex) {
            if (hypot(field.get_data()[pointIndex].get_point()[0], field.get_data()[pointIndex].get_point()[1]) > windingWindowRadialHeight) {
                field.get_mutable_data()[pointIndex].set_real(additionalField.get_data()[pointIndex].get_real());
                field.get_mutable_data()[pointIndex].set_imaginary(additionalField.get_data()[pointIndex].get_imaginary());
                quadratureField.get_mutable_data()[pointIndex].set_real(additionalQuadratureField.get_data()[pointIndex].get_real());
                quadratureField.get_mutable_data()[pointIndex].set_imaginary(additionalQuadratureField.get_data()[pointIndex].get_imaginary());
            }
        }
    }

    return {field, quadratureField, dA};
}

LeakageInductanceOutput LeakageInductance::calculate_leakage_inductance(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    // On an unprocessed magnetic (core shape/material still name strings, no
    // effective parameters) the field/reluctance path fails deep with a cryptic
    // "std::get: wrong index for variant". Surface it clearly instead.
    if (!magnetic.get_core().get_processed_description()) {
        throw CoreNotProcessedException(
            "Cannot calculate leakage inductance: the core has no processed description "
            "(effective parameters/shape unresolved). Run magnetic autocomplete / process the core first.");
    }

    // Magnetic shunts (MAS-RFC 0015, ABT #1176). A sheet across the window carries the leakage flux
    // through its own reluctance network, which the air-field integral below cannot see.
    if (magnetic.get_shunts() && !magnetic.get_shunts()->empty()) {
        MagneticShuntModel::check_supported_placements(magnetic);
        if (MagneticShuntModel::has_leakage_shunts(magnetic)) {
            if (ReluctanceNetwork::has_non_main_placement(magnetic)) {
                throw NotImplementedException("Leakage inductance with a magnetic shunt for windings placed on several columns");
            }
            auto shuntResult = calculate_shunt_leakage(magnetic, frequency, sourceIndex, destinationIndex, harmonicIndex);
            LeakageInductanceOutput shuntOutput;
            shuntOutput.set_method_used("Shunt");
            shuntOutput.set_origin(ResultOrigin::SIMULATION);
            DimensionWithTolerance shuntDimensionWithTolerance;
            shuntDimensionWithTolerance.set_nominal(shuntResult.leakageInductance);
            shuntOutput.set_leakage_inductance_per_winding({shuntDimensionWithTolerance});
            return shuntOutput;
        }
    }

    // Multi-column winding: the window-energy method below integrates the field of ONE
    // window revolved around the main column — meaningless for a winding pair sitting
    // on different columns, whose leakage flux closes through the other window and the
    // air outside the core. Use the per-column reluctance network's short-circuit
    // inductance (L_ss − L_sd²/L_dd, referred to the source): the magnetic-circuit
    // leakage between leg-separated windings.
    if (ReluctanceNetwork::has_non_main_placement(magnetic)) {
        auto columnIndexPerWinding = ReluctanceNetwork::resolve_winding_column_indexes(magnetic);
        if (columnIndexPerWinding[sourceIndex] != columnIndexPerWinding[destinationIndex]) {
            auto reluctanceModel = ReluctanceModel::factory();
            auto reluctanceOutput = reluctanceModel->get_core_reluctance(magnetic.get_mutable_core(), std::optional<OperatingPoint>(std::nullopt));
            ReluctanceNetwork reluctanceNetwork(magnetic.get_core(), reluctanceOutput.get_ungapped_core_reluctance().value(),
                                                reluctanceOutput.get_reluctance_per_gap().value_or(std::vector<AirGapReluctanceOutput>{}));
            auto inductanceMatrix = reluctanceNetwork.calculate_magnetizing_inductance_matrix(magnetic);
            double shortCircuitInductance = inductanceMatrix[sourceIndex][sourceIndex] -
                                            inductanceMatrix[sourceIndex][destinationIndex] * inductanceMatrix[sourceIndex][destinationIndex] /
                                                inductanceMatrix[destinationIndex][destinationIndex];
            LeakageInductanceOutput leakageInductanceOutput;
            leakageInductanceOutput.set_method_used("ReluctanceNetwork");
            leakageInductanceOutput.set_origin(ResultOrigin::SIMULATION);
            DimensionWithTolerance dimensionWithTolerance;
            dimensionWithTolerance.set_nominal(shortCircuitInductance);
            leakageInductanceOutput.set_leakage_inductance_per_winding({dimensionWithTolerance});
            return leakageInductanceOutput;
        }
    }

    // Toroids: exact ring-plane field with per-region images and the body-of-revolution 3-D correction
    // (see LeakageInductance.h). Ampere-turn balanced pair: 1 A in the source against the turns-ratio current
    // in the destination; the result is referred to the source.
    if (magnetic.get_mutable_core().get_shape_family() == CoreShapeFamily::T) {
        if (sourceIndex == destinationIndex) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Leakage inductance needs two different windings");
        }
        size_t numberWindings = magnetic.get_coil().get_functional_description().size();
        std::vector<double> currents(numberWindings, 0.0);
        currents[sourceIndex] = 1.0;
        currents[destinationIndex] = -static_cast<double>(magnetic.get_mutable_coil().get_number_turns(sourceIndex)) /
                                     static_cast<double>(magnetic.get_mutable_coil().get_number_turns(destinationIndex));
        auto toroidalEnergy = calculate_toroidal_leakage_energy(magnetic, currents);
        LeakageInductanceOutput toroidalOutput;
        toroidalOutput.set_method_used("Energy");
        toroidalOutput.set_origin(ResultOrigin::SIMULATION);
        DimensionWithTolerance toroidalDimensionWithTolerance;
        toroidalDimensionWithTolerance.set_nominal(2.0 * toroidalEnergy.energy);
        toroidalOutput.set_leakage_inductance_per_winding({toroidalDimensionWithTolerance});
        return toroidalOutput;
    }

    // RAII: any throw between the manual set/restore pair (several are right below) used
    // to leave fringing globally disabled for the rest of the process.
    SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing, &Settings::set_magnetic_field_include_fringing, false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw CoilNotProcessedException("Cannot calculate leakage inductance: bobbin description has not been processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Cannot calculate leakage inductance: bobbin column width is not defined");
    }


    OperatingPoint operatingPoint = create_leakage_operating_point(magnetic, sourceIndex, destinationIndex, frequency);

    // The leakage measurement requires the source and destination ampere-turns to OPPOSE.
    // The field model's default direction vector is [+1, -1, -1, ...] (winding 0 positive,
    // the rest negative), which only opposes correctly when the source is winding 0. For a
    // source/destination pair that does not involve winding 0 (e.g. between two secondaries)
    // the default makes both currents share a sign, so their MMFs ADD instead of opposing and
    // the "leakage" comes out as the additive energy. Pass an explicit direction vector that
    // opposes source and destination for any pair of indices.
    size_t numberWindingsForDirection = magnetic.get_coil().get_functional_description().size();
    std::vector<int8_t> currentDirectionPerWinding(numberWindingsForDirection, 0);
    currentDirectionPerWinding[sourceIndex] = 1;
    currentDirectionPerWinding[destinationIndex] = -1;

    auto magneticFieldResult = calculate_magnetic_field_phasor(operatingPoint, magnetic, sourceIndex, destinationIndex, harmonicIndex, currentDirectionPerWinding);
    ComplexField field = magneticFieldResult.inPhase;
    double dA = magneticFieldResult.dA;

    // Time-averaged energy of a phasor field: in-phase plus quadrature.
    double energy = integrate_leakage_energy(magnetic, field, dA) + integrate_leakage_energy(magnetic, magneticFieldResult.quadrature, dA);

    // The field model drives every turn with the PEAK amplitude of the current harmonic
    // (CoilMesher::generate_mesh_inducing_coil), so |H|^2 is a peak-squared field and the
    // integral is the peak stored energy W = L I_peak^2 / 2. Normalise by the same peak
    // amplitude. Dividing by I_rms^2 (as before ABT #1211) doubled every Energy-method result.
    double sourceCurrentPeak = harmonic_peak_current_at_field_frequency(operatingPoint.get_excitations_per_winding()[sourceIndex], field.get_frequency());
    double leakageInductance = 2.0 * energy / (sourceCurrentPeak * sourceCurrentPeak);
    LeakageInductanceOutput leakageInductanceOutput;

    leakageInductanceOutput.set_method_used("Energy");
    leakageInductanceOutput.set_origin(ResultOrigin::SIMULATION);
    DimensionWithTolerance dimensionWithTolerance;
    dimensionWithTolerance.set_nominal(leakageInductance);
    leakageInductanceOutput.set_leakage_inductance_per_winding({dimensionWithTolerance});

    return leakageInductanceOutput;
}

MagneticShuntLeakageResult LeakageInductance::calculate_shunt_leakage(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex,
                                                                     std::optional<std::vector<double>> relativePermeabilityPerShunt) {
    if (!magnetic.get_shunts() || magnetic.get_shunts()->empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Shunt leakage method called on a magnetic without shunts");
    }
    MagneticShuntModel::check_supported_placements(magnetic);
    auto shunts = magnetic.get_shunts().value();
    if (relativePermeabilityPerShunt && relativePermeabilityPerShunt->size() != shunts.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Shunt leakage method: " + std::to_string(relativePermeabilityPerShunt->size()) +
                                    " relative permeabilities given for " + std::to_string(shunts.size()) + " shunts");
    }

    // The windings' own leakage, by the air-field Energy method, without the sheets.
    Magnetic withoutShunts = magnetic;
    withoutShunts.set_shunts(std::nullopt);
    auto windingOutput = calculate_leakage_inductance(withoutShunts, frequency, sourceIndex, destinationIndex, harmonicIndex);
    if (windingOutput.get_method_used() != "Energy") {
        throw NotImplementedException("Shunt leakage method on top of the " + windingOutput.get_method_used() + " leakage method");
    }
    double windingLeakageInductance = windingOutput.get_leakage_inductance_per_winding()[0].get_nominal().value();

    return MagneticShuntModel::assemble_leakage(magnetic, windingLeakageInductance, frequency, sourceIndex, destinationIndex, relativePermeabilityPerShunt);
}

ComplexField LeakageInductance::calculate_leakage_magnetic_field(Magnetic magnetic, double frequency, size_t sourceIndex, size_t destinationIndex, size_t harmonicIndex) {
    // RAII: this function never restored the flag at all — one call permanently
    // disabled fringing for every later field computation in the process.
    SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing, &Settings::set_magnetic_field_include_fringing, false);

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

    leakageInductanceOutput.set_method_used(MagneticShuntModel::has_leakage_shunts(magnetic) ? "Shunt" : "Energy");
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

double LeakageInductance::integrate_leakage_energy(Magnetic& magnetic, ComplexField& field, double dA) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
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

    return energy;
}

double LeakageInductance::calculate_leakage_field_energy(Magnetic magnetic, const std::vector<double>& currentsRmsSigned, double frequency, size_t harmonicIndex) {
    if (MagneticShuntModel::has_leakage_shunts(magnetic)) {
        throw NotImplementedException("Leakage field energy (and the leakage inductance matrix) of a magnetic with a shunt in the window");
    }
    if (!magnetic.get_core().get_processed_description()) {
        throw CoreNotProcessedException(
            "Cannot calculate leakage field energy: the core has no processed description "
            "(effective parameters/shape unresolved). Run magnetic autocomplete / process the core first.");
    }
    size_t numberWindings = magnetic.get_coil().get_functional_description().size();
    if (currentsRmsSigned.size() != numberWindings) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot calculate leakage field energy: current vector size (" + std::to_string(currentsRmsSigned.size()) +
            ") does not match number of windings (" + std::to_string(numberWindings) + ")");
    }

    // RAII guard (see calculate_leakage_inductance above).
    SettingsGuard<bool> fringingGuard(settings, &Settings::get_magnetic_field_include_fringing, &Settings::set_magnetic_field_include_fringing, false);

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw CoilNotProcessedException("Cannot calculate leakage field energy: bobbin description has not been processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Cannot calculate leakage field energy: bobbin column width is not defined");
    }

    // Split the signed RMS current vector into magnitudes (operating point) and directions
    // (sign vector passed to the field model, which scales each winding's turn currents).
    std::vector<int8_t> directions;
    directions.reserve(numberWindings);
    for (auto current : currentsRmsSigned) {
        directions.push_back(current < 0 ? int8_t(-1) : int8_t(1));
    }

    OperatingPoint operatingPoint = create_excitation_operating_point(magnetic, currentsRmsSigned, frequency);

    if (magnetic.get_mutable_core().get_shape_family() == CoreShapeFamily::T) {
        // The same peak amplitudes the field model would drive the turns with (the harmonic at the field frequency).
        double fieldFrequency = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics()->get_frequencies()[harmonicIndex];
        std::vector<double> peakCurrents(numberWindings, 0.0);
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            if (currentsRmsSigned[windingIndex] != 0) {
                peakCurrents[windingIndex] = directions[windingIndex] * harmonic_peak_current_at_field_frequency(operatingPoint.get_excitations_per_winding()[windingIndex], fieldFrequency);
            }
        }
        return calculate_toroidal_leakage_energy(magnetic, peakCurrents).energy;
    }

    auto magneticFieldResult = calculate_magnetic_field_phasor(operatingPoint, magnetic, 0, 1, harmonicIndex, directions);
    ComplexField field = magneticFieldResult.inPhase;
    double dA = magneticFieldResult.dA;

    // Time-averaged energy of a phasor field: in-phase plus quadrature.
    double energy = integrate_leakage_energy(magnetic, field, dA) + integrate_leakage_energy(magnetic, magneticFieldResult.quadrature, dA);

    return energy;
}

std::vector<std::vector<double>> LeakageInductance::calculate_leakage_inductance_matrix(Magnetic magnetic, double frequency, size_t harmonicIndex) {
    size_t numberWindings = magnetic.get_coil().get_functional_description().size();
    if (numberWindings == 0) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS,
            "Cannot calculate leakage inductance matrix: no windings defined");
    }

    // Self-leakage energies W(e_a): excite a single winding with the unit reference current.
    std::vector<double> selfEnergy(numberWindings);
    for (size_t a = 0; a < numberWindings; ++a) {
        std::vector<double> currents(numberWindings, 0.0);
        currents[a] = 1.0;
        selfEnergy[a] = calculate_leakage_field_energy(magnetic, currents, frequency, harmonicIndex);
    }

    std::vector<std::vector<double>> leakageMatrix(numberWindings, std::vector<double>(numberWindings, 0.0));
    // The unit reference excitation is a sinusoid of 1 A peak; the field (and so W) is built from
    // the peak amplitude of its fundamental harmonic, a. With W = 1/2 i^T Λ i for peak currents:
    //   Λ_aa = 2·W(e_a)/a²,   Λ_ab = [W(e_a+e_b) − W(e_a) − W(e_b)]/a²
    // (ABT #1211: this used 4·W and 2·[...], the I_rms normalisation that doubled every value).
    std::vector<double> unitCurrents(numberWindings, 0.0);
    unitCurrents[0] = 1.0;
    auto unitOperatingPoint = create_excitation_operating_point(magnetic, unitCurrents, frequency);
    double referencePeak = harmonic_peak_current_at_field_frequency(unitOperatingPoint.get_excitations_per_winding()[0], frequency);
    double referencePeakSquared = referencePeak * referencePeak;
    for (size_t a = 0; a < numberWindings; ++a) {
        leakageMatrix[a][a] = 2.0 * selfEnergy[a] / referencePeakSquared;
    }
    // Off-diagonal via polarization.
    for (size_t a = 0; a < numberWindings; ++a) {
        for (size_t b = a + 1; b < numberWindings; ++b) {
            std::vector<double> currents(numberWindings, 0.0);
            currents[a] = 1.0;
            currents[b] = 1.0;
            double pairEnergy = calculate_leakage_field_energy(magnetic, currents, frequency, harmonicIndex);
            double mutualLeakage = (pairEnergy - selfEnergy[a] - selfEnergy[b]) / referencePeakSquared;
            leakageMatrix[a][b] = mutualLeakage;
            leakageMatrix[b][a] = mutualLeakage;
        }
    }

    return leakageMatrix;
}

double LeakageInductance::harmonic_peak_current_at_field_frequency(const OperatingPointExcitation& excitation, double fieldFrequency) {
    if (!excitation.get_current() || !excitation.get_current()->get_harmonics()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Leakage inductance: the source excitation has no current harmonics");
    }
    auto harmonics = excitation.get_current()->get_harmonics().value();
    auto& frequencies = harmonics.get_frequencies();
    auto& amplitudes = harmonics.get_amplitudes();
    for (size_t harmonicIndex = 0; harmonicIndex < frequencies.size(); ++harmonicIndex) {
        if (frequencies[harmonicIndex] > 0 && std::abs(frequencies[harmonicIndex] - fieldFrequency) <= 1e-9 * fieldFrequency) {
            if (!(amplitudes[harmonicIndex] > 0)) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Leakage inductance: the source current harmonic at " +
                                           std::to_string(fieldFrequency) + " Hz has no amplitude");
            }
            return amplitudes[harmonicIndex];
        }
    }
    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Leakage inductance: the source current has no harmonic at the field frequency " +
                               std::to_string(fieldFrequency) + " Hz");
}

OperatingPoint LeakageInductance::create_excitation_operating_point(Magnetic& magnetic, const std::vector<double>& currentsRmsSigned, double frequency) {
    size_t numberWindings = magnetic.get_coil().get_functional_description().size();
    std::vector<OperatingPointExcitation> excitationPerWinding;
    excitationPerWinding.reserve(numberWindings);

    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        // The field model takes current direction from the sign vector, so the operating
        // point carries the magnitude only. A zero-amplitude waveform still yields valid
        // harmonic frequency bins (read from winding 0 to set the field frequency).
        // Mirror create_leakage_operating_point's amplitude convention exactly so that an
        // entry of 1.0 here equals that method's "unit source" (peak_to_peak=2, rms=1/sqrt2).
        double magnitude = std::abs(currentsRmsSigned[windingIndex]);
        ProcessedWaveform processed;
        processed.set_peak_to_peak(SINUSOIDAL_PEAK_TO_PEAK * magnitude);
        processed.set_duty_cycle(SINUSOIDAL_DUTY_CYCLE);
        processed.set_offset(SINUSOIDAL_OFFSET);
        processed.set_rms(magnitude / std::sqrt(SINUSOIDAL_PEAK_TO_PEAK));
        processed.set_label(WaveformLabel::SINUSOIDAL);
        auto waveform = Inputs::create_waveform(processed, frequency);
        SignalDescriptor current;
        current.set_waveform(waveform);
        current.set_processed(processed);
        current.set_harmonics(Inputs::calculate_harmonics_data(waveform, frequency));
        OperatingPointExcitation excitation;
        excitation.set_current(current);
        excitationPerWinding.push_back(excitation);
    }

    OperatingPoint operatingPoint;
    operatingPoint.set_excitations_per_winding(excitationPerWinding);
    return operatingPoint;
}

OperatingPoint LeakageInductance::create_leakage_operating_point(Magnetic& magnetic, size_t sourceIndex, size_t destinationIndex, double frequency) {
    // Create source excitation (1 A RMS sinusoidal)
    ProcessedWaveform sourceProcessed;
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
    ProcessedWaveform destinationProcessed = sourceProcessed;
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
    ProcessedWaveform restProcessed = sourceProcessed;
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

// ---------------------------------------------------------------------------------------------------------
// Toroidal cores. The model is described in LeakageInductance.h.
// ---------------------------------------------------------------------------------------------------------

double LeakageInductance::calculate_ring_plane_region_energy_per_length(const std::vector<ToroidalLineCurrent>& conductors, double wallRadius, double imageFactor, bool interiorRegion) {
    if (!(wallRadius > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal leakage: the wall radius must be positive");
    }
    // A(z_j) = -(µ0 / 2π) Σ_k I_k [ln|z_j − z_k| + k·ln|z_j − R²/conj(z_k)| − k·ln|z_j|]. With
    // |z_j − R²/conj(z_k)| = |z_j·conj(z_k) − R²| / |z_k| the kernel is symmetric in (j, k). The self term uses the
    // conductor's geometric mean radius, which includes the energy inside the wire.
    //
    // Net current. The centre image makes the wall free of tangential field, which is right only when the region
    // carries no net current. A net current Σ I (a single winding, or any excitation whose ampere-turns do not
    // balance) links the core: the core field at the wall is the magnetizing H = Σ I_bore / (2π R), and the air
    // between the conductors and the wall carries that circulation too. Without it the kernel changes by ln(s)
    // when the geometry is scaled by s, so the energy of a net current would depend on the unit of length.
    // netCurrentTerm restores the circulation (exact for k = 1, and gauge invariant for any k; it vanishes when
    // Σ I = 0):
    //   bore:    + ln|z_j| + ln|z_k| − 3·ln R   (field: direct + wall images, no centre image)
    //   outside: − ln|z_j| − ln|z_k| + ln R     (field: direct + wall images + −(1 + k)·Σ I at the centre)
    // For a uniform ring of radius ρ this is the annulus energy µ0·I²/(4π)·ln(R/ρ) in the bore and
    // µ0·I²/(4π)·ln(ρ/R) outside.
    double vacuumPermeability = Constants().vacuumPermeability;
    double logWall = std::log(wallRadius);
    double sum = 0;
    for (size_t j = 0; j < conductors.size(); ++j) {
        std::complex<double> zj(conductors[j].x, conductors[j].y);
        double rj = std::abs(zj);
        if (interiorRegion ? !(rj < wallRadius) : !(rj > wallRadius)) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "Toroidal leakage: a turn crossing at radius " + std::to_string(rj) +
                                        " m is on the wrong side of the core wall at " + std::to_string(wallRadius) + " m");
        }
        if (!(conductors[j].geometricMeanRadius > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Toroidal leakage: a conductor has no positive geometric mean radius");
        }
        for (size_t k = 0; k < conductors.size(); ++k) {
            std::complex<double> zk(conductors[k].x, conductors[k].y);
            double rk = std::abs(zk);
            double direct = (j == k) ? std::log(conductors[j].geometricMeanRadius) : std::log(std::abs(zj - zk));
            double image = std::log(std::abs(zj * std::conj(zk) - wallRadius * wallRadius)) - std::log(rk) - std::log(rj);
            double netCurrentTerm = interiorRegion ? std::log(rj) + std::log(rk) - 3 * logWall : -std::log(rj) - std::log(rk) + logWall;
            sum += conductors[j].current * conductors[k].current * (direct + imageFactor * image + netCurrentTerm);
        }
    }
    return -0.5 * vacuumPermeability / (2 * std::numbers::pi) * sum;
}

std::array<double, 2> LeakageInductance::calculate_ring_plane_region_field(const std::vector<ToroidalLineCurrent>& conductors, double wallRadius, double imageFactor, bool interiorRegion, double x, double y) {
    // H of a z-directed line current I at c: (I / 2π)·(−(y − c_y), x − c_x) / ρ².
    auto lineField = [](double current, double cx, double cy, double px, double py) -> std::array<double, 2> {
        double dx = px - cx;
        double dy = py - cy;
        double rho2 = dx * dx + dy * dy;
        return {-current / (2 * std::numbers::pi) * dy / rho2, current / (2 * std::numbers::pi) * dx / rho2};
    };
    std::array<double, 2> field = {0, 0};
    for (auto& conductor : conductors) {
        double radius = conductor.geometricMeanRadius * std::exp(0.25);
        double dx = x - conductor.x;
        double dy = y - conductor.y;
        if (dx * dx + dy * dy < radius * radius) {
            // Uniform current density inside the wire: H = I·ρ / (2π r²).
            field[0] += -conductor.current / (2 * std::numbers::pi * radius * radius) * dy;
            field[1] += conductor.current / (2 * std::numbers::pi * radius * radius) * dx;
        }
        else {
            auto h = lineField(conductor.current, conductor.x, conductor.y, x, y);
            field[0] += h[0];
            field[1] += h[1];
        }
        double scale = wallRadius * wallRadius / (conductor.x * conductor.x + conductor.y * conductor.y);
        auto hImage = lineField(imageFactor * conductor.current, conductor.x * scale, conductor.y * scale, x, y);
        // Centre image (see calculate_ring_plane_region_energy_per_length): none in the bore, so the net current
        // keeps its circulation at the wall; −(1 + k)·I outside, so the wall carries the core's circulation and the
        // far field of the bore and outer crossings together vanishes.
        double centreCurrent = interiorRegion ? 0.0 : -(1 + imageFactor) * conductor.current;
        auto hCentre = lineField(centreCurrent, 0, 0, x, y);
        field[0] += hImage[0] + hCentre[0];
        field[1] += hImage[1] + hCentre[1];
    }
    return field;
}

double LeakageInductance::calculate_sheet_mode_energy_per_length(const std::vector<double>& angles, const std::vector<double>& currents, size_t mode) {
    if (mode == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal leakage: the MMF sheet harmonics start at m = 1");
    }
    if (angles.size() != currents.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal leakage: one angle per current is needed");
    }
    std::complex<double> harmonic = 0;
    for (size_t k = 0; k < angles.size(); ++k) {
        harmonic += currents[k] * std::exp(std::complex<double>(0, -static_cast<double>(mode) * angles[k]));
    }
    return Constants().vacuumPermeability * std::norm(harmonic) / (std::numbers::pi * static_cast<double>(mode));
}

namespace {
// Nodes from start towards end (either direction): first spacing firstStep, growing by ratio up to maximumStep.
std::vector<double> graded_nodes(double start, double end, double firstStep, double ratio, double maximumStep) {
    std::vector<double> nodes = {start};
    double direction = end > start ? 1.0 : -1.0;
    double length = std::abs(end - start);
    double position = 0;
    double step = firstStep;
    while (position < length) {
        position = std::min(length, position + step);
        nodes.push_back(start + direction * position);
        step = std::min(step * ratio, maximumStep);
    }
    if (nodes.size() > 2 && std::abs(nodes[nodes.size() - 1] - nodes[nodes.size() - 2]) < 0.3 * firstStep) {
        nodes.erase(nodes.end() - 2);
    }
    return nodes;
}

std::vector<double> merge_nodes(std::vector<double> nodes, double tolerance) {
    std::sort(nodes.begin(), nodes.end());
    std::vector<double> merged;
    for (auto node : nodes) {
        if (merged.empty() || node - merged.back() > tolerance) {
            merged.push_back(node);
        }
    }
    return merged;
}
} // namespace

double LeakageInductance::calculate_body_of_revolution_effective_height(double innerRadius, double outerRadius, double height, size_t mode, double refinement, double farFieldFactor) {
    if (!(innerRadius > 0) || !(outerRadius > innerRadius) || !(height > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal leakage: the conductor envelope needs 0 < inner radius < outer radius and a positive height");
    }
    if (mode == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal leakage: the body-of-revolution solve needs an azimuthal harmonic m >= 1");
    }
    if (!(refinement > 0) || !(farFieldFactor > 1)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal leakage: the body-of-revolution mesh needs refinement > 0 and farFieldFactor > 1");
    }
    // Quarter meridian plane r >= 0, z >= 0 (the problem is even in z). The 3-D potential is ψ(r, z)·cos(mφ) with
    // ψ = 1 on the envelope, ψ = 0 on the axis (m >= 1) and on the far boundary. Bilinear elements on a tensor grid
    // graded towards the envelope edges. The energy functional per element is ∫∫ (ψ_r² + ψ_z² + m²ψ²/r²) r dr dz;
    // ∫cos²(mφ)dφ = π and the lower half plane double it, so the 3-D energy is 2π times the quarter-plane value.
    double halfHeight = height / 2;
    double m = static_cast<double>(mode);
    double smallestFeature = std::min({outerRadius - innerRadius, halfHeight, innerRadius});
    double firstStep = smallestFeature / (40.0 * refinement);
    double ratio = 1.15;
    double farRadius = farFieldFactor * outerRadius;
    double interiorMaximumStep = smallestFeature / (6.0 * refinement);

    std::vector<double> rNodes;
    double middleRadius = (innerRadius + outerRadius) / 2;
    for (auto node : graded_nodes(innerRadius, 0, firstStep, ratio, farRadius)) rNodes.push_back(node);
    for (auto node : graded_nodes(innerRadius, middleRadius, firstStep, ratio, interiorMaximumStep)) rNodes.push_back(node);
    for (auto node : graded_nodes(outerRadius, middleRadius, firstStep, ratio, interiorMaximumStep)) rNodes.push_back(node);
    for (auto node : graded_nodes(outerRadius, farRadius, firstStep, ratio, farRadius)) rNodes.push_back(node);
    std::vector<double> zNodes;
    for (auto node : graded_nodes(halfHeight, 0, firstStep, ratio, farRadius)) zNodes.push_back(node);
    for (auto node : graded_nodes(halfHeight, farRadius, firstStep, ratio, farRadius)) zNodes.push_back(node);
    double tolerance = 1e-9 * smallestFeature;
    rNodes = merge_nodes(rNodes, tolerance);
    zNodes = merge_nodes(zNodes, tolerance);

    size_t nr = rNodes.size();
    size_t nz = zNodes.size();
    auto nodeIndex = [nz](size_t i, size_t j) { return i * nz + j; };
    auto insideEnvelope = [&](double r, double z) {
        return r > innerRadius + tolerance && r < outerRadius - tolerance && z < halfHeight - tolerance;
    };
    auto onEnvelope = [&](double r, double z) {
        return r >= innerRadius - tolerance && r <= outerRadius + tolerance && z <= halfHeight + tolerance;
    };

    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<bool> used(nr * nz, false);
    const double gauss = 1.0 / std::sqrt(3.0);
    for (size_t i = 0; i + 1 < nr; ++i) {
        double r0 = rNodes[i];
        double dr = rNodes[i + 1] - r0;
        for (size_t j = 0; j + 1 < nz; ++j) {
            double dz = zNodes[j + 1] - zNodes[j];
            if (insideEnvelope(r0 + dr / 2, zNodes[j] + dz / 2)) {
                continue;
            }
            double element[4][4] = {};
            for (double a : {-gauss, gauss}) {
                for (double b : {-gauss, gauss}) {
                    double xi = (a + 1) / 2;
                    double eta = (b + 1) / 2;
                    double r = r0 + xi * dr;
                    double weight = dr * dz / 4 * r;
                    double shape[4] = {(1 - xi) * (1 - eta), xi * (1 - eta), (1 - xi) * eta, xi * eta};
                    double dShapeDr[4] = {-(1 - eta) / dr, (1 - eta) / dr, -eta / dr, eta / dr};
                    double dShapeDz[4] = {-(1 - xi) / dz, -xi / dz, (1 - xi) / dz, xi / dz};
                    for (size_t p = 0; p < 4; ++p) {
                        for (size_t q = 0; q < 4; ++q) {
                            element[p][q] += weight * (dShapeDr[p] * dShapeDr[q] + dShapeDz[p] * dShapeDz[q] + m * m / (r * r) * shape[p] * shape[q]);
                        }
                    }
                }
            }
            size_t nodes[4] = {nodeIndex(i, j), nodeIndex(i + 1, j), nodeIndex(i, j + 1), nodeIndex(i + 1, j + 1)};
            for (size_t p = 0; p < 4; ++p) {
                used[nodes[p]] = true;
                for (size_t q = 0; q < 4; ++q) {
                    triplets.emplace_back(static_cast<int>(nodes[p]), static_cast<int>(nodes[q]), element[p][q]);
                }
            }
        }
    }
    size_t numberNodes = nr * nz;
    Eigen::SparseMatrix<double> stiffness(static_cast<int>(numberNodes), static_cast<int>(numberNodes));
    stiffness.setFromTriplets(triplets.begin(), triplets.end());

    // Dirichlet values. The free unknowns are the used nodes that are neither on the envelope nor on the axis or far box.
    Eigen::VectorXd potential = Eigen::VectorXd::Zero(static_cast<int>(numberNodes));
    std::vector<int> freeIndex(numberNodes, -1);
    int numberFree = 0;
    for (size_t i = 0; i < nr; ++i) {
        for (size_t j = 0; j < nz; ++j) {
            size_t n = nodeIndex(i, j);
            if (!used[n]) {
                continue;
            }
            if (onEnvelope(rNodes[i], zNodes[j])) {
                potential[static_cast<int>(n)] = 1.0;
            }
            else if (i == 0 || i == nr - 1 || j == nz - 1) {
                potential[static_cast<int>(n)] = 0.0;
            }
            else {
                freeIndex[n] = numberFree++;
            }
        }
    }
    std::vector<Eigen::Triplet<double>> freeTriplets;
    Eigen::VectorXd rightHandSide = Eigen::VectorXd::Zero(numberFree);
    for (int k = 0; k < stiffness.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness, k); it; ++it) {
            int row = freeIndex[static_cast<size_t>(it.row())];
            if (row < 0) {
                continue;
            }
            int column = freeIndex[static_cast<size_t>(it.col())];
            if (column >= 0) {
                freeTriplets.emplace_back(row, column, it.value());
            }
            else {
                rightHandSide[row] -= it.value() * potential[static_cast<int>(it.col())];
            }
        }
    }
    Eigen::SparseMatrix<double> freeStiffness(numberFree, numberFree);
    freeStiffness.setFromTriplets(freeTriplets.begin(), freeTriplets.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(freeStiffness);
    if (solver.info() != Eigen::Success) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Toroidal leakage: the body-of-revolution stiffness matrix could not be factorised");
    }
    Eigen::VectorXd freeSolution = solver.solve(rightHandSide);
    if (solver.info() != Eigen::Success) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Toroidal leakage: the body-of-revolution solve failed");
    }
    for (size_t n = 0; n < numberNodes; ++n) {
        if (freeIndex[n] >= 0) {
            potential[static_cast<int>(n)] = freeSolution[freeIndex[n]];
        }
    }
    double energyIntegral = 2 * std::numbers::pi * potential.dot(stiffness * potential);
    return energyIntegral / (2 * std::numbers::pi * m);
}

double LeakageInductance::cached_body_of_revolution_effective_height(double innerRadius, double outerRadius, double height, size_t mode) {
    // h_eff(m) depends only on the envelope; adviser loops evaluate the same core many times.
    static std::mutex cacheMutex;
    static std::map<std::tuple<long long, long long, long long, size_t>, double> cache;
    auto key = std::make_tuple(std::llround(innerRadius * 1e9), std::llround(outerRadius * 1e9), std::llround(height * 1e9), mode);
    {
        std::lock_guard<std::mutex> lock(cacheMutex);
        auto found = cache.find(key);
        if (found != cache.end()) {
            return found->second;
        }
    }
    double effectiveHeight = calculate_body_of_revolution_effective_height(innerRadius, outerRadius, height, mode);
    std::lock_guard<std::mutex> lock(cacheMutex);
    cache[key] = effectiveHeight;
    return effectiveHeight;
}

LeakageInductance::ToroidalRingPlaneConductors LeakageInductance::calculate_toroidal_ring_plane_conductors(Magnetic magnetic, const std::vector<double>& currentPerWinding) {
    auto& core = magnetic.get_mutable_core();
    if (core.get_shape_family() != CoreShapeFamily::T) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal leakage model called on a core that is not toroidal");
    }
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Toroidal leakage: the core has no processed description");
    }
    auto& coil = magnetic.get_mutable_coil();
    auto windings = coil.get_functional_description();
    if (currentPerWinding.size() != windings.size()) {
        throw InvalidInputException(ErrorCode::COIL_INVALID_TURNS, "Toroidal leakage: " + std::to_string(currentPerWinding.size()) +
                                    " currents given for " + std::to_string(windings.size()) + " windings");
    }
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Toroidal leakage: the coil has no turns description");
    }
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    if (dimensions.find("A") == dimensions.end() || dimensions.find("B") == dimensions.end() || dimensions.find("C") == dimensions.end()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal leakage: the toroid shape needs dimensions A, B and C");
    }
    double outerWallRadius = dimensions["A"] / 2;
    double innerWallRadius = dimensions["B"] / 2;
    // Core::get_number_stacks is MKF's single reading of the stack count (MAS: absent = one core).
    double coreHeight = dimensions["C"] * static_cast<double>(core.get_number_stacks());
    double corePermeability = core.get_initial_permeability(Defaults().ambientTemperature);
    if (!(corePermeability >= 1)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal leakage: the core initial permeability is below 1");
    }
    double imageFactor = (corePermeability - 1) / (corePermeability + 1);

    auto wires = coil.get_wires();
    std::vector<double> geometricMeanRadiusPerWinding;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        auto& wire = wires[windingIndex];
        double geometricMeanRadius;
        switch (wire.get_type()) {
            case WireType::ROUND:
            case WireType::LITZ:
                geometricMeanRadius = std::exp(-0.25) * wire.get_maximum_conducting_width() / 2;
                break;
            case WireType::RECTANGULAR:
            case WireType::FOIL:
            case WireType::PLANAR:
                // Geometric mean distance of a rectangle from itself, 0.2235·(a + b) (Grover).
                geometricMeanRadius = 0.2235 * (wire.get_maximum_conducting_width() + wire.get_maximum_conducting_height());
                break;
            default:
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Toroidal leakage: unsupported wire type");
        }
        if (!(geometricMeanRadius > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Toroidal leakage: winding " + windings[windingIndex].get_name() + " has no conducting dimensions");
        }
        geometricMeanRadiusPerWinding.push_back(geometricMeanRadius);
    }

    std::map<std::string, size_t> windingIndexByName;
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        windingIndexByName[windings[windingIndex].get_name()] = windingIndex;
    }

    ToroidalRingPlaneConductors conductors;
    conductors.innerWallRadius = innerWallRadius;
    conductors.outerWallRadius = outerWallRadius;
    conductors.coreHeight = coreHeight;
    conductors.imageFactor = imageFactor;
    const auto turns = coil.get_turns_description().value();
    for (auto& turn : turns) {
        auto found = windingIndexByName.find(turn.get_winding());
        if (found == windingIndexByName.end()) {
            throw CoilNotProcessedException("Toroidal leakage: turn " + turn.get_name() + " belongs to unknown winding " + turn.get_winding());
        }
        size_t windingIndex = found->second;
        double current = currentPerWinding[windingIndex] / static_cast<double>(windings[windingIndex].get_number_parallels());
        if (current == 0) {
            continue;
        }
        if (!turn.get_coordinate_system() || turn.get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            throw CoilNotProcessedException("Toroidal leakage: turn " + turn.get_name() + " coordinates are not cartesian");
        }
        if (!turn.get_additional_coordinates() || turn.get_additional_coordinates()->empty() || turn.get_additional_coordinates().value()[0].size() < 2) {
            throw CoilNotProcessedException("Toroidal leakage: turn " + turn.get_name() + " has no outer crossing (additional coordinates)");
        }
        auto inner = turn.get_coordinates();
        auto outer = turn.get_additional_coordinates().value()[0];
        double geometricMeanRadius = geometricMeanRadiusPerWinding[windingIndex];
        // The bore crossing and the outer crossing of one turn carry its current in opposite directions.
        conductors.bore.push_back({inner[0], inner[1], current, geometricMeanRadius});
        conductors.outer.push_back({outer[0], outer[1], -current, geometricMeanRadius});
        conductors.turnAngles.push_back(std::atan2(inner[1], inner[0]));
        conductors.turnCurrents.push_back(current);
    }
    return conductors;
}

std::array<double, 2> LeakageInductance::calculate_toroidal_ring_plane_field(const ToroidalRingPlaneConductors& conductors, double x, double y) {
    double radius = std::hypot(x, y);
    if (radius < conductors.innerWallRadius) {
        return calculate_ring_plane_region_field(conductors.bore, conductors.innerWallRadius, conductors.imageFactor, true, x, y);
    }
    if (radius > conductors.outerWallRadius) {
        return calculate_ring_plane_region_field(conductors.outer, conductors.outerWallRadius, conductors.imageFactor, false, x, y);
    }
    // Ferrite: H of a z-directed line current I at c is (I / 2π)·(−(y − c_y), x − c_x) / ρ².
    auto add_line = [&](std::array<double, 2>& field, double current, double cx, double cy) {
        double dx = x - cx;
        double dy = y - cy;
        double rho2 = dx * dx + dy * dy;
        field[0] += -current / (2 * std::numbers::pi) * dy / rho2;
        field[1] += current / (2 * std::numbers::pi) * dx / rho2;
    };
    std::array<double, 2> field = {0, 0};
    double transmission = 1 - conductors.imageFactor;
    double boreNetCurrent = 0;
    for (auto& conductor : conductors.bore) {
        add_line(field, transmission * conductor.current, conductor.x, conductor.y);
        boreNetCurrent += conductor.current;
    }
    for (auto& conductor : conductors.outer) {
        add_line(field, transmission * conductor.current, conductor.x, conductor.y);
    }
    add_line(field, conductors.imageFactor * boreNetCurrent, 0, 0);
    return field;
}

LeakageInductance::ToroidalLeakageEnergy LeakageInductance::calculate_toroidal_leakage_energy(Magnetic magnetic, const std::vector<double>& currentPerWinding) {
    auto conductors = calculate_toroidal_ring_plane_conductors(magnetic, currentPerWinding);
    if (conductors.bore.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Toroidal leakage: no turn carries current");
    }
    double innerWallRadius = conductors.innerWallRadius;
    double outerWallRadius = conductors.outerWallRadius;
    double coreHeight = conductors.coreHeight;
    double imageFactor = conductors.imageFactor;
    const auto& boreConductors = conductors.bore;
    const auto& outerConductors = conductors.outer;
    const auto& turnAngles = conductors.turnAngles;
    const auto& turnCurrents = conductors.turnCurrents;
    double boreRadiusSum = 0;
    double outerRadiusSum = 0;
    double currentWeightSum = 0;
    for (size_t index = 0; index < boreConductors.size(); ++index) {
        double weight = std::abs(boreConductors[index].current);
        boreRadiusSum += weight * std::hypot(boreConductors[index].x, boreConductors[index].y);
        outerRadiusSum += weight * std::hypot(outerConductors[index].x, outerConductors[index].y);
        currentWeightSum += weight;
    }

    ToroidalLeakageEnergy result;
    result.ringPlaneEnergyPerLength = calculate_ring_plane_region_energy_per_length(boreConductors, innerWallRadius, imageFactor, true) +
                                      calculate_ring_plane_region_energy_per_length(outerConductors, outerWallRadius, imageFactor, false);

    // Envelope through the conductor centres (current-weighted mean crossing radii). Its axial offset from the core
    // faces is the mean of the radial offsets from the bore and outer walls.
    result.envelopeInnerRadius = boreRadiusSum / currentWeightSum;
    result.envelopeOuterRadius = outerRadiusSum / currentWeightSum;
    double axialOffset = ((innerWallRadius - result.envelopeInnerRadius) + (result.envelopeOuterRadius - outerWallRadius)) / 2;
    result.envelopeHeight = coreHeight + 2 * axialOffset;
    result.extrusionLength = result.envelopeHeight + (result.envelopeOuterRadius - result.envelopeInnerRadius);

    result.threeDimensionalCorrection = 0;
    for (size_t mode = 1; mode <= TOROIDAL_LEAKAGE_NUMBER_MODES; ++mode) {
        double modeEnergyPerLength = calculate_sheet_mode_energy_per_length(turnAngles, turnCurrents, mode);
        if (modeEnergyPerLength == 0) {
            continue;
        }
        double effectiveHeight = cached_body_of_revolution_effective_height(result.envelopeInnerRadius, result.envelopeOuterRadius, result.envelopeHeight, mode);
        result.threeDimensionalCorrection += modeEnergyPerLength * (effectiveHeight - result.extrusionLength);
    }
    result.energy = result.ringPlaneEnergyPerLength * result.extrusionLength + result.threeDimensionalCorrection;
    return result;
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
