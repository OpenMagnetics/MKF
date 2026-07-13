#include "support/CoilMesher.h"
#include "physical_models/WindingOhmicLosses.h"
#include "constructive_models/Coil.h"
#include "support/Utils.h"

#include <cmath>
#include <limits>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"
#include "support/Logger.h"

using OpenMagnetics::linspace;

namespace OpenMagnetics {

std::shared_ptr<CoilMesherModel> CoilMesherModel::factory(CoilMesherModels modelName){
    if (modelName == CoilMesherModels::CENTER) {
        return std::make_shared<CoilMesherCenterModel>();
    }
    else if (modelName == CoilMesherModels::WANG) {
        return std::make_shared<CoilMesherWangModel>();
    }
    else
        throw ModelNotAvailableException("Unknown coil breaker mode, available options are: {WANG, CENTER}");
}


std::vector<size_t> CoilMesher::get_common_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold) {
    auto commonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);
    if (commonHarmonicIndexes.empty()) {
        // An empty selection means every winding's current has zero AC content
        // (a DC-only excitation, e.g. a choke bias point). That is a valid
        // operating point whose AC inducing field — and thus proximity losses —
        // is exactly zero, not an error. Mesh the fundamental with its zero
        // amplitude so both the inducing and induced meshes keep one consistent
        // harmonic and downstream consumers don't see COIL_NOT_PROCESSED.
        // Only degrade this way when a fundamental with a real frequency exists;
        // otherwise the excitation is genuinely malformed and the mesher throws.
        for (auto& excitation : operatingPoint.get_excitations_per_winding()) {
            auto current = excitation.get_current();
            if (!current || !current->get_harmonics()) {
                continue;
            }
            auto frequencies = current->get_harmonics().value().get_frequencies();
            if (frequencies.size() > 1 && frequencies[1] > 0) {
                return {1};
            }
        }
        return commonHarmonicIndexes;
    }
    if (commonHarmonicIndexes.size() > operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_amplitudes().size() * _quickModeForManyHarmonicsThreshold) {
        return get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold * 3);
    }
    else {
        return commonHarmonicIndexes;
    }
}

bool is_point_inside_turn(const Turn& turn, double pointX, double pointY, double factor = 1.0) {
    double distanceX = fabs(turn.get_coordinates()[0] - pointX) * factor;
    double distanceY = fabs(turn.get_coordinates()[1] - pointY) * factor;
    if (!turn.get_dimensions()) {
        throw CoilNotProcessedException("Turns is missing dimensions, which is needed for leakage inductance calculation");
    }
    if (!turn.get_cross_sectional_shape()) {
        if (distanceX < turn.get_dimensions().value()[0] / 2 && distanceY < turn.get_dimensions().value()[1] / 2) {
            return true;
        }
    }
    else if (turn.get_cross_sectional_shape().value() == TurnCrossSectionalShape::ROUND) {
        if (hypot(distanceX, distanceY) < turn.get_dimensions().value()[0] / 2) {
            return true;
        }
    }
    else {
        if (distanceX < turn.get_dimensions().value()[0] / 2 && distanceY < turn.get_dimensions().value()[1] / 2) {
            return true;
        }
    }
    return false;
}

bool is_inside_turns(std::vector<Turn> turns, double pointX, double pointY) {
    for (auto turn : turns) {
        // Check main coordinates
        if (is_point_inside_turn(turn, pointX, pointY, settings.get_coil_mesher_inside_turns_factor())) {
            return true;
        }
        // Check additional coordinates for toroidal turns (outer half)
        if (turn.get_additional_coordinates()) {
            auto additionalCoords = turn.get_additional_coordinates().value();
            for (const auto& coord : additionalCoords) {
                // Create a temporary turn with additional coordinates to check
                Turn tempTurn = turn;
                tempTurn.set_coordinates(coord);
                if (is_point_inside_turn(tempTurn, pointX, pointY, settings.get_coil_mesher_inside_turns_factor())) {
                    return true;
                }
            }
        }
    }

    return false;
}

bool is_far_from_turns(std::vector<Turn> turns, double pointX, double pointY) {
    for (auto turn : turns) {
        // Check distance to main coordinates
        double distanceX = fabs(turn.get_coordinates()[0] - pointX);
        double distanceY = fabs(turn.get_coordinates()[1] - pointY);
        if (!turn.get_dimensions()) {
            throw CoilNotProcessedException("Turns is missing dimensions, which is needed for leakage inductance calculation");
        }
        if (hypot(distanceX, distanceY) < std::max(turn.get_dimensions().value()[0], turn.get_dimensions().value()[1]) * 2) {
            return false;
        }
        // Check distance to additional coordinates for toroidal turns (outer half)
        if (turn.get_additional_coordinates()) {
            auto additionalCoords = turn.get_additional_coordinates().value();
            for (const auto& coord : additionalCoords) {
                double addDistanceX = fabs(coord[0] - pointX);
                double addDistanceY = fabs(coord[1] - pointY);
                if (hypot(addDistanceX, addDistanceY) < std::max(turn.get_dimensions().value()[0], turn.get_dimensions().value()[1]) * 2) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool is_passed_from_all_turns(std::vector<Turn> turns, double pointX, double pointY) {
    for (auto turn : turns) {
        if (pointX < turn.get_coordinates()[0]) {
            return false;
        }
    }
    bool noTurnsAbove = true;
    for (auto turn : turns) {
        if (pointY < turn.get_coordinates()[1]) {
            noTurnsAbove = false;
        }
    }

    bool noTurnsBelow = true;
    for (auto turn : turns) {
        if (pointY > turn.get_coordinates()[1]) {
            noTurnsBelow = false;
        }
    }

    return noTurnsBelow ^ noTurnsAbove;
}

std::pair<Field, double> CoilMesher::generate_mesh_induced_grid(Magnetic magnetic, double frequency, size_t numberPointsX, size_t numberPointsY, bool ignoreTurns, bool includeInsideTurns, bool meshAllWindows) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();

    std::vector<FieldPoint> points;
    auto extraDimension = Coil::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    // One (xPoints, yPoints) grid per meshed region.
    std::vector<std::pair<std::vector<double>, std::vector<double>>> regionGrids;
    double coreWidth = magnetic.get_mutable_core().get_width();
    double coreHeight = magnetic.get_mutable_core().get_height();
    double dA;

    auto bobbinWindingWindows = bobbin.get_processed_description().value().get_winding_windows();
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR &&
        meshAllWindows && bobbinWindingWindows.size() > 1) {
        // Multi-column cores: mesh every DISTINCT winding-window region. Windows
        // sharing one region (the same region wound from two different columns)
        // are meshed once; the regions of today's multi-window cores are congruent
        // (mirrors or duplicates of window 0), which keeps dA uniform.
        std::vector<std::vector<double>> meshedRegionCenters;
        std::optional<double> commonPixelArea;
        for (auto& windingWindow : bobbinWindingWindows) {
            if (!windingWindow.get_coordinates() || !windingWindow.get_width() || !windingWindow.get_height()) {
                continue;
            }
            auto regionCenter = windingWindow.get_coordinates().value();
            bool alreadyMeshed = false;
            for (auto& meshedCenter : meshedRegionCenters) {
                if (std::abs(meshedCenter[0] - regionCenter[0]) < 1e-9 && std::abs(meshedCenter[1] - regionCenter[1]) < 1e-9) {
                    alreadyMeshed = true;
                    break;
                }
            }
            if (alreadyMeshed) {
                continue;
            }
            meshedRegionCenters.push_back(regionCenter);
            double regionWidth = windingWindow.get_width().value();
            double regionHeight = windingWindow.get_height().value();
            double pixelXDimension = regionWidth / numberPointsX;
            double pixelYDimension = regionHeight / numberPointsY;
            double pixelArea = pixelXDimension * pixelYDimension;
            if (!commonPixelArea) {
                commonPixelArea = pixelArea;
            }
            else if (std::abs(pixelArea - commonPixelArea.value()) > commonPixelArea.value() * 1e-6) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                                            "Cannot mesh winding windows of different sizes with a uniform grid");
            }
            regionGrids.push_back({linspace(regionCenter[0] - regionWidth / 2 + pixelXDimension / 2,
                                            regionCenter[0] + regionWidth / 2 - pixelXDimension / 2, numberPointsX),
                                   linspace(regionCenter[1] - regionHeight / 2 + pixelYDimension / 2,
                                            regionCenter[1] + regionHeight / 2 - pixelYDimension / 2, numberPointsY)});
        }
        if (!commonPixelArea) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "No meshable winding window found");
        }
        dA = commonPixelArea.value();
    }
    else if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double bobbinWidthStart = bobbinWindingWindows[0].get_coordinates().value()[0] - bobbinWindingWindows[0].get_width().value() / 2;
        double bobbinWidth = bobbinWindingWindows[0].get_width().value();
        double coreColumnWidth = magnetic.get_mutable_core().get_columns()[0].get_width();
        double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();

        double totalWidthInGrid = bobbinWidthStart + bobbinWidth - coreColumnWidth / 2;
        double pixelXDimension = totalWidthInGrid / numberPointsX;
        double pixelYDimension = coreColumnHeight / numberPointsY;

        regionGrids.push_back({linspace(coreColumnWidth / 2 + pixelXDimension / 2, bobbinWidthStart + bobbinWidth - pixelXDimension / 2, numberPointsX),
                               linspace(-coreColumnHeight / 2 + pixelYDimension / 2, coreColumnHeight / 2 - pixelYDimension / 2, numberPointsY)});

        double dx = totalWidthInGrid / numberPointsX;
        double dy = coreColumnHeight / numberPointsY;
        dA = dx * dy;
    }
    else {
        // For toroidal cores with round winding windows, use standard Cartesian grid
        // This provides consistent behavior across all core types
        regionGrids.push_back({linspace(-coreWidth / 2 * extraDimension, coreWidth / 2 * extraDimension, numberPointsX),
                               linspace(-coreHeight / 2 * extraDimension, coreHeight / 2 * extraDimension, numberPointsY)});
        double dx = coreWidth * extraDimension / numberPointsX;
        double dy = coreHeight * extraDimension / numberPointsY;
        dA = dx * dy;
    }
    auto isPlanar = magnetic.get_wires()[0].get_type() == WireType::PLANAR;
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");

    }
    auto turns = coil.get_turns_description().value();
    auto windingOrientation = bobbin.get_winding_orientation();
    bool checkOnlyDistance = true;
    if (windingOrientation && windingOrientation.value() == WindingOrientation::CONTIGUOUS) {
        checkOnlyDistance = false;
    }
    for (auto& [bobbinPointsX, bobbinPointsY] : regionGrids) {
        for (size_t j = 0; j < bobbinPointsY.size(); ++j) {
            for (size_t i = 0; i < bobbinPointsX.size(); ++i) {
                if (!ignoreTurns) {
                    if (is_far_from_turns(turns, bobbinPointsX[i], bobbinPointsY[j]) && (checkOnlyDistance || is_passed_from_all_turns(turns, bobbinPointsX[i], bobbinPointsY[j]))) {
                        continue;
                    }
                }
                if (isPlanar) {
                    if (!ignoreTurns) {
                        // Planar are so thin and can be so close, that we need to remove he copper part to avoid having a much larger value
                        // TODO: Evaluate for other wires
                        if (is_inside_turns(turns, bobbinPointsX[i], bobbinPointsY[j])) {
                            continue;
                        }
                    }
                }
                else if (!includeInsideTurns) {
                    if (is_inside_turns(turns, bobbinPointsX[i], bobbinPointsY[j])) {
                        continue;
                    }
                }
                FieldPoint fieldPoint;
                fieldPoint.set_point(std::vector<double>{bobbinPointsX[i], bobbinPointsY[j]});
                points.push_back(fieldPoint);
            }
        }
    }
    Field inducedField;
    inducedField.set_data(points);
    inducedField.set_frequency(frequency);

    return {inducedField, dA};
}


std::vector<Field> CoilMesher::generate_mesh_inducing_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding, std::optional<CoilMesherModels> coilMesherModel) {
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");

    }
    auto wirePerWinding = coil.get_wires();
    auto turns = coil.get_turns_description().value();

    std::vector<int8_t> currentDirectionPerWinding;
    if (!customCurrentDirectionPerWinding) {
        currentDirectionPerWinding.push_back(1);
        for (size_t windingIndex = 1; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
            currentDirectionPerWinding.push_back(-1);
        }
    }
    else {
        currentDirectionPerWinding = customCurrentDirectionPerWinding.value();
    }

    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, defaults.ambientTemperature);
    auto currentDividerPerTurn = windingLossesOutput.get_current_divider_per_turn().value();

    if (!operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform() || 
        operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform()->get_data().size() == 0)
    {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Input has no waveform. TODO: get waveform from processed data");
    }

    std::vector<std::shared_ptr<CoilMesherModel>> breakdownModelPerWinding;

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        std::shared_ptr<CoilMesherModel> model;

        if (coilMesherModel) {
            model = CoilMesherModel::factory(coilMesherModel.value());
        }
        else {
            switch(coil.get_wire_type(windingIndex)) {
                case WireType::ROUND: {
                    model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                    break;
                }
                case WireType::LITZ: {
                    model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                    break;
                }
                case WireType::PLANAR: {
                    model = CoilMesherModel::factory(CoilMesherModels::WANG);
                    break;
                }
                case WireType::RECTANGULAR: {
                    model = CoilMesherModel::factory(CoilMesherModels::WANG);
                    break;
                }
                case WireType::FOIL: {
                    model = CoilMesherModel::factory(CoilMesherModels::WANG);
                    break;
                }
                default:
                    throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
            }
        }

        breakdownModelPerWinding.push_back(model);
    }

    auto commonHarmonicIndexes = get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    std::map<size_t, Field> tempFieldPerHarmonic;
    for (auto harmonicIndex : commonHarmonicIndexes) {
        Field field;
        for (auto excitation : operatingPoint.get_excitations_per_winding()) {
            if (harmonicIndex < excitation.get_current()->get_harmonics().value().get_frequencies().size()) {
                field.set_frequency(excitation.get_current()->get_harmonics().value().get_frequencies()[harmonicIndex]);
                break;
            }
        }
        if (field.get_frequency() == 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "0 frequency found in Coil Mesher");
        }
        tempFieldPerHarmonic[harmonicIndex] = field;
    }



    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = wirePerWinding[windingIndex];
        auto harmonics = operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics().value();

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->generate_mesh_inducing_turn(turn, wire, turnIndex, turn.get_length(), magnetic.get_core());

        for (auto harmonicIndex : commonHarmonicIndexes) {
            double harmonicCurrentPeak = 0;
            if (harmonicIndex < harmonics.get_amplitudes().size()) {
                harmonicCurrentPeak = harmonics.get_amplitudes()[harmonicIndex];  // Because a harmonic is always sinusoidal
            }
            auto harmonicCurrentPeakInTurn = harmonicCurrentPeak * currentDividerPerTurn[turnIndex];
            if (std::isnan(harmonicCurrentPeakInTurn)) {
                throw NaNResultException("NaN found in harmonicCurrentPeakInTurn value");
            }
            harmonicCurrentPeakInTurn *= currentDirectionPerWinding[windingIndex];
            for (auto& fieldPoint : fieldPoints) {
                FieldPoint fieldPointThisHarmonic(fieldPoint);
                fieldPointThisHarmonic.set_value(fieldPoint.get_value() * harmonicCurrentPeakInTurn);
                tempFieldPerHarmonic[harmonicIndex].get_mutable_data().push_back(fieldPointThisHarmonic);
            }
        }
    }
    std::vector<Field> fieldPerHarmonic;
    for (size_t harmonicIndex = 0; harmonicIndex < tempFieldPerHarmonic.size(); ++harmonicIndex){
        if (tempFieldPerHarmonic[harmonicIndex].get_data().size() > 0) {
            fieldPerHarmonic.push_back(tempFieldPerHarmonic[harmonicIndex]);
        }
    }

    // Guard the harmonic[0] access — if every tempFieldPerHarmonic
    // entry had zero data the vector is empty and the loop below
    // SEGVs. Triggered by isolated_buck whose primary-loop excitation
    // produced no usable harmonic content with the picked core.
    if (fieldPerHarmonic.empty()) {
        throw CoilNotProcessedException(
            "generate_mesh_inducing_coil: no harmonic produced inducing-field data "
            "(every tempFieldPerHarmonic[i].get_data() was empty) — the picked "
            "core / operating-point combination does not excite the mesh");
    }

    for (auto& inducingFieldPoint : fieldPerHarmonic[0].get_data()) {

        if (std::isnan(inducingFieldPoint.get_value())) {
            throw NaNResultException("NaN found in inducingFieldPoint value");
        }
    }

    return fieldPerHarmonic;
}

std::vector<Field> CoilMesher::generate_mesh_induced_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold) {
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");

    }
    auto wirePerWinding = coil.get_wires();
    auto turns = coil.get_turns_description().value();

    std::vector<std::shared_ptr<CoilMesherModel>> breakdownModelPerWinding;

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex)
    {
        std::shared_ptr<CoilMesherModel> model;

        switch(coil.get_wire_type(windingIndex)) {
            case WireType::ROUND: {
                model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                break;
            }
            case WireType::LITZ: {
                model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                break;
            }
            case WireType::PLANAR: {
                model = CoilMesherModel::factory(CoilMesherModels::WANG);
                break;
            }
            case WireType::RECTANGULAR: {
                model = CoilMesherModel::factory(CoilMesherModels::WANG);
                break;
            }
            case WireType::FOIL: {
                model = CoilMesherModel::factory(CoilMesherModels::WANG);
                break;
            }
            default:
                throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Unknown type of wire");
        }

        breakdownModelPerWinding.push_back(model);
    }

    auto commonHarmonicIndexes = get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    std::map<size_t, Field> tempFieldPerHarmonic;
    for (auto harmonicIndex : commonHarmonicIndexes) {
        Field field;
        for (auto excitation : operatingPoint.get_excitations_per_winding()) {
            if (harmonicIndex < excitation.get_current()->get_harmonics().value().get_frequencies().size()) {
                field.set_frequency(excitation.get_current()->get_harmonics().value().get_frequencies()[harmonicIndex]);
                break;
            }
        }
        tempFieldPerHarmonic[harmonicIndex] = field;
    }

    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = wirePerWinding[windingIndex];
        auto harmonics = operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics().value();

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->generate_mesh_induced_turn(turn, wire, turnIndex);

        for (auto harmonicIndex : commonHarmonicIndexes) {
            for (auto& fieldPoint : fieldPoints) {
                tempFieldPerHarmonic[harmonicIndex].get_mutable_data().push_back(fieldPoint);
            }
        }
    }

    std::vector<Field> fieldPerHarmonic;
    for (size_t harmonicIndex = 0; harmonicIndex < tempFieldPerHarmonic.size(); ++harmonicIndex){
        if (tempFieldPerHarmonic[harmonicIndex].get_data().size() > 0) {
            fieldPerHarmonic.push_back(tempFieldPerHarmonic[harmonicIndex]);
        }
    }
    return fieldPerHarmonic;
}

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_inducing_turn(Turn turn, [[maybe_unused]] Wire wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, Core core) {
    auto mirroringDimension = settings.get_magnetic_field_mirroring_dimension();
    std::vector<FieldPoint> fieldPoints;

    int M = mirroringDimension;
    int N = mirroringDimension;

    double corePermeability = core.get_initial_permeability(defaults.ambientTemperature);
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core is not processed");
    }

    auto processedDescription = core.get_processed_description().value();
    auto coreFamily = core.get_shape_family();

    // Mirror against the winding window that CONTAINS the turn (multi-column winding
    // support). Single-window cores keep the historical frame: window 0's inner edge
    // coincides with the main column surface, so the expressions below reduce exactly
    // to the old columnWidth/2-based ones.
    auto windingWindows = processedDescription.get_winding_windows();
    WindingWindowElement windingWindow = windingWindows[0];
    if (windingWindows.size() > 1 && coreFamily != CoreShapeFamily::T) {
        double turnX = turn.get_coordinates()[0];
        double bestDistance = std::numeric_limits<double>::max();
        for (auto& candidate : windingWindows) {
            if (!candidate.get_coordinates() || !candidate.get_width()) {
                continue;
            }
            double distance = std::abs(turnX - candidate.get_coordinates().value()[0]);
            if (distance < bestDistance) {
                bestDistance = distance;
                windingWindow = candidate;
            }
        }
    }

    if (coreFamily != CoreShapeFamily::T) {
        double A = windingWindow.get_width().value();
        double B = windingWindow.get_height().value();

        double windowLeftEdgeX;
        double windowBottomY;
        if (windingWindows.size() > 1 && windingWindow.get_coordinates()) {
            windowLeftEdgeX = windingWindow.get_coordinates().value()[0] - A / 2;
            windowBottomY = windingWindow.get_coordinates().value()[1] - B / 2;
        }
        else {
            // Historical single-window frame: window bounded left by the main column
            // surface, vertically centered. Kept by construction so single-window
            // magnetics are bit-identical to the pre-multi-column behavior.
            double coreColumnWidth = core.get_columns()[0].get_width();
            windowLeftEdgeX = coreColumnWidth / 2;
            windowBottomY = -B / 2;
        }

        double turnA = turn.get_coordinates()[0] - windowLeftEdgeX;
        double turnB = turn.get_coordinates()[1] - windowBottomY;

        for (int m = -M; m <= M; ++m)
        {
            for (int n = -N; n <= N; ++n)
            {
                FieldPoint mirroredFieldPoint;
                double currentMultiplier = (corePermeability - std::max(fabs(m), fabs(n))) / (corePermeability + std::max(fabs(m), fabs(n)));
                mirroredFieldPoint.set_value(currentMultiplier);  // Will be multiplied later
                if (turnLength) {
                    mirroredFieldPoint.set_turn_length(turnLength.value());
                }
                if (turnIndex) {
                    mirroredFieldPoint.set_turn_index(turnIndex.value());
                }
                double a;
                double b;
                if (m % 2 == 0) {
                    a = m * A + turnA;
                }
                else {
                    a = m * A + A - turnA;
                }
                if (n % 2 == 0) {
                    b = n * B + turnB;
                }
                else {
                    b = n * B + B - turnB;
                }
                mirroredFieldPoint.set_point(std::vector<double>{a + windowLeftEdgeX, b + windowBottomY});
                fieldPoints.push_back(mirroredFieldPoint);
            }
        }
    }
    else {
        // Toroidal cores: use Kelvin transformation (circle inversion) to enforce
        // high-μ boundary conditions at inner and outer core walls.
        // Based on Mühlethaler (2026) "Winding Loss Model for Toroidal Core Magnetics
        // Using the Mirroring Method", extending the classical method of images with
        // Kelvin inversion for circular boundaries.
        //
        // For a conductor at position z_k inside a circle of radius R:
        //   z_image = R² / conj(z_k) = (R²/r_k) * exp(j*θ_k)
        // Image current has SAME sign (magnetostatics, not electrostatics).

        auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        double R_inner = dimensions["B"] / 2.0;  // Inner boundary (winding window wall)
        double R_outer = dimensions["A"] / 2.0;  // Outer boundary (core outer wall)
        double attenuation = (corePermeability - 1.0) / (corePermeability + 1.0);

        // Collect all conductor coordinates (main + outer half for toroidal turns)
        std::vector<std::vector<double>> coordinatesToProcess;
        coordinatesToProcess.push_back({turn.get_coordinates()[0], turn.get_coordinates()[1]});

        if (turn.get_additional_coordinates()) {
            auto additionalCoords = turn.get_additional_coordinates().value();
            for (const auto& coord : additionalCoords) {
                if (coord.size() >= 2) {
                    coordinatesToProcess.push_back({coord[0], coord[1]});
                }
            }
        }

        if (!turn.get_rotation()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal cores should have rotation in the turn, even if it is 0");
        }
        if (!turn.get_coordinate_system()) {
            throw CoilNotProcessedException("Turn is missing coordinate system");
        }
        if (turn.get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            throw CoilNotProcessedException("CoilMesher: Turn coordinates are not in cartesian");
        }

        for (const auto& coords : coordinatesToProcess) {
            double x = coords[0];
            double y = coords[1];
            double r_k_sq = x * x + y * y;

            // Real conductor
            FieldPoint realPoint;
            realPoint.set_value(1.0);
            realPoint.set_rotation(turn.get_rotation().value());
            if (turnLength) {
                realPoint.set_turn_length(turnLength.value());
            }
            if (turnIndex) {
                realPoint.set_turn_index(turnIndex.value());
            }
            realPoint.set_point({x, y});
            fieldPoints.push_back(realPoint);

            // Kelvin images (only if mirroring is enabled and conductor is not at origin)
            if (mirroringDimension >= 1 && r_k_sq > 1e-24) {
                // Inner boundary image: z_img = R_inner² / conj(z_k)
                // In Cartesian: scale = R_inner² / |z_k|²
                double scale_inner = (R_inner * R_inner) / r_k_sq;
                FieldPoint innerImage;
                innerImage.set_value(attenuation);
                innerImage.set_rotation(turn.get_rotation().value());
                if (turnLength) {
                    innerImage.set_turn_length(turnLength.value());
                }
                if (turnIndex) {
                    innerImage.set_turn_index(turnIndex.value());
                }
                innerImage.set_point({x * scale_inner, y * scale_inner});
                fieldPoints.push_back(innerImage);

                // Outer boundary image: z_img = R_outer² / conj(z_k)
                double scale_outer = (R_outer * R_outer) / r_k_sq;
                FieldPoint outerImage;
                outerImage.set_value(attenuation);
                outerImage.set_rotation(turn.get_rotation().value());
                if (turnLength) {
                    outerImage.set_turn_length(turnLength.value());
                }
                if (turnIndex) {
                    outerImage.set_turn_index(turnIndex.value());
                }
                outerImage.set_point({x * scale_outer, y * scale_outer});
                fieldPoints.push_back(outerImage);
            }
        }
    }

    return fieldPoints;
}

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_induced_turn(Turn turn, [[maybe_unused]] Wire wire, std::optional<size_t> turnIndex) {
    std::vector<FieldPoint> fieldPoints;
    FieldPoint fieldPoint;

    fieldPoint.set_point(turn.get_coordinates());

    fieldPoint.set_value(0);
    if (turnIndex) {
        fieldPoint.set_turn_index(turnIndex.value());
    }
    fieldPoint.set_label("center");
    fieldPoints.push_back(fieldPoint);

    return fieldPoints;
}

std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_induced_turn(Turn turn, Wire wire, std::optional<size_t> turnIndex) {
    std::vector<FieldPoint> fieldPoints;
    FieldPoint fieldPoint;
    fieldPoint.set_value(0);
    if (turnIndex) {
        fieldPoint.set_turn_index(turnIndex.value());
    }

    fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    fieldPoint.set_label("right");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    fieldPoint.set_label("left");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("top");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("bottom");
    fieldPoints.push_back(fieldPoint);

    // Width-resolved samples across the wide face, consumed by the Wang proximity
    // model's perpendicular-field integral (Roshen 2007 gap fringing + thin-strip
    // screening bridge). The gap-fringing field varies strongly (up to ~50x) across
    // a wide flat conductor, so the 4 lumped surface points cannot represent
    // integral(Hperp(x)^2 dx); these samples make the superposed total field
    // (proximity + fringing, including their cross term) available along the width.
    // Only generated when fringing is enabled: without a gap field the lumped
    // mesh is sufficient and much cheaper.
    if (settings.get_magnetic_field_include_fringing()) {
        bool wideAlongY = (wire.get_type() == WireType::FOIL);
        double wideDimension = wideAlongY ? wire.get_maximum_conducting_height() : wire.get_maximum_conducting_width();
        double thinDimension = wideAlongY ? wire.get_maximum_conducting_width() : wire.get_maximum_conducting_height();
        // Self-scaling sample count: enough to resolve the near-gap field decay on
        // wide traces, cheap on near-square conductors.
        size_t numberSamples = std::min(size_t(32), std::max(size_t(8), size_t(std::round(wideDimension / thinDimension))));
        double sampleStep = wideDimension / double(numberSamples);
        for (size_t sampleIndex = 0; sampleIndex < numberSamples; ++sampleIndex) {
            double offset = -wideDimension / 2 + (sampleIndex + 0.5) * sampleStep;
            if (wideAlongY) {
                fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + offset});
            }
            else {
                fieldPoint.set_point({turn.get_coordinates()[0] + offset, turn.get_coordinates()[1]});
            }
            fieldPoint.set_label("widthsample");
            fieldPoints.push_back(fieldPoint);
        }
    }

    return fieldPoints;
}


std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_inducing_turn(Turn turn, Wire wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, Core core) {
    std::vector<FieldPoint> fieldPoints;
    FieldPoint fieldPoint;
    fieldPoint.set_value(1);
    double c;
    double h;

    WindingWindowElement windingWindow = core.get_processed_description().value().get_winding_windows()[0]; // Hardcoded
    auto bobbinColumnShape = core.get_processed_description().value().get_winding_windows()[0].get_shape();

    if (bobbinColumnShape == WindingWindowShape::ROUND) {
        throw NotImplementedException("Wang Mesher model not implemented yet for toroidal cores");
    }
    // if (wire.get_type() != WireType::PLANAR) {
    //     throw std::runtime_error("Wang model only supported for Planar");
    // }

    if (wire.get_type() == WireType::FOIL) {
        c = wire.get_maximum_conducting_height();
        h = wire.get_maximum_conducting_width();
    }
    else {
        c = wire.get_maximum_conducting_width();
        h = wire.get_maximum_conducting_height();
    }
    double k = c / h;
    double lambda = std::min(0.99, 0.01 * k + 0.66);
    double w = lambda * h;

    if (turnIndex) {
        fieldPoint.set_turn_index(turnIndex.value());
    }
    if (turnLength) {
        fieldPoint.set_turn_length(turnLength.value());
    }

    if (wire.get_type() == WireType::FOIL) {
        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2 - w});
        fieldPoint.set_label("top");
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2 + w});
        fieldPoint.set_label("bottom");
        fieldPoints.push_back(fieldPoint);
    }
    else if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::PLANAR) {
        fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2 - w, turn.get_coordinates()[1]});
        fieldPoint.set_label("right");
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2 + w, turn.get_coordinates()[1]});
        fieldPoint.set_label("left");
        fieldPoints.push_back(fieldPoint);
    }
    else {
        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2 - w});
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2 + w});
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2 - w, turn.get_coordinates()[1]});
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2 + w, turn.get_coordinates()[1]});
        fieldPoints.push_back(fieldPoint);
    }

    return fieldPoints;
}

} // namespace OpenMagnetics

