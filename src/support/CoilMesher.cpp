#include "support/CoilMesher.h"
#include "physical_models/AxisymmetricWindow.h"
#include "support/MasAccess.h"
#include "physical_models/Reluctance.h"
#include "physical_models/ReluctanceNetwork.h"
#include "physical_models/WindingOhmicLosses.h"
#include "processors/Inputs.h"
#include "processors/WaveformProcessor.h"
#include "constructive_models/Coil.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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
    if (!operatingPoint.get_excitations_per_winding()[0].get_current() ||
        !operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Primary excitation current is missing harmonics in Coil Mesher");
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

bool is_inside_turns(const std::vector<Turn>& turns, double pointX, double pointY) {
    for (const auto& turn : turns) {
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

bool is_far_from_turns(const std::vector<Turn>& turns, double pointX, double pointY) {
    for (const auto& turn : turns) {
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

bool is_passed_from_all_turns(const std::vector<Turn>& turns, double pointX, double pointY) {
    for (const auto& turn : turns) {
        if (pointX < turn.get_coordinates()[0]) {
            return false;
        }
    }
    bool noTurnsAbove = true;
    for (const auto& turn : turns) {
        if (pointY < turn.get_coordinates()[1]) {
            noTurnsAbove = false;
        }
    }

    bool noTurnsBelow = true;
    for (const auto& turn : turns) {
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


std::pair<Field, double> CoilMesher::generate_mesh_core_winding_window_grid(Magnetic magnetic, double frequency, size_t numberPointsX, size_t numberPointsY) {
    if (numberPointsX == 0 || numberPointsY == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Core winding window grid needs at least one point per axis");
    }
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("needed to grid the winding window");
    }
    if (core.get_shape_family() == CoreShapeFamily::T) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core winding window grid: toroidal windows are not rectangular");
    }
    auto windingWindows = core.get_processed_description()->get_winding_windows();
    if (windingWindows.size() != 1) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core winding window grid: expected a single winding window, the core has " + std::to_string(windingWindows.size()));
    }
    if (!windingWindows[0].get_width() || !windingWindows[0].get_height()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core winding window grid: the winding window has no width or height");
    }
    // Same frame as CoilMesherCenterModel's single-window image lattice: the window starts at the main column
    // face and is centred vertically.
    double windowLeftEdgeX = core.get_columns()[0].get_width() / 2;
    double windowWidth = windingWindows[0].get_width().value();
    double windowHeight = windingWindows[0].get_height().value();
    double pixelXDimension = windowWidth / static_cast<double>(numberPointsX);
    double pixelYDimension = windowHeight / static_cast<double>(numberPointsY);

    std::vector<FieldPoint> points;
    points.reserve(numberPointsX * numberPointsY);
    for (size_t j = 0; j < numberPointsY; ++j) {
        double y = -windowHeight / 2 + (static_cast<double>(j) + 0.5) * pixelYDimension;
        for (size_t i = 0; i < numberPointsX; ++i) {
            FieldPoint fieldPoint;
            fieldPoint.set_point(std::vector<double>{windowLeftEdgeX + (static_cast<double>(i) + 0.5) * pixelXDimension, y});
            points.push_back(fieldPoint);
        }
    }
    Field inducedField;
    inducedField.set_data(points);
    inducedField.set_frequency(frequency);
    return {inducedField, pixelXDimension * pixelYDimension};
}

size_t CoilMesher::get_reference_winding_index(Coil coil) {
    auto functionalDescription = coil.get_functional_description();
    for (size_t windingIndex = 0; windingIndex < functionalDescription.size(); ++windingIndex) {
        if (functionalDescription[windingIndex].get_isolation_side() == IsolationSide::PRIMARY) {
            return windingIndex;
        }
    }
    throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
        "MAS excitation convention: the coil has no winding with isolationSide \"primary\", so there is no "
        "reference winding to refer the winding currents to");
}

std::vector<int8_t> CoilMesher::calculate_current_direction_per_winding(Coil coil) {
    // MAS excitation convention (2026-09-24): primary-side windings are passive (+ into the
    // dot), all others source (+ out of the dot). Direction by ISOLATION SIDE, not by index:
    // a second primary-side winding (push-pull Primary Half 2, a forward's demagnetisation
    // winding) is +1, not -1.
    std::vector<int8_t> currentDirectionPerWinding;
    for (auto& winding : coil.get_functional_description()) {
        currentDirectionPerWinding.push_back(winding.get_isolation_side() == IsolationSide::PRIMARY ? int8_t(1) : int8_t(-1));
    }
    return currentDirectionPerWinding;
}

std::complex<double> CoilMesher::calculate_harmonic_phasor(const Waveform& waveform, double excitationFrequency, double harmonicFrequency) {
    // The same Fourier coefficients the MAS amplitudes come from (ABT #1880), never a second,
    // resampled DFT: that one added its own resampling error and, at noise-level harmonics,
    // compared noise against noise.
    return WaveformProcessor::calculate_harmonic_phasor(waveform, excitationFrequency, harmonicFrequency, settings.get_inputs_number_points_sampled_waveforms());
}

std::vector<std::vector<double>> CoilMesher::calculate_current_phase_per_winding(Coil coil, OperatingPoint operatingPoint, const std::vector<size_t>& harmonicIndexes, std::vector<std::optional<double>>* gaugePhasePerHarmonic, std::vector<std::vector<double>>* currentAmplitudePerHarmonic) {
    size_t numberWindings = coil.get_functional_description().size();
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    if (excitations.size() < numberWindings) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "The operating point has " + std::to_string(excitations.size()) +
                                    " excitations for a coil with " + std::to_string(numberWindings) + " windings");
    }
    size_t referenceWindingIndex = get_reference_winding_index(coil);
    // Each winding's Fourier series, computed once, on first use (ABT #1880).
    std::vector<std::optional<WaveformProcessor::FourierSeries>> seriesPerWinding(numberWindings);
    // And its largest harmonic, the scale is_zero_harmonic reads, taken with it.
    std::vector<double> largestHarmonicAmplitudePerWinding(numberWindings, 0.0);

    std::vector<std::vector<double>> phasePerHarmonicPerWinding;
    for (auto harmonicIndex : harmonicIndexes) {
        std::vector<double> phasePerWinding(numberWindings, 0.0);
        std::vector<double> amplitudePerWinding(numberWindings, 0.0);
        std::vector<bool> excitedPerWinding(numberWindings, false);
        for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
            // Read in place: the current (its waveform included) and its harmonics were copied per harmonic.
            const auto& current = mas_access::current_ref(excitations[windingIndex]);
            if (!current || !mas_access::harmonics_ref(*current)) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current for winding " + std::to_string(windingIndex) + " is missing harmonics in Coil Mesher");
            }
            const auto& harmonics = mas_access::harmonics_ref(*current).value();
            if (harmonicIndex >= harmonics.get_amplitudes().size() || harmonics.get_amplitudes()[harmonicIndex] == 0) {
                continue;  // no current at this harmonic: its phase is irrelevant
            }
            double amplitude = harmonics.get_amplitudes()[harmonicIndex];
            double harmonicFrequency = harmonics.get_frequencies()[harmonicIndex];
            if (!mas_access::waveform_ref(*current)) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Current for winding " + std::to_string(windingIndex) +
                                            " has a harmonic at " + std::to_string(harmonicFrequency) +
                                            " Hz but no waveform to take its phase from");
            }
            if (!seriesPerWinding[windingIndex]) {
                auto waveform = current->get_waveform().value();
                if (!waveform.get_time()) {
                    // As MKF lists its harmonics: Inputs::calculate_harmonics_data(waveform, sampled, f).
                    waveform = Inputs::calculate_sampled_waveform(waveform, excitations[windingIndex].get_frequency());
                }
                seriesPerWinding[windingIndex] = WaveformProcessor::calculate_fourier_series(waveform, excitations[windingIndex].get_frequency(),
                                                                                             settings.get_inputs_number_points_sampled_waveforms());
                largestHarmonicAmplitudePerWinding[windingIndex] = WaveformProcessor::calculate_largest_harmonic_amplitude(seriesPerWinding[windingIndex].value());
            }
            auto phasor = WaveformProcessor::calculate_harmonic_phasor(seriesPerWinding[windingIndex].value(), excitations[windingIndex].get_frequency(), harmonicFrequency);
            double phasorAmplitude = std::abs(phasor);
            if (WaveformProcessor::is_zero_harmonic(largestHarmonicAmplitudePerWinding[windingIndex], amplitude, phasor)) {
                continue;  // zero in the list and in the waveform (a symmetric wave's even harmonics): no current, no phase (ABT #1880)
            }
            if (!(phasorAmplitude >= WaveformProcessor::harmonicConsistencyTolerance * amplitude)) {
                std::ostringstream message;
                message << std::setprecision(6) << "MAS harmonics list a component of " << amplitude << " A at " << harmonicFrequency
                        << " Hz that the current waveform of winding " << windingIndex << " does not contain (its Fourier amplitude there is "
                        << phasorAmplitude << " A, its largest harmonic "
                        << WaveformProcessor::calculate_largest_harmonic_amplitude(seriesPerWinding[windingIndex].value()) << " A): its phase is undefined";
                throw InvalidInputException(ErrorCode::INVALID_INPUT, message.str());
            }
            phasePerWinding[windingIndex] = std::arg(phasor);
            amplitudePerWinding[windingIndex] = amplitude;
            excitedPerWinding[windingIndex] = true;
        }

        std::optional<size_t> gaugeWindingIndex;
        if (excitedPerWinding[referenceWindingIndex]) {
            gaugeWindingIndex = referenceWindingIndex;
        }
        else {
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                if (excitedPerWinding[windingIndex]) {
                    gaugeWindingIndex = windingIndex;
                    break;
                }
            }
        }
        std::optional<double> gaugePhaseThisHarmonic;
        if (gaugeWindingIndex) {
            double gaugePhase = phasePerWinding[gaugeWindingIndex.value()];
            gaugePhaseThisHarmonic = gaugePhase;
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                if (excitedPerWinding[windingIndex]) {
                    phasePerWinding[windingIndex] -= gaugePhase;
                }
            }
        }
        phasePerHarmonicPerWinding.push_back(phasePerWinding);
        if (gaugePhasePerHarmonic) {
            gaugePhasePerHarmonic->push_back(gaugePhaseThisHarmonic);
        }
        if (currentAmplitudePerHarmonic) {
            currentAmplitudePerHarmonic->push_back(amplitudePerWinding);
        }
    }
    return phasePerHarmonicPerWinding;
}

std::vector<Field> CoilMesher::generate_mesh_inducing_coil(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding, std::optional<CoilMesherModels> coilMesherModel) {
    return generate_mesh_inducing_coil_phasors(magnetic, operatingPoint, windingLossesHarmonicAmplitudeThreshold, customCurrentDirectionPerWinding, coilMesherModel).fieldPerHarmonic;
}

InducingCoilMesh CoilMesher::generate_mesh_inducing_coil_phasors(Magnetic magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding, std::optional<CoilMesherModels> coilMesherModel, bool unitCurrentForSilentParallels) {
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Winding does not have turns description");

    }
    auto wirePerWinding = coil.get_wires();
    auto turns = coil.get_turns_description().value();

    std::vector<int8_t> currentDirectionPerWinding;
    if (!customCurrentDirectionPerWinding) {
        currentDirectionPerWinding = calculate_current_direction_per_winding(coil);
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
    std::map<size_t, std::vector<double>> amplitudePerWindingPerHarmonic;
    for (auto harmonicIndex : commonHarmonicIndexes) {
        amplitudePerWindingPerHarmonic[harmonicIndex].assign(coil.get_functional_description().size(), 0.0);
        Field field;
        for (const auto& excitation : operatingPoint.get_excitations_per_winding()) {
            const auto& frequencies = mas_access::harmonics_ref(*(mas_access::current_ref(excitation))).value().get_frequencies();
            if (harmonicIndex < frequencies.size()) {
                field.set_frequency(frequencies[harmonicIndex]);
                break;
            }
        }
        if (field.get_frequency() == 0) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "0 frequency found in Coil Mesher");
        }
        tempFieldPerHarmonic[harmonicIndex] = field;
    }



    // Each winding's phase and the current it carries per common harmonic: its MAS amplitude, 0 where it carries none
    // (ABT #1880: a harmonic zero both in the list and in the waveform contributes nothing, whatever its phase).
    std::vector<std::optional<double>> gaugePhasePerCommonHarmonic;
    std::vector<std::vector<double>> currentAmplitudePerCommonHarmonic;
    auto phasePerCommonHarmonic = calculate_current_phase_per_winding(coil, operatingPoint, commonHarmonicIndexes, &gaugePhasePerCommonHarmonic, &currentAmplitudePerCommonHarmonic);

    // Hoisted: the MAS getter copies the whole Core (cached material datasets included);
    // fetching it per turn dominated this loop's cost.
    auto meshCore = magnetic.get_core();
    // One window frame for the whole coil (AxisymmetricWindow::applies), as the field step and the split decide it.
    const bool axisymmetricWindow = OpenMagnetics::Settings::GetInstance().get_magnetic_field_mirroring_neutral_cells() &&
                                    AxisymmetricWindow::applies(meshCore, coil);
    for (auto& model : breakdownModelPerWinding) {
        if (auto centerModel = std::dynamic_pointer_cast<CoilMesherCenterModel>(model)) {
            centerModel->axisymmetricWindow = axisymmetricWindow;
        }
    }
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        // Read in place: the turn and the winding's current were copied per turn; the mesher models only read
        // the wire (its type), so they get the winding's wire itself.
        const auto& turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto& wire = wirePerWinding[windingIndex];
        const auto& excitationCurrent = mas_access::current_ref(operatingPoint.get_excitations_per_winding()[windingIndex]);
        if (!excitationCurrent || !excitationCurrent->get_harmonics()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Current for winding " + std::to_string(windingIndex) + " is missing harmonics in Coil Mesher");
        }

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->generate_mesh_inducing_turn(turn, wire, turnIndex, turn.get_length(), meshCore);
        // MAS classes declare a destructor, so they have no move constructor: every reallocation of a harmonic's
        // point list copies each point. Room for every turn meshed like the first one is made once.
        if (turnIndex == 0) {
            for (auto& [harmonicIndex, field] : tempFieldPerHarmonic) {
                field.get_mutable_data().reserve(turns.size() * fieldPoints.size());
            }
        }

        for (size_t commonPosition = 0; commonPosition < commonHarmonicIndexes.size(); ++commonPosition) {
            auto harmonicIndex = commonHarmonicIndexes[commonPosition];
            double harmonicCurrentPeak = currentAmplitudePerCommonHarmonic[commonPosition][windingIndex];  // Because a harmonic is always sinusoidal
            if (unitCurrentForSilentParallels && harmonicCurrentPeak == 0 && coil.get_number_parallels(windingIndex) > 1) {
                harmonicCurrentPeak = 1;
            }
            amplitudePerWindingPerHarmonic[harmonicIndex][windingIndex] = harmonicCurrentPeak;
            auto harmonicCurrentPeakInTurn = harmonicCurrentPeak * currentDividerPerTurn[turnIndex];
            if (std::isnan(harmonicCurrentPeakInTurn)) {
                throw NaNResultException("NaN found in harmonicCurrentPeakInTurn value");
            }
            harmonicCurrentPeakInTurn *= currentDirectionPerWinding[windingIndex];
            auto& harmonicData = tempFieldPerHarmonic[harmonicIndex].get_mutable_data();
            for (auto& fieldPoint : fieldPoints) {
                FieldPoint fieldPointThisHarmonic(fieldPoint);
                fieldPointThisHarmonic.set_value(fieldPoint.get_value() * harmonicCurrentPeakInTurn);
                harmonicData.push_back(std::move(fieldPointThisHarmonic));
            }
        }
    }
    std::vector<Field> fieldPerHarmonic;
    std::vector<std::vector<double>> phasePerHarmonicPerWinding;
    std::vector<std::optional<double>> gaugePhasePerHarmonic;
    std::vector<std::vector<double>> amplitudePerHarmonicPerWinding;
    std::vector<std::vector<double>> currentAmplitudePerHarmonicPerWinding;
    // In ascending harmonic order. A Field has no move constructor (see above): its points are handed over by swapping
    // the point lists, into a list with room for every harmonic.
    fieldPerHarmonic.reserve(tempFieldPerHarmonic.size());
    for (auto& [harmonicIndex, meshedField] : tempFieldPerHarmonic) {
        if (meshedField.get_data().size() > 0) {
            fieldPerHarmonic.emplace_back();
            fieldPerHarmonic.back().set_frequency(meshedField.get_frequency());
            fieldPerHarmonic.back().get_mutable_data().swap(meshedField.get_mutable_data());
            auto commonPosition = std::find(commonHarmonicIndexes.begin(), commonHarmonicIndexes.end(), harmonicIndex);
            if (commonPosition == commonHarmonicIndexes.end()) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "generate_mesh_inducing_coil: meshed harmonic " + std::to_string(harmonicIndex) +
                                           " is not one of the common harmonics");
            }
            phasePerHarmonicPerWinding.push_back(phasePerCommonHarmonic[static_cast<size_t>(commonPosition - commonHarmonicIndexes.begin())]);
            gaugePhasePerHarmonic.push_back(gaugePhasePerCommonHarmonic[static_cast<size_t>(commonPosition - commonHarmonicIndexes.begin())]);
            amplitudePerHarmonicPerWinding.push_back(amplitudePerWindingPerHarmonic.at(harmonicIndex));
            currentAmplitudePerHarmonicPerWinding.push_back(currentAmplitudePerCommonHarmonic[static_cast<size_t>(commonPosition - commonHarmonicIndexes.begin())]);
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

    return {std::move(fieldPerHarmonic), std::move(phasePerHarmonicPerWinding), std::move(gaugePhasePerHarmonic), std::move(amplitudePerHarmonicPerWinding),
            std::move(currentAmplitudePerHarmonicPerWinding)};
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
        // Read in place: each excitation (its waveforms included) was copied per harmonic.
        for (const auto& excitation : operatingPoint.get_excitations_per_winding()) {
            const auto& frequencies = mas_access::harmonics_ref(*mas_access::current_ref(excitation)).value().get_frequencies();
            if (harmonicIndex < frequencies.size()) {
                field.set_frequency(frequencies[harmonicIndex]);
                break;
            }
        }
        tempFieldPerHarmonic[harmonicIndex] = field;
    }

    // Hoisted like the inducing loop: one Core copy for the whole mesh, not one per turn.
    auto meshCore = magnetic.get_core();
    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        // Read in place, as in the inducing loop: the turn, its wire and its winding's current were copied per turn.
        const auto& turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto& wire = wirePerWinding[windingIndex];
        const auto& excitationCurrent = mas_access::current_ref(operatingPoint.get_excitations_per_winding()[windingIndex]);
        if (!excitationCurrent || !excitationCurrent->get_harmonics()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Current for winding " + std::to_string(windingIndex) + " is missing harmonics in Coil Mesher");
        }

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->generate_mesh_induced_turn(turn, wire, turnIndex, &meshCore);
        // A FieldPoint has no move constructor (MAS classes declare a destructor): room for every turn meshed like the
        // first one is made once, so the lists are not copied as they grow.
        if (turnIndex == 0) {
            for (auto& [harmonicIndex, field] : tempFieldPerHarmonic) {
                field.get_mutable_data().reserve(turns.size() * fieldPoints.size());
            }
        }

        for (auto harmonicIndex : commonHarmonicIndexes) {
            auto& harmonicData = tempFieldPerHarmonic[harmonicIndex].get_mutable_data();
            for (auto& fieldPoint : fieldPoints) {
                harmonicData.push_back(fieldPoint);
            }
        }
    }

    // In ascending harmonic order; the point lists are handed over by swapping them (a Field has no move constructor).
    std::vector<Field> fieldPerHarmonic;
    fieldPerHarmonic.reserve(tempFieldPerHarmonic.size());
    for (auto& [harmonicIndex, meshedField] : tempFieldPerHarmonic) {
        if (meshedField.get_data().size() > 0) {
            fieldPerHarmonic.emplace_back();
            fieldPerHarmonic.back().set_frequency(meshedField.get_frequency());
            fieldPerHarmonic.back().get_mutable_data().swap(meshedField.get_mutable_data());
        }
    }
    return fieldPerHarmonic;
}

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_inducing_turn(const Turn& turn, [[maybe_unused]] Wire& wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, Core& core) {
    auto mirroringDimension = settings.get_magnetic_field_mirroring_dimension();
    std::vector<FieldPoint> fieldPoints;

    int M = mirroringDimension;
    int N = mirroringDimension;
    bool completeCells = settings.get_magnetic_field_mirroring_complete_cells();

    double corePermeability = core.get_initial_permeability(defaults.ambientTemperature);
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException();
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

    if (!axisymmetricWindow) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Coil mesher: the window frame (planar or axisymmetric) was not decided for this coil");
    }
    if (coreFamily != CoreShapeFamily::T) {
        double A = windingWindow.get_width().value();
        double B = windingWindow.get_height().value();

        // Multi-window cores (ABT #2267): the window regions between x-adjacent legs, the same frames the neutral window's
        // wall sheets are imaged in. The core's per-column window entries carry no coordinates, and the regions they
        // stand for are these: each crossing is imaged in the region that contains it.
        std::vector<WindowRegion> regions;
        if (windingWindows.size() > 1) {
            regions = CoilMesher::window_regions(core);
        }
        auto region_index_containing = [&](const std::vector<double>& point) -> std::optional<size_t> {
            for (size_t regionIndex = 0; regionIndex < regions.size(); ++regionIndex) {
                const auto& region = regions[regionIndex];
                if (point[0] >= region.x0 && point[0] <= region.x1 && point[1] >= region.y0 && point[1] <= region.y1) {
                    return regionIndex;
                }
            }
            return std::nullopt;
        };
        auto region_containing = [&](const std::vector<double>& point) -> std::optional<WindowRegion> {
            auto index = region_index_containing(point);
            if (!index) {
                return std::nullopt;
            }
            return regions[index.value()];
        };

        double windowLeftEdgeX;
        double windowBottomY;
        bool outsideTheCore = false;
        if (!regions.empty()) {
            auto region = region_containing(turn.get_coordinates());
            if (!region) {
                // A conductor of the parallel split that is a lateral turn's return crossing outside the core (ABT #2267):
                // a lone filament in open air, as its additional crossing is when the whole turn is meshed.
                outsideTheCore = true;
                region = regions.front();
            }
            A = region->x1 - region->x0;
            B = region->y1 - region->y0;
            windowLeftEdgeX = region->x0;
            windowBottomY = region->y0;
        }
        else if (windingWindows.size() > 1 && windingWindow.get_coordinates()) {
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

        // On a multi-window core every point is labelled with the window it belongs to (window_label): it acts only there.
        auto appendImageLattice = [&](const std::vector<double>& crossingPoint, double frameA, double frameB,
                                      double frameLeftEdgeX, double frameBottomY, double signMultiplier, std::optional<size_t> window = std::nullopt) {
            double crossingA = crossingPoint[0] - frameLeftEdgeX;
            double crossingB = crossingPoint[1] - frameBottomY;
            // Image index m places the conductor at m*A + a (m even) or m*A + A - a (m odd), |m| reflections
            // away. The symmetric range [-M, M] cuts the 2A x 2B mirror cell in two, so the lattice keeps a net
            // dipole and its sum does not converge to the enclosed-window field: on side-by-side windings
            // (field crossing the whole window) the leakage energy read 0.68-0.85 of OMFEM 2D (ABT #1240).
            // Complete cells, m in [-2Kx-1, 2Kx], pair every image with its mirror, so the remainder decays as
            // a quadrupole lattice. That remainder is set by the lattice's physical extent, so the short side of
            // the window gets proportionally more cells: Kx = M ceil(B/A), Ky = M ceil(A/B). On a two-chamber
            // ETD 49 (B/A = 3.5) Kx = Ky = 1 read 0.86 of the converged lattice; Kx = 4 is within 1 %.
            int cellsX = M;
            int cellsY = N;
            if (completeCells) {
                cellsX = M * static_cast<int>(std::max(1.0, std::ceil(frameB / frameA)));
                cellsY = N * static_cast<int>(std::max(1.0, std::ceil(frameA / frameB)));
                // A neutral lattice (MMF_BOUNDARY) has no dipole remainder to damp: square cells converge (ABT #1304).
                if (settings.get_magnetic_field_mirroring_neutral_cells()) {
                    cellsX = M;
                    cellsY = N;
                }
            }
            int mMinimum = completeCells && cellsX > 0 ? -2 * cellsX - 1 : -cellsX;
            int mMaximum = completeCells && cellsX > 0 ? 2 * cellsX : cellsX;
            int nMinimum = completeCells && cellsY > 0 ? -2 * cellsY - 1 : -cellsY;
            int nMaximum = completeCells && cellsY > 0 ? 2 * cellsY : cellsY;
            // A round-post window on the neutral path is a solid of revolution (AxisymmetricWindow): a turn is a ring,
            // imaged exactly in the two plates only. The cylinders (post face, outer wall) have no ring image; the
            // field model meets their condition with its wall correction instead.
            if (axisymmetricWindow.value()) {
                mMinimum = 0;
                mMaximum = 0;
                nMinimum = -AxisymmetricWindow::planeImages;
                nMaximum = AxisymmetricWindow::planeImages;
            }
            for (int m = mMinimum; m <= mMaximum; ++m)
            {
                for (int n = nMinimum; n <= nMaximum; ++n)
                {
                    FieldPoint mirroredFieldPoint;
                    double currentMultiplier;
                    if (completeCells) {
                        // A line current reflected in a permeable half-space has the image I (mu - 1) / (mu + 1);
                        // the image at (m, n) is |m| + |n| reflections away and carries that factor as many times.
                        // The legacy factor below turns negative once max(|m|, |n|) exceeds mu, which the larger
                        // complete-cell lattices reach on low-permeability cores.
                        currentMultiplier = std::pow((corePermeability - 1.0) / (corePermeability + 1.0), std::abs(m) + std::abs(n));
                    }
                    else {
                        currentMultiplier = (corePermeability - std::max(fabs(m), fabs(n))) / (corePermeability + std::max(fabs(m), fabs(n)));
                    }
                    mirroredFieldPoint.set_value(signMultiplier * currentMultiplier);  // Will be multiplied later
                    // ABT #1304: an image is the core's response, not the conductor. A window-boundary fringing model
                    // (MMF_BOUNDARY) keeps a turn's own images when it evaluates that turn, so it must tell them apart.
                    // Only the Wang field model reads inducing labels, and it refuses image points outright.
                    if (m != 0 || n != 0) {
                        mirroredFieldPoint.set_label(window_label(true, window));
                    }
                    else if (window) {
                        mirroredFieldPoint.set_label(window_label(false, window));
                    }
                    if (turnLength) {
                        mirroredFieldPoint.set_turn_length(turnLength.value());
                    }
                    if (turnIndex) {
                        mirroredFieldPoint.set_turn_index(turnIndex.value());
                    }
                    double a;
                    double b;
                    if (m % 2 == 0) {
                        a = m * frameA + crossingA;
                    }
                    else {
                        a = m * frameA + frameA - crossingA;
                    }
                    if (n % 2 == 0) {
                        b = n * frameB + crossingB;
                    }
                    else {
                        b = n * frameB + frameB - crossingB;
                    }
                    mirroredFieldPoint.set_point(std::vector<double>{a + frameLeftEdgeX, b + frameBottomY});
                    fieldPoints.push_back(mirroredFieldPoint);
                }
            }
        };

        if (outsideTheCore) {
            FieldPoint conductorPoint;
            conductorPoint.set_value(1.0);  // Will be multiplied later
            if (turnLength) {
                conductorPoint.set_turn_length(turnLength.value());
            }
            if (turnIndex) {
                conductorPoint.set_turn_index(turnIndex.value());
            }
            conductorPoint.set_point(std::vector<double>{turn.get_coordinates()[0], turn.get_coordinates()[1]});
            fieldPoints.push_back(conductorPoint);
        }
        else {
            appendImageLattice(turn.get_coordinates(), A, B, windowLeftEdgeX, windowBottomY, 1.0,
                               regions.empty() ? std::nullopt : region_index_containing(turn.get_coordinates()));
        }

        // Multi-column winding: the turn's SECOND cross-section crossing carries the
        // same current with the opposite out-of-plane direction. A crossing inside
        // another window gets THAT window's mirror images (the main winding's
        // annulus crossings seen in the far window); a crossing outside the core (a
        // lateral winding's return conductor) sits in open air and gets a single
        // unmirrored filament.
        if (turn.get_additional_coordinates()) {
            auto additionalCoordinatesList = turn.get_additional_coordinates().value();
            for (auto& additionalCoordinates : additionalCoordinatesList) {
                if (additionalCoordinates.size() < 2) {
                    continue;
                }
                auto containingRegion = region_containing(additionalCoordinates);
                if (containingRegion) {
                    appendImageLattice(additionalCoordinates, containingRegion->x1 - containingRegion->x0, containingRegion->y1 - containingRegion->y0,
                                       containingRegion->x0, containingRegion->y0, -1.0, region_index_containing(additionalCoordinates));
                }
                else {
                    FieldPoint returnConductorPoint;
                    returnConductorPoint.set_value(-1.0);  // Will be multiplied later
                    if (turnLength) {
                        returnConductorPoint.set_turn_length(turnLength.value());
                    }
                    if (turnIndex) {
                        returnConductorPoint.set_turn_index(turnIndex.value());
                    }
                    returnConductorPoint.set_point(std::vector<double>{additionalCoordinates[0], additionalCoordinates[1]});
                    fieldPoints.push_back(returnConductorPoint);
                }
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

        // Collect all conductor coordinates (main + outer half for toroidal turns).
        //
        // A toroidal turn pierces this plane TWICE and the two crossings carry current in
        // OPPOSITE directions: the wire runs through the bore, over the top, back down
        // OUTSIDE the core and under the bottom, so the bore crossing and the outer
        // crossing are the go and return halves of the same loop. The direction is
        // therefore part of the geometry, not a free choice (ABT #320).
        std::vector<std::pair<std::vector<double>, double>> coordinatesToProcess;
        coordinatesToProcess.push_back({{turn.get_coordinates()[0], turn.get_coordinates()[1]}, 1.0});

        if (turn.get_additional_coordinates()) {
            auto additionalCoords = turn.get_additional_coordinates().value();
            for (const auto& coord : additionalCoords) {
                if (coord.size() >= 2) {
                    coordinatesToProcess.push_back({{coord[0], coord[1]}, -1.0});
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

        // The ring is two separate regions for the field: the bore (r < R_inner) and the outside (r > R_outer). With
        // the ferrite at mu -> inf, each region sees its own conductors plus the core's reply to them, and nothing of
        // the other region's conductors:
        //  - a bore crossing is seen in the bore with its Kelvin image in R_inner (same sign, attenuated). Outside the
        //    core it appears as a line current of its full value at the centre (the ring carries its MMF around).
        //  - an outside crossing is seen outside with its Kelvin image in R_outer (same sign, attenuated) plus a
        //    centre line of -attenuation times its current, so the outside field has zero net image current.
        // Hence outside the core the circulation is the enclosed bore current and the far field of a complete winding
        // is zero, as for a real toroid. The label carries the region; MagneticField applies each source only in its
        // region. Without mirroring the conductors are bare filaments in free space, seen everywhere.
        auto emit = [&](double pointX, double pointY, double value, std::optional<std::string> label) {
            FieldPoint point;
            point.set_value(value);
            point.set_rotation(turn.get_rotation().value());
            if (turnLength) {
                point.set_turn_length(turnLength.value());
            }
            if (turnIndex) {
                point.set_turn_index(turnIndex.value());
            }
            point.set_point({pointX, pointY});
            if (label) {
                point.set_label(label.value());
            }
            fieldPoints.push_back(point);
        };
        for (const auto& [coords, direction] : coordinatesToProcess) {
            double x = coords[0];
            double y = coords[1];
            double r_k_sq = x * x + y * y;
            if (mirroringDimension < 1) {
                emit(x, y, direction, std::nullopt);
                continue;
            }
            if (r_k_sq < R_inner * R_inner) {
                emit(x, y, direction, std::string(toroidBoreConductorLabel));
                // A crossing at the very centre has its image at infinity, which adds nothing.
                if (r_k_sq > 0) {
                    double scale = (R_inner * R_inner) / r_k_sq;
                    emit(x * scale, y * scale, direction * attenuation, std::string(toroidBoreReplyLabel));
                }
                emit(0, 0, direction, std::string(toroidExteriorReplyLabel));
            }
            else if (r_k_sq > R_outer * R_outer) {
                emit(x, y, direction, std::string(toroidExteriorConductorLabel));
                double scale = (R_outer * R_outer) / r_k_sq;
                emit(x * scale, y * scale, direction * attenuation, std::string(toroidExteriorReplyLabel));
                emit(0, 0, -direction * attenuation, std::string(toroidExteriorReplyLabel));
            }
            else {
                throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION, "CoilMesher: a crossing of turn " + turn.get_name() +
                                            " lies inside the toroid's ferrite (r = " + std::to_string(std::sqrt(r_k_sq)) + " m, between " +
                                            std::to_string(R_inner) + " and " + std::to_string(R_outer) + " m)");
            }
        }
    }

    return fieldPoints;
}

// ABT #227.2: a turn with a genuine second plane crossing (additional_coordinates --
// the out-of-core return half of a rectangular multi-column lateral turn) sees a
// materially different proximity field at each crossing. Sampling only the primary
// crossing and billing the WHOLE turn length (both halves) at that single density
// overestimates the out-of-core half's loss. The models sample both crossings and let
// WindingProximityEffectLosses weight each by its own length share instead.
// Toroidal turns too: the bore crossing and the outside crossing sit in different regions of
// the ring's field (see generate_mesh_inducing_turn) and see very different fields (the bore
// is crowded, the outside is open), so each books its own half of the turn. Shared by the
// Center AND Wang models (the Wang half of the ticket was originally missed, leaving
// foil/rectangular/planar wires on single-crossing).
// Returns the secondary crossing's {x, y} when the split applies.
static std::optional<std::vector<double>> lateral_multi_column_secondary_crossing(const Turn& turn, const Core* core) {
    if (core == nullptr) {
        return std::nullopt;
    }
    // Bind once: the MAS getter returns the nested vectors BY VALUE.
    const auto additionalCoordinates = turn.get_additional_coordinates();
    if (!additionalCoordinates || additionalCoordinates->empty() || (*additionalCoordinates)[0].size() < 2) {
        return std::nullopt;
    }
    return std::vector<double>{(*additionalCoordinates)[0][0], (*additionalCoordinates)[0][1]};
}

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_induced_turn(const Turn& turn, [[maybe_unused]] Wire& wire, std::optional<size_t> turnIndex, const Core* core) {
    std::vector<FieldPoint> fieldPoints;

    auto secondaryCrossing = lateral_multi_column_secondary_crossing(turn, core);

    FieldPoint fieldPoint;
    fieldPoint.set_point(turn.get_coordinates());
    fieldPoint.set_value(0);
    if (turnIndex) {
        fieldPoint.set_turn_index(turnIndex.value());
    }
    fieldPoint.set_label("center");
    if (secondaryCrossing) {
        fieldPoint.set_turn_length(turn.get_length() / 2);
    }
    fieldPoints.push_back(fieldPoint);

    if (secondaryCrossing) {
        FieldPoint secondaryFieldPoint;
        secondaryFieldPoint.set_point(secondaryCrossing.value());
        secondaryFieldPoint.set_value(0);
        if (turnIndex) {
            secondaryFieldPoint.set_turn_index(turnIndex.value());
        }
        secondaryFieldPoint.set_label("center_secondary");
        secondaryFieldPoint.set_turn_length(turn.get_length() / 2);
        fieldPoints.push_back(secondaryFieldPoint);
    }

    return fieldPoints;
}

std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_induced_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex, const Core* core) {
    std::vector<FieldPoint> fieldPoints;

    // A toroid's turns are laid radially: the turn's rotation is its polar angle, and its width and height are read in
    // that rotated frame (WindingProximityEffectLosses reads the field in the same frame). Elsewhere the cluster is
    // axis-aligned.
    const bool toroidal = core != nullptr && core->get_shape_family() == CoreShapeFamily::T;
    double frameAngle = 0;
    if (toroidal) {
        if (!turn.get_rotation()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Toroidal cores should have rotation in the turn, even if it is 0");
        }
        frameAngle = turn.get_rotation().value() / 180 * std::numbers::pi;
    }
    const double frameCos = std::cos(frameAngle);
    const double frameSin = std::sin(frameAngle);

    // One full surface cluster per crossing, centred on {centerX, centerY}.
    auto emitCluster = [&](double centerX, double centerY) {
        auto place = [&](double offsetX, double offsetY) {
            return std::vector<double>{centerX + offsetX * frameCos - offsetY * frameSin, centerY + offsetX * frameSin + offsetY * frameCos};
        };
        FieldPoint fieldPoint;
        fieldPoint.set_value(0);
        if (turnIndex) {
            fieldPoint.set_turn_index(turnIndex.value());
        }

        fieldPoint.set_point(place(wire.get_maximum_conducting_width() / 2, 0));
        fieldPoint.set_label("right");
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point(place(-wire.get_maximum_conducting_width() / 2, 0));
        fieldPoint.set_label("left");
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point(place(0, wire.get_maximum_conducting_height() / 2));
        fieldPoint.set_label("top");
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point(place(0, -wire.get_maximum_conducting_height() / 2));
        fieldPoint.set_label("bottom");
        fieldPoints.push_back(fieldPoint);

        // Width-resolved samples across the wide face, consumed by the Wang proximity
        // model's perpendicular-field integral (Roshen 2007 gap fringing + thin-strip
        // screening bridge). The gap-fringing field varies strongly (up to ~50x) across
        // a wide flat conductor, so the 4 lumped surface points cannot represent
        // integral(Hperp(x)^2 dx); these samples make the superposed total field
        // (proximity + fringing, including their cross term) available along the width.
        // Generated when fringing is enabled, and always on toroids: there the bore walls
        // and the neighbouring crossings make the field vary across a thick bar's width
        // with no gap at all. Elsewhere, without a gap field, the lumped mesh is
        // sufficient and much cheaper.
        if (settings.get_magnetic_field_include_fringing() || toroidal) {
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
                    fieldPoint.set_point(place(0, offset));
                }
                else {
                    fieldPoint.set_point(place(offset, 0));
                }
                fieldPoint.set_label("widthsample");
                fieldPoints.push_back(fieldPoint);
            }
        }
    };

    emitCluster(turn.get_coordinates()[0], turn.get_coordinates()[1]);

    // ABT #227.2 parity with CoilMesherCenterModel: the factory routes PLANAR/RECTANGULAR/
    // FOIL wires here, and those are exactly the wire types lateral multi-column windings
    // use — without this second cluster their whole loop length was billed at the
    // in-window crossing's field density. WindingProximityEffectLosses classifies points
    // by nearest crossing and weights each half by its own length share.
    auto secondaryCrossing = lateral_multi_column_secondary_crossing(turn, core);
    if (secondaryCrossing) {
        emitCluster(secondaryCrossing.value()[0], secondaryCrossing.value()[1]);
    }

    return fieldPoints;
}


std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_inducing_turn(const Turn& turn, Wire& wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, [[maybe_unused]] Core& core) {
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


namespace {
void append_gap_mouth(std::vector<WindowWallSheet>& sheets, double wallX, double yStart, double yEnd, double wallStart, double wallEnd, double current,
                      std::optional<size_t> winding, size_t window);
}  // namespace

std::vector<WindowWallSheet> CoilMesher::generate_window_wall_sheet_sources(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint) {
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException();
    }
    if (core.get_shape_family() == CoreShapeFamily::T) {
        throw NotImplementedException("Neutral window field: toroids have no rectangular window boundary");
    }
    auto processedDescription = core.get_processed_description().value();
    auto windingWindows = processedDescription.get_winding_windows();
    if (windingWindows.size() != 1) {
        return generate_multi_window_wall_sheet_sources(magnetic, operatingPoint);
    }
    auto windingWindow = windingWindows[0];
    if (!windingWindow.get_width() || !windingWindow.get_height()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Neutral window field: the winding window has no width/height");
    }
    auto columns = core.get_columns();
    size_t mainColumnIndex = core.get_main_column_index();

    // The frame CoilMesherCenterModel mirrors in: the window lies between the main column face and the lateral
    // column face, vertically centred.
    const double A = windingWindow.get_width().value();
    const double B = windingWindow.get_height().value();
    const double x0 = columns[mainColumnIndex].get_width() / 2;
    const double y0 = -B / 2;
    const double x1 = x0 + A;
    const double y1 = y0 + B;

    // MMF split of the magnetic circuit, exactly as MKF composes it for the inductance: the ferrite, the main
    // column's gaps in series, and the lateral branches in parallel.
    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctanceOutput = reluctanceModel->get_core_reluctance(core, operatingPoint);
    const double totalReluctance = reluctanceOutput.get_core_reluctance();
    if (!reluctanceOutput.get_ungapped_core_reluctance()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the reluctance model returned no ferrite reluctance");
    }
    const double ferriteReluctance = reluctanceOutput.get_ungapped_core_reluctance().value();
    if (!(totalReluctance > 0) || !(ferriteReluctance >= 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: non-positive core reluctance");
    }

    std::vector<WindowWallSheet> originals;
    auto gapping = core.get_functional_description().get_gapping();
    std::vector<AirGapReluctanceOutput> reluctancePerGap;
    if (!gapping.empty()) {
        if (!reluctanceOutput.get_reluctance_per_gap() || reluctanceOutput.get_reluctance_per_gap()->size() != gapping.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the reluctance model did not return one reluctance per gap");
        }
        reluctancePerGap = reluctanceOutput.get_reluctance_per_gap().value();
    }

    // Which lateral limb this window's outer wall (x1) stands for. An E-type core has its lateral limbs displaced in x
    // and the window faces the one at x > 0. The EP family's single back leg is displaced only in z (x = 0): it is the
    // whole return path of the window, and the 2-D section represents it by the outer wall x1, as the CENTER mesher's
    // image frame already does (ABT #1759: EP 17 with a lateral residual gap at x = 0 used to throw here). A core mixing
    // both kinds of lateral limb has no such single wall and is refused.
    bool anyLateralInX = false;
    bool anyLateralInZOnly = false;
    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        if (columnIndex == mainColumnIndex) {
            continue;
        }
        if (columns[columnIndex].get_coordinates()[0] == 0) {
            anyLateralInZOnly = true;
        }
        else {
            anyLateralInX = true;
        }
    }
    if (anyLateralInX && anyLateralInZOnly) {
        throw NotImplementedException("Neutral window field on a core with lateral limbs displaced both in x and only in z");
    }
    auto facesThisWindow = [&](double lateralX) {
        return anyLateralInZOnly ? lateralX == 0 : lateralX > 0;
    };

    // Lateral branches: their parallel reluctance carries the lateral MMF; on THIS window's side it drops across
    // the lateral gaps of the limb the window faces, shared among them by reluctance (series within that limb), so
    // the window's sources stay exactly neutral however MKF groups them.
    double inverseLateralReluctance = 0;
    double lateralReluctanceThisSide = 0;
    for (size_t gapIndex = 0; gapIndex < gapping.size(); ++gapIndex) {
        if (!gapping[gapIndex].get_coordinates()) {
            throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Neutral window field: a gap has no coordinates");
        }
        auto coordinates = gapping[gapIndex].get_coordinates().value();
        if (static_cast<size_t>(core.find_closest_column_index_by_coordinates(coordinates)) != mainColumnIndex) {
            inverseLateralReluctance += 1 / reluctancePerGap[gapIndex].get_reluctance();
            if (facesThisWindow(coordinates[0])) {
                lateralReluctanceThisSide += reluctancePerGap[gapIndex].get_reluctance();
            }
        }
    }
    const double lateralMmfFraction = inverseLateralReluctance > 0 ? (1 / inverseLateralReluctance) / totalReluctance : 0.0;
    if (lateralMmfFraction > 0 && !(lateralReluctanceThisSide > 0)) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Neutral window field: the core has lateral gaps but none on the lateral limb this window faces");
    }

    // A gap's slot mouth on the window boundary: its axial extent, clipped to the wall it opens onto, carrying the gap's
    // drop with the slot's own distribution (slot_mouth_boundaries).
    auto gapSheet = [&](double wallX, double gapY, double gapLength, double current) {
        double yStart = std::max(y0, gapY - gapLength / 2);
        double yEnd = std::min(y1, gapY + gapLength / 2);
        if (!(yEnd > yStart)) {
            throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Neutral window field: a gap at y = " + std::to_string(gapY) +
                               " opens nowhere onto the winding window (y in [" + std::to_string(y0) + ", " + std::to_string(y1) + "])");
        }
        append_gap_mouth(originals, wallX, yStart, yEnd, y0, y1, current, std::nullopt, 0);
    };

    double gapMmfFractionThisWindow = 0;
    for (size_t gapIndex = 0; gapIndex < gapping.size(); ++gapIndex) {
        auto coordinates = gapping[gapIndex].get_coordinates().value();
        double gapLength = gapping[gapIndex].get_length();
        double reluctance = reluctancePerGap[gapIndex].get_reluctance();
        if (static_cast<size_t>(core.find_closest_column_index_by_coordinates(coordinates)) == mainColumnIndex) {
            double fraction = reluctance / totalReluctance;
            gapSheet(x0, coordinates[1], gapLength, -fraction);
            gapMmfFractionThisWindow += fraction;
        }
        else if (facesThisWindow(coordinates[0])) {
            double fraction = lateralMmfFraction * reluctance / lateralReluctanceThisSide;
            gapSheet(x1, coordinates[1], gapLength, -fraction);
            gapMmfFractionThisWindow += fraction;
        }
    }

    // Ferrite: what Ampere leaves for the walls once the gaps are placed. Distributed by length x tangential H, with
    // H proportional to the flux through each wall's limb over that limb's section: the main column carries all of
    // it, this side's lateral limb and the two plates its area share.
    const double ferriteMmfFraction = 1 - gapMmfFractionThisWindow;
    if (ferriteMmfFraction < -1e-9) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the gaps' MMF exceeds the winding's");
    }
    if (ferriteMmfFraction > 1e-12) {
        double lateralAreaTotal = 0;
        double lateralAreaThisSide = 0;
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            if (columnIndex == mainColumnIndex) {
                continue;
            }
            lateralAreaTotal += columns[columnIndex].get_area();
            if (facesThisWindow(columns[columnIndex].get_coordinates()[0])) {
                lateralAreaThisSide += columns[columnIndex].get_area();
            }
        }
        const double mainArea = columns[mainColumnIndex].get_area();
        const double plateThickness = (processedDescription.get_height() - B) / 2;
        const double plateArea = processedDescription.get_depth() * plateThickness;
        if (!(mainArea > 0) || !(lateralAreaThisSide > 0) || !(lateralAreaTotal > 0) || !(plateArea > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Neutral window field: the core's limb sections are not all positive");
        }
        const double sideFluxShare = lateralAreaThisSide / lateralAreaTotal;
        const double weightMain = B / mainArea;
        const double weightLateral = B * sideFluxShare / lateralAreaThisSide;
        const double weightPlate = A * sideFluxShare / plateArea;
        const double weightTotal = weightMain + weightLateral + 2 * weightPlate;
        // The whole wall carries the sheet; where a gap sits on it, the gap's own sheet adds the gap's drop.
        originals.push_back({{x0, y0}, {x0, y1}, -ferriteMmfFraction * weightMain / weightTotal});
        originals.push_back({{x1, y0}, {x1, y1}, -ferriteMmfFraction * weightLateral / weightTotal});
        originals.push_back({{x0, y0}, {x1, y0}, -ferriteMmfFraction * weightPlate / weightTotal});
        originals.push_back({{x0, y1}, {x1, y1}, -ferriteMmfFraction * weightPlate / weightTotal});
    }

    double sum = 0;
    for (auto& sheet : originals) {
        sum += sheet.current;
    }
    if (std::abs(sum + 1) > 1e-9) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the wall sources return " +
                                   std::to_string(-sum) + " of the winding MMF instead of all of it");
    }
    return originals;
}

// The potential across a gap's slot mouth, from the Schwarz-Christoffel map of the slot onto the upper half t-plane
// (unit slot width, the leg corner at z = 0 = z(1), the slot's far end at t = 0, the window at t -> infinity):
//  - a slot between two leg corners, opening into the half-space: dz/dt = A sqrt(t - 1) sqrt(t + 1) / t;
//  - a slot between a leg corner and a plate face flush with the window's end wall: dz/dt = A sqrt(t - 1) / t;
// A = -j / pi. The irons are the real axis (the leg on t > 0, the other side on t < 0), so the magnetic scalar potential
// is arg(t) / pi, 0 on the leg and 1 on the other side, and the mouth (the wall line x = 0 across the slot) is where the
// equipotential ray arg t = pi psi crosses Re z = 0. Uniform pieces of equal current between those crossings reproduce
// the exact window field to 0.9 % (two corners) and 0.6 % (one corner) at a quarter of the slot width from the mouth,
// 0.2 % at one width; a uniform sheet over the whole mouth is 9 % and 14 % off there (checked against the exact field of
// the same map).
const std::vector<double>& CoilMesher::slot_mouth_boundaries(bool atWindowCorner) {
    static const auto compute = [](bool corner) {
        static_assert(slotMouthPieces % 2 == 0, "the two-corner mouth is computed on its half");
        // Gauss-Legendre nodes and weights on [-1, 1].
        constexpr size_t order = 64;
        std::vector<double> nodes(order);
        std::vector<double> weights(order);
        for (size_t i = 0; i < order; ++i) {
            double x = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(order) + 0.5));
            double derivative = 0;
            for (int iteration = 0; iteration < 100; ++iteration) {
                double p0 = 1;
                double p1 = x;
                for (size_t k = 2; k <= order; ++k) {
                    const double p2 = ((2.0 * static_cast<double>(k) - 1) * x * p1 - (static_cast<double>(k) - 1) * p0) / static_cast<double>(k);
                    p0 = p1;
                    p1 = p2;
                }
                derivative = static_cast<double>(order) * (x * p1 - p0) / (x * x - 1);
                const double step = p1 / derivative;
                x -= step;
                if (std::abs(step) < 1e-15) {
                    break;
                }
            }
            nodes[i] = x;
            weights[i] = 2 / ((1 - x * x) * derivative * derivative);
        }
        const std::complex<double> j(0, 1);
        const std::complex<double> scale = -j / std::numbers::pi;
        auto slope = [&](std::complex<double> t) {
            std::complex<double> value = scale * std::sqrt(t - 1.0) / t;
            return corner ? value : value * std::sqrt(t + 1.0);
        };
        auto integrate = [&](std::complex<double> from, std::complex<double> to) {
            std::complex<double> sum = 0;
            for (size_t i = 0; i < order; ++i) {
                sum += weights[i] * slope(from + (to - from) * (nodes[i] + 1.0) / 2.0);
            }
            return sum * (to - from) / 2.0;
        };
        // z(t) along a path inside the upper half plane: up from the corner, across, down to t.
        auto z_of = [&](std::complex<double> t) {
            const double height = t.imag() + 1;
            const std::complex<double> up(1, height);
            const std::complex<double> across(t.real(), height);
            return integrate(1.0, up) + integrate(up, across) + integrate(across, t);
        };
        auto crossing = [&](double potential) {
            const double angle = std::numbers::pi * potential;
            double low = -25;
            double high = 25;
            if (!(z_of(std::polar(std::exp(low), angle)).real() < 0) || !(z_of(std::polar(std::exp(high), angle)).real() > 0)) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Slot mouth: the equipotential " + std::to_string(potential) + " does not cross the mouth");
            }
            for (int iteration = 0; iteration < 90; ++iteration) {
                const double middle = (low + high) / 2;
                if (middle == low || middle == high) {
                    break;  // adjacent doubles: every further step would leave low and high as they are
                }
                if (z_of(std::polar(std::exp(middle), angle)).real() < 0) {
                    low = middle;
                }
                else {
                    high = middle;
                }
            }
            return z_of(std::polar(std::exp((low + high) / 2), angle)).imag();
        };
        std::vector<double> boundaries(slotMouthPieces + 1);
        boundaries.front() = 0;
        boundaries.back() = 1;
        if (corner) {
            for (size_t k = 1; k < slotMouthPieces; ++k) {
                boundaries[k] = crossing(static_cast<double>(k) / static_cast<double>(slotMouthPieces));
            }
        }
        else {
            // Symmetric about the middle of the mouth, where the potential is 1/2.
            for (size_t k = 1; k < slotMouthPieces / 2; ++k) {
                boundaries[k] = crossing(static_cast<double>(k) / static_cast<double>(slotMouthPieces));
                boundaries[slotMouthPieces - k] = 1 - boundaries[k];
            }
            boundaries[slotMouthPieces / 2] = 0.5;
        }
        for (size_t k = 0; k < slotMouthPieces; ++k) {
            if (!(boundaries[k + 1] > boundaries[k])) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Slot mouth: the pieces' boundaries are not increasing");
            }
        }
        return boundaries;
    };
    static const std::vector<double> twoCorners = compute(false);
    static const std::vector<double> oneCorner = compute(true);
    return atWindowCorner ? oneCorner : twoCorners;
}

void CoilMesher::append_slot_mouth(std::vector<WindowWallSheet>& sheets, double wallX, double yStart, double yEnd, double current, bool atWindowCorner,
                                   bool cornerAtStart, std::optional<size_t> winding, size_t window) {
    const auto& boundaries = slot_mouth_boundaries(atWindowCorner);
    const double length = yEnd - yStart;
    auto position = [&](double fraction) {
        return cornerAtStart ? yStart + fraction * length : yEnd - fraction * length;
    };
    for (size_t k = 0; k < slotMouthPieces; ++k) {
        sheets.push_back({{wallX, position(boundaries[k])}, {wallX, position(boundaries[k + 1])}, current / static_cast<double>(slotMouthPieces), winding, window});
    }
}

namespace {
// A gap's mouth on a wall running from wallStart to wallEnd (y): the slot between two leg corners when it is inside the
// wall, the slot between a leg corner and the window's end wall when it reaches one end (the gap sits at the window's
// corner, as a piece-and-plate gap against its plate). A mouth reaching both ends has neither shape and is refused.
void append_gap_mouth(std::vector<WindowWallSheet>& sheets, double wallX, double yStart, double yEnd, double wallStart, double wallEnd, double current,
                      std::optional<size_t> winding, size_t window) {
    const double tolerance = 1e-9 * (wallEnd - wallStart);
    const bool reachesStart = yStart <= wallStart + tolerance;
    const bool reachesEnd = yEnd >= wallEnd - tolerance;
    if (reachesStart && reachesEnd) {
        throw NotImplementedException("Neutral window field: a gap mouth spanning the whole window wall (y in [" + std::to_string(wallStart) + ", " +
                                      std::to_string(wallEnd) + "]) has no slot shape");
    }
    if (!reachesStart && !reachesEnd) {
        CoilMesher::append_slot_mouth(sheets, wallX, yStart, yEnd, current, false, true, winding, window);
    }
    else {
        // The leg corner is the end away from the window's end wall.
        CoilMesher::append_slot_mouth(sheets, wallX, yStart, yEnd, current, true, reachesEnd, winding, window);
    }
}
}  // namespace

std::vector<WindowRegion> CoilMesher::window_regions(Core core) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException();
    }
    auto windingWindows = core.get_processed_description()->get_winding_windows();
    if (!windingWindows[0].get_height()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Neutral window field: the winding window has no height");
    }
    const double B = windingWindows[0].get_height().value();
    auto columns = core.get_columns();
    std::vector<size_t> ordered(columns.size());
    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        ordered[columnIndex] = columnIndex;
        if (columnIndex != core.get_main_column_index() && columns[columnIndex].get_coordinates()[0] == 0) {
            throw NotImplementedException("Neutral window field on a multi-window core with a leg displaced only in z");
        }
    }
    std::sort(ordered.begin(), ordered.end(), [&](size_t left, size_t right) {
        return columns[left].get_coordinates()[0] < columns[right].get_coordinates()[0];
    });
    std::vector<WindowRegion> regions;
    for (size_t position = 0; position + 1 < ordered.size(); ++position) {
        const auto& left = columns[ordered[position]];
        const auto& right = columns[ordered[position + 1]];
        WindowRegion region{left.get_coordinates()[0] + left.get_width() / 2, -B / 2, right.get_coordinates()[0] - right.get_width() / 2, B / 2};
        if (!(region.x1 > region.x0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Neutral window field: two adjacent legs overlap");
        }
        regions.push_back(region);
    }
    return regions;
}

std::vector<WindowWallSheet> CoilMesher::generate_multi_window_wall_sheet_sources(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint) {
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException();
    }
    const auto regions = window_regions(core);
    auto reluctanceModel = ReluctanceModel::factory();
    auto reluctanceOutput = reluctanceModel->get_core_reluctance(core, operatingPoint);
    if (!reluctanceOutput.get_ungapped_core_reluctance()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the reluctance model returned no ferrite reluctance");
    }
    auto gapping = core.get_functional_description().get_gapping();
    std::vector<AirGapReluctanceOutput> reluctancePerGap;
    if (!gapping.empty()) {
        if (!reluctanceOutput.get_reluctance_per_gap() || reluctanceOutput.get_reluctance_per_gap()->size() != gapping.size()) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the reluctance model did not return one reluctance per gap");
        }
        reluctancePerGap = reluctanceOutput.get_reluctance_per_gap().value();
    }
    ReluctanceNetwork network(core, reluctanceOutput.get_ungapped_core_reluctance().value(), reluctancePerGap);
    const auto& ordered = network.get_ordered_column_indexes();
    const auto& legReluctances = network.get_column_reluctances();
    const auto& yokeReluctances = network.get_yoke_segment_reluctances();
    if (regions.size() + 1 != ordered.size() || yokeReluctances.size() != regions.size()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Neutral window field: the window regions do not match the reluctance network's loops");
    }
    // Each leg's gaps, and its ferrite share (the leg's branch reluctance without them).
    std::vector<std::vector<size_t>> gapsPerColumn(core.get_columns().size());
    std::vector<double> ferriteReluctancePerColumn = legReluctances;
    for (size_t gapIndex = 0; gapIndex < gapping.size(); ++gapIndex) {
        if (!gapping[gapIndex].get_coordinates()) {
            throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Neutral window field: a gap has no coordinates");
        }
        size_t columnIndex = static_cast<size_t>(core.find_closest_column_index_by_coordinates(gapping[gapIndex].get_coordinates().value()));
        gapsPerColumn[columnIndex].push_back(gapIndex);
        ferriteReluctancePerColumn[columnIndex] -= reluctancePerGap[gapIndex].get_reluctance();
    }

    const size_t numberWindings = magnetic.get_coil().get_functional_description().size();
    std::vector<WindowWallSheet> sources;
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        // The leg fluxes this winding drives with one ampere (its polarity and turns included, ReluctanceNetwork).
        std::vector<double> currents(numberWindings, 0.0);
        currents[windingIndex] = 1;
        const auto legFlux = network.calculate_column_magnetic_fluxes(magnetic, currents);
        // Mesh loop k circulates in the window between ordered legs k and k + 1 (ReluctanceNetwork's mesh solve): the
        // drops round it are R_k phi_k on leg k, -R_{k+1} phi_{k+1} on leg k + 1 and Y_k sum_{j <= k} phi_j on the
        // yoke segment, and they add up to the current the winding encloses in that window.
        double cumulativeFlux = 0;
        for (size_t windowIndex = 0; windowIndex < regions.size(); ++windowIndex) {
            const auto& region = regions[windowIndex];
            const size_t leftColumn = ordered[windowIndex];
            const size_t rightColumn = ordered[windowIndex + 1];
            cumulativeFlux += legFlux[leftColumn];
            auto add = [&](std::vector<double> start, std::vector<double> end, double drop) {
                if (drop != 0) {
                    sources.push_back({start, end, -drop, windingIndex, windowIndex});
                }
            };
            auto gapSheets = [&](size_t columnIndex, double wallX, double sign) {
                for (size_t gapIndex : gapsPerColumn[columnIndex]) {
                    if (gapping[gapIndex].get_type() == GapType::RESIDUAL) {
                        // A residual gap's drop is spread along its leg's face with the ferrite (it opens onto no slot).
                        continue;
                    }
                    const double gapY = gapping[gapIndex].get_coordinates().value()[1];
                    const double gapLength = gapping[gapIndex].get_length();
                    const double yStart = std::max(region.y0, gapY - gapLength / 2);
                    const double yEnd = std::min(region.y1, gapY + gapLength / 2);
                    if (!(yEnd > yStart)) {
                        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Neutral window field: a gap at y = " + std::to_string(gapY) +
                                           " opens nowhere onto its winding window");
                    }
                    const double drop = sign * reluctancePerGap[gapIndex].get_reluctance() * legFlux[columnIndex];
                    if (drop != 0) {
                        append_gap_mouth(sources, wallX, yStart, yEnd, region.y0, region.y1, -drop, windingIndex, windowIndex);
                    }
                }
            };
            auto ferriteReluctance = [&](size_t columnIndex) {
                double reluctance = ferriteReluctancePerColumn[columnIndex];
                for (size_t gapIndex : gapsPerColumn[columnIndex]) {
                    if (gapping[gapIndex].get_type() == GapType::RESIDUAL) {
                        reluctance += reluctancePerGap[gapIndex].get_reluctance();
                    }
                }
                return reluctance;
            };
            gapSheets(leftColumn, region.x0, 1.0);
            gapSheets(rightColumn, region.x1, -1.0);
            add({region.x0, region.y0}, {region.x0, region.y1}, ferriteReluctance(leftColumn) * legFlux[leftColumn]);
            add({region.x1, region.y0}, {region.x1, region.y1}, -ferriteReluctance(rightColumn) * legFlux[rightColumn]);
            const double yokeDrop = yokeReluctances[windowIndex] * cumulativeFlux;
            add({region.x0, region.y0}, {region.x1, region.y0}, yokeDrop / 2);
            add({region.x0, region.y1}, {region.x1, region.y1}, yokeDrop / 2);
        }
    }
    return sources;
}

std::vector<WindowWallSheet> CoilMesher::generate_window_wall_sheets(Magnetic magnetic, std::optional<OperatingPoint> operatingPoint) {
    auto originals = generate_window_wall_sheet_sources(magnetic, operatingPoint);
    auto core = magnetic.get_core();
    if (core.get_processed_description()->get_winding_windows().size() != 1) {
        // Each window's sheets in that window's complete-cell lattice, as CoilMesherCenterModel images the conductors
        // of the window that contains them.
        const auto regions = window_regions(core);
        const int M = settings.get_magnetic_field_mirroring_dimension();
        if (M < 1) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Neutral window field: needs the image lattice (mirroring dimension >= 1)");
        }
        const double corePermeability = core.get_initial_permeability(defaults.ambientTemperature);
        std::vector<WindowWallSheet> sheets;
        for (const auto& source : originals) {
            const auto& region = regions[source.window];
            const double A = region.x1 - region.x0;
            const double B = region.y1 - region.y0;
            const bool neutralCells = settings.get_magnetic_field_mirroring_neutral_cells();
            const int cellsX = neutralCells ? M : M * static_cast<int>(std::max(1.0, std::ceil(B / A)));
            const int cellsY = neutralCells ? M : M * static_cast<int>(std::max(1.0, std::ceil(A / B)));
            auto imagePosition = [&](const std::vector<double>& point, int m, int n) {
                double a = point[0] - region.x0;
                double b = point[1] - region.y0;
                double imageA = (m % 2 == 0) ? m * A + a : m * A + A - a;
                double imageB = (n % 2 == 0) ? n * B + b : n * B + B - b;
                return std::vector<double>{imageA + region.x0, imageB + region.y0};
            };
            for (int m = -2 * cellsX - 1; m <= 2 * cellsX; ++m) {
                for (int n = -2 * cellsY - 1; n <= 2 * cellsY; ++n) {
                    double imageFactor = std::pow((corePermeability - 1.0) / (corePermeability + 1.0), std::abs(m) + std::abs(n));
                    sheets.push_back({imagePosition(source.start, m, n), imagePosition(source.end, m, n), source.current * imageFactor, source.winding, source.window});
                }
            }
        }
        return sheets;
    }
    auto windingWindow = core.get_processed_description().value().get_winding_windows()[0];
    const double A = windingWindow.get_width().value();
    const double B = windingWindow.get_height().value();
    const double x0 = core.get_columns()[core.get_main_column_index()].get_width() / 2;
    const double y0 = -B / 2;

    // The complete-cell lattice of #1240, the same index ranges and positions CoilMesherCenterModel uses for the
    // conductors, so conductors and wall sources are imaged identically.
    const int M = settings.get_magnetic_field_mirroring_dimension();
    if (M < 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Neutral window field: needs the image lattice (mirroring dimension >= 1)");
    }
    const double corePermeability = core.get_initial_permeability(defaults.ambientTemperature);
    // Same extent as the conductors' lattice (CoilMesherCenterModel), so both are imaged identically.
    const bool neutralCells = settings.get_magnetic_field_mirroring_neutral_cells();
    const int cellsX = neutralCells ? M : M * static_cast<int>(std::max(1.0, std::ceil(B / A)));
    const int cellsY = neutralCells ? M : M * static_cast<int>(std::max(1.0, std::ceil(A / B)));
    auto imagePosition = [&](const std::vector<double>& point, int m, int n) {
        double a = point[0] - x0;
        double b = point[1] - y0;
        double imageA = (m % 2 == 0) ? m * A + a : m * A + A - a;
        double imageB = (n % 2 == 0) ? n * B + b : n * B + B - b;
        return std::vector<double>{imageA + x0, imageB + y0};
    };
    std::vector<WindowWallSheet> sheets;
    for (int m = -2 * cellsX - 1; m <= 2 * cellsX; ++m) {
        for (int n = -2 * cellsY - 1; n <= 2 * cellsY; ++n) {
            double imageFactor = std::pow((corePermeability - 1.0) / (corePermeability + 1.0), std::abs(m) + std::abs(n));
            for (auto& sheet : originals) {
                sheets.push_back({imagePosition(sheet.start, m, n), imagePosition(sheet.end, m, n), sheet.current * imageFactor});
            }
        }
    }
    return sheets;
}


} // namespace OpenMagnetics

