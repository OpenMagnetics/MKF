#include "CoilMesher.h"
#include "WindingOhmicLosses.h"
#include "CoilWrapper.h"
#include "Utils.h"
#include "Defaults.h"
#include "Settings.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <matplot/matplot.h>
#include <streambuf>
#include <vector>


namespace OpenMagnetics {

std::shared_ptr<CoilMesherModel> CoilMesherModel::factory(CoilMesherModels modelName){
    if (modelName == CoilMesherModels::CENTER) {
        return std::make_shared<CoilMesherCenterModel>();
    }
    else if (modelName == CoilMesherModels::WANG) {
        return std::make_shared<CoilMesherWangModel>();
    }
    else
        throw std::runtime_error("Unknown coil breaker mode, available options are: {WANG, CENTER}");
}


std::vector<size_t> CoilMesher::get_common_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold) {
    auto commonHarmonicIndexes = get_main_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    if (commonHarmonicIndexes.size() > operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_amplitudes().size() * _quickModeForManyHarmonicsThreshold) {
        return get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold * 3);
    }
    else {
        return commonHarmonicIndexes;
    }
}


std::pair<Field, double> CoilMesher::generate_mesh_induced_grid(MagneticWrapper magnetic, double frequency, size_t numberPointsX, size_t numberPointsY) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();

    std::vector<FieldPoint> points;
    auto extraDimension = CoilWrapper::calculate_external_proportion_for_wires_in_toroidal_cores(magnetic.get_core(), magnetic.get_coil());
    auto bobbinWindingWindowShape = bobbin.get_winding_window_shape();
    std::vector<double> bobbinPointsX;
    std::vector<double> bobbinPointsY;
    double coreWidth = magnetic.get_mutable_core().get_width();
    double coreHeight = magnetic.get_mutable_core().get_height();
    double dA;

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        double bobbinWidthStart = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_coordinates().value()[0] - magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value() / 2;
        double bobbinWidth = magnetic.get_mutable_coil().resolve_bobbin().get_processed_description().value().get_winding_windows()[0].get_width().value();
        double coreColumnWidth = magnetic.get_mutable_core().get_columns()[0].get_width();
        double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();

        bobbinPointsX = matplot::linspace(coreColumnWidth / 2, bobbinWidthStart + bobbinWidth, numberPointsX);
        bobbinPointsY = matplot::linspace(-coreColumnHeight / 2, coreColumnHeight / 2, numberPointsY);

        double dx = (bobbinWidthStart + bobbinWidth - coreColumnWidth / 2) / numberPointsX;
        double dy = coreColumnHeight / numberPointsY;
        dA = dx * dy;
    }
    else {
        auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();
        bobbinPointsX = matplot::linspace(-coreWidth / 2 * extraDimension, coreWidth / 2 * extraDimension, numberPointsX);
        bobbinPointsY = matplot::linspace(-coreHeight / 2 * extraDimension, coreHeight / 2 * extraDimension, numberPointsY);
        double dx = coreWidth * extraDimension / numberPointsX;
        double dy = coreHeight * extraDimension / numberPointsY;
        dA = dx * dy;
    }

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

    return {inducedField, dA};
}


std::vector<Field> CoilMesher::generate_mesh_inducing_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold, std::optional<std::vector<int8_t>> customCurrentDirectionPerWinding) {
    auto defaults = Defaults();
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding does not have turns description");

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
        throw std::runtime_error("Input has no waveform. TODO: get waveform from processed data");
    }

    std::vector<std::shared_ptr<CoilMesherModel>> breakdownModelPerWinding;

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
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
            case WireType::PLANAR:
            case WireType::RECTANGULAR:
            {
                model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                break;
            }
            case WireType::FOIL: {
                model = CoilMesherModel::factory(CoilMesherModels::CENTER);
                break;
            }
            default:
                throw std::runtime_error("Unknown type of wire");
        }

        breakdownModelPerWinding.push_back(model);
    }

    auto commonHarmonicIndexes = get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    std::vector<Field> tempFieldPerHarmonic;
    for (size_t harmonicIndex = 0; harmonicIndex < operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_amplitudes().size(); ++harmonicIndex){
        Field field;
        field.set_frequency(operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_frequencies()[harmonicIndex]);
        tempFieldPerHarmonic.push_back(field);
    }


    for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
        auto turn = turns[turnIndex];
        int windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = wirePerWinding[windingIndex];
        auto harmonics = operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics().value();

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->generate_mesh_inducing_turn(turn, wire, turnIndex, turn.get_length(), magnetic.get_core());

        for (auto harmonicIndex : commonHarmonicIndexes) {
            auto harmonicCurrentPeak = harmonics.get_amplitudes()[harmonicIndex];  // Because a harmonic is always sinusoidal
            auto harmonicCurrentPeakInTurn = harmonicCurrentPeak * currentDividerPerTurn[turnIndex];
            if (std::isnan(harmonicCurrentPeakInTurn)) {
                throw std::runtime_error("NaN found in harmonicCurrentPeakInTurn value");
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

    for (auto& inducingFieldPoint : fieldPerHarmonic[0].get_data()) {
        if (std::isnan(inducingFieldPoint.get_value())) {
            std::cout << "inducingFieldPoint.get_value(): " << inducingFieldPoint.get_value() << std::endl;
            throw std::runtime_error("NaN found in inducingFieldPoint value");
        }
    }

    return fieldPerHarmonic;
}

std::vector<Field> CoilMesher::generate_mesh_induced_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold) {
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding does not have turns description");

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
            case WireType::PLANAR:
            case WireType::RECTANGULAR:
            {
                model = CoilMesherModel::factory(CoilMesherModels::WANG);
                break;
            }
            case WireType::FOIL: {
                model = CoilMesherModel::factory(CoilMesherModels::WANG);
                break;
            }
            default:
                throw std::runtime_error("Unknown type of wire");
        }

        breakdownModelPerWinding.push_back(model);
    }

    auto commonHarmonicIndexes = get_common_harmonic_indexes(operatingPoint, windingLossesHarmonicAmplitudeThreshold);

    std::vector<Field> tempFieldPerHarmonic;
    for (size_t harmonicIndex = 0; harmonicIndex < operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_amplitudes().size(); ++harmonicIndex){
        Field field;
        field.set_frequency(operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value().get_frequencies()[harmonicIndex]);
        tempFieldPerHarmonic.push_back(field);
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

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_inducing_turn(Turn turn, [[maybe_unused]] WireWrapper wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    auto mirroringDimension = settings->get_magnetic_field_mirroring_dimension();
    std::vector<FieldPoint> fieldPoints;

    int M = mirroringDimension;
    int N = mirroringDimension;

    double corePermeability = core.get_initial_permeability(Defaults().ambientTemperature);
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }

    auto processedDescription = core.get_processed_description().value();
    auto coreFamily = core.get_shape_family();

    WindingWindowElement windingWindow = processedDescription.get_winding_windows()[0]; // Hardcoded

    if (coreFamily != CoreShapeFamily::T) {
        double A = windingWindow.get_width().value();
        double B = windingWindow.get_height().value();
        double coreColumnWidth = core.get_columns()[0].get_width();

        double turn_a = turn.get_coordinates()[0] - coreColumnWidth / 2;
        double turn_b = turn.get_coordinates()[1] + B / 2;

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
                    a = m * A + turn_a;
                }
                else {
                    a = m * A + A - turn_a;
                }
                if (n % 2 == 0) {
                    b = n * B + turn_b;
                }
                else {
                    b = n * B + B - turn_b;
                }
                mirroredFieldPoint.set_point(std::vector<double>{a + coreColumnWidth / 2, b - B / 2});
                fieldPoints.push_back(mirroredFieldPoint);
            }
        }
    }
    else {
        auto mainColumn = processedDescription.get_columns()[0];

        FieldPoint fieldPoint;
        fieldPoint.set_value(1);  // Will be multiplied later
        if (!turn.get_rotation()) {
            throw std::runtime_error("Toroidal cores should have rotation in the turn, even if it is 0");
        }

        fieldPoint.set_rotation(turn.get_rotation().value());
        if (turnLength) {
            fieldPoint.set_turn_length(turnLength.value());
        }
        if (turnIndex) {
            fieldPoint.set_turn_index(turnIndex.value());
        }

        if (!turn.get_coordinate_system()) {
            throw std::runtime_error("Turn is missing coordinate system");
        }
        if (turn.get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            throw std::runtime_error("CoilMesher: Turn coordinates are not in cartesian");
        }

        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1]});
        fieldPoints.push_back(fieldPoint);
    }

    return fieldPoints;
}

std::vector<FieldPoint> CoilMesherCenterModel::generate_mesh_induced_turn(Turn turn, [[maybe_unused]] WireWrapper wire, std::optional<size_t> turnIndex) {
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

std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_induced_turn(Turn turn, WireWrapper wire, std::optional<size_t> turnIndex) {
    std::vector<FieldPoint> fieldPoints;
    FieldPoint fieldPoint;
    fieldPoint.set_value(0);
    if (turnIndex) {
        fieldPoint.set_turn_index(turnIndex.value());
    }

    // fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    // fieldPoint.set_label("right");
    // fieldPoints.push_back(fieldPoint);

    // fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    // fieldPoint.set_label("left");
    // fieldPoints.push_back(fieldPoint);

    // fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2});
    // fieldPoint.set_label("top");
    // fieldPoints.push_back(fieldPoint);

    // fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2});
    // fieldPoint.set_label("bottom");
    // fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("top right");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("top center");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1] + wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("top center");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    fieldPoint.set_label("center right");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1]});
    fieldPoint.set_label("center center");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1]});
    fieldPoint.set_label("center center");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("bottom right");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("bottom center");
    fieldPoints.push_back(fieldPoint);

    fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2, turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2});
    fieldPoint.set_label("bottom center");
    fieldPoints.push_back(fieldPoint);


    return fieldPoints;
}


std::vector<FieldPoint> CoilMesherWangModel::generate_mesh_inducing_turn(Turn turn, WireWrapper wire, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core) {
    std::vector<FieldPoint> fieldPoints;
    FieldPoint fieldPoint;
    fieldPoint.set_value(0.5);
    double c;
    double h;

    WindingWindowElement windingWindow = core.get_processed_description().value().get_winding_windows()[0]; // Hardcoded
    auto bobbinColumnShape = core.get_processed_description().value().get_winding_windows()[0].get_shape();

    if (bobbinColumnShape == WindingWindowShape::ROUND) {
        throw std::runtime_error("Wang Mesher model not implemented yet for toroidal cores");
    }

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
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0], turn.get_coordinates()[1] - wire.get_maximum_conducting_height() / 2 + w});
        fieldPoints.push_back(fieldPoint);
    }
    else if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::PLANAR) {
        fieldPoint.set_point({turn.get_coordinates()[0] + wire.get_maximum_conducting_width() / 2 - w, turn.get_coordinates()[1]});
        fieldPoints.push_back(fieldPoint);

        fieldPoint.set_point({turn.get_coordinates()[0] - wire.get_maximum_conducting_width() / 2 + w, turn.get_coordinates()[1]});
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

