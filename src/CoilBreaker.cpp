#include "CoilBreaker.h"
#include "WindingOhmicLosses.h"
#include "Defaults.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include <../tests/TestingUtils.h>

namespace OpenMagnetics {

std::shared_ptr<CoilBreakerModel>  CoilBreakerModel::factory(CoilBreakerModels modelName){
    if (modelName == CoilBreakerModels::CENTER) {
        return std::make_shared<CoilBreakerCenterModel>();
    }
    // else if (modelName == CoilBreakerModels::ALBACH) {
    //     return std::make_shared<WindingSkinEffectLossesAlbachModel>();
    // }
    // else if (modelName == CoilBreakerModels::PAYNE) {
    //     return std::make_shared<WindingSkinEffectLossesPayneModel>();
    // }
    else
        throw std::runtime_error("Unknown coil breaker mode, available options are: {SQUARE, RECTANGULAR, CENTER, WANG}");
}

WindingWindowCurrentFieldOutput CoilBreaker::breakdown_coil(MagneticWrapper magnetic, OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold) {
    auto defaults = Defaults();
    auto coil = magnetic.get_coil();
    if (!coil.get_turns_description()) {
        throw std::runtime_error("Winding does not have turns description");

    }
    auto wirePerWinding = coil.get_wires();
    auto turns = coil.get_turns_description().value();


    auto windingLossesOutput = WindingOhmicLosses::calculate_ohmic_losses(coil, operatingPoint, defaults.ambientTemperature);
    auto currentDividerPerTurn = windingLossesOutput.get_current_divider_per_turn().value();

    if (!operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform() || 
        operatingPoint.get_excitations_per_winding()[0].get_current()->get_waveform()->get_data().size() == 0)
    {
        throw std::runtime_error("Input has no waveform. TODO: get waveform from processed data");
    }

    std::vector<std::shared_ptr<CoilBreakerModel>> breakdownModelPerWinding;

    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex)
    {
        std::shared_ptr<CoilBreakerModel> model;

        switch(coil.get_wire_type(windingIndex)) {
            case WireType::ROUND: {
                model = CoilBreakerModel::factory(CoilBreakerModels::CENTER);
                break;
            }
            case WireType::LITZ: {
                model = CoilBreakerModel::factory(CoilBreakerModels::CENTER);
                break;
            }
            case WireType::RECTANGULAR: {
                model = CoilBreakerModel::factory(CoilBreakerModels::CENTER);
                break;
            }
            case WireType::FOIL: {
                model = CoilBreakerModel::factory(CoilBreakerModels::CENTER);
                break;
            }
            default:
                throw std::runtime_error("Unknown type of wire");
        }

        std::cout << "79 _mirroringDimension: " << _mirroringDimension << std::endl;
        model->set_mirroring_dimension(_mirroringDimension);
        breakdownModelPerWinding.push_back(model);
    }



    std::vector<double> maximumHarmonicAmplitudeTimesRootFrequencyPerWinding(coil.get_functional_description().size(), 0);
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto harmonics = operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics().value();
        for (size_t harmonicIndex = 1; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            maximumHarmonicAmplitudeTimesRootFrequencyPerWinding[windingIndex] = std::max(harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex]), maximumHarmonicAmplitudeTimesRootFrequencyPerWinding[windingIndex]);
        }
    }

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

        auto fieldPoints = breakdownModelPerWinding[windingIndex]->breakdown_turn(turn, wire, 0, turnIndex, turn.get_length(), magnetic.get_core());

        for (size_t harmonicIndex = 0; harmonicIndex < harmonics.get_amplitudes().size(); ++harmonicIndex) {
            if ((harmonics.get_amplitudes()[harmonicIndex] * sqrt(harmonics.get_frequencies()[harmonicIndex])) < maximumHarmonicAmplitudeTimesRootFrequencyPerWinding[windingIndex] * windingLossesHarmonicAmplitudeThreshold) {
                continue;
            }

            auto harmonicCurrentPeak = harmonics.get_amplitudes()[harmonicIndex];  // Because a harmonic is always sinusoidal
            auto harmonicFrequency = harmonics.get_frequencies()[harmonicIndex];
            auto harmonicCurrentPeakInTurn = harmonicCurrentPeak * currentDividerPerTurn[turnIndex];
            if (windingIndex > 0) {
                harmonicCurrentPeakInTurn *= -1;
            }
            for (auto& fieldPoint : fieldPoints) {
                fieldPoint.set_value(harmonicCurrentPeakInTurn);
                tempFieldPerHarmonic[harmonicIndex].get_mutable_data().push_back(fieldPoint);
            }
            // std::move(fieldPoints.begin(), fieldPoints.end(), std::back_inserter(tempFieldPerHarmonic[harmonicIndex].get_mutable_data()));
        }
    }
    std::vector<Field> fieldPerHarmonic;
    for (size_t harmonicIndex = 0; harmonicIndex < tempFieldPerHarmonic.size(); ++harmonicIndex){
        if (tempFieldPerHarmonic[harmonicIndex].get_data().size() > 0) {
            fieldPerHarmonic.push_back(tempFieldPerHarmonic[harmonicIndex]);
        }
    }

    WindingWindowCurrentFieldOutput windingWindowCurrentFieldOutput;

    windingWindowCurrentFieldOutput.set_origin(ResultOrigin::SIMULATION);
    windingWindowCurrentFieldOutput.set_method_used("AnalyticalModels");
    windingWindowCurrentFieldOutput.set_field_per_frequency(fieldPerHarmonic);
    return windingWindowCurrentFieldOutput;
}

std::vector<FieldPoint> CoilBreakerCenterModel::breakdown_turn(Turn turn, WireWrapper wire, double currentPeak, std::optional<size_t> turnIndex, std::optional<double> turnLength, CoreWrapper core) {
    // FieldPoint fieldPoint;
    // std::cout << "turn.get_coordinates()[0]: " << turn.get_coordinates()[0] << std::endl;
    // std::cout << "turn.get_coordinates()[1]: " << turn.get_coordinates()[1] << std::endl;
    // fieldPoint.set_point(turn.get_coordinates());

    // fieldPoint.set_value(currentPeak);
    // if (turnIndex) {
    //     fieldPoint.set_turn_index(turnIndex.value());
    // }
    // if (turnLength) {
    //     fieldPoint.set_turn_length(turnLength.value());
    // }
    // fieldPoints.push_back(fieldPoint);
    std::vector<FieldPoint> fieldPoints;

    std::cout << "158 _mirroringDimension: " << _mirroringDimension << std::endl;
    int M = _mirroringDimension;
    int N = _mirroringDimension;

    WindingWindowElement windingWindow = core.get_processed_description().value().get_winding_windows()[0]; // Hardcoded
    double A = windingWindow.get_width().value();
    double B = windingWindow.get_height().value();
    double coreColumnWidth = core.get_columns()[0].get_width();

    // double turn_a = turn.get_coordinates()[0];
    // double turn_b = turn.get_coordinates()[1];
    double turn_a = turn.get_coordinates()[0] - coreColumnWidth / 2;
    double turn_b = turn.get_coordinates()[1] + B / 2;

    for (int m = -M; m <= M; ++m)
    {
        for (int n = -N; n <= N; ++n)
        {
            FieldPoint mirroredFieldPoint;
            mirroredFieldPoint.set_value(currentPeak);
            if (turnIndex) {
                mirroredFieldPoint.set_turn_index(turnIndex.value());
            }
            if (turnLength) {
                mirroredFieldPoint.set_turn_length(turnLength.value());
            }
            double a;
            double b;
            if (m % 2 == 0) {
                a = m * A + turn_a;
            }
            else {
                a = m * (1 + A) - turn_a;
            }
            if (n % 2 == 0) {
                b = n * B + turn_b;
            }
            else {
                b = n * (1 - B) - turn_b;
            }

            // if (m == 0 && n == 0) {
            //     if (turnIndex) {
            //         mirroredFieldPoint.set_turn_index(turnIndex.value());
            //     }
            // }
            // else {
            //     mirroredFieldPoint.set_turn_index(-1);
            // }
            mirroredFieldPoint.set_point(std::vector<double>{a + coreColumnWidth / 2, b - B / 2});
            // mirroredFieldPoint.set_point(std::vector<double>{a, b});
            fieldPoints.push_back(mirroredFieldPoint);
        }
    }

    std::cout << "fieldPoints.size(): " << fieldPoints.size() << std::endl;
    // std::cout << "A: " << A << std::endl;
    // std::cout << "B: " << B << std::endl;
    // std::cout << "coreColumnWidth / 2: " << coreColumnWidth / 2 << std::endl;

    // for (auto& point : fieldPoints) {
    //     std::cout << "*********************" << std::endl;
    //     OpenMagneticsTesting::print(point.get_point());
    // }
    // std::cout << "-----------------------------------------" << std::endl;

    return fieldPoints;
}



} // namespace OpenMagnetics

