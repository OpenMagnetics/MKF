#include "support/Painter.h"
#include "converter_models/CurrentTransformer.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    double maximumError = 0.1;

    TEST_CASE("Test_Current_Transformer", "[converter-model][current-transformer-topology][smoke-test]") {
        json currentTransformerInputsJson;

        currentTransformerInputsJson["diodeVoltageDrop"] = 0.7;
        currentTransformerInputsJson["frequency"] = 150000;
        currentTransformerInputsJson["burdenResistor"] = 2;
        currentTransformerInputsJson["maximumDutyCycle"] = 0.9;
        currentTransformerInputsJson["maximumPrimaryCurrentPeak"] = 120;
        currentTransformerInputsJson["waveformLabel"] = "Sinusoidal";
        currentTransformerInputsJson["ambientTemperature"] = 25;
        double turnsRatio = 0.01;

        OpenMagnetics::CurrentTransformer currentTransformerInputs(currentTransformerInputsJson);

        auto inputs = currentTransformerInputs.process(turnsRatio);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset() == 0);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * maximumError * 0.01));

        REQUIRE_THAT(double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]), Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_peak().value(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * maximumError));

        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        REQUIRE_THAT(0, Catch::Matchers::WithinAbs(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * turnsRatio * maximumError * 0.01));
        REQUIRE(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() > 0);
    }

}  // namespace
