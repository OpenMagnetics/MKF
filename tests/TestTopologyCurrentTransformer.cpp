#include "support/Painter.h"
#include "converter_models/CurrentTransformer.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

SUITE(TopologyCurrentTransformer) {
    double maximumError = 0.1;

    TEST(Test_Current_Transformer) {
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
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding().size() == 2);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        CHECK_EQUAL(inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), 0);
        CHECK_CLOSE(0, inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_offset(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * maximumError * 0.01);

        CHECK_CLOSE(double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]), inputs.get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_peak().value(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * maximumError);

        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_voltage()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_label() == WaveformLabel::SINUSOIDAL);
        CHECK_CLOSE(0, inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_offset(), double(currentTransformerInputsJson["maximumPrimaryCurrentPeak"]) * turnsRatio * maximumError * 0.01);
        CHECK(inputs.get_operating_points()[0].get_excitations_per_winding()[1].get_current()->get_processed()->get_peak().value() > 0);
    }
}