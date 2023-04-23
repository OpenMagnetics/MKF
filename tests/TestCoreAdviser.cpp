#include "CoreAdviser.h"
#include "InputsWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>


SUITE(CoreAdviser) {
    void prepare_test_parameters(double dcCurrent, double ambientTemperature, double frequency, std::vector<double> turnsRatios,
                                 double desiredMagnetizingInductance, OpenMagnetics::InputsWrapper& inputs,
                                 double peakToPeak = 20, double dutyCycle = 0.5) {
        inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, OpenMagnetics::WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent);
    }

    TEST(test) {
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 20000;
        double desiredMagnetizingInductance = 1e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        OpenMagnetics::OperationPoint operationPoint;
        OpenMagnetics::CoreAdviser coreAdviser;
        auto core = coreAdviser.get_advised_core(inputs);
    }
}