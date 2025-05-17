#include "physical_models/MagneticEnergy.h"
#include "TestingUtils.h"
#include "support/Utils.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <typeinfo>
#include <vector>

using namespace MAS;
using namespace OpenMagnetics;

using json = nlohmann::json;

SUITE(MagneticEnergy) {
    double max_error = 0.05;
    void prepare_test_parameters(double dcCurrent, double ambientTemperature, double frequency,
                                 double desiredMagnetizingInductance, std::vector<CoreGap> gapping,
                                 std::string coreShape, std::string coreMaterial, Core& core,
                                  OpenMagnetics::Inputs& inputs, double peakToPeak = 20, int numberStacks = 1) {
        double dutyCycle = 0.5;

        inputs = OpenMagnetics::Inputs::create_quick_operating_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent);

        core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);
    }

    TEST(Test_Magnetic_Energy_Iron_Powder_Core) {
        settings->reset();
        clear_databases();

        double ambientTemperature = 25;
        double frequency = 100000;
        std::string coreShape = "ETD 49";
        std::string coreMaterial = "XFlux 60";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        OpenMagnetics::Inputs inputs;
        Core core;

        prepare_test_parameters(0, ambientTemperature, frequency, -1, gapping, coreShape,
                                coreMaterial, core, inputs);
        auto operatingPoint = inputs.get_operating_point(0);

        MagneticEnergy magneticEnergy(
            std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));

        double expectedValue = 1.34;

        double totalMagneticEnergy = magneticEnergy.calculate_core_maximum_magnetic_energy(core, operatingPoint);
        CHECK_CLOSE(expectedValue, totalMagneticEnergy, max_error * expectedValue);
    }

    TEST(Test_Magnetic_Energy_Ferrite_Core) {
        settings->reset();
        clear_databases();
        double ambientTemperature = 25;
        double frequency = 100000;
        std::string coreShape = "ETD 49";
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        OpenMagnetics::Inputs inputs;
        Core core;

        prepare_test_parameters(0, ambientTemperature, frequency, -1, gapping, coreShape,
                                coreMaterial, core, inputs);
        auto operatingPoint = inputs.get_operating_point(0);

        MagneticEnergy magneticEnergy(
            std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));

        double expectedValue = 0.124;

        double totalMagneticEnergy = magneticEnergy.calculate_core_maximum_magnetic_energy(core, operatingPoint);
        CHECK_CLOSE(expectedValue, totalMagneticEnergy, max_error * expectedValue);
    }

    TEST(Test_Magnetic_Energy_Gap) {
        settings->reset();
        clear_databases();
        int numberStacks = 1;
        double magneticFluxDensitySaturation = 0.42;
        std::string coreShape = "ETD 49";
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);

        Core core = OpenMagneticsTesting::get_quick_core(coreShape, gapping, numberStacks, coreMaterial);

        MagneticEnergy magneticEnergy(
            std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));

        std::vector<double> expectedValues = {0.07, 0.045, 0.045};

        for (size_t i = 0; i < core.get_functional_description().get_gapping().size(); ++i)
        {
            double gapMagneticEnergy = magneticEnergy.get_gap_maximum_magnetic_energy(core.get_functional_description().get_gapping()[i], magneticFluxDensitySaturation);
            CHECK_CLOSE(expectedValues[i], gapMagneticEnergy, max_error * expectedValues[i]);
        }
    }

    TEST(Test_Magnetic_Energy_Input) {
        settings->reset();
        clear_databases();
        double voltagePeakToPeak = 1000;
        double desiredMagnetizingInductance = 0.0002;
        double ambientTemperature = 25;
        double frequency = 100000;
        std::string coreShape = "ETD 49";
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_spacer_gap(0.003);
        OpenMagnetics::Inputs inputs;
        Core core;

        prepare_test_parameters(0, ambientTemperature, frequency, desiredMagnetizingInductance, gapping, coreShape,
                                coreMaterial, core, inputs, voltagePeakToPeak);

        MagneticEnergy magneticEnergy(
            std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));

        double expectedValue = 0.0016;

        double requiredMagneticEnergy = magneticEnergy.calculate_required_magnetic_energy(inputs).get_nominal().value();
        CHECK_CLOSE(expectedValue, requiredMagneticEnergy, max_error * expectedValue);
    }

}