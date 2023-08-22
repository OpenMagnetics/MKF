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
        inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, OpenMagnetics::WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent, turnsRatios);
    }

    TEST(Test_All_Cores) {
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        // std::cout << coreAdviser.read_log() << std::endl;

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "T 18/9.0/7.1 - Kool Mu Hf 40 - Ungapped");
    }

    TEST(Test_All_Cores_Two_Chosen_Ones) {
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(true);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, 2);

        CHECK(masMagnetics.size() == 2);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "T 18/9.0/7.1 - Kool Mu Hf 40 - Ungapped");
        CHECK(masMagnetics[1].get_magnetic().get_core().get_name() == "EP 20 - 3C91 - Gapped 0.605 mm");
    }

    TEST(Test_No_Toroids_High_Power) {
        double voltagePeakToPeak = 6000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "E 70/33/32 - 95 - Distributed gapped 1.36 mm");
        CHECK(masMagnetics[0].get_magnetic().get_core().get_functional_description().get_number_stacks() == 3);
    }

    TEST(Test_No_Toroids_High_Power_High_Frequency) {
        double voltagePeakToPeak = 600000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 500000;
        double desiredMagnetizingInductance = 10e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 0;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "E 114/46/35 - XFlux 26 - Ungapped");
        CHECK(masMagnetics[0].get_magnetic().get_core().get_functional_description().get_number_stacks() == 2);
    }

    TEST(Test_No_Toroids_Low_Power) {
        double voltagePeakToPeak = 60;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "EFD 10/5/3 - 3C95 - Gapped 0.13999999999999999 mm");
        CHECK(masMagnetics[0].get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
    }

    TEST(Test_No_Toroids_Low_Power_Low_Losses) {
        double voltagePeakToPeak = 60;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0.1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "ER 48/18/18 - 3C94 - Gapped 1.0 mm");
        CHECK(masMagnetics[0].get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
    }

    TEST(Test_No_Toroids_Redo_Culling) {
        double voltagePeakToPeak = 6000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "E 42/21/15 - Kool Mu Hf 40 - Ungapped");
        CHECK(masMagnetics[0].get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
    }

    TEST(Test_No_Toroids_Two_Windings) {
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {0.1};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::WINDING_WINDOW_AREA] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_LOSSES] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::CORE_TEMPERATURE] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, 2);

        CHECK(masMagnetics.size() == 2);
        CHECK(masMagnetics[0].get_magnetic().get_core().get_name() == "PQ 26/20 - 3C95 - Gapped 0.365 mm");
        CHECK(masMagnetics[1].get_magnetic().get_core().get_name() == "PQ 26/20 - 3C94 - Gapped 0.361 mm");
    }

}