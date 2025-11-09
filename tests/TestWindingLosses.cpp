#include "support/Settings.h"
#include "support/Painter.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingLosses.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"
#include "processors/Inputs.h"
#include "processors/MagneticSimulator.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <cmath>

using namespace MAS;
using namespace OpenMagnetics;

auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");

SUITE(WindingLossesRound) {
    auto settings = Settings::GetInstance();
    double maximumError = 0.15;

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Stacked) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "E 42/21/20";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);
        int64_t numberStacks = 2;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true,
                                                         numberStacks);

        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0049}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Tendency) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});

        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "ETD 34/17/11";

        Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);

        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto ohmicLosses100kHz = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(), ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0], ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0] * maximumError);
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses());
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);

        // auto scaledOperatingPoint = OpenMagnetics::Inputs::scale_time_to_frequency(inputs.get_operating_point(0), frequency * 10);
        OperatingPoint scaledOperatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::Inputs::scale_time_to_frequency(scaledOperatingPoint, frequency * 10);
        scaledOperatingPoint = OpenMagnetics::Inputs::process_operating_point(scaledOperatingPoint, frequency * 10);
        auto ohmicLosses1MHz = WindingLosses().calculate_losses(magnetic, scaledOperatingPoint, temperature);
        CHECK_CLOSE(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses() * maximumError);
        CHECK(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1] > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);
        CHECK(ohmicLosses1MHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses());
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal) {
        // Test to evaluate skin effect losses, as no fringing or proximity are present
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Round_Sinusoidal.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.0022348},
            {10000, 0.002238},
            {20000, 0.0022476},
            {30000, 0.0022633},
            {40000, 0.0022546},
            {50000, 0.0023109},
            {60000, 0.0023419},
            {70000, 0.0023769},
            {80000, 0.0024153},
            {90000, 0.0024569},
            {100000, 0.0025011},
            {200000, 0.0030259},
            {300000, 0.0035737},
            {400000, 0.0040654},
            {500000, 0.0044916},
            {600000, 0.0048621},
            {700000, 0.0051882},
            {800000, 0.0054789},
            {900000, 0.0057414},
            {1000000, 0.0059805},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_One_Turn_Round_Sinusoidal_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal) {
        // Test to evaluate proximity effect losses, as there is no fringing and the wire is small enough to avoid skin
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.17371},
            {10000, 0.17372},
            {20000, 0.17373},
            {30000, 0.17374},
            {40000, 0.17375},
            {50000, 0.17378},
            {60000, 0.1738},
            {70000, 0.17384},
            {80000, 0.17387},
            {90000, 0.17391},
            {100000, 0.17396},
            {200000, 0.1747},
            {300000, 0.17593},
            {400000, 0.17764},
            {500000, 0.17983},
            {600000, 0.18248},
            {700000, 0.1856},
            {800000, 0.18916},
            {900000, 0.19315},
            {1000000, 0.19755},
            {3000000, 0.34496},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 167.89},
            {10000, 169.24},
            {20000, 174.77},
            {30000, 183.33},
            {40000, 194.12},
            {50000, 206.33},
            {60000, 219.3},
            {70000, 232.5},
            {80000, 245.61},
            {90000, 258.43},
            {100000, 270.86},
            {200000, 376.14},
            {300000, 460.7},
            {400000, 532.27},
            {500000, 594.6},
            {600000, 649.64},
            {700000, 699.9},
            {800000, 746.3},
            {900000, 789.66},
            {1000000, 830.49},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing_Far) {
        // Worst error in this one
        double maximumError = 0.4;
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing_Far.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 204.23},
            {10000, 204.61},
            {20000, 205.73},
            {30000, 207.52},
            {40000, 209.9},
            {50000, 212.74},
            {60000, 215.94},
            {70000, 219.41},
            {80000, 223.07},
            {90000, 226.85},
            {100000, 230.71},
            {200000, 269.05},
            {300000, 303.53},
            {400000, 333.71},
            {500000, 360.06},
            {600000, 383.12},
            {700000, 403.36},
            {800000, 421.2},
            {900000, 436.95},
            {1000000, 450.91},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_One_Turn_Round_Sinusoidal_Fringing_Far_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Eight_Turns_Round_Sinusoidal_Rectangular_Column) {
        // Test to evaluate proximity effect losses, as there is no fringing and the wire is small enough to avoid skin
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Eight_Turns_Round_Sinusoidal_Rectangular_Column.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.1194420289},
            {10000, 0.1194431089},
            {20000, 0.1194463487},
            {30000, 0.1194517488},
            {40000, 0.1194593093},
            {50000, 0.1194690305},
            {60000, 0.1194809123},
            {70000, 0.1194949548},
            {80000, 0.1195111578},
            {90000, 0.119529521},
            {100000, 0.1195500418},
            {200000, 0.1198743283},
            {300000, 0.120416096},
            {400000, 0.1211778912},
            {500000, 0.1221636039},
            {600000, 0.1233738382},
            {700000, 0.1248007533},
            {800000, 0.1264334757},
            {900000, 0.1282616081},
            {1000000, 0.1302673627},
        });



        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Eight_Turns_Round_Sinusoidal_Rectangular_Column_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_With_DC) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 4.2;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.03782},
                                                                        {25000, 0.03784},
                                                                        {50000, 0.03791},
                                                                        {100000, 0.03811},
                                                                        {200000, 0.03871},
                                                                        {250000, 0.03903},
                                                                        {500000, 0.04035}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);

            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);

            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_Interleaving) {
        // Test to evaluate proximity effect losses, as there is no fringing and the wire is small enough to avoid skin
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_Interleaving.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.14837},
            {10000, 0.14837},
            {20000, 0.14837},
            {30000, 0.14837},
            {40000, 0.14838},
            {50000, 0.14838},
            {60000, 0.14839},
            {70000, 0.14839},
            {80000, 0.1484},
            {90000, 0.14841},
            {100000, 0.14842},
            {200000, 0.14856},
            {300000, 0.14881},
            {400000, 0.14941},
            {500000, 0.14957},
            {600000, 0.1501},
            {700000, 0.15071},
            {800000, 0.15141},
            {900000, 0.1522},
            {1000000, 0.15307},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_Interleaving_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving) {
        // Test to evaluate proximity effect losses, as there is no fringing and the wire is small enough to avoid skin
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.13843},
            {10000, 0.13843},
            {20000, 0.13843},
            {30000, 0.13844},
            {40000, 0.13844},
            {50000, 0.13846},
            {60000, 0.13846},
            {70000, 0.13847},
            {80000, 0.13848},
            {90000, 0.13849},
            {100000, 0.1385},
            {200000, 0.13873},
            {300000, 0.1391},
            {400000, 0.13963},
            {500000, 0.14029},
            {600000, 0.1411},
            {700000, 0.14206},
            {800000, 0.14314},
            {900000, 0.14437},
            {1000000, 0.14572},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving_2) {
        // Test to evaluate proximity effect losses, as there is no fringing and the wire is small enough to avoid skin
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving_2.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();
        MagnetizingInductance magnetizingInductanceModel("ZHANG");

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.48177},
            {10000, 0.48177},
            {20000, 0.48178},
            {30000, 0.48179},
            {40000, 0.48181},
            {50000, 0.48183},
            {60000, 0.48186},
            {70000, 0.4819},
            {80000, 0.48194},
            {90000, 0.48198},
            {100000, 0.48203},
            {200000, 0.48284},
            {300000, 0.48418},
            {400000, 0.48605},
            {500000, 0.48845},
            {600000, 0.49138},
            {700000, 0.49483},
            {800000, 0.49879},
            {900000, 0.50325},
            {1000000, 0.50821},
            {3000000, 0.69729},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            OpenMagnetics::Mas auxMas;
            auxMas.set_magnetic(magnetic);
            auxMas.set_inputs(inputs);
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Twelve_Turns_Round_Sinusoidal_No_Interleaving_2_" + std::to_string(frequency) +".json");
            OpenMagnetics::to_file(outFile, auxMas);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty_With_DC) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::TRIANGULAR;
        double offset = 4.2;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.03783},
                                                                        {100000, 0.03811},
                                                                        {500000, 0.040374}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
}

SUITE(WindingLossesLitz) {
    auto settings = Settings::GetInstance();
    double maximumError = 0.15;
    bool plot = false;

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.0000762);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000873);
        wire.set_number_conductors(60);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.003374},
                                                                        {25000, 0.003371},
                                                                        {50000, 0.00336},
                                                                        {100000, 0.003387},
                                                                        {200000, 0.003415},
                                                                        {250000, 0.003435},
                                                                        {500000, 0.003629}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Strands) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00004);
        strandOuterDiameter.set_nominal(0.000049);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.001576);
        wire.set_number_conductors(600);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.001118},
                                                                        {25000, 0.00113},
                                                                        {50000, 0.001117},
                                                                        {100000, 0.001113},
                                                                        {200000, 0.00117},
                                                                        {250000, 0.001139},
                                                                        {500000, 0.001205}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Triangular_With_DC_Many_Strands) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00004);
        strandOuterDiameter.set_nominal(0.000049);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.001576);
        wire.set_number_conductors(600);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::TRIANGULAR;
        double offset = 10;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{500000, 0.113}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Few_Strands) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.0001);
        strandOuterDiameter.set_nominal(0.00011);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000551);
        wire.set_number_conductors(12);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.008411},
                                                                        {25000, 0.008412},
                                                                        {50000, 0.008416},
                                                                        {100000, 0.008430},
                                                                        {200000, 0.008489},
                                                                        {250000, 0.008433},
                                                                        {500000, 0.008800}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Many_Strands) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00002);
        strandOuterDiameter.set_nominal(0.000029);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.004384);
        wire.set_number_conductors(20000);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0001572},
                                                                        {25000, 0.0001578},
                                                                        {50000, 0.0001586},
                                                                        {100000, 0.000159},
                                                                        {200000, 0.0001516},
                                                                        {250000, 0.0001547},
                                                                        {500000, 0.00016}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Litz_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.0000762);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000873);
        wire.set_number_conductors(60);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.03375},
                                                                        {25000, 0.03371},
                                                                        {50000, 0.03388},
                                                                        {100000, 0.03399},
                                                                        {200000, 0.03479},
                                                                        {250000, 0.0353},
                                                                        {500000, 0.04004}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            if (plot) {
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_Ten_Turns_Litz_Sinusoidal_" + std::to_string(expectedValue) +".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile, true);
                settings->set_painter_mode(PainterModes::CONTOUR);
                settings->set_painter_logarithmic_scale(false);
                settings->set_painter_include_fringing(false);
                settings->set_painter_number_points_x(200);
                settings->set_painter_number_points_y(200);
                settings->set_painter_maximum_value_colorbar(std::nullopt);
                settings->set_painter_minimum_value_colorbar(std::nullopt);
                painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
                painter.paint_core(magnetic);
                painter.paint_bobbin(magnetic);
                painter.paint_coil_turns(magnetic);
                painter.export_svg();
            }
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Thirty_Turns_Litz_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({30});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "P 26/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        WireRound strand;
        OpenMagnetics::Wire wire;
        DimensionWithTolerance strandConductingDiameter;
        DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.0000762);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000873);
        wire.set_number_conductors(60);
        wire.set_type(WireType::LITZ);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.1133},
                                                                        {25000, 0.1137},
                                                                        {50000, 0.1146},
                                                                        {100000, 0.1186},
                                                                        {200000, 0.1343},
                                                                        {250000, 0.1460},
                                                                        {500000, 0.2442}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
}

SUITE(WindingLossesRectangular) {
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    auto settings = Settings::GetInstance();
    double maximumError = 0.2;

    TEST(Test_Winding_Losses_One_Turn_Rectangular_Sinusoidal_No_Fringing) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Rectangular_Sinusoidal_No_Fringing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 0.0004333},
            {10000, 0.00045385},
            {20000, 0.00048219},
            {30000, 0.00050866},
            {40000, 0.00053534},
            {50000, 0.00056317},
            {60000, 0.0005923},
            {70000, 0.00062263},
            {80000, 0.00065399},
            {90000, 0.0006862},
            {100000, 0.00071907},
            {200000, 0.0010607},
            {300000, 0.0013908},
            {400000, 0.0016968},
            {500000, 0.0019801},
            {600000, 0.0022432},
            {700000, 0.0024883},
            {800000, 0.002717},
            {900000, 0.0029309},
            {1000000, 0.0031313},
            {3000000, 0.005539},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            settings->set_magnetic_field_include_fringing(false);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Ungapped_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.00011974 + 0.0003736},
                                                                        {100000, 0.00087387 + 0.0003736},
                                                                        {200000, 0.0012227 + 0.0003736},
                                                                        {300000, 0.001487 + 0.0003736},
                                                                        {400000, 0.0017156 + 0.0003736},
                                                                        {500000, 0.001927 + 0.0003736},
                                                                        {600000, 0.0021237 + 0.0003736},
                                                                        {700000, 0.0023129 + 0.0003736},
                                                                        {800000, 0.0024954 + 0.0003736},
                                                                        {900000, 0.0026719 + 0.0003736},
                                                                        {1000000, 0.002843 + 0.0003736}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Ungapped_Sinusoidal_7_Amps) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 14;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0058551 * 5},
                                                                        {100000, 0.03363 + 0.0058551 * 5},
                                                                        {200000, 0.06063 + 0.0058551 * 5},
                                                                        {300000, 0.073721 + 0.0058551 * 5},
                                                                        {400000, 0.084989 + 0.0058551 * 5},
                                                                        {500000, 0.095376 + 0.0058551 * 5},
                                                                        {600000, 0.10525 + 0.0058551 * 5},
                                                                        {700000, 0.11477 + 0.0058551 * 5},
                                                                        {800000, 0.12399 + 0.0058551 * 5},
                                                                        {900000, 0.13296 + 0.0058551 * 5},
                                                                        {1000000, 0.141661 + 0.0058551 * 5}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Gapped_Sinusoidal_7_Amps) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 14;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0058551 * 5},
                                                                        {100000, 0.03363 + 0.0058551 * 5},
                                                                        {200000, 0.06063 + 0.0058551 * 5},
                                                                        {300000, 0.073721 + 0.0058551 * 5},
                                                                        {400000, 0.084989 + 0.0058551 * 5},
                                                                        {500000, 0.095376 + 0.0058551 * 5},
                                                                        {600000, 0.10525 + 0.0058551 * 5},
                                                                        {700000, 0.11477 + 0.0058551 * 5},
                                                                        {800000, 0.12399 + 0.0058551 * 5},
                                                                        {900000, 0.13296 + 0.0058551 * 5},
                                                                        {1000000, 0.141661 + 0.0058551 * 5}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Seven_Turns_Rectangular_Ungapped_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({7});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 27/17";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_height(0.0007676);
        wire.set_nominal_value_outer_width(0.003838);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);
        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);

        auto bobbin = coil.resolve_bobbin();
        auto bobbinProcessedDescription = bobbin.get_processed_description().value();
        double bobbinColumnThickness = 0.001;
        bobbinProcessedDescription.set_wall_thickness(bobbinColumnThickness);
        bobbinProcessedDescription.set_column_thickness(bobbinColumnThickness);

        WindingWindowElement windingWindowElement;

        auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
        auto coreCentralColumn = core.get_processed_description()->get_columns()[0];
        windingWindowElement.set_width(coreWindingWindow.get_width().value() - bobbinColumnThickness);
        windingWindowElement.set_height(coreWindingWindow.get_height().value() - bobbinColumnThickness * 2);
        windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + bobbinColumnThickness + windingWindowElement.get_width().value() / 2, 0, 0}));
        bobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + bobbinColumnThickness);
        bobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + bobbinColumnThickness);
        bobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));

        bobbin.set_processed_description(bobbinProcessedDescription);
        coil.set_bobbin(bobbin);
        coil.wind();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.00019313 + 0.00019313 * 6},
                                                                        {100000, 0.0025456},
                                                                        {200000, 0.00342804},
                                                                        {300000, 0.00419621},
                                                                        {400000, 0.00484483},
                                                                        {500000, 0.0054165},
                                                                        {600000, 0.0059334},
                                                                        {700000, 0.00640876},
                                                                        {800000, 0.00685123},
                                                                        {900000, 0.00726681},
                                                                        {1000000, 0.00765988}
                                                                        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
            // auto outFile = outputFilePath;
            // outFile.append("Test_Winding_Losses_Seven_Turns_Rectangular_Ungapped_Sinusoidal_" + std::to_string(expectedValue) +".svg");
            // std::filesystem::remove(outFile);
            // Painter painter(outFile, PainterModes::QUIVER);
            // painter.set_number_points_x(100);
            // painter.set_number_points_y(100);
            // painter.set_logarithmic_scale(false);
            // painter.set_fringing_effect(true);
            // painter.set_maximum_scale_value(std::nullopt);
            // painter.set_minimum_scale_value(std::nullopt);
            // painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            // painter.paint_core(magnetic);
            // painter.paint_bobbin(magnetic);
            // painter.paint_coil_turns(magnetic);
            // painter.export_svg();

        }
        settings->reset();
    }
}

SUITE(WindingLossesFoil) {
    auto settings = Settings::GetInstance();
    double maximumError = 0.3;
    TEST(Test_Winding_Losses_One_Turn_Foil_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.02);
        wire.set_nominal_value_outer_width(0.00011);
        wire.set_nominal_value_outer_height(0.02);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.00037759},
                                                                        {25000, 0.00039609},
                                                                        {50000, 0.0004512},
                                                                        {100000, 0.00058852},
                                                                        {200000, 0.00061066},
                                                                        {250000, 0.00061684},
                                                                        {500000, 0.00064413}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Foil_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.02);
        wire.set_nominal_value_outer_width(0.00018);
        wire.set_nominal_value_outer_height(0.02);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0041536},
                                                                        {25000, 0.0039304},
                                                                        {50000, 0.00440837},
                                                                        {100000, 0.00498607},
                                                                        {200000, 0.00565498},
                                                                        {250000, 0.00589129},
                                                                        {500000, 0.00670688}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Short_Turns_Foil_Sinusoidal) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.007);
        wire.set_nominal_value_outer_width(0.00015);
        wire.set_nominal_value_outer_height(0.007);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::FOIL);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         false);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0090968},
                                                                        {25000, 0.0096968},
                                                                        {50000, 0.0104159},
                                                                        {100000, 0.0115721},
                                                                        {200000, 0.0130536},
                                                                        {250000, 0.013589}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);

        }
        settings->reset();
    }
}

SUITE(WindingLossesToroidalCores) {
    auto settings = Settings::GetInstance();
    double maximumError = 0.15;
    bool plot = false;

    TEST(Test_Winding_Losses_Toroidal_Core_One_Turn_Round_Tendency) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::vector<double> turnsRatios;

        auto label = WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

        Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset,
                                                                                         turnsRatios);

        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::SPREAD;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        auto ohmicLosses100kHz = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);


        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(), ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0], ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0] * maximumError);
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses());
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);

        OperatingPoint scaledOperatingPoint = inputs.get_operating_point(0);
        OpenMagnetics::Inputs::scale_time_to_frequency(scaledOperatingPoint, frequency * 10);
        scaledOperatingPoint = OpenMagnetics::Inputs::process_operating_point(scaledOperatingPoint, frequency * 10);
        auto ohmicLosses1MHz =    WindingLosses().calculate_losses(magnetic, scaledOperatingPoint, temperature);
        CHECK_CLOSE(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses() * maximumError);
        CHECK(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1] > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);
        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_proximity_effect_losses().value().get_losses_per_harmonic()[1],
            ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_proximity_effect_losses().value().get_losses_per_harmonic()[1],
            ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_proximity_effect_losses().value().get_losses_per_harmonic()[1] * maximumError);
        CHECK(ohmicLosses1MHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses());
        settings->reset();

        if (plot) {
            auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Toroidal_Core_One_Turn_Round_Tendency.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings->set_painter_mode(PainterModes::CONTOUR);
            settings->set_painter_mode(PainterModes::QUIVER);
            settings->set_painter_logarithmic_scale(false);
            settings->set_painter_include_fringing(true);
            settings->set_painter_number_points_x(50);
            settings->set_painter_number_points_y(50);
            settings->set_painter_maximum_value_colorbar(std::nullopt);
            settings->set_painter_minimum_value_colorbar(std::nullopt);
            painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
            painter.paint_core(magnetic);
            // painter.paint_coil_sections(magnetic);
            painter.paint_coil_turns(magnetic);
            painter.export_svg();
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Toroidal_Core) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.002032},
                                                                      {25000, 0.002053},
                                                                      {50000, 0.002121},
                                                                      {100000, 0.002355},
                                                                      {200000, 0.002987},
                                                                      {250000, 0.003293},
                                                                      {500000, 0.004466}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Toroidal_Core) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::ROUND);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 100e-6;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0213484},
                                                                      {25000, 0.0217045},
                                                                      {50000, 0.0227473},
                                                                      {100000, 0.0265887},
                                                                      {200000, 0.0409422},
                                                                      {250000, 0.0460606},
                                                                      {500000, 0.0658813}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto windingLosses = WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0058);
        wire.set_nominal_value_conducting_height(0.002);
        wire.set_nominal_value_outer_width(0.0059);
        wire.set_nominal_value_outer_height(0.002076);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 9.94277e-05},
                                                                      {25000, 0.000216944},
                                                                      {50000, 0.00030592},
                                                                      {100000, 0.000432524},
                                                                      {200000, 0.000611667},
                                                                      {250000, 0.000683864},
                                                                      {500000, 0.000967129}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);

            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_One_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire" + std::to_string(expectedValue) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile);
                OpenMagnetics::Magnetic magnetic;
                magnetic.set_core(core);
                magnetic.set_coil(coil);
                // settings->set_painter_mode(PainterModes::CONTOUR);
                settings->set_painter_mode(PainterModes::QUIVER);
                settings->set_painter_logarithmic_scale(false);
                settings->set_painter_include_fringing(true);
                settings->set_painter_number_points_x(50);
                settings->set_painter_number_points_y(50);
                settings->set_painter_maximum_value_colorbar(std::nullopt);
                settings->set_painter_minimum_value_colorbar(std::nullopt);
                painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
                painter.paint_core(magnetic);
                // painter.paint_coil_sections(magnetic);
                painter.paint_coil_turns(magnetic);
                painter.export_svg();
            }

            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire) {
        settings->reset();
        clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::CENTERED;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        OpenMagnetics::Wire wire;
        wire.set_nominal_value_conducting_width(0.0058);
        wire.set_nominal_value_conducting_height(0.002);
        wire.set_nominal_value_outer_width(0.0059);
        wire.set_nominal_value_outer_height(0.002076);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(WireType::RECTANGULAR);
        wires.push_back(wire);

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires);

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 9.94277e-04},
                                                                      {25000, 0.00216944},
                                                                      {50000, 0.0030592},
                                                                      {100000, 0.00432524},
                                                                      {200000, 0.00611667},
                                                                      {250000, 0.00683864},
                                                                      {500000, 0.00967129}});

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {

            Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::Inputs::create_quick_operating_point_only_current(expectedValue,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);

            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_Ten_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire" + std::to_string(expectedValue) + ".svg");
                std::filesystem::remove(outFile);
                Painter painter(outFile);
                OpenMagnetics::Magnetic magnetic;
                magnetic.set_core(core);
                magnetic.set_coil(coil);
                // settings->set_painter_mode(PainterModes::CONTOUR);
                settings->set_painter_mode(PainterModes::QUIVER);
                settings->set_painter_logarithmic_scale(false);
                settings->set_painter_include_fringing(true);
                settings->set_painter_number_points_x(50);
                settings->set_painter_number_points_y(50);
                settings->set_painter_maximum_value_colorbar(std::nullopt);
                settings->set_painter_minimum_value_colorbar(std::nullopt);
                painter.paint_magnetic_field(inputs.get_operating_point(0), magnetic);
                painter.paint_core(magnetic);
                // painter.paint_coil_sections(magnetic);
                painter.paint_coil_turns(magnetic);
                painter.export_svg();
            }

            auto ohmicLosses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

}

SUITE(WindingLossesPlanar) {
    auto settings = Settings::GetInstance();
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double maximumError = 0.3;

    TEST(Test_Winding_Losses_One_Turn_Planar_Sinusoidal_No_Fringing) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Planar_Sinusoidal_No_Fringing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 87.383},
            {10000, 87.385},
            {20000, 87.389},
            {30000, 87.395},
            {40000, 87.403},
            {50000, 87.413},
            {60000, 87.423},
            {70000, 87.435},
            {80000, 87.446},
            {90000, 87.458},
            {100000, 87.470},
            {200000, 87.577},
            {300000, 87.660},
            {400000, 87.723},
            {500000, 87.771},
            {600000, 87.808},
            {700000, 87.838},
            {800000, 87.862},
            {900000, 87.882},
            {1000000, 87.898},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            settings->set_magnetic_field_include_fringing(false);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Planar_Sinusoidal_Fringing) {
        // Not sure about that many losses due to fringing losses in a small piece
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_One_Turn_Planar_Sinusoidal_Fringing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 87.648},
            // {10000, 275.05},
            // {20000, 356.31},
            // {30000, 410.13},
            // {40000, 4542.13},
            // {50000, 487.38},
            // {60000, 518.12},
            // {70000, 545.53},
            // {80000, 570.38},
            // {90000, 593.2},
            // {100000, 614.37},
            // {200000, 778.35},
            // {300000, 904.07},
            // {400000, 1011.6},
            // {500000, 1106.8},
            // {600000, 1192.3},
            // {700000, 1269.8},
            // {800000, 1340.6},
            // {900000, 1405.9},
            // {1000000, 1466.4},
            // {1500000, 1708.1},
        });

        settings->set_magnetic_field_include_fringing(true);
        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 5.8488},
            {10000, 13.251},
            {20000, 15.197},
            {30000, 16.110},
            {40000, 16.717},
            {50000, 17.2},
            {60000, 17.619},
            {70000, 18.000},
            {80000, 18.354},
            {90000, 18.686},
            {100000, 19.002},
            {200000, 21.636},
            {300000, 23.821},
            // {400000, 25.829},
            // {500000, 27.704},
            // {600000, 29.416},
            // {700000, 30.925},
            // {800000, 32.231},
            // {900000, 33.353},
            // {1000000, 34.308},
        });

        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            // settings->set_magnetic_field_include_fringing(false);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Close) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Close.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 5.53},
            {10000, 117.63},
            {20000, 167.38},
            {30000, 200.47},
            {40000, 224.01},
            {50000, 241.59},
            {60000, 255.19},
            {70000, 266.0},
            {80000, 274.77},
            {90000, 282.02},
            {100000, 288.1},
            {200000, 318.41},
            {300000, 329.73},
            {400000, 336.19},
            {500000, 340.86},
            {600000, 344.67},
            {700000, 347.9},
            {800000, 350.69},
            {900000, 353.12},
            {1000000, 355.25},
        });

        settings->set_magnetic_field_include_fringing(true);
        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Far) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Far.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 5.8408},
            {10000, 78.113},
            {20000, 105.33},
            {30000, 122.53},
            {40000, 134.58},
            {50000, 143.52},
            {60000, 150.44},
            {70000, 155.94},
            {80000, 160.43},
            {90000, 164.17},
            {100000, 167.33},
            {200000, 183.69},
            {300000, 189.99},
            {400000, 193.37},
            {500000, 195.62},
            {600000, 197.29},
            {700000, 198.61},
            {800000, 199.68},
            {900000, 200.56},
            {1000000, 201.29},
        });

        settings->set_magnetic_field_include_fringing(true);
        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing_Interleaving) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing_Interleaving.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        settings->reset();
        clear_databases();

        double temperature = 22;
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        std::vector<std::pair<double, double>> expectedWindingLosses({
            {1, 38.429},
            {10000, 40.235},
            {20000, 40.602},
            {30000, 40.841},
            {40000, 41.032},
            {50000, 41.199},
            {60000, 41.352},
            {70000, 41.494},
            {80000, 41.629},
            {90000, 41.756},
            {100000, 41.878},
            {200000, 42.915},
            {300000, 43.733},
            {400000, 44.417},
            {500000, 45.019},
            {600000, 45.555},
            {700000, 46.007},
            {800000, 46.374},
            {900000, 46.665},
            {1000000, 46.876},
        });

        settings->set_magnetic_field_include_fringing(true);
        for (auto& [frequency, expectedValue] : expectedWindingLosses) {
            OperatingPoint operatingPoint = inputs.get_operating_point(0);
            OpenMagnetics::Inputs::scale_time_to_frequency(operatingPoint, frequency, true);
            double magnetizingInductance = OpenMagnetics::resolve_dimensional_values(magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil(), &operatingPoint).get_magnetizing_inductance());
            operatingPoint = OpenMagnetics::Inputs::process_operating_point(operatingPoint, magnetizingInductance);

            WindingLosses windingLosses;
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, operatingPoint, temperature);
            CHECK_CLOSE(expectedValue, ohmicLosses.get_winding_losses(), expectedValue * maximumError);
        }
        settings->reset();
    }
}


SUITE(WindingLossesResistanceMatrix) {
    TEST(Test_Resistance_Matrix) {
        std::vector<int64_t> numberTurns = {80, 8, 6};
        std::vector<int64_t> numberParallels = {1, 2, 6};
        std::vector<double> turnsRatios = {16, 13};
        std::string shapeName = "ER 28";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = WindingOrientation::OVERLAPPING;
        auto layersOrientation = WindingOrientation::OVERLAPPING;
        auto turnsAlignment = CoilAlignment::SPREAD;
        auto sectionsAlignment = CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::Wire> wires;
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::Wire wire = find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         shapeName,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment,
                                                         wires,
                                                         true);

        coil.wind({0, 1, 2}, 1);

        double temperature = 20;
        double frequency = 123456;
        int64_t numberStacks = 1;
        std::string coreMaterial = "3C95";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0000008);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto resistanceMatrixAtFrequency = WindingLosses().calculate_resistance_matrix(magnetic, temperature, frequency);

        CHECK(resistanceMatrixAtFrequency.get_magnitude().size() == magnetic.get_coil().get_functional_description().size());
        for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
            CHECK(resistanceMatrixAtFrequency.get_magnitude()[windingIndex].size() == magnetic.get_coil().get_functional_description().size());
            CHECK(resolve_dimensional_values(resistanceMatrixAtFrequency.get_magnitude()[windingIndex][windingIndex]) > 0);
        }

    }
}

SUITE(WindingLossesWeb) {
    TEST(Test_Winding_Losses_Web_0) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/negative_losses.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        auto losses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), 25);

        json mierda;
        to_json(mierda, losses);
        // std::cout << mierda << std::endl;

        CHECK(losses.get_dc_resistance_per_winding().value()[0] > 0);
    }

    TEST(Test_Winding_Losses_Web_1) {
        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/slow_simulation.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        auto losses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), 25);

        json mierda;
        to_json(mierda, losses);
        // std::cout << mierda << std::endl;

        CHECK(losses.get_dc_resistance_per_winding().value()[0] > 0);
    }

    TEST(Test_Winding_Losses_Web_2) {
        settings->set_magnetic_field_include_fringing(false);

        std::string file_path = __FILE__;
        auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/huge_losses.json");
        auto mas = OpenMagneticsTesting::mas_loader(path);

        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        auto losses = WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), 25);

        CHECK(losses.get_winding_losses() < 2);
        settings->reset();
    }

    TEST(Test_Winding_Losses_Web_3) {
        settings->set_magnetic_field_include_fringing(false);
        settings->set_magnetic_field_mirroring_dimension(3);

        std::string file_path = __FILE__;
        OpenMagnetics::Mas mas1;
        OpenMagnetics::Mas mas2;
        OpenMagnetics::Mas mas3;
        {
            auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/planar_proximity_losses_1.json");
            mas1 = OpenMagneticsTesting::mas_loader(path);
        }
        {
            auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/planar_proximity_losses_2.json");
            mas2 = OpenMagneticsTesting::mas_loader(path);
        }
        {
            auto path = file_path.substr(0, file_path.rfind("/")).append("/testData/planar_proximity_losses_3.json");
            mas3 = OpenMagneticsTesting::mas_loader(path);
        } 

        auto magnetic1 = mas1.get_magnetic();
        auto inputs1 = mas1.get_inputs();

        auto losses1 = WindingLosses().calculate_losses(magnetic1, inputs1.get_operating_point(0), 25);
        // auto lossesPerTurns1 = losses1.get_winding_losses_per_turn().value();
        // for (auto losses : lossesPerTurns1) {
        //     std::cout << "losses.get_proximity_effect_losses(): " << losses.get_proximity_effect_losses()->get_losses_per_harmonic()[1] << std::endl;
        // }

        auto magnetic2 = mas2.get_magnetic();
        auto inputs2 = mas2.get_inputs();

        auto losses2 = WindingLosses().calculate_losses(magnetic2, inputs2.get_operating_point(0), 25);
        // auto lossesPerTurns2 = losses2.get_winding_losses_per_turn().value();
        // for (auto losses : lossesPerTurns2) {
        //     std::cout << "losses.get_proximity_effect_losses(): " << losses.get_proximity_effect_losses()->get_losses_per_harmonic()[1] << std::endl;
        // }

        auto magnetic3 = mas3.get_magnetic();
        auto inputs3 = mas3.get_inputs();

        auto losses3 = WindingLosses().calculate_losses(magnetic3, inputs3.get_operating_point(0), 25);
        // auto lossesPerTurns3 = losses3.get_winding_losses_per_turn().value();
        // for (auto losses : lossesPerTurns3) {
        //     std::cout << "losses.get_proximity_effect_losses(): " << losses.get_proximity_effect_losses()->get_losses_per_harmonic()[1] << std::endl;
        // }

        settings->set_painter_include_fringing(false);
        // CHECK(losses.get_winding_losses() <  2);
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Web_3_1.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            painter.paint_magnetic_field(inputs1.get_operating_point(0), magnetic1);
            painter.paint_core(magnetic1);
            painter.paint_bobbin(magnetic1);
            painter.paint_coil_turns(magnetic1);
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Web_3_2.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            painter.paint_magnetic_field(inputs2.get_operating_point(0), magnetic2);
            painter.paint_core(magnetic2);
            painter.paint_bobbin(magnetic2);
            painter.paint_coil_turns(magnetic2);
            painter.export_svg();
        }
        {
            auto outFile = outputFilePath;
            outFile.append("Test_Winding_Losses_Web_3_3.svg");
            std::filesystem::remove(outFile);
            Painter painter(outFile, true);
            painter.paint_magnetic_field(inputs3.get_operating_point(0), magnetic3);
            painter.paint_core(magnetic3);
            painter.paint_bobbin(magnetic3);
            painter.paint_coil_turns(magnetic3);
            painter.export_svg();
        }
        settings->reset();
    }
}