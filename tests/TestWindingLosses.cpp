#include "WindingLosses.h"
#include "Utils.h"
#include "CoreWrapper.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <cmath>

SUITE(WindingLosses) {
    double maximumError = 0.05;

    TEST(Test_Winding_Losses_One_Turn_Round_Tendency) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;

        OpenMagnetics::Processed processed;
        processed.set_label(label);
        processed.set_offset(offset);
        processed.set_peak_to_peak(peakToPeak);
        processed.set_duty_cycle(dutyCycle);
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                         magnetizingInductance,
                                                                                         temperature,
                                                                                         label,
                                                                                         peakToPeak,
                                                                                         dutyCycle,
                                                                                         offset);


        double bobbinHeight = 0.0035;
        double bobbinWidth = 0.0031;
        std::vector<double> bobbinCenterCoodinates({0.01, 0, 0});
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                         numberParallels,
                                                         bobbinHeight,
                                                         bobbinWidth,
                                                         bobbinCenterCoodinates,
                                                         interleavingLevel,
                                                         windingOrientation,
                                                         layersOrientation,
                                                         turnsAlignment,
                                                         sectionsAlignment);


        auto ohmicLosses100kHz = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(), ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0], ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0] * maximumError);
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses());
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);

        auto scaledOperationPoint = OpenMagnetics::InputsWrapper::scale_time_to_frequency(inputs.get_operating_point(0), frequency * 10);
        scaledOperationPoint = OpenMagnetics::InputsWrapper::process_operating_point(scaledOperationPoint, frequency * 10);
        auto ohmicLosses1MHz = OpenMagnetics::WindingLosses::calculate_losses(coil, scaledOperationPoint, temperature);
        CHECK_CLOSE(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses() * maximumError);
        CHECK(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1] > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);
        CHECK(ohmicLosses1MHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses());
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::ROUND);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.002031},
                                                                        {25000, 0.002054},
                                                                        {50000, 0.002120},
                                                                        {100000, 0.002356},
                                                                        {200000, 0.002988},
                                                                        {250000, 0.003294},
                                                                        {500000, 0.004467}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_With_DC) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::ROUND);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 4.2;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.037826952089878296},
                                                                        {25000, 0.037848903327517246},
                                                                        {50000, 0.037910312818917893},
                                                                        {100000, 0.038119608801605503},
                                                                        {200000, 0.038714247929474387},
                                                                        {250000, 0.039032903977111739},
                                                                        {500000, 0.040351988490429848}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);

            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);

            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::ROUND);
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

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.00205},
                                                                        {100000, 0.00234},
                                                                        {500000, 0.00459}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty_With_DC) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::ROUND);
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

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 4.2;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.03784},
                                                                        {100000, 0.03812},
                                                                        {500000, 0.04037354}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_90_Duty) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_diameter(0.00071);
        wire.set_nominal_value_outer_diameter(0.000762);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::ROUND);
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

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.9;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.00209},
                                                                        {100000, 0.00261},
                                                                        {500000, 0.00512}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.0000762);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000873);
        wire.set_number_conductors(60);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0033734},
                                                                        {25000, 0.0033741},
                                                                        {50000, 0.003376},
                                                                        {100000, 0.0033837},
                                                                        {200000, 0.0034145},
                                                                        {250000, 0.0034375},
                                                                        {500000, 0.0036299}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Strands) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00004);
        strandOuterDiameter.set_nominal(0.000049);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.001576);
        wire.set_number_conductors(600);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0011128},
                                                                        {25000, 0.001113},
                                                                        {50000, 0.0011137},
                                                                        {100000, 0.0011163},
                                                                        {200000, 0.001127},
                                                                        {250000, 0.0011349},
                                                                        {500000, 0.0012015}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Triangular_With_DC_Many_Strands) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00004);
        strandOuterDiameter.set_nominal(0.000049);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.001576);
        wire.set_number_conductors(600);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 10;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{500000, 0.11111}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Few_Strands) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.0001);
        strandOuterDiameter.set_nominal(0.00011);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000551);
        wire.set_number_conductors(12);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0084111},
                                                                        {25000, 0.0084123},
                                                                        {50000, 0.008416},
                                                                        {100000, 0.0084307},
                                                                        {200000, 0.0084895},
                                                                        {250000, 0.0084336},
                                                                        {500000, 0.0088009}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Many_Strands) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00002);
        strandOuterDiameter.set_nominal(0.000029);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.004384);
        wire.set_number_conductors(20000);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.00015792},
                                                                        {25000, 0.00015798},
                                                                        {50000, 0.00015816},
                                                                        {100000, 0.0001589},
                                                                        {200000, 0.00016186},
                                                                        {250000, 0.00016407},
                                                                        {500000, 0.00018254}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Rectangular_Sinusoidal) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "EQ 13";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.002);
        wire.set_nominal_value_conducting_height(0.0008);
        wire.set_nominal_value_outer_width(0.00223);
        wire.set_nominal_value_outer_height(0.00103);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.000238331},
                                                                        {25000, 0.000308836},
                                                                        {50000, 0.000379671},
                                                                        {100000, 0.000482108},
                                                                        {200000, 0.000639395},
                                                                        {250000, 0.000704667},
                                                                        {500000, 0.000969564}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Foil_Sinusoidal) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({1});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.02);
        wire.set_nominal_value_outer_width(0.0001);
        wire.set_nominal_value_outer_height(0.02);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::FOIL);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.000238331},
                                                                        {25000, 0.000308836},
                                                                        {50000, 0.000379671},
                                                                        {100000, 0.000482108},
                                                                        {200000, 0.000639395},
                                                                        {250000, 0.000704667},
                                                                        {500000, 0.000969564}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }

    TEST(Test_Winding_Losses_One_Turn_Thin_Rectangular_Sinusoidal) {

        double temperature = 20;
        std::vector<uint64_t> numberTurns({20});
        std::vector<uint64_t> numberParallels({1});
        std::string shapeName = "ETD 44/22/15";
        uint64_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0069);
        wire.set_nominal_value_conducting_height(0.0007);
        wire.set_nominal_value_outer_width(0.00709);
        wire.set_nominal_value_outer_height(0.00093);
        wire.set_number_conductors(1);
        wire.set_material("copper");
        wire.set_type(OpenMagnetics::WireType::RECTANGULAR);
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0052798},
                                                                        {25000, 0.057152},
                                                                        {50000, 0.10666},
                                                                        {100000, 0.1748},
                                                                        {200000, 0.24893},
                                                                        {250000, 0.27669},
                                                                        {500000, 0.39098}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset);


            auto ohmicLosses = OpenMagnetics::WindingLosses::calculate_losses(coil, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
    }
}