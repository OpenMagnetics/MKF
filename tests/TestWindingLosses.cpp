#include "Settings.h"
#include "Painter.h"
#include "WindingLosses.h"
#include "Utils.h"
#include "CoreWrapper.h"
#include "InputsWrapper.h"
#include "MagneticSimulator.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>
#include <cmath>


auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");

SUITE(WindingLossesRound) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.15;

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Stacked) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "E 42/21/20";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0049}});

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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }


    TEST(Test_Winding_Losses_One_Turn_Round_Tendency) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "ETD 34/17/11";

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

        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

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

        auto ohmicLosses100kHz = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(), ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0], ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0] * maximumError);
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses());
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);

        // auto scaledOperationPoint = OpenMagnetics::InputsWrapper::scale_time_to_frequency(inputs.get_operating_point(0), frequency * 10);
        OpenMagnetics::OperatingPoint scaledOperationPoint = inputs.get_operating_point(0);
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(scaledOperationPoint, frequency * 10);
        scaledOperationPoint = OpenMagnetics::InputsWrapper::process_operating_point(scaledOperationPoint, frequency * 10);
        auto ohmicLosses1MHz = OpenMagnetics::WindingLosses().calculate_losses(magnetic, scaledOperationPoint, temperature);
        CHECK_CLOSE(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(),
                    ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses() * maximumError);
        CHECK(ohmicLosses1MHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1] > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);
        CHECK(ohmicLosses1MHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses());
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.05e-3);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 100e-6;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.02044},
                                                                      {25000, 0.02171},
                                                                      {50000, 0.02512},
                                                                      {100000, 0.03373},
                                                                      {200000, 0.098},
                                                                      {250000, 0.115},
                                                                      {500000, 0.182}});

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


            auto windingLosses = OpenMagnetics::WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Interleaving) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({20, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.01e-3);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.095454},
                                                                      {25000, 0.088692},
                                                                      {50000, 0.099},
                                                                      {100000, 0.132},
                                                                      {200000, 0.26954},
                                                                      {250000, 0.30431},
                                                                      {500000, 0.51}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            std::vector<double> turnsRatios = {1};
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset,
                                                                                                  turnsRatios);


            auto windingLosses = OpenMagnetics::WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);

        }
        settings->reset();
    }


    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Interleaving_Rectangular_Column) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({20, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "E 34/14/9";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.01e-3);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.095454},
                                                                      {25000, 0.103},
                                                                      {50000, 0.113},
                                                                      {100000, 0.153},
                                                                      {200000, 0.344},
                                                                      {250000, 0.399},
                                                                      {500000, 0.61}});

        for (auto& testPoint : expectedWindingLosses) {

            OpenMagnetics::Processed processed;
            processed.set_label(label);
            processed.set_offset(offset);
            processed.set_peak_to_peak(peakToPeak);
            processed.set_duty_cycle(dutyCycle);
            std::vector<double> turnsRatios = {1};
            auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(testPoint.first,
                                                                                                  magnetizingInductance,
                                                                                                  temperature,
                                                                                                  label,
                                                                                                  peakToPeak,
                                                                                                  dutyCycle,
                                                                                                  offset,
                                                                                                  turnsRatios);


            auto windingLosses = OpenMagnetics::WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);

        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_With_DC) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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

            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);

            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.00204},
                                                                        {100000, 0.00233},
                                                                        {500000, 0.00458}});

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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty_With_DC) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 4.2;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.03783},
                                                                        {100000, 0.03811},
                                                                        {500000, 0.040374}});

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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_90_Duty) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.9;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{25000, 0.00208},
                                                                        {100000, 0.00262},
                                                                        {500000, 0.00513}});

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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
}

SUITE(WindingLossesLitz) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.15;
    bool plot = false;

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Strands) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Triangular_With_DC_Many_Strands) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 10;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{500000, 0.113}});

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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Few_Strands) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Many_Strands) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Litz_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
            if (plot) {
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_Ten_Turns_Litz_Sinusoidal_" + std::to_string(testPoint.first) +".svg");
                std::filesystem::remove(outFile);
                OpenMagnetics::Painter painter(outFile, true);
                settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
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
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({30});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "P 26/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
}

SUITE(WindingLossesRectangular) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.2;

    TEST(Test_Winding_Losses_One_Turn_Rectangular_Sinusoidal_No_Fringing) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 44/22/15";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_width(0.0039);
        wire.set_nominal_value_outer_height(0.0007676);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);
        auto turns = coil.get_turns_description().value();
        turns[0].set_length(0.1);
        coil.set_turns_description(turns);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0002985},
                                                                        {100000, 0.00062751},
                                                                        {200000, 0.00086631},
                                                                        {300000, 0.001052},
                                                                        {400000, 0.0012088},
                                                                        {500000, 0.0013523},
                                                                        {600000, 0.0015184},
                                                                        {700000, 0.0016552},
                                                                        {800000, 0.0017905},
                                                                        {900000, 0.0019245},
                                                                        {1000000, 0.0020574}});

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


            OpenMagnetics::WindingLosses windingLosses;
            settings->set_magnetic_field_include_fringing(false);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Ungapped_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            OpenMagnetics::WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Ungapped_Sinusoidal_7_Amps) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            OpenMagnetics::WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Five_Turns_Rectangular_Gapped_Sinusoidal_7_Amps) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({5});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 20/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.004);
        wire.set_nominal_value_conducting_height(0.0009);
        wire.set_nominal_value_outer_height(0.000909);
        wire.set_nominal_value_outer_width(0.00404);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            OpenMagnetics::WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }
    TEST(Test_Winding_Losses_Seven_Turns_Rectangular_Ungapped_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({7});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 27/17";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0038);
        wire.set_nominal_value_conducting_height(0.00076);
        wire.set_nominal_value_outer_height(0.0007676);
        wire.set_nominal_value_outer_width(0.003838);
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
        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_residual_gap();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);

        auto bobbin = coil.resolve_bobbin();
        auto bobbinProcessedDescription = bobbin.get_processed_description().value();
        double bobbinColumnThickness = 0.001;
        bobbinProcessedDescription.set_wall_thickness(bobbinColumnThickness);
        bobbinProcessedDescription.set_column_thickness(bobbinColumnThickness);

        OpenMagnetics::WindingWindowElement windingWindowElement;

        auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];
        auto coreCentralColumn = core.get_processed_description()->get_columns()[0];
        windingWindowElement.set_width(coreWindingWindow.get_width().value() - bobbinColumnThickness);
        windingWindowElement.set_height(coreWindingWindow.get_height().value() - bobbinColumnThickness * 2);
        windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + bobbinColumnThickness + windingWindowElement.get_width().value() / 2, 0, 0}));
        bobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + bobbinColumnThickness);
        bobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + bobbinColumnThickness);
        bobbinProcessedDescription.set_winding_windows(std::vector<OpenMagnetics::WindingWindowElement>({windingWindowElement}));

        bobbin.set_processed_description(bobbinProcessedDescription);
        coil.set_bobbin(bobbin);
        coil.wind();

        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            OpenMagnetics::WindingLosses windingLosses;
            settings->set_magnetic_field_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
            // auto outFile = outputFilePath;
            // outFile.append("Test_Winding_Losses_Seven_Turns_Rectangular_Ungapped_Sinusoidal_" + std::to_string(testPoint.first) +".svg");
            // std::filesystem::remove(outFile);
            // OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
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
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.3;
    TEST(Test_Winding_Losses_One_Turn_Foil_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.02);
        wire.set_nominal_value_outer_width(0.00011);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Foil_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.02);
        wire.set_nominal_value_outer_width(0.00018);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Short_Turns_Foil_Sinusoidal) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.007);
        wire.set_nominal_value_outer_width(0.00015);
        wire.set_nominal_value_outer_height(0.007);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_ground_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);

        }
        settings->reset();
    }
}

SUITE(WindingLossesToroidalCores) {
    auto settings = OpenMagnetics::Settings::GetInstance();
    double maximumError = 0.15;
    bool plot = false;

    TEST(Test_Winding_Losses_Toroidal_Core_One_Turn_Round_Tendency) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::vector<double> turnsRatios;

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 0;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double frequency = 100000;
        double magnetizingInductance = 1e-3;
        std::string shapeName = "T 20/10/7";

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
                                                                                         offset,
                                                                                         turnsRatios);

        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

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
        auto ohmicLosses100kHz = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);


        CHECK_CLOSE(ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses(), ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0], ohmicLosses100kHz.get_dc_resistance_per_turn().value()[0] * maximumError);
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_ohmic_losses().value().get_losses());
        CHECK(ohmicLosses100kHz.get_winding_losses() > ohmicLosses100kHz.get_winding_losses_per_winding().value()[0].get_skin_effect_losses().value().get_losses_per_harmonic()[1]);

        OpenMagnetics::OperatingPoint scaledOperationPoint = inputs.get_operating_point(0);
        OpenMagnetics::InputsWrapper::scale_time_to_frequency(scaledOperationPoint, frequency * 10);
        scaledOperationPoint = OpenMagnetics::InputsWrapper::process_operating_point(scaledOperationPoint, frequency * 10);
        auto ohmicLosses1MHz = OpenMagnetics::   WindingLosses().calculate_losses(magnetic, scaledOperationPoint, temperature);
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
            OpenMagnetics::Painter painter(outFile);
            OpenMagnetics::MagneticWrapper magnetic;
            magnetic.set_core(core);
            magnetic.set_coil(coil);
            // settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
            settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
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
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Toroidal_Core) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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


            auto windingLosses = OpenMagnetics::WindingLosses();
            settings->set_magnetic_field_mirroring_dimension(1);
            settings->set_magnetic_field_include_fringing(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0058);
        wire.set_nominal_value_conducting_height(0.002);
        wire.set_nominal_value_outer_width(0.0059);
        wire.set_nominal_value_outer_height(0.002076);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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

            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_One_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire" + std::to_string(testPoint.first) + ".svg");
                std::filesystem::remove(outFile);
                OpenMagnetics::Painter painter(outFile);
                OpenMagnetics::MagneticWrapper magnetic;
                magnetic.set_core(core);
                magnetic.set_coil(coil);
                // settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
                settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
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

            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
        }
        settings->reset();
    }

    TEST(Test_Winding_Losses_Ten_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire) {
        settings->reset();
        OpenMagnetics::clear_databases();

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "T 40/24/16";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0058);
        wire.set_nominal_value_conducting_height(0.002);
        wire.set_nominal_value_outer_width(0.0059);
        wire.set_nominal_value_outer_height(0.002076);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
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

            if (plot) {
                auto outputFilePath = std::filesystem::path{ __FILE__ }.parent_path().append("..").append("output");
                auto outFile = outputFilePath;
                outFile.append("Test_Winding_Losses_Ten_Turn_Round_Sinusoidal_Toroidal_Core_Rectangular_Wire" + std::to_string(testPoint.first) + ".svg");
                std::filesystem::remove(outFile);
                OpenMagnetics::Painter painter(outFile);
                OpenMagnetics::MagneticWrapper magnetic;
                magnetic.set_core(core);
                magnetic.set_coil(coil);
                // settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::CONTOUR);
                settings->set_painter_mode(OpenMagnetics::Painter::PainterModes::QUIVER);
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

            auto ohmicLosses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
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
        auto windingOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto layersOrientation = OpenMagnetics::WindingOrientation::OVERLAPPING;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round T21A01TXXX-1");
            wires.push_back(wire);
        }
        {
            OpenMagnetics::WireWrapper wire = OpenMagnetics::find_wire_by_name("Round 0.25 - FIW 6");
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

        auto resistanceMatrixAtFrequency = OpenMagnetics::WindingLosses().calculate_resistance_matrix(magnetic, temperature, frequency);

        CHECK(resistanceMatrixAtFrequency.get_magnitude().size() == magnetic.get_coil().get_functional_description().size());
        for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
            CHECK(resistanceMatrixAtFrequency.get_magnitude()[windingIndex].size() == magnetic.get_coil().get_functional_description().size());
            CHECK(OpenMagnetics::resolve_dimensional_values(resistanceMatrixAtFrequency.get_magnitude()[windingIndex][windingIndex]) > 0);
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

        auto losses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, inputs.get_operating_point(0), 25);

        json mierda;
        to_json(mierda, losses);
        // std::cout << mierda << std::endl;

        CHECK(losses.get_dc_resistance_per_winding().value()[0] > 0);

    }
}