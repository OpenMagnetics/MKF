#include "Painter.h"
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


std::string filePath = __FILE__;
auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");

SUITE(WindingLossesRound) {
    double maximumError = 0.15;

    TEST(Test_Winding_Losses_One_Turn_Round_Tendency) {

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
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.05e-3);
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
                                                                      {200000, 0.05962},
                                                                      {250000, 0.06861},
                                                                      {500000, 0.103}});

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
            windingLosses.set_mirroring_dimension(1);
            windingLosses.set_fringing_effect(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);
            // auto outFile = outputFilePath;
            // outFile.append("Test_Winding_Losses_Ten_Turns_Round_Sinusoidal" + std::to_string(testPoint.first) +".svg");
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
    }

    TEST(Test_Winding_Losses_Ten_Turns_Round_Sinusoidal_Interleaving) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({20, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 2;
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.01e-3);
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
                                                                      {50000, 0.11888},
                                                                      {100000, 0.132},
                                                                      {200000, 0.26954},
                                                                      {250000, 0.30431},
                                                                      {500000, 0.43865}});

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
            windingLosses.set_mirroring_dimension(1);
            windingLosses.set_fringing_effect(true);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);

        }
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Sinusoidal_With_DC) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_50_Duty_With_DC) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Round_Triangular_90_Duty) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
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
    }
}
SUITE(WindingLossesLitz) {
    double maximumError = 0.15;

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Strands) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Triangular_With_DC_Many_Strands) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double offset = 10;
        double peakToPeak = 2 * 1.73205;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{500000, 0.11112}});

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
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Few_Strands) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

    TEST(Test_Winding_Losses_One_Turn_Litz_Sinusoidal_Many_Many_Strands) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
                                                                        {200000, 0.0001616},
                                                                        {250000, 0.0001647},
                                                                        {500000, 0.0001824}});

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
    }

    TEST(Test_Winding_Losses_Ten_Turns_Litz_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
        }
    }

    TEST(Test_Winding_Losses_Thirty_Turns_Litz_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({30});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "P 26/16";
        uint8_t interleavingLevel = 1;
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
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
    }

}
SUITE(WindingLossesRectangular) {
    double maximumError = 0.2;

    TEST(Test_Winding_Losses_Twenty_Turns_Thin_Rectangular_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({20});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 44/22/15";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.00709);
        wire.set_nominal_value_conducting_height(0.0007);
        wire.set_nominal_value_outer_width(0.0072);
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

        int64_t numberStacks = 1;
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0052798},
                                                                        {25000, 0.0249616},
                                                                        {50000, 0.10666},
                                                                        {100000, 0.1748},
                                                                        {200000, 0.24893},
                                                                        {250000, 0.27669},
                                                                        {500000, 0.45924}});

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
    }
    TEST(Test_Winding_Losses_Seven_Turns_Rectangular_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({7});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 27/17";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.00228);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0010804},
                                                                        {100000, 0.052998},
                                                                        {200000, 0.075264},
                                                                        {300000, 0.091907},
                                                                        {400000, 0.10598},
                                                                        {500000, 0.11843},
                                                                        {600000, 0.12968},
                                                                        {700000, 0.14002},
                                                                        {800000, 0.14963},
                                                                        {900000, 0.15864},
                                                                        {1000000, 0.16715}});

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
            // auto outFile = outputFilePath;
            // outFile.append("Test_Winding_Losses_Seven_Turns_Rectangular_Sinusoidal" + std::to_string(testPoint.first) +".svg");
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
    }
    TEST(Test_Winding_Losses_Seven_Turns_Rectangular_Ungapped_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({7});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "PQ 27/17";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.0010804},
                                                                        {100000, 0.02606},
                                                                        {200000, 0.035916},
                                                                        {300000, 0.04906},
                                                                        {400000, 0.061195},
                                                                        {500000, 0.072619},
                                                                        {600000, 0.08396},
                                                                        {700000, 0.09393},
                                                                        {800000, 0.10399},
                                                                        {900000, 0.11373},
                                                                        {1000000, 0.12319}});

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
            windingLosses.set_mirroring_dimension(2);
            auto ohmicLosses = windingLosses.calculate_losses(magnetic, inputs.get_operating_point(0), temperature);
            CHECK_CLOSE(testPoint.second, ohmicLosses.get_winding_losses(), testPoint.second * maximumError);

        }
    }
}
SUITE(WindingLossesFoil) {
    double maximumError = 0.3;
    TEST(Test_Winding_Losses_One_Turn_Foil_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({1});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
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
                                                                        {25000, 0.00049609},
                                                                        {50000, 0.00055119},
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
    }

    TEST(Test_Winding_Losses_Ten_Turns_Foil_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
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
                                                                        {25000, 0.0099166},
                                                                        {50000, 0.013048},
                                                                        {100000, 0.017562},
                                                                        {200000, 0.024423},
                                                                        {250000, 0.027482},
                                                                        {500000, 0.042515}});

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
    }

    TEST(Test_Winding_Losses_Ten_Short_Turns_Foil_Sinusoidal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({10});
        std::vector<int64_t> numberParallels({1});
        std::string shapeName = "ETD 34/17/11";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_conducting_width(0.0001);
        wire.set_nominal_value_conducting_height(0.007);
        wire.set_nominal_value_outer_width(0.00011);
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double peakToPeak = 2 * 1.4142;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;

        std::vector<std::pair<double, double>> expectedWindingLosses({{0.01, 0.011868},
                                                                        {25000, 0.015571},
                                                                        {50000, 0.018916},
                                                                        {100000, 0.021},
                                                                        {200000, 0.022},
                                                                        {250000, 0.024},
                                                                        {500000, 0.03}});

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
    }

}