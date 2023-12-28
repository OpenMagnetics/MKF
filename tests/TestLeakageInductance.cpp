#include "LeakageInductance.h"
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


SUITE(LeakageInductance) {
    double maximumError = 0.1;
    std::string filePath = __FILE__;
    auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");


    TEST(Test_Leakage_Inductance_E_0) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({69, 69});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
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
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000352);
        wire.set_number_conductors(25);
        wire.set_type(OpenMagnetics::WireType::LITZ);
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
        auto gapping = OpenMagneticsTesting::get_grinded_gap(2e-5);
        auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto label = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double offset = 0;
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 6.7e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_E_1) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({64, 20});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
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
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.000352);
        wire.set_number_conductors(25);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(0.001056);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 13e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_E_2) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({16, 6});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 42/33/20";
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
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(370);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(370) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.0001);
        strandOuterDiameter.set_nominal(0.00011);
        wire.set_strand(strand);
        wire.set_number_conductors(666);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(666) * 0.00011);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 4e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_E_3) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({36, 26});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "E 65/32/27";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(650);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(650) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.000073);
        wire.set_strand(strand);
        wire.set_number_conductors(700);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(700) * 0.000073);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 8e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_Parallels_No_Interleaving) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({24, 6});
        std::vector<int64_t> numberParallels({2, 4});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 32/25";
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
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(75);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(75) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_strand(strand);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(225) * 0.000055);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 7e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_Parallels_Interleaving) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({24, 6});
        std::vector<int64_t> numberParallels({2, 4});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 32/25";
        uint8_t interleavingLevel = 2;
        auto windingOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(75);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(75) * 0.000055);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_strand(strand);
        wire.set_number_conductors(225);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(225) * 0.000055);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 2e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_ETD_0) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({60, 59});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "ETD 39";
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
        strandOuterDiameter.set_nominal(0.000073);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(135);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(135) * 0.000073);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        wire.set_number_conductors(69);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(69) * 0.000073);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 40e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_PQ_26_0) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({27, 3});
        std::vector<int64_t> numberParallels({1, 1});
        std::vector<double> turnsRatios({double(numberTurns[0]) / numberTurns[1]});
        std::string shapeName = "PQ 26/25";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.000071);
        strandOuterDiameter.set_nominal(0.000073);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_number_conductors(60);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(60) * 0.000073);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);

        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        wire.set_strand(strand);
        wire.set_number_conductors(370);
        wire.set_nominal_value_outer_diameter(1.28 * sqrt(370) * 0.000055);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 100e-6;


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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_PQ_40_Horizontal) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({20, 2});
        std::vector<int64_t> numberParallels({1, 3});
        std::vector<double> turnsRatios({10});
        std::string shapeName = "PQ 40/40";
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
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.0021);
        wire.set_number_conductors(800);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_nominal_value_outer_diameter(0.002469);
        wire.set_number_conductors(1000);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 9.9e-6;

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

        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

    TEST(Test_Leakage_Inductance_PQ_40_Vertical) {

        double temperature = 20;
        std::vector<int64_t> numberTurns({20, 2});
        std::vector<int64_t> numberParallels({1, 3});
        std::vector<double> turnsRatios({10});
        std::string shapeName = "PQ 40/40";
        uint8_t interleavingLevel = 1;
        auto windingOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        auto turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        auto sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        std::vector<OpenMagnetics::WireWrapper> wires;
        OpenMagnetics::WireRound strand;
        OpenMagnetics::WireWrapper wire;
        OpenMagnetics::DimensionWithTolerance strandConductingDiameter;
        OpenMagnetics::DimensionWithTolerance strandOuterDiameter;
        strandConductingDiameter.set_nominal(0.00005);
        strandOuterDiameter.set_nominal(0.000055);
        strand.set_conducting_diameter(strandConductingDiameter);
        strand.set_outer_diameter(strandOuterDiameter);
        strand.set_number_conductors(1);
        strand.set_material("copper");
        strand.set_type(OpenMagnetics::WireType::ROUND);

        wire.set_strand(strand);
        wire.set_nominal_value_outer_diameter(0.0021);
        wire.set_number_conductors(800);
        wire.set_type(OpenMagnetics::WireType::LITZ);
        wires.push_back(wire);
        wire.set_nominal_value_outer_diameter(0.002469);
        wire.set_number_conductors(1000);
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
        double frequency = 100000;
        double peakToPeak = 2;
        double dutyCycle = 0.5;
        double magnetizingInductance = 1e-3;
        double expectedLeakageInductance = 40e-6;

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


        auto leakageInductance = OpenMagnetics::LeakageInductance().calculate_leakage_inductance(inputs.get_operating_point(0), magnetic).get_leakage_inductance_per_winding()[0].get_nominal().value();
        CHECK_CLOSE(expectedLeakageInductance, leakageInductance, expectedLeakageInductance * maximumError);
    }

}
